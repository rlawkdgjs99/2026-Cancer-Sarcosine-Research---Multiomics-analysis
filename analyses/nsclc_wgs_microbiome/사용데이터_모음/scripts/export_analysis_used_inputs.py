#!/usr/bin/env python3
"""Export the exact NSCLC WGS R-vs-NR inputs used by the frozen workflow.

The four protected HGMT cohort folders are read only. New CSV files are written
only below ``사용데이터_모음``. The exporter mirrors the cohort scripts:

* exactly one ``selected_project_*.txt`` per cohort;
* Assay type == WGS;
* Phenotype name == Carcinoma, Non-Small-Cell Lung;
* response_group (preferred) or response parsed as R/NR;
* species rows containing ``|s__`` and not ``|t__``;
* all KEGG-KO rows for selected Run IDs, zero-filled after widening;
* saved per-cohort sarcosine scores copied into one separate, tidy table.

Existing CSV outputs are never overwritten unless ``--force`` is explicitly
provided. Even with ``--force``, paths outside the new output folder are never
opened for writing.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "사용데이터_모음"
NSCLC_LABEL = "Carcinoma, Non-Small-Cell Lung"


@dataclass(frozen=True)
class CohortSpec:
    folder: str
    role: str


COHORTS: Dict[str, CohortSpec] = {
    "PRJNA751792": CohortSpec("NSCLC_PRJNA751792", "Discovery"),
    "PRJNA1023797": CohortSpec("NSCLC_PRJNA1023797", "Discovery"),
    "PRJEB22863": CohortSpec("NSCLC_RCC_PRJEB22863", "Discovery"),
    "PRJEB26531": CohortSpec("NSCLC_PRJEB26531", "External_validation"),
}

DEGRADATION_KOS = (
    "K00301", "K00302", "K00303", "K00304", "K00305", "K00306",
    "K00314", "K18897",
)
PRODUCTION_KOS = ("K08688", "K08687", "K00315")
BIDIRECTIONAL_KOS = ("K21833", "K21834")


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block_size):
            digest.update(chunk)
    return digest.hexdigest()


def relative_to_root(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


def atomic_csv_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    handle = tmp.open("w", newline="", encoding="utf-8")
    return handle, csv.writer(handle), tmp


def finalize_atomic(handle, tmp: Path, path: Path) -> None:
    handle.flush()
    os.fsync(handle.fileno())
    handle.close()
    os.replace(tmp, path)


def write_dict_rows_atomic(path: Path, rows: Sequence[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    handle, _, tmp = atomic_csv_writer(path)
    try:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        finalize_atomic(handle, tmp, path)
    except Exception:
        handle.close()
        if tmp.exists():
            tmp.unlink()
        raise


def parse_response(description: str) -> Optional[str]:
    text = description or ""
    match = re.search(r"response_group:[^;]*", text)
    if match is None:
        match = re.search(r"response:[^;]*", text)
    if match is None:
        return None
    value = re.sub(r"^response(?:_group)?:\s*", "", match.group(0)).strip()
    return value if value in {"R", "NR"} else None


def unique_source(cohort_dir: Path, pattern: str) -> Path:
    files = sorted(cohort_dir.glob(pattern))
    if len(files) != 1:
        raise ValueError(
            f"Expected exactly one {pattern} in {cohort_dir}; found {len(files)}"
        )
    return files[0]


def load_selected_metadata(cohort: str, spec: CohortSpec):
    cohort_dir = ROOT / spec.folder
    source = unique_source(cohort_dir, "selected_project_*.txt")
    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            preamble = next(reader)
            headers = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Incomplete metadata file: {source}") from exc
        if not preamble or not headers:
            raise ValueError(f"Invalid metadata preamble/header: {source}")
        rows = []
        for line_number, fields in enumerate(reader, start=3):
            if not fields or not any(fields):
                continue
            if len(fields) < len(headers):
                raise ValueError(
                    f"Metadata row {line_number} has {len(fields)} fields, "
                    f"expected at least {len(headers)}: {source}"
                )
            rows.append(dict(zip(headers, fields[: len(headers)])))

    required = {"Run ID", "Sample name", "Assay type", "Phenotype name", "Sample description"}
    missing = required.difference(headers)
    if missing:
        raise ValueError(f"Missing metadata columns {sorted(missing)} in {source}")

    selected = []
    for row in rows:
        group = parse_response(row["Sample description"])
        if (
            row["Assay type"] == "WGS"
            and row["Phenotype name"] == NSCLC_LABEL
            and group is not None
        ):
            out = dict(row)
            out["Analysis Group"] = group
            out["Cohort"] = cohort
            out["Cohort Role"] = spec.role
            selected.append(out)

    run_ids = [row["Run ID"] for row in selected]
    if not selected:
        raise ValueError(f"No WGS NSCLC R/NR records selected for {cohort}")
    if len(run_ids) != len(set(run_ids)):
        duplicates = [run for run, count in Counter(run_ids).items() if count > 1]
        raise ValueError(f"Duplicate selected Run IDs for {cohort}: {duplicates}")
    return source, headers, selected, run_ids


def iter_json_array(path: Path, chunk_size: int = 8 * 1024 * 1024) -> Iterator[dict]:
    """Stream a top-level JSON array while preserving numeric token text."""
    decoder = json.JSONDecoder(parse_float=str, parse_int=str)
    with path.open("r", encoding="utf-8") as handle:
        buffer = ""
        position = 0
        started = False
        eof = False
        while True:
            if not eof and len(buffer) - position < chunk_size:
                if position:
                    buffer = buffer[position:]
                    position = 0
                chunk = handle.read(chunk_size)
                if chunk:
                    buffer += chunk
                else:
                    eof = True

            if not started:
                while position < len(buffer) and buffer[position].isspace():
                    position += 1
                if position >= len(buffer):
                    if eof:
                        raise ValueError(f"Empty JSON file: {path}")
                    continue
                if buffer[position] != "[":
                    raise ValueError(f"Expected top-level JSON array in {path}")
                position += 1
                started = True
                continue

            while position < len(buffer) and (
                buffer[position].isspace() or buffer[position] == ","
            ):
                position += 1
            if position < len(buffer) and buffer[position] == "]":
                return
            if position >= len(buffer):
                if eof:
                    raise ValueError(f"Unexpected end of JSON array in {path}")
                continue
            try:
                item, end = decoder.raw_decode(buffer, position)
            except json.JSONDecodeError:
                if eof:
                    raise
                if position:
                    buffer = buffer[position:]
                    position = 0
                chunk = handle.read(chunk_size)
                if chunk:
                    buffer += chunk
                else:
                    eof = True
                continue
            if not isinstance(item, dict):
                raise ValueError(f"Expected JSON objects in {path}")
            yield item
            position = end


def scan_grouped_table(
    source: Path,
    allowed_runs: set[str],
    include_row: Callable[[dict], bool],
) -> Tuple[List[str], List[str], int]:
    sample_order: List[str] = []
    feature_order: List[str] = []
    seen_samples: set[str] = set()
    seen_features: set[str] = set()
    current_run = None
    current_features: set[str] = set()
    retained = 0

    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Taxa", "Run ID", "Abundance"]:
            raise ValueError(f"Unexpected species header in {source}: {reader.fieldnames}")
        for row in reader:
            run_id = row["Run ID"]
            if not include_row(row):
                continue
            feature = row["Taxa"]
            # The frozen R read_species() widens the complete Bacteria file first
            # and only then subsets rows to the selected NSCLC R/NR Run IDs.
            # Preserve even features that are all-zero in the selected samples.
            if feature not in seen_features:
                seen_features.add(feature)
                feature_order.append(feature)
            if run_id not in allowed_runs:
                continue
            if run_id != current_run:
                if run_id in seen_samples:
                    raise ValueError(f"Non-contiguous Run ID {run_id} in {source}")
                seen_samples.add(run_id)
                sample_order.append(run_id)
                current_run = run_id
                current_features = set()
            if feature in current_features:
                raise ValueError(f"Duplicate species cell: {run_id}, {feature}")
            current_features.add(feature)
            retained += 1
    return sample_order, feature_order, retained


def scan_grouped_json(source: Path, allowed_runs: set[str]):
    sample_order: List[str] = []
    feature_order: List[str] = []
    seen_samples: set[str] = set()
    seen_features: set[str] = set()
    current_run = None
    current_features: set[str] = set()
    retained = 0
    for row in iter_json_array(source):
        run_id = row.get("run_id")
        if run_id not in allowed_runs:
            continue
        feature = row.get("ko")
        abundance = row.get("abundance")
        if not isinstance(feature, str) or abundance is None:
            raise ValueError(f"Malformed KO row in {source}: {row}")
        if run_id != current_run:
            if run_id in seen_samples:
                raise ValueError(f"Non-contiguous Run ID {run_id} in {source}")
            seen_samples.add(run_id)
            sample_order.append(run_id)
            current_run = run_id
            current_features = set()
        if feature in current_features:
            raise ValueError(f"Duplicate KO cell: {run_id}, {feature}")
        current_features.add(feature)
        if feature not in seen_features:
            seen_features.add(feature)
            feature_order.append(feature)
        retained += 1
    return sample_order, feature_order, retained


def require_source_order(selected_runs: Sequence[str], source_runs: Sequence[str], label: str) -> None:
    expected = [run for run in selected_runs if run in set(source_runs)]
    if expected != list(source_runs):
        raise ValueError(f"{label} source order does not follow selected metadata order")


def export_species_matrix(
    source: Path,
    output: Path,
    allowed_runs: set[str],
    sample_order: Sequence[str],
    feature_order: Sequence[str],
) -> None:
    include_row = lambda row: "|s__" in row["Taxa"] and "|t__" not in row["Taxa"]
    feature_index = {feature: index for index, feature in enumerate(feature_order)}
    expected_index = 0
    current_run = None
    values: Dict[int, str] = {}
    handle, writer, tmp = atomic_csv_writer(output)
    try:
        writer.writerow(["Run_ID", *feature_order])

        def flush() -> None:
            nonlocal expected_index, current_run, values
            if current_run is None:
                return
            if expected_index >= len(sample_order) or current_run != sample_order[expected_index]:
                raise ValueError(f"Species run-order mismatch at {current_run}")
            row_values = ["0"] * len(feature_order)
            for index, value in values.items():
                row_values[index] = value
            writer.writerow([current_run, *row_values])
            expected_index += 1

        with source.open("r", encoding="utf-8", newline="") as source_handle:
            for row in csv.DictReader(source_handle, delimiter="\t"):
                run_id = row["Run ID"]
                if run_id not in allowed_runs or not include_row(row):
                    continue
                if run_id != current_run:
                    flush()
                    current_run = run_id
                    values = {}
                index = feature_index[row["Taxa"]]
                if index in values:
                    raise ValueError(f"Duplicate species cell while exporting: {run_id}, {row['Taxa']}")
                values[index] = row["Abundance"]
        flush()
        if expected_index != len(sample_order):
            raise ValueError(f"Species rows written {expected_index}, expected {len(sample_order)}")
        finalize_atomic(handle, tmp, output)
    except Exception:
        handle.close()
        if tmp.exists():
            tmp.unlink()
        raise


def export_ko_matrix(
    source: Path,
    output: Path,
    allowed_runs: set[str],
    sample_order: Sequence[str],
    feature_order: Sequence[str],
) -> None:
    feature_index = {feature: index for index, feature in enumerate(feature_order)}
    expected_index = 0
    current_run = None
    values: Dict[int, str] = {}
    handle, writer, tmp = atomic_csv_writer(output)
    try:
        writer.writerow(["Run_ID", *feature_order])

        def flush() -> None:
            nonlocal expected_index, current_run, values
            if current_run is None:
                return
            if expected_index >= len(sample_order) or current_run != sample_order[expected_index]:
                raise ValueError(f"KO run-order mismatch at {current_run}")
            row_values = ["0"] * len(feature_order)
            for index, value in values.items():
                row_values[index] = value
            writer.writerow([current_run, *row_values])
            expected_index += 1

        for row in iter_json_array(source):
            run_id = row["run_id"]
            if run_id not in allowed_runs:
                continue
            if run_id != current_run:
                flush()
                current_run = run_id
                values = {}
            index = feature_index[row["ko"]]
            if index in values:
                raise ValueError(f"Duplicate KO cell while exporting: {run_id}, {row['ko']}")
            values[index] = row["abundance"]
        flush()
        if expected_index != len(sample_order):
            raise ValueError(f"KO rows written {expected_index}, expected {len(sample_order)}")
        finalize_atomic(handle, tmp, output)
    except Exception:
        handle.close()
        if tmp.exists():
            tmp.unlink()
        raise


def export_metadata(
    output: Path,
    headers: Sequence[str],
    selected_rows: Sequence[dict],
    species_samples: set[str],
    ko_samples: set[str],
) -> None:
    extras = [
        "Analysis Group", "Cohort", "Cohort Role",
        "Included in species matrix", "Included in KO matrix",
    ]
    handle, writer, tmp = atomic_csv_writer(output)
    try:
        writer.writerow([*headers, *extras])
        for row in selected_rows:
            run_id = row["Run ID"]
            writer.writerow(
                [row.get(header, "") for header in headers]
                + [
                    row["Analysis Group"], row["Cohort"], row["Cohort Role"],
                    "TRUE" if run_id in species_samples else "FALSE",
                    "TRUE" if run_id in ko_samples else "FALSE",
                ]
            )
        finalize_atomic(handle, tmp, output)
    except Exception:
        handle.close()
        if tmp.exists():
            tmp.unlink()
        raise


def output_paths() -> List[Path]:
    paths: List[Path] = []
    for cohort in COHORTS:
        folder = OUTPUT / cohort
        paths.extend(
            [
                folder / f"{cohort}_patient_metadata_WGS_NSCLC_R_NR.csv",
                folder / f"{cohort}_species_relative_abundance_matrix.csv",
                folder / f"{cohort}_KEGG_KO_relative_abundance_matrix.csv",
            ]
        )
    paths.extend(
        [
            OUTPUT / "Derived_Sarcosine_Functional_Scores" / "NSCLC_WGS_sarcosine_functional_scores_all_cohorts.csv",
            OUTPUT / "SOURCE_MANIFEST.csv",
            OUTPUT / "VALIDATION_SUMMARY.csv",
        ]
    )
    return paths


def ensure_output_policy(force: bool) -> None:
    existing = [path for path in output_paths() if path.exists()]
    if existing and not force:
        joined = "\n".join(str(path) for path in existing)
        raise FileExistsError(
            "Existing output files would be replaced. Re-run with --force only after review:\n"
            + joined
        )


def export_combined_scores(cohort_state: Dict[str, dict]) -> Tuple[Path, List[Path], int]:
    output = (
        OUTPUT / "Derived_Sarcosine_Functional_Scores"
        / "NSCLC_WGS_sarcosine_functional_scores_all_cohorts.csv"
    )
    source_files: List[Path] = []
    total = 0
    handle, writer, tmp = atomic_csv_writer(output)
    try:
        writer.writerow(
            [
                "Run_ID", "Cohort", "Cohort_Role", "Analysis_Group",
                "Degradation_sum", "Production_sum",
                "Log2_Production_Degradation_Ratio",
            ]
        )
        for cohort, spec in COHORTS.items():
            state = cohort_state[cohort]
            score_source = ROOT / spec.folder / "results" / f"sarcosine_scores_{cohort}.csv"
            source_files.append(score_source)
            if not score_source.exists():
                raise FileNotFoundError(f"Missing frozen score output: {score_source}")
            with score_source.open("r", encoding="utf-8", newline="") as source_handle:
                reader = csv.DictReader(source_handle)
                required = {
                    "RunID", "group", "degradation", "production", "prod_deg_log2ratio"
                }
                if required.difference(reader.fieldnames or []):
                    raise ValueError(f"Unexpected score schema: {score_source}")
                rows = list(reader)
            runs = [row["RunID"] for row in rows]
            if runs != state["selected_runs"]:
                raise ValueError(f"Frozen score Run-ID order mismatch for {cohort}")
            groups = {row["Run ID"]: row["Analysis Group"] for row in state["selected_rows"]}
            for row in rows:
                if row["group"] != groups[row["RunID"]]:
                    raise ValueError(f"Frozen score group mismatch: {cohort}, {row['RunID']}")
                writer.writerow(
                    [
                        row["RunID"], cohort, spec.role, row["group"],
                        row["degradation"], row["production"], row["prod_deg_log2ratio"],
                    ]
                )
                total += 1
        finalize_atomic(handle, tmp, output)
    except Exception:
        handle.close()
        if tmp.exists():
            tmp.unlink()
        raise
    return output, source_files, total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Atomically replace only existing CSVs under the new output folder")
    args = parser.parse_args()
    ensure_output_policy(args.force)
    OUTPUT.mkdir(parents=True, exist_ok=True)

    source_hash_before: Dict[Path, str] = {}
    manifest_rows: List[dict] = []
    validation_rows: List[dict] = []
    cohort_state: Dict[str, dict] = {}

    for cohort, spec in COHORTS.items():
        cohort_dir = ROOT / spec.folder
        metadata_source, metadata_headers, selected_rows, selected_runs = load_selected_metadata(cohort, spec)
        allowed_runs = set(selected_runs)
        bacteria_source = unique_source(cohort_dir, "Bacteria_*.txt")
        ko_source = cohort_dir / "KO_relative_abundance.tsv"
        analysis_script = cohort_dir / "R_scripts" / "run_analysis.R"
        score_script = cohort_dir / "R_scripts" / "run_sarcosine.R"
        for source in (metadata_source, bacteria_source, ko_source, analysis_script, score_script):
            source_hash_before[source] = sha256_file(source)

        species_filter = lambda row: "|s__" in row["Taxa"] and "|t__" not in row["Taxa"]
        species_runs, species_features, species_records = scan_grouped_table(
            bacteria_source, allowed_runs, species_filter
        )
        ko_runs, ko_features, ko_records = scan_grouped_json(ko_source, allowed_runs)
        require_source_order(selected_runs, species_runs, f"{cohort} species")
        require_source_order(selected_runs, ko_runs, f"{cohort} KO")

        output_dir = OUTPUT / cohort
        metadata_output = output_dir / f"{cohort}_patient_metadata_WGS_NSCLC_R_NR.csv"
        species_output = output_dir / f"{cohort}_species_relative_abundance_matrix.csv"
        ko_output = output_dir / f"{cohort}_KEGG_KO_relative_abundance_matrix.csv"
        export_metadata(
            metadata_output, metadata_headers, selected_rows,
            set(species_runs), set(ko_runs),
        )
        export_species_matrix(
            bacteria_source, species_output, allowed_runs, species_runs, species_features
        )
        export_ko_matrix(ko_source, ko_output, allowed_runs, ko_runs, ko_features)

        counts = Counter(row["Analysis Group"] for row in selected_rows)
        missing_species = [run for run in selected_runs if run not in set(species_runs)]
        missing_ko = [run for run in selected_runs if run not in set(ko_runs)]
        validation_rows.append(
            {
                "Cohort": cohort,
                "Cohort_role": spec.role,
                "R_metadata_n": counts["R"],
                "NR_metadata_n": counts["NR"],
                "Metadata_n": len(selected_rows),
                "Unique_sample_names_within_cohort_n": len({row["Sample name"] for row in selected_rows}),
                "Sample_names_also_seen_in_other_cohorts_n": 0,
                "Global_unique_sample_names_across_all_cohorts_n": 0,
                "Species_matrix_samples_n": len(species_runs),
                "Species_features_n": len(species_features),
                "Species_source_records_n": species_records,
                "KO_matrix_samples_n": len(ko_runs),
                "KO_features_n": len(ko_features),
                "KO_source_records_n": ko_records,
                "Missing_from_species_matrix": ";".join(missing_species),
                "Missing_from_KO_matrix": ";".join(missing_ko),
                "Selection_rule": "WGS; NSCLC phenotype; parsed R/NR",
            }
        )

        output_records = [
            ("patient_metadata", metadata_output, metadata_source, len(selected_rows), len(metadata_headers) + 5),
            ("species_matrix", species_output, bacteria_source, len(species_runs), len(species_features) + 1),
            ("KEGG_KO_matrix", ko_output, ko_source, len(ko_runs), len(ko_features) + 1),
        ]
        for data_type, output_file, source_file, rows, columns in output_records:
            manifest_rows.append(
                {
                    "Cohort": cohort,
                    "Cohort_role": spec.role,
                    "Data_type": data_type,
                    "Output_file": relative_to_root(output_file),
                    "Data_rows_excluding_header": rows,
                    "Columns": columns,
                    "Source_file": relative_to_root(source_file),
                    "Source_SHA256": source_hash_before[source_file],
                    "Selection_script": relative_to_root(analysis_script),
                    "Selection_script_SHA256": source_hash_before[analysis_script],
                    "Score_script": relative_to_root(score_script),
                    "Score_script_SHA256": source_hash_before[score_script],
                    "Output_SHA256": sha256_file(output_file),
                }
            )

        cohort_state[cohort] = {
            "selected_rows": selected_rows,
            "selected_runs": selected_runs,
            "ko_source": ko_source,
            "ko_output": ko_output,
        }
        print(
            f"{cohort}: metadata={len(selected_rows)} (R={counts['R']}, NR={counts['NR']}), "
            f"species={len(species_runs)}x{len(species_features)}, "
            f"KO={len(ko_runs)}x{len(ko_features)}",
            flush=True,
        )

    sample_name_cohorts: Dict[str, set[str]] = {}
    for cohort, state in cohort_state.items():
        for row in state["selected_rows"]:
            sample_name_cohorts.setdefault(row["Sample name"], set()).add(cohort)
    global_unique_sample_names = len(sample_name_cohorts)
    for validation_row in validation_rows:
        cohort = validation_row["Cohort"]
        cohort_names = {
            row["Sample name"] for row in cohort_state[cohort]["selected_rows"]
        }
        validation_row["Sample_names_also_seen_in_other_cohorts_n"] = sum(
            len(sample_name_cohorts[name]) > 1 for name in cohort_names
        )
        validation_row["Global_unique_sample_names_across_all_cohorts_n"] = (
            global_unique_sample_names
        )

    score_output, score_sources, score_rows = export_combined_scores(cohort_state)
    for path in score_sources:
        source_hash_before[path] = sha256_file(path)
    manifest_rows.append(
        {
            "Cohort": "ALL",
            "Cohort_role": "3 discovery + 1 external validation",
            "Data_type": "derived_sarcosine_functional_scores",
            "Output_file": relative_to_root(score_output),
            "Data_rows_excluding_header": score_rows,
            "Columns": 7,
            "Source_file": ";".join(relative_to_root(path) for path in score_sources),
            "Source_SHA256": ";".join(source_hash_before[path] for path in score_sources),
            "Selection_script": "per-cohort run_sarcosine.R",
            "Selection_script_SHA256": "see per-cohort manifest rows",
            "Score_script": "per-cohort run_sarcosine.R",
            "Score_script_SHA256": "see per-cohort manifest rows",
            "Output_SHA256": sha256_file(score_output),
        }
    )

    write_dict_rows_atomic(OUTPUT / "VALIDATION_SUMMARY.csv", validation_rows)
    write_dict_rows_atomic(OUTPUT / "SOURCE_MANIFEST.csv", manifest_rows)

    changed = []
    for source, before in source_hash_before.items():
        after = sha256_file(source)
        if before != after:
            changed.append(relative_to_root(source))
    if changed:
        raise RuntimeError(f"Protected source files changed during export: {changed}")
    print(f"Combined score rows: {score_rows}; protected source hashes unchanged", flush=True)


if __name__ == "__main__":
    main()
