#!/usr/bin/env python3
"""Export the exact per-cohort CRC WGS inputs used by integrated_analysis_pooled.R.

The exporter deliberately mirrors the frozen R workflow:
  * most recent selected_project_*.txt by descending filename;
  * Health plus the first non-Health/non-adenoma phenotype;
  * WGS only;
  * species rows containing ``|s__`` but not ``|t__``;
  * KO rows restricted to the selected Run IDs;
  * zero-filled wide matrices in source first-occurrence order.

Raw source files are read only. Output is written atomically into cohort subfolders.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[2]   # module root (this script sits in <module>/사용데이터_모음/scripts/)
INTEGRATED = ROOT / "Healthy_vs_Cancer_4_CRC_cohorts_integrated"
OUTPUT = ROOT / "사용데이터_모음"

COHORTS = {
    "PRJEB6070": "PRJEB6070_CRC_AdenomatousPolyps",
    "PRJNA429097": "PRJNA429097_CRC",
    "PRJEB10878": "PRJEB10878_CRC",
    "PRJEB27928": "PRJEB27928_CRC",
}


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block_size):
            digest.update(chunk)
    return digest.hexdigest()


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


def load_selected_metadata(cohort: str, cohort_dir: Path):
    files = sorted(cohort_dir.glob("selected_project_*.txt"), reverse=True)
    if not files:
        raise FileNotFoundError(f"No selected_project_*.txt in {cohort_dir}")
    source = files[0]

    with source.open("r", encoding="utf-8", newline="") as handle:
        first_line = handle.readline()
        if not first_line:
            raise ValueError(f"Empty metadata file: {source}")
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        headers = list(reader.fieldnames or [])

    phenotype_column = next(
        (name for name in headers if name.replace(" ", ".") == "Phenotype.name"),
        None,
    )
    if phenotype_column is None:
        raise ValueError(f"Phenotype.name column not found in {source}")
    if "Run ID" not in headers or "Assay type" not in headers:
        raise ValueError(f"Required metadata columns missing in {source}")

    phenotype_order: List[str] = []
    for row in rows:
        value = row.get(phenotype_column, "")
        if value and value not in phenotype_order:
            phenotype_order.append(value)
    cancer_candidates = [
        value
        for value in phenotype_order
        if value not in {"Health", "Adenomatous Polyps"}
    ]
    if not cancer_candidates:
        raise ValueError(f"No cancer phenotype detected in {source}")
    cancer_label = cancer_candidates[0]

    selected = [
        row
        for row in rows
        if row.get(phenotype_column) in {"Health", cancer_label}
        and row.get("Assay type") == "WGS"
    ]
    run_ids = [row["Run ID"] for row in selected]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError(f"Duplicate Run IDs in selected metadata for {cohort}")

    for row in selected:
        row["Analysis Group"] = (
            "Healthy" if row[phenotype_column] == "Health" else "Cancer"
        )
        row["Cohort"] = cohort

    return source, headers, phenotype_column, cancer_label, selected, run_ids


def scan_grouped_table(
    source: Path,
    delimiter: str,
    allowed_runs: set[str],
    run_column: str,
    feature_column: str,
    include_row,
) -> Tuple[List[str], List[str], int, int]:
    """Return sample order, feature order, retained rows, and duplicate cells."""
    sample_order: List[str] = []
    feature_order: List[str] = []
    seen_samples: set[str] = set()
    seen_features: set[str] = set()
    current_run = None
    current_features: set[str] = set()
    retained = 0
    duplicates = 0

    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            run_id = row[run_column]
            if run_id not in allowed_runs or not include_row(row):
                continue
            feature = row[feature_column]
            if run_id != current_run:
                if run_id in seen_samples:
                    raise ValueError(
                        f"Non-contiguous repeated Run ID {run_id} in {source}"
                    )
                seen_samples.add(run_id)
                sample_order.append(run_id)
                current_run = run_id
                current_features = set()
            if feature in current_features:
                duplicates += 1
            current_features.add(feature)
            if feature not in seen_features:
                seen_features.add(feature)
                feature_order.append(feature)
            retained += 1

    return sample_order, feature_order, retained, duplicates


def export_grouped_wide_csv(
    source: Path,
    delimiter: str,
    output: Path,
    allowed_runs: set[str],
    run_column: str,
    feature_column: str,
    value_column: str,
    include_row,
    sample_order: Sequence[str],
    feature_order: Sequence[str],
) -> int:
    feature_index = {feature: index for index, feature in enumerate(feature_order)}
    expected_sample_index = 0
    current_run = None
    values: Dict[int, str] = {}
    written = 0

    handle, writer, tmp = atomic_csv_writer(output)
    try:
        writer.writerow(["Run_ID", *feature_order])

        def flush_row() -> None:
            nonlocal expected_sample_index, written, values, current_run
            if current_run is None:
                return
            if expected_sample_index >= len(sample_order):
                raise ValueError(f"Unexpected extra sample {current_run} in {source}")
            expected = sample_order[expected_sample_index]
            if current_run != expected:
                raise ValueError(
                    f"Sample-order mismatch in {source}: expected {expected}, got {current_run}"
                )
            row_values = ["0"] * len(feature_order)
            for index, value in values.items():
                row_values[index] = value
            writer.writerow([current_run, *row_values])
            expected_sample_index += 1
            written += 1

        with source.open("r", encoding="utf-8", newline="") as source_handle:
            reader = csv.DictReader(source_handle, delimiter=delimiter)
            for row in reader:
                run_id = row[run_column]
                if run_id not in allowed_runs or not include_row(row):
                    continue
                if run_id != current_run:
                    flush_row()
                    current_run = run_id
                    values = {}
                feature = row[feature_column]
                index = feature_index[feature]
                if index in values:
                    raise ValueError(
                        f"Duplicate ({run_id}, {feature}) cell in {source}"
                    )
                values[index] = row[value_column]
        flush_row()
        if written != len(sample_order):
            raise ValueError(
                f"Expected {len(sample_order)} samples but wrote {written} from {source}"
            )
        finalize_atomic(handle, tmp, output)
    except Exception:
        handle.close()
        if tmp.exists():
            tmp.unlink()
        raise
    return written


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


def scan_grouped_json(
    source: Path, allowed_runs: set[str]
) -> Tuple[List[str], List[str], int, int]:
    sample_order: List[str] = []
    feature_order: List[str] = []
    seen_samples: set[str] = set()
    seen_features: set[str] = set()
    current_run = None
    current_features: set[str] = set()
    retained = 0
    duplicates = 0

    for row in iter_json_array(source):
        run_id = row["run_id"]
        if run_id not in allowed_runs:
            continue
        feature = row["ko"]
        if run_id != current_run:
            if run_id in seen_samples:
                raise ValueError(f"Non-contiguous repeated Run ID {run_id} in {source}")
            seen_samples.add(run_id)
            sample_order.append(run_id)
            current_run = run_id
            current_features = set()
        if feature in current_features:
            duplicates += 1
        current_features.add(feature)
        if feature not in seen_features:
            seen_features.add(feature)
            feature_order.append(feature)
        retained += 1
    return sample_order, feature_order, retained, duplicates


def export_grouped_json_wide_csv(
    source: Path,
    output: Path,
    allowed_runs: set[str],
    sample_order: Sequence[str],
    feature_order: Sequence[str],
) -> int:
    feature_index = {feature: index for index, feature in enumerate(feature_order)}
    expected_sample_index = 0
    current_run = None
    values: Dict[int, str] = {}
    written = 0

    handle, writer, tmp = atomic_csv_writer(output)
    try:
        writer.writerow(["Run_ID", *feature_order])

        def flush_row() -> None:
            nonlocal expected_sample_index, written, values, current_run
            if current_run is None:
                return
            expected = sample_order[expected_sample_index]
            if current_run != expected:
                raise ValueError(
                    f"KO sample-order mismatch: expected {expected}, got {current_run}"
                )
            row_values = ["0"] * len(feature_order)
            for index, value in values.items():
                row_values[index] = value
            writer.writerow([current_run, *row_values])
            expected_sample_index += 1
            written += 1

        for row in iter_json_array(source):
            run_id = row["run_id"]
            if run_id not in allowed_runs:
                continue
            if run_id != current_run:
                flush_row()
                current_run = run_id
                values = {}
            feature = row["ko"]
            index = feature_index[feature]
            if index in values:
                raise ValueError(f"Duplicate ({run_id}, {feature}) cell in {source}")
            values[index] = row["abundance"]
        flush_row()
        if written != len(sample_order):
            raise ValueError(
                f"Expected {len(sample_order)} KO samples but wrote {written}"
            )
        finalize_atomic(handle, tmp, output)
    except Exception:
        handle.close()
        if tmp.exists():
            tmp.unlink()
        raise
    return written


def export_metadata(
    output: Path,
    headers: Sequence[str],
    selected_rows: Sequence[dict],
    species_samples: set[str],
    ko_samples: set[str],
) -> None:
    extra = [
        "Analysis Group",
        "Cohort",
        "Included in species matrix",
        "Included in KO matrix",
    ]
    handle, writer, tmp = atomic_csv_writer(output)
    try:
        writer.writerow([*headers, *extra])
        for row in selected_rows:
            run_id = row["Run ID"]
            writer.writerow(
                [row.get(header, "") for header in headers]
                + [
                    row["Analysis Group"],
                    row["Cohort"],
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


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    validation_rows = []

    for cohort, relative_dir in COHORTS.items():
        cohort_dir = INTEGRATED / relative_dir
        output_dir = OUTPUT / cohort
        output_dir.mkdir(parents=True, exist_ok=True)

        (
            metadata_source,
            metadata_headers,
            phenotype_column,
            cancer_label,
            selected_rows,
            selected_runs,
        ) = load_selected_metadata(cohort, cohort_dir)
        allowed_runs = set(selected_runs)

        bacteria_files = sorted(cohort_dir.glob("Bacteria_*.txt"))
        if len(bacteria_files) != 1:
            raise ValueError(
                f"Expected one Bacteria_*.txt in {cohort_dir}; found {len(bacteria_files)}"
            )
        bacteria_source = bacteria_files[0]
        species_filter = lambda row: "|s__" in row["Taxa"] and "|t__" not in row["Taxa"]
        (
            species_sample_order,
            species_feature_order,
            species_retained_rows,
            species_duplicates,
        ) = scan_grouped_table(
            bacteria_source,
            "\t",
            allowed_runs,
            "Run ID",
            "Taxa",
            species_filter,
        )
        if species_duplicates:
            raise ValueError(f"{cohort}: {species_duplicates} duplicate species cells")

        ko_source = cohort_dir / "KO_relative_abundance.tsv"
        (
            ko_sample_order,
            ko_feature_order,
            ko_retained_rows,
            ko_duplicates,
        ) = scan_grouped_json(ko_source, allowed_runs)
        if ko_duplicates:
            raise ValueError(f"{cohort}: {ko_duplicates} duplicate KO cells")

        metadata_output = output_dir / f"{cohort}_patient_metadata_WGS_Healthy_Cancer.csv"
        species_output = output_dir / f"{cohort}_species_relative_abundance_matrix.csv"
        ko_output = output_dir / f"{cohort}_KEGG_KO_relative_abundance_matrix.csv"

        export_metadata(
            metadata_output,
            metadata_headers,
            selected_rows,
            set(species_sample_order),
            set(ko_sample_order),
        )
        export_grouped_wide_csv(
            bacteria_source,
            "\t",
            species_output,
            allowed_runs,
            "Run ID",
            "Taxa",
            "Abundance",
            species_filter,
            species_sample_order,
            species_feature_order,
        )
        export_grouped_json_wide_csv(
            ko_source,
            ko_output,
            allowed_runs,
            ko_sample_order,
            ko_feature_order,
        )

        group_counts = Counter(row["Analysis Group"] for row in selected_rows)
        missing_species = [run for run in selected_runs if run not in set(species_sample_order)]
        missing_ko = [run for run in selected_runs if run not in set(ko_sample_order)]
        validation_rows.append(
            {
                "Cohort": cohort,
                "Healthy_metadata_n": group_counts["Healthy"],
                "Cancer_metadata_n": group_counts["Cancer"],
                "Metadata_n": len(selected_rows),
                "Species_matrix_samples_n": len(species_sample_order),
                "Species_features_n": len(species_feature_order),
                "Species_nonzero_rows_n": species_retained_rows,
                "KO_matrix_samples_n": len(ko_sample_order),
                "KO_features_n": len(ko_feature_order),
                "KO_nonzero_rows_n": ko_retained_rows,
                "Missing_from_species_matrix": ";".join(missing_species),
                "Missing_from_KO_matrix": ";".join(missing_ko),
                "Cancer_phenotype_label": cancer_label,
            }
        )

        outputs = [
            ("patient_metadata", metadata_output, metadata_source, len(selected_rows), len(metadata_headers) + 4),
            ("species_matrix", species_output, bacteria_source, len(species_sample_order), len(species_feature_order) + 1),
            ("KEGG_KO_matrix", ko_output, ko_source, len(ko_sample_order), len(ko_feature_order) + 1),
        ]
        for data_type, output_file, source_file, rows, columns in outputs:
            manifest_rows.append(
                {
                    "Cohort": cohort,
                    "Data_type": data_type,
                    "Output_file": str(output_file.relative_to(OUTPUT)),
                    "Data_rows_excluding_header": rows,
                    "Columns": columns,
                    "Source_file": str(source_file),
                    "Source_SHA256": sha256_file(source_file),
                    "Output_SHA256": sha256_file(output_file),
                }
            )
        print(
            f"{cohort}: metadata={len(selected_rows)}, "
            f"species={len(species_sample_order)}x{len(species_feature_order)}, "
            f"KO={len(ko_sample_order)}x{len(ko_feature_order)}",
            flush=True,
        )

    validation_path = OUTPUT / "VALIDATION_SUMMARY.csv"
    with validation_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(validation_rows[0]))
        writer.writeheader()
        writer.writerows(validation_rows)

    manifest_path = OUTPUT / "SOURCE_MANIFEST.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
