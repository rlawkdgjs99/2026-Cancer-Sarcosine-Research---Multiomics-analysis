#!/usr/bin/env python3
"""Independent full-cell validation of the exported NSCLC WGS CSV package."""

from __future__ import annotations

import csv
import importlib.util
import math
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Sequence, Tuple


HERE = Path(__file__).resolve().parent
EXPORTER_PATH = HERE / "export_analysis_used_inputs.py"
SPEC = importlib.util.spec_from_file_location("nsclc_exporter", EXPORTER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot import {EXPORTER_PATH}")
exporter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = exporter
SPEC.loader.exec_module(exporter)


def write_dict_rows_atomic(path: Path, rows: Sequence[dict]) -> None:
    exporter.write_dict_rows_atomic(path, rows)


def grouped_species_rows(
    source: Path, allowed_runs: set[str]
) -> Iterator[Tuple[str, Dict[str, str]]]:
    current_run = None
    values: Dict[str, str] = {}
    with source.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run_id = row["Run ID"]
            if (
                run_id not in allowed_runs
                or "|s__" not in row["Taxa"]
                or "|t__" in row["Taxa"]
            ):
                continue
            if run_id != current_run:
                if current_run is not None:
                    yield current_run, values
                current_run = run_id
                values = {}
            feature = row["Taxa"]
            if feature in values:
                raise ValueError(f"Duplicate source species cell: {run_id}, {feature}")
            values[feature] = row["Abundance"]
    if current_run is not None:
        yield current_run, values


def grouped_ko_rows(
    source: Path, allowed_runs: set[str]
) -> Iterator[Tuple[str, Dict[str, str]]]:
    current_run = None
    values: Dict[str, str] = {}
    for row in exporter.iter_json_array(source):
        run_id = row["run_id"]
        if run_id not in allowed_runs:
            continue
        if run_id != current_run:
            if current_run is not None:
                yield current_run, values
            current_run = run_id
            values = {}
        feature = row["ko"]
        if feature in values:
            raise ValueError(f"Duplicate source KO cell: {run_id}, {feature}")
        values[feature] = row["abundance"]
    if current_run is not None:
        yield current_run, values


def validate_matrix(
    output_path: Path,
    source_groups: Iterator[Tuple[str, Dict[str, str]]],
) -> Tuple[int, int, int, float, float]:
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if not header or header[0] != "Run_ID":
            raise ValueError(f"Invalid matrix header: {output_path}")
        features = header[1:]
        if len(features) != len(set(features)):
            raise ValueError(f"Duplicate feature columns: {output_path}")

        rows_checked = 0
        cells_checked = 0
        min_row_sum = math.inf
        max_row_sum = -math.inf
        for output_row, (source_run, source_values) in zip(reader, source_groups):
            if len(output_row) != len(header):
                raise ValueError(
                    f"Column-count mismatch for {output_row[0]} in {output_path}"
                )
            if output_row[0] != source_run:
                raise ValueError(
                    f"Run-order mismatch in {output_path}: {output_row[0]} != {source_run}"
                )
            row_sum = 0.0
            for feature, observed in zip(features, output_row[1:]):
                expected = source_values.get(feature, "0")
                if observed != expected:
                    raise ValueError(
                        f"Value mismatch in {output_path}: {source_run}, {feature}: "
                        f"{observed} != {expected}"
                    )
                value = float(observed)
                if not math.isfinite(value) or value < 0:
                    raise ValueError(
                        f"Invalid abundance in {output_path}: {source_run}, {feature}, {observed}"
                    )
                row_sum += value
                cells_checked += 1
            min_row_sum = min(min_row_sum, row_sum)
            max_row_sum = max(max_row_sum, row_sum)
            rows_checked += 1

        try:
            extra_output = next(reader)
        except StopIteration:
            extra_output = None
        if extra_output is not None:
            raise ValueError(f"Extra output row: {extra_output[0]} in {output_path}")
        try:
            extra_source = next(source_groups)
        except StopIteration:
            extra_source = None
        if extra_source is not None:
            raise ValueError(f"Missing output row: {extra_source[0]} in {output_path}")
    return rows_checked, len(features), cells_checked, min_row_sum, max_row_sum


def validate_metadata(
    output_path: Path,
    original_headers: Sequence[str],
    selected_rows: Sequence[dict],
    species_runs: set[str],
    ko_runs: set[str],
) -> Tuple[int, int]:
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        output_rows = list(reader)
    expected_headers = list(original_headers) + [
        "Analysis Group", "Cohort", "Cohort Role",
        "Included in species matrix", "Included in KO matrix",
    ]
    if reader.fieldnames != expected_headers:
        raise ValueError(f"Metadata schema mismatch in {output_path}")
    if len(output_rows) != len(selected_rows):
        raise ValueError(f"Metadata row-count mismatch in {output_path}")
    if len({row["Run ID"] for row in output_rows}) != len(output_rows):
        raise ValueError(f"Duplicate metadata Run IDs in {output_path}")

    for observed, expected in zip(output_rows, selected_rows):
        run_id = expected["Run ID"]
        for header in original_headers:
            if observed[header] != expected[header]:
                raise ValueError(f"Metadata mismatch: {run_id}, {header}")
        for header in ("Analysis Group", "Cohort", "Cohort Role"):
            if observed[header] != expected[header]:
                raise ValueError(f"Derived metadata mismatch: {run_id}, {header}")
        if observed["Included in species matrix"] != (
            "TRUE" if run_id in species_runs else "FALSE"
        ):
            raise ValueError(f"Species inclusion flag mismatch: {run_id}")
        if observed["Included in KO matrix"] != (
            "TRUE" if run_id in ko_runs else "FALSE"
        ):
            raise ValueError(f"KO inclusion flag mismatch: {run_id}")
    return len(output_rows), len(expected_headers)


def close_enough(observed: float, expected: float) -> bool:
    return math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-15)


def parse_optional_float(value: str):
    return None if value == "" else float(value)


def validate_scores(cohort_state: Dict[str, dict]) -> List[dict]:
    combined_path = (
        exporter.OUTPUT / "Derived_Sarcosine_Functional_Scores"
        / "NSCLC_WGS_sarcosine_functional_scores_all_cohorts.csv"
    )
    with combined_path.open("r", encoding="utf-8", newline="") as handle:
        combined = list(csv.DictReader(handle))
    combined_by_key = {(row["Cohort"], row["Run_ID"]): row for row in combined}
    if len(combined_by_key) != len(combined):
        raise ValueError("Duplicate Cohort/Run_ID keys in combined score table")

    report_rows: List[dict] = []
    expected_total = 0
    for cohort, spec in exporter.COHORTS.items():
        state = cohort_state[cohort]
        frozen_path = (
            exporter.ROOT / spec.folder / "results" / f"sarcosine_scores_{cohort}.csv"
        )
        with frozen_path.open("r", encoding="utf-8", newline="") as handle:
            frozen_rows = list(csv.DictReader(handle))
        frozen_by_run = {row["RunID"]: row for row in frozen_rows}
        if len(frozen_by_run) != len(frozen_rows):
            raise ValueError(f"Duplicate RunID in frozen score table: {cohort}")

        max_deg_diff = 0.0
        max_prod_diff = 0.0
        max_ratio_diff = 0.0
        ratio_na_n = 0
        ko_path = state["ko_output"]
        with ko_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            index = {feature: position for position, feature in enumerate(header)}
            for row in reader:
                run_id = row[0]
                frozen = frozen_by_run.get(run_id)
                combined_row = combined_by_key.get((cohort, run_id))
                if frozen is None or combined_row is None:
                    raise ValueError(f"Missing score record: {cohort}, {run_id}")
                if combined_row["Analysis_Group"] != frozen["group"]:
                    raise ValueError(f"Score group mismatch: {cohort}, {run_id}")
                if combined_row["Cohort_Role"] != spec.role:
                    raise ValueError(f"Score cohort-role mismatch: {cohort}, {run_id}")
                if (
                    combined_row["Degradation_sum"] != frozen["degradation"]
                    or combined_row["Production_sum"] != frozen["production"]
                    or combined_row["Log2_Production_Degradation_Ratio"]
                    != frozen["prod_deg_log2ratio"]
                ):
                    raise ValueError(f"Combined score is not an exact frozen copy: {cohort}, {run_id}")

                def ko_value(ko: str) -> float:
                    return float(row[index[ko]]) if ko in index else 0.0

                degradation = sum(ko_value(ko) for ko in exporter.DEGRADATION_KOS)
                production = sum(ko_value(ko) for ko in exporter.PRODUCTION_KOS)
                ratio = (
                    math.log2(production / degradation)
                    if production > 0 and degradation > 0
                    else None
                )
                frozen_deg = float(frozen["degradation"])
                frozen_prod = float(frozen["production"])
                frozen_ratio = parse_optional_float(frozen["prod_deg_log2ratio"])
                if not close_enough(degradation, frozen_deg):
                    raise ValueError(f"Degradation mismatch: {cohort}, {run_id}")
                if not close_enough(production, frozen_prod):
                    raise ValueError(f"Production mismatch: {cohort}, {run_id}")
                if ratio is None:
                    ratio_na_n += 1
                    if frozen_ratio is not None:
                        raise ValueError(f"Ratio missingness mismatch: {cohort}, {run_id}")
                else:
                    if frozen_ratio is None or not close_enough(ratio, frozen_ratio):
                        raise ValueError(f"Ratio mismatch: {cohort}, {run_id}")
                    max_ratio_diff = max(max_ratio_diff, abs(ratio - frozen_ratio))
                max_deg_diff = max(max_deg_diff, abs(degradation - frozen_deg))
                max_prod_diff = max(max_prod_diff, abs(production - frozen_prod))

        metadata_groups = {
            row["Run ID"]: row["Analysis Group"] for row in state["selected_rows"]
        }
        if set(frozen_by_run) != set(metadata_groups):
            raise ValueError(f"Frozen score/metadata Run-ID set mismatch: {cohort}")
        for run_id, row in frozen_by_run.items():
            if row["group"] != metadata_groups[run_id]:
                raise ValueError(f"Frozen score/metadata group mismatch: {cohort}, {run_id}")

        expected_total += len(frozen_rows)
        report_rows.append(
            {
                "Cohort": cohort,
                "Cohort_role": spec.role,
                "Samples_checked": len(frozen_rows),
                "R_n": sum(row["group"] == "R" for row in frozen_rows),
                "NR_n": sum(row["group"] == "NR" for row in frozen_rows),
                "Ratio_NA_n": ratio_na_n,
                "Max_abs_difference_degradation": format(max_deg_diff, ".17g"),
                "Max_abs_difference_production": format(max_prod_diff, ".17g"),
                "Max_abs_difference_log2_ratio": format(max_ratio_diff, ".17g"),
                "Frozen_score_file_SHA256": exporter.sha256_file(frozen_path),
                "Status": "PASS",
            }
        )
    if len(combined) != expected_total:
        raise ValueError(
            f"Combined score row count {len(combined)} != expected {expected_total}"
        )
    return report_rows


def validate_source_immutability() -> List[dict]:
    manifest_path = exporter.OUTPUT / "SOURCE_MANIFEST.csv"
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        manifest = list(csv.DictReader(handle))
    expected_hashes: Dict[str, str] = {}
    for row in manifest:
        source_paths = row["Source_file"].split(";")
        source_hashes = row["Source_SHA256"].split(";")
        if len(source_paths) != len(source_hashes):
            raise ValueError("Source path/hash cardinality mismatch in manifest")
        for source_path, source_hash in zip(source_paths, source_hashes):
            if source_path in expected_hashes and expected_hashes[source_path] != source_hash:
                raise ValueError(f"Conflicting manifest hashes for {source_path}")
            expected_hashes[source_path] = source_hash
        for path_field, hash_field in (
            ("Selection_script", "Selection_script_SHA256"),
            ("Score_script", "Score_script_SHA256"),
        ):
            path_text = row[path_field]
            hash_text = row[hash_field]
            candidate = exporter.ROOT / path_text
            if candidate.is_file() and len(hash_text) == 64:
                expected_hashes[path_text] = hash_text

    checks = []
    for relative_path in sorted(expected_hashes):
        path = exporter.ROOT / relative_path
        current = exporter.sha256_file(path)
        expected = expected_hashes[relative_path]
        status = "PASS" if current == expected else "FAIL"
        checks.append(
            {
                "Source_file": relative_path,
                "Manifest_SHA256": expected,
                "Current_SHA256": current,
                "Status": status,
            }
        )
        if status != "PASS":
            raise ValueError(f"Protected source hash changed: {relative_path}")
    return checks


def validate_cross_cohort_sample_names(cohort_state: Dict[str, dict]) -> None:
    sample_name_cohorts: Dict[str, set[str]] = {}
    for cohort, state in cohort_state.items():
        names = [row["Sample name"] for row in state["selected_rows"]]
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate Sample name within cohort: {cohort}")
        for name in names:
            sample_name_cohorts.setdefault(name, set()).add(cohort)

    summary_path = exporter.OUTPUT / "VALIDATION_SUMMARY.csv"
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        summary = {row["Cohort"]: row for row in csv.DictReader(handle)}
    global_unique = len(sample_name_cohorts)
    for cohort, state in cohort_state.items():
        names = {row["Sample name"] for row in state["selected_rows"]}
        shared = sum(len(sample_name_cohorts[name]) > 1 for name in names)
        row = summary[cohort]
        if int(row["Unique_sample_names_within_cohort_n"]) != len(names):
            raise ValueError(f"Within-cohort Sample-name count mismatch: {cohort}")
        if int(row["Sample_names_also_seen_in_other_cohorts_n"]) != shared:
            raise ValueError(f"Cross-cohort Sample-name reuse mismatch: {cohort}")
        if int(row["Global_unique_sample_names_across_all_cohorts_n"]) != global_unique:
            raise ValueError(f"Global Sample-name count mismatch: {cohort}")


def main() -> None:
    cell_reports: List[dict] = []
    cohort_state: Dict[str, dict] = {}
    for cohort, spec in exporter.COHORTS.items():
        cohort_dir = exporter.ROOT / spec.folder
        metadata_source, headers, selected_rows, selected_runs = exporter.load_selected_metadata(
            cohort, spec
        )
        allowed = set(selected_runs)
        bacteria_source = exporter.unique_source(cohort_dir, "Bacteria_*.txt")
        ko_source = cohort_dir / "KO_relative_abundance.tsv"
        species_runs, _, _ = exporter.scan_grouped_table(
            bacteria_source,
            allowed,
            lambda row: "|s__" in row["Taxa"] and "|t__" not in row["Taxa"],
        )
        ko_runs, _, _ = exporter.scan_grouped_json(ko_source, allowed)

        output_dir = exporter.OUTPUT / cohort
        metadata_path = output_dir / f"{cohort}_patient_metadata_WGS_NSCLC_R_NR.csv"
        species_path = output_dir / f"{cohort}_species_relative_abundance_matrix.csv"
        ko_path = output_dir / f"{cohort}_KEGG_KO_relative_abundance_matrix.csv"
        metadata_rows, metadata_columns = validate_metadata(
            metadata_path, headers, selected_rows, set(species_runs), set(ko_runs)
        )
        species_stats = validate_matrix(
            species_path, grouped_species_rows(bacteria_source, allowed)
        )
        ko_stats = validate_matrix(ko_path, grouped_ko_rows(ko_source, allowed))

        if not (99.0 <= species_stats[3] <= 101.0 and 99.0 <= species_stats[4] <= 101.0):
            raise ValueError(f"Species TSS row sums outside expected percentage scale: {cohort}")
        if not (99.0 <= ko_stats[3] <= 101.0 and 99.0 <= ko_stats[4] <= 101.0):
            raise ValueError(f"KO TSS row sums outside expected percentage scale: {cohort}")

        cell_reports.append(
            {
                "Cohort": cohort,
                "Metadata_rows_checked": metadata_rows,
                "Metadata_columns_checked": metadata_columns,
                "Species_rows_checked": species_stats[0],
                "Species_feature_columns_checked": species_stats[1],
                "Species_matrix_cells_checked": species_stats[2],
                "Species_row_sum_min": format(species_stats[3], ".17g"),
                "Species_row_sum_max": format(species_stats[4], ".17g"),
                "KO_rows_checked": ko_stats[0],
                "KO_feature_columns_checked": ko_stats[1],
                "KO_matrix_cells_checked": ko_stats[2],
                "KO_row_sum_min": format(ko_stats[3], ".17g"),
                "KO_row_sum_max": format(ko_stats[4], ".17g"),
                "Status": "PASS",
            }
        )
        cohort_state[cohort] = {
            "selected_rows": selected_rows,
            "selected_runs": selected_runs,
            "ko_output": ko_path,
        }
        print(
            f"{cohort}: PASS; metadata={metadata_rows}, "
            f"species cells={species_stats[2]}, KO cells={ko_stats[2]}",
            flush=True,
        )

    validate_cross_cohort_sample_names(cohort_state)
    score_reports = validate_scores(cohort_state)
    source_checks = validate_source_immutability()
    write_dict_rows_atomic(
        exporter.OUTPUT / "POST_EXPORT_CELL_LEVEL_VALIDATION.csv", cell_reports
    )
    write_dict_rows_atomic(
        exporter.OUTPUT / "DERIVED_SCORE_VALIDATION.csv", score_reports
    )
    write_dict_rows_atomic(
        exporter.OUTPUT / "SOURCE_IMMUTABILITY_CHECK.csv", source_checks
    )
    print(
        f"Derived scores: PASS ({sum(row['Samples_checked'] for row in score_reports)} samples); "
        f"protected sources: PASS ({len(source_checks)} files)",
        flush=True,
    )


if __name__ == "__main__":
    main()
