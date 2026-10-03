#!/usr/bin/env python3
"""Export the frozen CRC-WGS sarcosine functional scores as one tidy CSV.

The calculation is copied from the frozen integrated R analysis:
  Production = K00315 + K08688 (K00552 is globally absent)
  Degradation = K00301 + K00302 + K00303 + K00305 + K00306
                (K00304 is globally absent)
  Ratio = log2((Production + 1e-8) / (Degradation + 1e-8))

KO columns absent from an individual cohort matrix are treated as zero, matching the
union-matrix zero filling in integrated_analysis_pooled.R.
"""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
CRC_ROOT = PACKAGE_DIR.parent
WORKSPACE_ROOT = PACKAGE_DIR.parents[2]

OUTPUT_DIR = PACKAGE_DIR / "Derived_Sarcosine_Functional_Scores"
OUTPUT_CSV = OUTPUT_DIR / "CRC_WGS_sarcosine_functional_scores_all_cohorts.csv"
VALIDATION_CSV = OUTPUT_DIR / "CRC_WGS_sarcosine_functional_scores_validation.csv"

COHORTS = (
    ("PRJEB6070", "SupFig5a.csv"),
    ("PRJEB10878", "SupFig5b.csv"),
    ("PRJEB27928", "SupFig5c.csv"),
    ("PRJNA429097", "SupFig5d.csv"),
)

DEGRADATION_KOS = ("K00301", "K00302", "K00303", "K00305", "K00306")
PRODUCTION_KOS = ("K00315", "K08688")
PSEUDOCOUNT = 1e-8

ARCHIVE_DIR = (
    WORKSPACE_ROOT
    / "각 FIgure 패널데이터_Archive"
    / "26.07.23"
    / "Figure_Panel_Source_Data_26.07.23 2"
)
STATS_CSV = (
    CRC_ROOT
    / "Healthy_vs_Cancer_4_CRC_cohorts_integrated"
    / "results_integrated"
    / "cross_cohort"
    / "Fig1i_pathway_bubble_matrix_provenance_26.08.19"
    / "Fig1i_cross_cohort_pathway_bubble_matrix_statistics.csv"
)


def parse_bool(value: str) -> bool:
    return value.strip().upper() == "TRUE"


def load_group_map(metadata_csv: Path) -> dict[str, str]:
    group_map: dict[str, str] = {}
    with metadata_csv.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if not parse_bool(row["Included in KO matrix"]):
                continue
            run_id = row["Run ID"]
            group = row["Analysis Group"]
            if run_id in group_map:
                raise ValueError(f"Duplicate KO-included Run ID in metadata: {run_id}")
            if group not in {"Healthy", "Cancer"}:
                raise ValueError(f"Unexpected Analysis Group for {run_id}: {group}")
            group_map[run_id] = group
    return group_map


def calculate_cohort(cohort: str) -> tuple[list[dict[str, object]], set[str]]:
    cohort_dir = PACKAGE_DIR / cohort
    metadata_csv = cohort_dir / f"{cohort}_patient_metadata_WGS_Healthy_Cancer.csv"
    ko_csv = cohort_dir / f"{cohort}_KEGG_KO_relative_abundance_matrix.csv"
    group_map = load_group_map(metadata_csv)

    records: list[dict[str, object]] = []
    seen_runs: set[str] = set()
    with ko_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        column_index = {name: idx for idx, name in enumerate(header)}
        if "Run_ID" not in column_index:
            raise ValueError(f"Run_ID column missing from {ko_csv}")

        present_kos = set(column_index).intersection(DEGRADATION_KOS + PRODUCTION_KOS)
        degradation_indices = [column_index[ko] for ko in DEGRADATION_KOS if ko in column_index]
        production_indices = [column_index[ko] for ko in PRODUCTION_KOS if ko in column_index]

        for row in reader:
            run_id = row[column_index["Run_ID"]]
            if run_id in seen_runs:
                raise ValueError(f"Duplicate Run ID in KO matrix: {cohort}/{run_id}")
            seen_runs.add(run_id)
            if run_id not in group_map:
                raise ValueError(f"KO Run ID absent from included metadata: {cohort}/{run_id}")

            degradation = sum(float(row[idx]) for idx in degradation_indices)
            production = sum(float(row[idx]) for idx in production_indices)
            ratio = math.log2(
                (production + PSEUDOCOUNT) / (degradation + PSEUDOCOUNT)
            )
            records.append(
                {
                    "Run_ID": run_id,
                    "Cohort": cohort,
                    "Analysis_Group": group_map[run_id],
                    "Production_sum": production,
                    "Degradation_sum": degradation,
                    "Log2_Production_Degradation_Ratio": ratio,
                }
            )

    if seen_runs != set(group_map):
        missing = sorted(set(group_map) - seen_runs)
        raise ValueError(
            f"Included metadata Run IDs absent from KO matrix for {cohort}: {missing[:10]}"
        )
    return records, present_kos


def load_archived_values(archive_csv: Path) -> dict[tuple[str, str], tuple[str, float]]:
    values: dict[tuple[str, str], tuple[str, float]] = {}
    with archive_csv.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            key = (row["Run.ID"], row["Metric"])
            if key in values:
                raise ValueError(f"Duplicate archived metric row: {archive_csv.name}/{key}")
            values[key] = (row["Group"], float(row["Value"]))
    return values


def nearly_equal(left: float, right: float, tolerance: float = 2e-14) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def median_by_group(records: list[dict[str, object]], metric: str, group: str) -> float:
    values = [float(row[metric]) for row in records if row["Analysis_Group"] == group]
    return statistics.median(values)


def zeros_by_group(records: list[dict[str, object]], metric: str, group: str) -> int:
    return sum(
        float(row[metric]) == 0.0
        for row in records
        if row["Analysis_Group"] == group
    )


def validate_against_stats(
    all_records: list[dict[str, object]],
) -> dict[tuple[str, str], str]:
    records_by_cohort: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in all_records:
        records_by_cohort[str(row["Cohort"])].append(row)

    metric_map = {
        "Production_sum": "Production_sum",
        "Degradation_sum": "Degradation_sum",
        "Log2_Prod_Deg_Ratio": "Log2_Production_Degradation_Ratio",
    }
    result: dict[tuple[str, str], str] = {}
    with STATS_CSV.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            cohort = row["Cohort"]
            frozen_metric = row["Metric"]
            output_metric = metric_map[frozen_metric]
            cohort_records = records_by_cohort[cohort]

            checks = [
                sum(r["Analysis_Group"] == "Healthy" for r in cohort_records)
                == int(row["n_Healthy"]),
                sum(r["Analysis_Group"] == "Cancer" for r in cohort_records)
                == int(row["n_Cancer"]),
                zeros_by_group(cohort_records, output_metric, "Healthy")
                == int(row["zero_Healthy"]),
                zeros_by_group(cohort_records, output_metric, "Cancer")
                == int(row["zero_Cancer"]),
                nearly_equal(
                    median_by_group(cohort_records, output_metric, "Healthy"),
                    float(row["median_Healthy"]),
                ),
                nearly_equal(
                    median_by_group(cohort_records, output_metric, "Cancer"),
                    float(row["median_Cancer"]),
                ),
            ]
            result[(cohort, frozen_metric)] = "PASS" if all(checks) else "FAIL"
    return result


def main() -> None:
    all_records: list[dict[str, object]] = []
    cohort_present_kos: dict[str, set[str]] = {}
    for cohort, _ in COHORTS:
        records, present_kos = calculate_cohort(cohort)
        all_records.extend(records)
        cohort_present_kos[cohort] = present_kos

    if len(all_records) != 1647:
        raise ValueError(f"Expected 1,647 KO-profile samples; found {len(all_records):,}")
    run_ids = [str(row["Run_ID"]) for row in all_records]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("Run IDs are not unique across the combined output")

    stats_status = validate_against_stats(all_records)
    if not stats_status or any(status != "PASS" for status in stats_status.values()):
        raise ValueError(f"Frozen-statistics validation failed: {stats_status}")

    output_metric_map = {
        "Production_sum": "Production_sum",
        "Degradation_sum": "Degradation_sum",
        "Log2_Production_Degradation_Ratio": "Log2_Prod_Deg_Ratio",
    }
    validation_rows: list[dict[str, object]] = []

    for cohort, archive_name in COHORTS:
        cohort_records = [row for row in all_records if row["Cohort"] == cohort]
        archived = load_archived_values(ARCHIVE_DIR / archive_name)
        mismatch_count = 0
        max_abs_difference = 0.0

        for record in cohort_records:
            for output_metric, archived_metric in output_metric_map.items():
                key = (str(record["Run_ID"]), archived_metric)
                if key not in archived:
                    raise ValueError(f"Archived value missing: {archive_name}/{key}")
                archived_group, archived_value = archived[key]
                if archived_group != record["Analysis_Group"]:
                    mismatch_count += 1
                calculated = float(record[output_metric])
                difference = abs(calculated - archived_value)
                max_abs_difference = max(max_abs_difference, difference)
                if not nearly_equal(calculated, archived_value):
                    mismatch_count += 1

        expected_archive_rows = len(cohort_records) * 3
        if len(archived) != expected_archive_rows:
            raise ValueError(
                f"Unexpected archived row count for {archive_name}: "
                f"{len(archived)} != {expected_archive_rows}"
            )
        if mismatch_count:
            raise ValueError(
                f"Archived per-sample validation failed for {cohort}: "
                f"{mismatch_count} mismatches"
            )

        group_counts = defaultdict(int)
        for row in cohort_records:
            group_counts[str(row["Analysis_Group"])] += 1
        present = cohort_present_kos[cohort]
        validation_rows.append(
            {
                "Cohort": cohort,
                "n_Healthy": group_counts["Healthy"],
                "n_Cancer": group_counts["Cancer"],
                "n_Total": len(cohort_records),
                "Production_KOs_present": ";".join(
                    ko for ko in PRODUCTION_KOS if ko in present
                ),
                "Degradation_KOs_present": ";".join(
                    ko for ko in DEGRADATION_KOS if ko in present
                ),
                "Archived_source": archive_name,
                "Archived_value_rows_checked": expected_archive_rows,
                "Max_absolute_difference": max_abs_difference,
                "Per_sample_validation": "PASS",
                "Frozen_summary_validation": "PASS",
            }
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_fields = (
        "Run_ID",
        "Cohort",
        "Analysis_Group",
        "Production_sum",
        "Degradation_sum",
        "Log2_Production_Degradation_Ratio",
    )
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_records)

    validation_fields = tuple(validation_rows[0])
    with VALIDATION_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=validation_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(validation_rows)

    print(f"Wrote {OUTPUT_CSV} ({len(all_records):,} rows)")
    print(f"Wrote {VALIDATION_CSV} ({len(validation_rows)} cohort checks)")


if __name__ == "__main__":
    main()
