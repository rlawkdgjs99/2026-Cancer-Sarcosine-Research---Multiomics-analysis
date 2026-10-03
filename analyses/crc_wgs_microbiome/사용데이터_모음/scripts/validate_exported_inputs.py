#!/usr/bin/env python3
"""Independent cell-level validation for the exported CRC WGS CSV inputs."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve().parent
EXPORTER_PATH = HERE / "export_analysis_used_inputs.py"
SPEC = importlib.util.spec_from_file_location("crc_exporter", EXPORTER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot import {EXPORTER_PATH}")
exporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exporter)


def grouped_table_rows(source, allowed_runs, include_row):
    current_run = None
    values = {}
    with source.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run_id = row["Run ID"]
            if run_id not in allowed_runs or not include_row(row):
                continue
            if run_id != current_run:
                if current_run is not None:
                    yield current_run, values
                current_run = run_id
                values = {}
            feature = row["Taxa"]
            if feature in values:
                raise ValueError(f"Duplicate species cell: {run_id}, {feature}")
            values[feature] = row["Abundance"]
    if current_run is not None:
        yield current_run, values


def grouped_json_rows(source, allowed_runs):
    current_run = None
    values = {}
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
            raise ValueError(f"Duplicate KO cell: {run_id}, {feature}")
        values[feature] = row["abundance"]
    if current_run is not None:
        yield current_run, values


def validate_matrix(output_path, source_groups):
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if not header or header[0] != "Run_ID":
            raise ValueError(f"Invalid matrix header in {output_path}")
        features = header[1:]
        if len(features) != len(set(features)):
            raise ValueError(f"Duplicate feature columns in {output_path}")
        rows_checked = 0
        cells_checked = 0
        for output_row, (source_run, source_values) in zip(reader, source_groups):
            if len(output_row) != len(header):
                raise ValueError(
                    f"Column-count mismatch for {output_row[0]} in {output_path}"
                )
            if output_row[0] != source_run:
                raise ValueError(
                    f"Run-order mismatch in {output_path}: "
                    f"{output_row[0]} != {source_run}"
                )
            for feature, observed in zip(features, output_row[1:]):
                expected = source_values.get(feature, "0")
                if observed != expected:
                    raise ValueError(
                        f"Value mismatch in {output_path}: "
                        f"{source_run}, {feature}: {observed} != {expected}"
                    )
                cells_checked += 1
            rows_checked += 1
        try:
            extra_output = next(reader)
        except StopIteration:
            extra_output = None
        if extra_output is not None:
            raise ValueError(f"Extra output row in {output_path}: {extra_output[0]}")
        try:
            extra_source = next(source_groups)
        except StopIteration:
            extra_source = None
        if extra_source is not None:
            raise ValueError(f"Missing output row in {output_path}: {extra_source[0]}")
    return rows_checked, len(features), cells_checked


def validate_metadata(output_path, original_headers, selected_rows):
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if len(rows) != len(selected_rows):
        raise ValueError(f"Metadata row mismatch in {output_path}")
    if len({row["Run ID"] for row in rows}) != len(rows):
        raise ValueError(f"Duplicate metadata Run IDs in {output_path}")
    for observed, expected in zip(rows, selected_rows):
        for header in original_headers:
            if observed.get(header, "") != expected.get(header, ""):
                raise ValueError(
                    f"Metadata mismatch in {output_path}: {expected['Run ID']}, {header}"
                )
        if observed["Analysis Group"] != expected["Analysis Group"]:
            raise ValueError(f"Analysis Group mismatch in {output_path}")
        if observed["Cohort"] != expected["Cohort"]:
            raise ValueError(f"Cohort mismatch in {output_path}")
    return len(rows), len(reader.fieldnames or [])


def main():
    report_rows = []
    for cohort, relative_dir in exporter.COHORTS.items():
        cohort_dir = exporter.INTEGRATED / relative_dir
        output_dir = exporter.OUTPUT / cohort
        (
            _,
            metadata_headers,
            _,
            _,
            selected_rows,
            selected_runs,
        ) = exporter.load_selected_metadata(cohort, cohort_dir)
        allowed = set(selected_runs)

        metadata_path = output_dir / f"{cohort}_patient_metadata_WGS_Healthy_Cancer.csv"
        metadata_rows, metadata_columns = validate_metadata(
            metadata_path, metadata_headers, selected_rows
        )

        bacteria_source = next(cohort_dir.glob("Bacteria_*.txt"))
        species_filter = lambda row: "|s__" in row["Taxa"] and "|t__" not in row["Taxa"]
        species_path = output_dir / f"{cohort}_species_relative_abundance_matrix.csv"
        species_rows, species_features, species_cells = validate_matrix(
            species_path,
            grouped_table_rows(bacteria_source, allowed, species_filter),
        )

        ko_source = cohort_dir / "KO_relative_abundance.tsv"
        ko_path = output_dir / f"{cohort}_KEGG_KO_relative_abundance_matrix.csv"
        ko_rows, ko_features, ko_cells = validate_matrix(
            ko_path,
            grouped_json_rows(ko_source, allowed),
        )

        report_rows.append(
            {
                "Cohort": cohort,
                "Metadata_rows_checked": metadata_rows,
                "Metadata_columns_checked": metadata_columns,
                "Species_rows_checked": species_rows,
                "Species_feature_columns_checked": species_features,
                "Species_matrix_cells_checked": species_cells,
                "KO_rows_checked": ko_rows,
                "KO_feature_columns_checked": ko_features,
                "KO_matrix_cells_checked": ko_cells,
                "Status": "PASS",
            }
        )
        print(
            f"{cohort}: PASS; metadata={metadata_rows}, "
            f"species cells={species_cells}, KO cells={ko_cells}",
            flush=True,
        )

    report_path = exporter.OUTPUT / "POST_EXPORT_CELL_LEVEL_VALIDATION.csv"
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(report_rows[0]))
        writer.writeheader()
        writer.writerows(report_rows)


if __name__ == "__main__":
    main()
