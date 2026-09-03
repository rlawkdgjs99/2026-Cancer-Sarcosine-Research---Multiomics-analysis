#!/usr/bin/env python3
"""Preflight audit for matched sarcosine metabolomics and bulk RNA-seq.

Inputs are the immutable public Zenodo matrices plus Supplementary Table 1.
Outputs are descriptive/QC artifacts only; no pathway inference is performed.
"""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ModuleNotFoundError:
    matplotlib = None
    plt = None
    HAS_MATPLOTLIB = False


PROJECT = Path(__file__).resolve().parents[1]
RNA_PATH = PROJECT / "RNAseq_Data" / "bulkRNA_matrix_TPM.csv"
LC_PATH = PROJECT / "Metabolomics_Data" / "LC_matrix.csv"
GC_PATH = PROJECT / "Metabolomics_Data" / "GC_matrix.csv"
SUP_TABLES = PROJECT / "논문_SupData" / "41588_2024_1662_MOESM4_ESM.xlsx"
OUT_DIR = PROJECT / "results" / "00_input_audit_26.09.02"

EXPECTED_MD5 = {
    RNA_PATH: "90fc32db521cf1b4ed2142ca429f09b2",
    LC_PATH: "00eb5c2817d1e1fe8e7c4b305f9e6338",
    GC_PATH: "589086409d8a6375f4a06d3f0389dcb7",
}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sample_columns(frame: pd.DataFrame) -> list[str]:
    return [str(column) for column in frame.columns if str(column).endswith(("_T", "_N"))]


def summarize(values: pd.Series) -> dict:
    numeric = pd.to_numeric(values, errors="coerce")
    return {
        "n": int(numeric.notna().sum()),
        "missing": int(numeric.isna().sum()),
        "nonpositive": int((numeric <= 0).sum()),
        "min": float(numeric.min()),
        "q1": float(numeric.quantile(0.25)),
        "median": float(numeric.median()),
        "q3": float(numeric.quantile(0.75)),
        "max": float(numeric.max()),
        "mean": float(numeric.mean()),
        "sd": float(numeric.std(ddof=1)),
    }


def native_scalar(value: object) -> object:
    """Convert pandas/numpy scalars to JSON-native values."""
    if pd.isna(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    observed_md5 = {str(path.relative_to(PROJECT)): md5(path) for path in EXPECTED_MD5}
    for path, expected in EXPECTED_MD5.items():
        observed = observed_md5[str(path.relative_to(PROJECT))]
        if observed != expected:
            raise ValueError(f"MD5 mismatch for {path}: {observed} != {expected}")

    rna = pd.read_csv(RNA_PATH, index_col=0)
    # The public LC file contains a few legacy single-byte characters in
    # annotation text; latin-1 is lossless for those bytes and sample headers.
    lc = pd.read_csv(LC_PATH, encoding="latin-1")
    gc = pd.read_csv(GC_PATH, encoding="latin-1")
    clinical = pd.read_excel(SUP_TABLES, sheet_name="Table S1", header=1, index_col=0)
    clinical.index = clinical.index.astype(str)

    rna_samples = rna.columns.astype(str).tolist()
    lc_samples = sample_columns(lc)
    gc_samples = sample_columns(gc)
    clinical_samples = clinical.index.tolist()

    if set(rna_samples) != set(lc_samples) or set(rna_samples) != set(gc_samples):
        raise ValueError("RNA, LC-MS and GC-MS sample sets are not identical")
    tumor_samples = sorted(sample for sample in rna_samples if sample.endswith("_T"))
    normal_samples = sorted(sample for sample in rna_samples if sample.endswith("_N"))
    if set(tumor_samples) != set(clinical_samples):
        raise ValueError("Tumor sample IDs do not match Supplementary Table 1")

    sarcosine_rows = gc.loc[gc["Metabolite name"].astype(str).str.casefold() == "sarcosine"]
    if len(sarcosine_rows) != 1:
        raise ValueError(f"Expected one exact Sarcosine row, found {len(sarcosine_rows)}")
    sarcosine_row = sarcosine_rows.iloc[0]
    sarcosine = pd.to_numeric(sarcosine_row[rna_samples], errors="coerce")
    sarcosine.index = rna_samples
    if sarcosine.isna().any() or (sarcosine <= 0).any():
        raise ValueError("Sarcosine contains missing or non-positive values")

    tumor_sarcosine = sarcosine.loc[tumor_samples].sort_index()
    normal_sarcosine = sarcosine.loc[normal_samples].sort_index()
    tumor_median = float(tumor_sarcosine.median())
    high = tumor_sarcosine > tumor_median
    low = tumor_sarcosine <= tumor_median
    if int(high.sum()) + int(low.sum()) != len(tumor_sarcosine):
        raise AssertionError("High/Low partition does not cover all tumors")

    group_map = pd.DataFrame(
        {
            "sample_id": tumor_sarcosine.index,
            "patient_id": tumor_sarcosine.index.str.replace(r"_T$", "", regex=True),
            "sarcosine_normalized_intensity": tumor_sarcosine.values,
            "log2_sarcosine_intensity": np.log2(tumor_sarcosine.values),
            "tumor_median_cutpoint": tumor_median,
            "sarcosine_group": np.where(high.values, "High", "Low"),
        }
    ).set_index("sample_id")
    group_map = group_map.join(clinical, how="left")
    group_map.reset_index().to_csv(OUT_DIR / "tumor_sarcosine_group_map.csv", index=False)

    normal_patient_to_sample = {sample.removesuffix("_N"): sample for sample in normal_samples}
    tumor_patient_to_sample = {sample.removesuffix("_T"): sample for sample in tumor_samples}
    paired_patients = sorted(set(normal_patient_to_sample) & set(tumor_patient_to_sample))
    paired = pd.DataFrame(
        {
            "patient_id": paired_patients,
            "normal_sample": [normal_patient_to_sample[x] for x in paired_patients],
            "tumor_sample": [tumor_patient_to_sample[x] for x in paired_patients],
        }
    )
    paired["normal_sarcosine"] = [sarcosine.loc[x] for x in paired["normal_sample"]]
    paired["tumor_sarcosine"] = [sarcosine.loc[x] for x in paired["tumor_sample"]]
    paired["log2_tumor_normal_ratio"] = np.log2(
        paired["tumor_sarcosine"] / paired["normal_sarcosine"]
    )
    paired.to_csv(OUT_DIR / "paired_tumor_normal_sarcosine.csv", index=False)
    clinical_balance = []
    for variable in ["sex", "age", "pT", "pN", "pM", "stage", "grade", "immune subtype", "batch"]:
        table = pd.crosstab(group_map["sarcosine_group"], group_map[variable], dropna=False)
        for group_name, row in table.iterrows():
            for level, count in row.items():
                clinical_balance.append(
                    {
                        "variable": variable,
                        "level": str(level),
                        "sarcosine_group": group_name,
                        "n": int(count),
                    }
                )
    pd.DataFrame(clinical_balance).to_csv(OUT_DIR / "sarcosine_group_clinical_counts.csv", index=False)

    genes_of_interest = [
        "SARDH", "PIPOX", "DMGDH", "GNMT", "GLDC", "GCSH", "AMT", "CD28",
        "IFNG", "IL12RB1", "IL12RB2", "STAT4", "NFKB1", "RELA", "CD3D", "CD8A",
    ]
    gene_presence = pd.DataFrame(
        {
            "gene": genes_of_interest,
            "present_in_rna_matrix": [gene in rna.index for gene in genes_of_interest],
            "duplicate_row_count": [int((rna.index == gene).sum()) for gene in genes_of_interest],
        }
    )
    gene_presence.to_csv(OUT_DIR / "rna_gene_presence_check.csv", index=False)

    # A value range of 0--17 is incompatible with untransformed TPM for common
    # highly expressed genes.  We report it descriptively and do not assign an
    # undocumented exact transformation formula.
    rna_numeric = rna.to_numpy(dtype=float)
    rna_summary = {
        "features": int(rna.shape[0]),
        "samples": int(rna.shape[1]),
        "duplicate_gene_rows": int(rna.index.duplicated().sum()),
        "missing_values": int(np.isnan(rna_numeric).sum()),
        "min": float(np.nanmin(rna_numeric)),
        "max": float(np.nanmax(rna_numeric)),
        "zero_fraction": float(np.mean(rna_numeric == 0)),
        "observed_scale_interpretation": (
            "already transformed/log-scale-like; exact formula is not documented in the supplied matrix"
        ),
    }

    audit = {
        "input_md5_verified": observed_md5,
        "matrix_dimensions": {
            "RNA_features_by_samples": list(rna.shape),
            "LC_rows_by_total_columns": list(lc.shape),
            "GC_rows_by_total_columns": list(gc.shape),
            "clinical_tumors_by_variables": list(clinical.shape),
        },
        "sample_alignment": {
            "identical_sample_sets_RNA_LC_GC": True,
            "RNA_LC_same_order": rna_samples == lc_samples,
            "RNA_GC_same_order": rna_samples == gc_samples,
            "LC_GC_same_order": lc_samples == gc_samples,
            "total_samples": len(rna_samples),
            "tumors": len(tumor_samples),
            "NATs": len(normal_samples),
            "paired_tumor_NAT_patients": len(paired_patients),
            "tumor_IDs_match_clinical_table": True,
        },
        "sarcosine_annotation": {
            key: native_scalar(sarcosine_row[key])
            for key in [
                "Alignment ID", "Average Rt(min)", "Quant mass", "Metabolite name",
                "KEGG", "Total score", "HMDB", "Super Class", "Class", "Sub Class",
            ]
        },
        "sarcosine_values": {
            "all_samples": summarize(sarcosine),
            "tumors": summarize(tumor_sarcosine),
            "NATs": summarize(normal_sarcosine),
            "tumor_median_cutpoint": tumor_median,
            "tumor_high_n": int(high.sum()),
            "tumor_low_n": int(low.sum()),
            "values_equal_to_median": int((tumor_sarcosine == tumor_median).sum()),
        },
        "paired_tumor_NAT_descriptive": {
            "n_pairs": len(paired),
            "median_log2_tumor_normal_ratio": float(paired["log2_tumor_normal_ratio"].median()),
            "mean_log2_tumor_normal_ratio": float(paired["log2_tumor_normal_ratio"].mean()),
            "tumor_higher_pairs": int((paired["log2_tumor_normal_ratio"] > 0).sum()),
            "tumor_lower_pairs": int((paired["log2_tumor_normal_ratio"] < 0).sum()),
        },
        "rna_matrix": rna_summary,
        "clinical_variables": clinical.columns.tolist(),
        "software": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__ if HAS_MATPLOTLIB else "not installed; QC figure skipped",
        },
    }
    (OUT_DIR / "matched_multiomics_preflight.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # QC figure: descriptive display only, with no causal interpretation.
    if HAS_MATPLOTLIB:
        colors = {"NAT": "#1F9E89", "Tumor": "#C83E3E", "Low": "#1F9E89", "High": "#C83E3E"}
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)

        rng = np.random.default_rng(260902)
        for x, (label, values) in enumerate([("NAT", normal_sarcosine), ("Tumor", tumor_sarcosine)]):
            axes[0].boxplot(
                np.log2(values), positions=[x], widths=0.45, patch_artist=True,
                boxprops={"facecolor": colors[label], "alpha": 0.8},
                medianprops={"color": "black"}, whiskerprops={"color": "black"},
                capprops={"color": "black"}, flierprops={"marker": ""},
            )
            axes[0].scatter(
                x + rng.uniform(-0.12, 0.12, len(values)), np.log2(values), s=13,
                color=colors[label], alpha=0.55, edgecolor="none",
            )
        axes[0].set_xticks([0, 1], ["NAT\n(n=50)", "Tumor\n(n=100)"])
        axes[0].set_ylabel("log2 normalized GC-MS intensity")
        axes[0].set_title("Sarcosine distribution")

        for _, row in paired.iterrows():
            axes[1].plot([0, 1], np.log2([row["normal_sarcosine"], row["tumor_sarcosine"]]), color="#777777", alpha=0.28, lw=0.7)
        axes[1].scatter(np.zeros(len(paired)), np.log2(paired["normal_sarcosine"]), s=15, color=colors["NAT"], zorder=3)
        axes[1].scatter(np.ones(len(paired)), np.log2(paired["tumor_sarcosine"]), s=15, color=colors["Tumor"], zorder=3)
        axes[1].set_xticks([0, 1], ["NAT", "Matched tumor"])
        axes[1].set_ylabel("log2 normalized GC-MS intensity")
        axes[1].set_title("Paired samples (n=50)")

        for x, group_name in enumerate(["Low", "High"]):
            values = group_map.loc[group_map["sarcosine_group"] == group_name, "sarcosine_normalized_intensity"]
            axes[2].boxplot(
                np.log2(values), positions=[x], widths=0.45, patch_artist=True,
                boxprops={"facecolor": colors[group_name], "alpha": 0.8},
                medianprops={"color": "black"}, whiskerprops={"color": "black"},
                capprops={"color": "black"}, flierprops={"marker": ""},
            )
            axes[2].scatter(
                x + rng.uniform(-0.12, 0.12, len(values)), np.log2(values), s=13,
                color=colors[group_name], alpha=0.55, edgecolor="none",
            )
        axes[2].set_xticks([0, 1], ["Low\n(n=50)", "High\n(n=50)"])
        axes[2].set_ylabel("log2 normalized GC-MS intensity")
        axes[2].set_title("Tumor median split")

        for axis in axes:
            axis.spines[["top", "right"]].set_visible(False)
            axis.grid(axis="y", color="#E8E8E8", lw=0.7)
        fig.suptitle("TJ-RCC tissue sarcosine preflight", fontsize=15, fontweight="bold")
        fig.savefig(OUT_DIR / "Fig_Sarcosine_Matched_Data_Preflight.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
