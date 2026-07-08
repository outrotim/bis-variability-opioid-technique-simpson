#!/usr/bin/env python3
"""Reproduce core Gamma GLM estimates from a local analysis dataset.

No patient-level data are distributed with this repository. Users must obtain
source data under the applicable data terms and reconstruct the analysis
dataset locally before running this script.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.genmod import families
from statsmodels.genmod.generalized_linear_model import GLM


MODEL2_COVARIATES = [
    "CV_BIS",
    "age",
    "male",
    "bmi",
    "asa_num",
    "opdur_min",
    "TWA_BIS",
]


def clean_json(value: Any) -> Any:
    if isinstance(value, (np.floating, float)):
        value = float(value)
        if math.isnan(value) or math.isinf(value):
            return None
        return value
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, dict):
        return {key: clean_json(val) for key, val in value.items()}
    if isinstance(value, list):
        return [clean_json(val) for val in value]
    return value


def load_dataset(analysis_dataset: Path, mac_source: Path | None = None) -> pd.DataFrame:
    df = pd.read_csv(analysis_dataset)
    if "TWA_MAC" not in df.columns:
        if mac_source is None:
            raise ValueError("TWA_MAC is missing; provide --mac-source or include TWA_MAC in the analysis dataset.")
        mac = pd.read_csv(mac_source, usecols=["caseid", "TWA_MAC"])
        df = df.merge(mac, on="caseid", how="left")
    return df


def add_model_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["male"] = (out["sex"] == "M").astype(int)
    out["asa_num"] = out["asa"].astype(float)
    if "opdur_min" not in out.columns:
        if {"opend", "opstart"}.issubset(out.columns):
            out["opdur_min"] = (out["opend"] - out["opstart"]) / 60.0
        else:
            raise ValueError("Missing opdur_min or opstart/opend columns.")
    out["tiva"] = (out["TWA_MAC"] < 0.3).astype(float)
    return out


def fit_gamma(data: pd.DataFrame, outcome: str, covariates: list[str], term: str = "CV_BIS") -> dict[str, Any]:
    subset = data.dropna(subset=[outcome] + covariates).copy()
    subset = subset[subset[outcome] > 0]
    if subset.empty:
        raise ValueError(f"No analyzable rows for outcome {outcome}.")

    x = sm.add_constant(subset[covariates].astype(float), has_constant="add")
    y = subset[outcome].astype(float)
    fit = GLM(y, x, family=families.Gamma(link=families.links.Log())).fit()

    beta = float(fit.params[term])
    se = float(fit.bse[term])
    return {
        "n": int(len(subset)),
        "outcome": outcome,
        "term": term,
        "beta": beta,
        "se": se,
        "exp_beta": float(np.exp(beta)),
        "exp_ci_low": float(np.exp(beta - 1.96 * se)),
        "exp_ci_high": float(np.exp(beta + 1.96 * se)),
        "p": float(fit.pvalues[term]),
        "aic": float(fit.aic),
    }


def fit_interaction(data: pd.DataFrame) -> dict[str, Any]:
    work = data.dropna(subset=MODEL2_COVARIATES + ["TWA_MAC", "ome_per_kg_hr"]).copy()
    work = work[work["ome_per_kg_hr"] > 0]
    work["CV_BIS_x_TWA_MAC"] = work["CV_BIS"] * work["TWA_MAC"]
    covariates = MODEL2_COVARIATES + ["TWA_MAC", "CV_BIS_x_TWA_MAC"]
    x = sm.add_constant(work[covariates].astype(float), has_constant="add")
    y = work["ome_per_kg_hr"].astype(float)
    fit = GLM(y, x, family=families.Gamma(link=families.links.Log())).fit()

    terms: dict[str, Any] = {}
    for term in ["CV_BIS", "TWA_MAC", "CV_BIS_x_TWA_MAC"]:
        beta = float(fit.params[term])
        se = float(fit.bse[term])
        terms[term] = {
            "beta": beta,
            "se": se,
            "exp_beta": float(np.exp(beta)),
            "exp_ci_low": float(np.exp(beta - 1.96 * se)),
            "exp_ci_high": float(np.exp(beta + 1.96 * se)),
            "p": float(fit.pvalues[term]),
        }
    return {"n": int(len(work)), "aic": float(fit.aic), "terms": terms}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dataset", type=Path, required=True, help="Local reconstructed analysis CSV.")
    parser.add_argument("--mac-source", type=Path, default=None, help="Optional CSV with caseid and TWA_MAC.")
    parser.add_argument("--out", type=Path, default=Path("outputs/reanalysis.json"), help="Output JSON path.")
    args = parser.parse_args()

    df = add_model_columns(load_dataset(args.analysis_dataset, args.mac_source))
    primary = df[(df["ome_per_kg_hr"] > 0) & df["asa"].notna()].copy()
    mac_available = primary[primary["TWA_MAC"].notna()].copy()
    tiva = mac_available[mac_available["TWA_MAC"] < 0.3].copy()
    volatile = mac_available[mac_available["TWA_MAC"] >= 0.3].copy()

    results: dict[str, Any] = {
        "cohort_counts": {
            "processed_rows": int(len(df)),
            "primary_complete_case": int(len(primary)),
            "mac_available": int(len(mac_available)),
            "tiva_mac_lt_0_3": int(len(tiva)),
            "volatile_supplemented_mac_ge_0_3": int(len(volatile)),
            "missing_mac_in_primary": int(primary["TWA_MAC"].isna().sum()),
        },
        "core_models": {
            "pooled_model2": fit_gamma(primary, "ome_per_kg_hr", MODEL2_COVARIATES),
            "mac_available_pooled_model2": fit_gamma(mac_available, "ome_per_kg_hr", MODEL2_COVARIATES),
            "tiva_model2": fit_gamma(tiva, "ome_per_kg_hr", MODEL2_COVARIATES),
            "volatile_supplemented_model2": fit_gamma(volatile, "ome_per_kg_hr", MODEL2_COVARIATES),
            "continuous_mac_interaction": fit_interaction(mac_available),
        },
    }

    if "rftn_rate_mcg_kg_hr" in mac_available.columns:
        rftn = mac_available[mac_available["rftn_rate_mcg_kg_hr"] > 0].copy()
        results["raw_remifentanil_models"] = {
            "pooled": fit_gamma(rftn, "rftn_rate_mcg_kg_hr", MODEL2_COVARIATES),
            "tiva": fit_gamma(rftn[rftn["TWA_MAC"] < 0.3], "rftn_rate_mcg_kg_hr", MODEL2_COVARIATES),
            "volatile_supplemented": fit_gamma(rftn[rftn["TWA_MAC"] >= 0.3], "rftn_rate_mcg_kg_hr", MODEL2_COVARIATES),
        }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(clean_json(results), indent=2) + "\n", encoding="utf-8")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
