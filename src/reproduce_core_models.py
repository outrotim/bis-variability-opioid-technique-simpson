#!/usr/bin/env python3
"""Reproduce the core regimen-aware models from a locally rebuilt dataset.

No patient-level data are distributed with this repository. Users must obtain
VitalDB under its applicable terms and recreate the locked analysis variables
before running this script.
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


CORE_COVARIATES = [
    "CV_BIS",
    "age",
    "male",
    "bmi",
    "asa_num",
    "opdur_min",
    "TWA_BIS",
]


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    return value


def load_dataset(analysis_dataset: Path, mac_source: Path | None = None) -> pd.DataFrame:
    data = pd.read_csv(analysis_dataset)
    if "TWA_MAC" not in data.columns:
        if mac_source is None:
            raise ValueError(
                "TWA_MAC is missing; provide --mac-source or include it in the analysis dataset."
            )
        mac = pd.read_csv(mac_source, usecols=["caseid", "TWA_MAC"])
        data = data.merge(mac, on="caseid", how="left", validate="one_to_one")
    return data


def add_model_columns(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    if "male" not in out.columns:
        if "sex" not in out.columns:
            raise ValueError("The dataset must contain male or sex.")
        out["male"] = out["sex"].astype(str).str.upper().eq("M").astype(int)
    if "asa_num" not in out.columns:
        if "asa" not in out.columns:
            raise ValueError("The dataset must contain asa_num or asa.")
        out["asa_num"] = pd.to_numeric(out["asa"], errors="coerce")
    if "opdur_min" not in out.columns:
        if {"opstart", "opend"}.issubset(out.columns):
            out["opdur_min"] = (out["opend"] - out["opstart"]) / 60.0
        else:
            raise ValueError("The dataset must contain opdur_min or opstart/opend.")
    return out


def resolve_outcome(data: pd.DataFrame, preferred: str, fallback: str) -> str:
    if preferred in data.columns:
        return preferred
    if fallback in data.columns:
        return fallback
    raise ValueError(f"Missing outcome column: {preferred} or {fallback}.")


def fit_gamma(
    data: pd.DataFrame,
    outcome: str,
    covariates: list[str] | None = None,
) -> dict[str, Any]:
    covariates = covariates or CORE_COVARIATES
    subset = data.dropna(subset=[outcome] + covariates).copy()
    subset = subset[subset[outcome] > 0]
    if subset.empty:
        raise ValueError(f"No analyzable positive observations for {outcome}.")
    design = sm.add_constant(subset[covariates].astype(float), has_constant="add")
    fit = GLM(
        subset[outcome].astype(float),
        design,
        family=families.Gamma(link=families.links.Log()),
    ).fit(cov_type="HC3")
    beta = float(fit.params["CV_BIS"])
    se = float(fit.bse["CV_BIS"])
    return {
        "n": int(len(subset)),
        "outcome": outcome,
        "term": "CV_BIS",
        "covariance": "HC3",
        "beta": beta,
        "se": se,
        "ratio": float(np.exp(beta)),
        "ci_low": float(np.exp(beta - 1.96 * se)),
        "ci_high": float(np.exp(beta + 1.96 * se)),
        "p": float(fit.pvalues["CV_BIS"]),
    }


def fit_outcome_suite(
    primary: pd.DataFrame,
    mac_available: pd.DataFrame,
    outcome: str,
) -> dict[str, Any]:
    tiva = mac_available[mac_available["TWA_MAC"] < 0.3]
    volatile = mac_available[mac_available["TWA_MAC"] >= 0.3]
    context = mac_available.copy()
    context["volatile_regimen"] = (context["TWA_MAC"] >= 0.3).astype(int)
    return {
        "pooled_all": fit_gamma(primary, outcome),
        "pooled_mac_available": fit_gamma(mac_available, outcome),
        "technique_adjusted": fit_gamma(
            context,
            outcome,
            CORE_COVARIATES + ["volatile_regimen"],
        ),
        "tiva": fit_gamma(tiva, outcome),
        "volatile_supplemented": fit_gamma(volatile, outcome),
    }


def run_core_models(data: pd.DataFrame) -> dict[str, Any]:
    work = add_model_columns(data)
    ome = resolve_outcome(work, "ome_reset_safe_p99", "ome_per_kg_hr")
    remifentanil = resolve_outcome(
        work,
        "rftn_rate_reset_safe_p99",
        "rftn_rate_mcg_kg_hr",
    )
    eligible = work["eligible"].fillna(False).astype(bool) if "eligible" in work else True
    primary = work[eligible & work["asa_num"].notna() & work[ome].gt(0)].copy()
    mac_available = primary[primary["TWA_MAC"].notna()].copy()
    remif_primary = primary[primary[remifentanil].gt(0)].copy()
    remif_mac = remif_primary[remif_primary["TWA_MAC"].notna()].copy()

    return {
        "cohort_counts": {
            "primary_analysis": int(len(primary)),
            "mac_available": int(len(mac_available)),
            "tiva": int((mac_available["TWA_MAC"] < 0.3).sum()),
            "volatile_supplemented": int((mac_available["TWA_MAC"] >= 0.3).sum()),
            "remifentanil_positive": int(len(remif_primary)),
        },
        "ome_models": fit_outcome_suite(primary, mac_available, ome),
        "remifentanil_only_models": fit_outcome_suite(
            remif_primary,
            remif_mac,
            remifentanil,
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dataset", type=Path, required=True)
    parser.add_argument("--mac-source", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=Path("outputs/reanalysis.json"))
    args = parser.parse_args()

    results = run_core_models(load_dataset(args.analysis_dataset, args.mac_source))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(clean_json(results), indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
