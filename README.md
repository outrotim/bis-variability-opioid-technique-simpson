# BIS Variability and Anesthetic-Regimen-Specific Opioid Administration

This repository contains the minimal public code and aggregate estimates for:

**Bispectral index variability reveals anesthetic-regimen-specific patterns of intraoperative opioid administration: a retrospective cohort study**

The study asks whether routine BIS variability has a consistent medication-management meaning across anesthetic regimens. The public package supports inspection of the model formulas, local re-estimation after an independent VitalDB reconstruction, and regeneration of aggregate figures. It intentionally excludes patient-level and event-window data.

## Repository Contents

| Path | Purpose |
|---|---|
| `README.md` | Scope, data availability, use, licenses, and citation placeholder. |
| `LICENSE` | MIT license for code. |
| `LICENSE-CC-BY-4.0.md` | CC BY 4.0 notice for aggregate estimates and non-code parameters. |
| `requirements.txt` | Minimal Python dependencies. |
| `aggregate_estimates.json` | Aggregate cohort flow, preprocessing parameters, quality-audit counts, and reported estimates. |
| `src/reproduce_core_models.py` | Fits pooled, technique-adjusted, and technique-stratified Gamma GLMs with HC3 covariance from a locally reconstructed dataset. |
| `src/generate_public_figures.py` | Regenerates the aggregate study-flow and forest figures. |

## Data Availability

The source data are available from the [VitalDB repository](https://doi.org/10.13026/czw8-9p62) and should be obtained directly under its applicable terms. This repository does not redistribute:

- patient-level source or derived data;
- high-resolution physiologic or infusion-pump tracks;
- five-minute event-window records;
- intermediate feature tables;
- local paths, credentials, or internal workflow metadata.

`aggregate_estimates.json` contains only non-identifiable parameters, counts, and estimates reported in the manuscript or supplement.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Reproduce Core Models

The model script expects a local dataset reconstructed under the manuscript specification. The dataset is not included.

```bash
python src/reproduce_core_models.py \
  --analysis-dataset path/to/local_rebuilt_analysis.csv \
  --out outputs/reanalysis.json
```

If `TWA_MAC` is stored separately, add `--mac-source path/to/local_mac.csv`. The required analysis fields are:

`CV_BIS`, `age`, `male` or `sex`, `bmi`, `asa_num` or `asa`, `opdur_min` or `opstart`/`opend`, `TWA_BIS`, and `TWA_MAC`.

The preferred locked outcomes are `ome_reset_safe_p99` and `rftn_rate_reset_safe_p99`; the aliases `ome_per_kg_hr` and `rftn_rate_mcg_kg_hr` are also accepted. When present, `eligible` applies the strict signal-quality mask.

## Regenerate Aggregate Figures

```bash
python src/generate_public_figures.py --outdir outputs
```

This creates:

- `outputs/figure1_flow_public.png/.pdf/.svg`
- `outputs/figure2_forest_public.png/.pdf/.svg`

The aggregate forest plot recreates the reported pooled, technique-adjusted, and technique-stratified estimates. The patient-level distribution panel in the manuscript cannot be regenerated without locally reconstructing the source data and is therefore not included here.

## Interpretation

The observed information is anesthetic-regimen specific. The positive pooled association largely reflects separation between regimen-level operating states, whereas adjusted and within-regimen estimates are modestly inverse or null. These results support regimen-aware interpretation of BIS dynamics in monitoring and perioperative algorithm development; they do not define a causal opioid-titration rule or bedside cutoff.

## Licenses

- Code in `src/` is available under the MIT License.
- Aggregate estimates and non-code parameters are available under CC BY 4.0.
- VitalDB data are not redistributed and remain subject to their source terms.

## Citation

Please cite the manuscript after publication:

> Citation to be added after publication.
