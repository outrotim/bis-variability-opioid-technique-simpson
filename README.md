# BIS Variability, Anesthetic Technique, and Intraoperative Opioid Administration

This repository provides the minimal public code and aggregate estimates for the manuscript:

**Anesthetic Technique Drives a Simpson's Paradox in the Association Between Bispectral Index Variability and Intraoperative Opioid Administration: A Retrospective Cohort Study**

The repository is intentionally small. It supports transparency for the model specification, core aggregate estimates, and public figure regeneration without redistributing patient-level data.

## Repository Contents

| Path | Purpose |
|---|---|
| `README.md` | Repository scope, data availability, usage, licenses, and citation placeholder. |
| `LICENSE` | MIT license for code. |
| `LICENSE-CC-BY-4.0.md` | CC BY 4.0 license notice for aggregate estimates and non-code parameters. |
| `requirements.txt` | Minimal Python dependencies. |
| `aggregate_estimates.json` | Aggregate-only cohort counts and model estimates reported in the manuscript. |
| `src/reproduce_core_models.py` | Recomputes the core Gamma GLM models from a user-provided local analysis dataset. |
| `src/generate_public_figures.py` | Regenerates public aggregate versions of the study flow diagram and forest plot. |

## Data Availability

The source database used for the study is the public VitalDB database. Users should obtain data directly from VitalDB under its applicable terms and reconstruct any analysis dataset locally.

This repository does **not** redistribute:

- patient-level case data;
- derived patient-level analysis datasets;
- high-resolution physiologic waveform files;
- intermediate per-case feature tables;
- local audit outputs containing file paths or internal workflow metadata.

The file `aggregate_estimates.json` contains only aggregate cohort counts and model estimates that are already reported or directly support reported manuscript tables and figures.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Reproduce Core Models

The model script expects a locally reconstructed analysis dataset. The dataset is not included in this repository.

```bash
python src/reproduce_core_models.py \
  --analysis-dataset path/to/final_analysis_dataset.csv \
  --mac-source path/to/round5_sensitivity_data.csv \
  --out outputs/reanalysis.json
```

If `TWA_MAC` is already present in the analysis dataset, `--mac-source` can be omitted.

Expected core columns include:

`CV_BIS`, `ome_per_kg_hr`, `age`, `sex`, `bmi`, `asa`, `opdur_min` or `opstart`/`opend`, `TWA_BIS`, and `TWA_MAC` for technique-aware analyses.

For raw remifentanil-rate sensitivity, include `rftn_rate_mcg_kg_hr`.

## Regenerate Public Figures

The public figure script uses only `aggregate_estimates.json`.

```bash
python src/generate_public_figures.py --outdir outputs
```

This regenerates:

- `outputs/figure1_flow_public.png/.pdf/.svg`
- `outputs/figure2_forest_public.png/.pdf/.svg`

The public Figure 2 is an aggregate forest-plot version. It does not recreate patient-level scatter plots or restricted cubic spline curves because those require patient-level data.

## Interpretation Caveats

- The study is observational. Estimates should be interpreted as associations, not causal effects.
- BIS variability should not be interpreted as a stand-alone nociception or opioid-requirement signal.
- The key result is technique-dependent: pooled estimates differ from within-technique estimates.
- External users who reconstruct the dataset may see minor numerical differences if preprocessing, signal-quality filtering, or opioid conversion choices differ.

## Licenses

- Code in `src/` is released under the MIT License.
- Aggregate estimates and non-code parameters in `aggregate_estimates.json` are released under CC BY 4.0.
- The VitalDB source data are not redistributed and remain subject to VitalDB's own terms.

## Citation

Please cite the manuscript after publication:

> [Citation to be added after publication]
