# Data dictionary

## Inputs and source notes

| File | Content |
|---|---|
| `data/age_control.csv` | Thirty-one measured or tuned age controls. Measured controls correspond to Supplementary Table 2; MIS/MISP/TP tuning controls correspond to Supplementary Table 3. |
| `data/mtd_events.csv` | Identifiers and original and event-free depths for 11 MTD events. |
| `data/environment_proxy_metrics.csv` | Processed 0-236 ka RSL and EASM series used by the environmental analysis. |
| `data/processed/bacon_model.rds` | Archived Bacon model used by the default analysis. |
| `data/source_notes.md` | Input-field definitions, units and data provenance. |

## Outputs

| File | Content | Supplementary mapping |
|---|---|---|
| `results/bacon/event_age_posterior_draws_wide.csv` | Joint event-age posterior draws | - |
| `results/event_age_summary.csv` | Event-age posterior summary | Supplementary Table 5 |
| `results/bacon/age_depth_model_1cm.csv` | One-centimetre age-depth summary | Supplementary Data 2 |
| `results/bacon/ZJ05_LSR_step_10ka_pseudoPosterior.csv` | 10 ka-window linear sedimentation-rate summary | - |
| `results/bacon/stage_sedimentation_rate_posterior_summary.csv` | Stage-scale sedimentation-rate summary | - |
| `results/bacon/mtd_linked_sedimentation_rate_posterior_summary.csv` | Event- and group-linked sedimentation-rate summary | - |
| `results/categorized_timeseries.csv` | Environmental time series and category labels | - |
| `results/univariate_enrichment_summary.csv` | Univariate enrichment and incidence summary | Supplementary Table 6 |
| `results/bivariate_enrichment_summary.csv` | Bivariate enrichment and incidence summary | Supplementary Data 3 |
| `results/univariate_draws_ALL.csv`; `results/midRSL_bivariate_draws_ALL.csv` | Posterior enrichment draws | - |
| `results/posterior_mean_phase_shift_null_summary.csv`; `results/posterior_mean_phase_shift_null_long.csv` | Phase-shift null summary and phase-level results | Supplementary Data 5 |
| `results/circular_shift_robustness_univariate.csv`; `results/circular_shift_robustness_midRSL_bivariate.csv`; `results/circular_shift_by_offset_univariate.csv`; `results/circular_shift_by_offset_midRSL_bivariate.csv` | Circular-shift summaries and results by offset | Supplementary Data 4 |

## Rate sources

Bacon sedimentation-rate outputs are calculated from the selected Bacon model
as `(depth_right - depth_left) / (age_right - age_left)`. The RSL and EASM
rate-of-change columns are processed inputs read from
`data/environment_proxy_metrics.csv`; their definitions and units are given in
`data/source_notes.md`.
