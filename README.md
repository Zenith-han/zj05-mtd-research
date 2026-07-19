# ZJ05 MTD analysis package

This package contains the R code, inputs and results for the 236 ka ZJ05 MTD
analysis.

## Requirements

- R 4.5.1
- rbacon 3.5.2

Install rbacon if needed:

```r
install.packages("rbacon")
```

## Run

Run commands from the package root. The default analysis uses
`data/processed/bacon_model.rds` and writes 17 CSV files under `results/`:

```r
source("code/run_all.R")
```

| Task | Command | Output |
|---|---|---|
| Run the archived model explicitly | `Rscript code/run_all.R --model-source=archived` | `results/` |
| Build a new Bacon model only | `Rscript code/run_all.R --build-bacon-model` | `results/rebuilt_runs/<run-id>/` |
| Build a new model and run the analysis | `Rscript code/run_all.R --build-bacon-model --model-source=rebuilt` | `results/rebuilt_runs/<run-id>/` |
| Rerun an existing rebuilt model | `Rscript code/run_all.R --model-source=rebuilt --run-id=<id>` | `results/rebuilt_runs/<id>/` |

New builds use a timestamp run ID unless `--run-id=<id>` is supplied.

## Files

```text
code/            R scripts
data/            inputs and source notes
data/processed/  archived Bacon model
results/         archived-model results
```

See `DATA_DICTIONARY.md` for the input, output and supplementary-material
mapping. Source data and units are documented in `data/source_notes.md`.
