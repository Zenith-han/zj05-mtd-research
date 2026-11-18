# ZJ05 MTD analysis package

Version 1.0.2 accompanies *Submarine landslide deposits concentrate at intermediate sea levels on a monsoon-influenced margin*. The analysis covers 11 mass-transport deposits (MTDs) in core ZJ05 over the past 237 kyr.

## Reproduce the published results

Download or clone the repository. From the package root, run the following command.

```sh
Rscript --vanilla code/run_all.R --model-source=archived
```

The archived workflow uses base R and reads the compressed model directly. It writes 17 CSV files under `results/`, including the eight summary tables distributed here. The other nine files contain the age–depth curve, 10-kyr sedimentation rates, the classified environmental series, posterior draws and shift-level results.

## Inputs and methods

| Input | Contents |
| --- | --- |
| `data/age_control.csv` | 31 age controls. Depth is in centimetres on the event-free depth scale. Ages, errors and reservoir corrections are in years. The ten AMS radiocarbon ages use Marine20 with a regional correction of −149 ± 40 years. |
| `data/mtd_events.csv` | Event identifiers, top and base depths in metres below seafloor, and event-free depths in centimetres. |
| `data/environment_proxy_metrics.csv` | Relative sea-level proxy (RSL), East Asian summer monsoon index (EASM-PC1) and their rates of change on a 0.2-kyr grid. |
| `data/processed/bacon_model.rds` | Archived Bacon chronology with 72,000 posterior age–depth curves. The model was fitted with R 4.5.3 and rbacon 3.5.2, using seed 20264622 and a normal measurement-error likelihood. |

The RSL proxy follows Spratt and Lisiecki (2016). EASM-PC1 combines the loess grain-size record of Yang and Ding (2014), the magnetic-susceptibility record of Sun et al. (2006) and the sign-adjusted speleothem oxygen-isotope record of Cheng et al. (2016). These sources are cited in the article. The input table supplies the continuous proxy series and their rates in metres per kyr (RSL) or index units per kyr (EASM). Rates use the difference between the younger and older proxy values 0.6 kyr on either side of each grid point, divided by 1.2 kyr. The edge grid points inherit the nearest available rate. The code computes variability as the sample standard deviation in a centered 10-kyr window, using the available values at the boundaries.

The environmental analysis samples 5,000 complete sets of event ages from the archived model. The eight distributed tables report event ages, sedimentation rates by stage and event group, univariate and bivariate enrichment and incidence, full-cycle phase-shift tests and bounded-shift sensitivity. Incidence is expressed in recorded events per kyr. The age controls correspond to Supplementary Tables 2 and 3, event ages to Supplementary Table 5 and univariate enrichment to Supplementary Table 6. The age–depth curve corresponds to Supplementary Data 2, bivariate results to Supplementary Data 3, bounded shifts to Supplementary Data 4 and full-cycle shifts to Supplementary Data 5.

## Fit a new chronology

With the R package `rbacon` installed, run `Rscript code/run_all.R --build-bacon-model --model-source=rebuilt`. The new fit and its analysis outputs are written to a separate directory under `results/rebuilt_runs/`. MCMC draws can vary across software versions.

## Licensing

Custom code is distributed under the MIT License (see `LICENSE`). Author-generated processed data and analysis outputs are distributed under CC BY 4.0 (see `LICENSE-DATA.md`). The article identifies the GMGS datasets subject to access restrictions and the original sources of the published environmental records.
