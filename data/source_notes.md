# Source Notes

The processed metric file used by the analysis is:

```text
data/environment_proxy_metrics.csv
```

The palaeoenvironmental records used to construct or contextualize these metrics include:

- Pleistocene relative sea-level stack: Spratt and Lisiecki (2016).
- ODP Site 1146 planktonic foraminiferal oxygen-isotope record:
  Caballero-Gill et al. (2012).
- Chinese loess grain-size stack: Yang and Ding (2014).
- Chinese loess magnetic-susceptibility stack: Sun et al. (2006).
- Chinese speleothem oxygen-isotope stack: Cheng et al. (2016).

The metric file contains the following core fields:

```text
age_ka,RSL_State_Value,EASM_Intensity_Value,RSL_RateofChange_dka_h1.2,EASM_RateofChange_dka_h1.2
```

All five fields are required. The two rate-of-change columns are processed
metrics calculated as manuscript Equation (1):

```text
[X(t - 0.6 ka) - X(t + 0.6 ka)] / 1.2 ka
```

`RSL_RateofChange_dka_h1.2` is in m/ka, and
`EASM_RateofChange_dka_h1.2` is in EASM index units/ka. Script
`04_build_environment_timeseries.R` reads these CSV columns directly and
stops if either is missing; the public workflow does not recalculate them.

Script `04_build_environment_timeseries.R` derives both variability series
from the continuous state/intensity metrics using the symmetric 51-point
window defined in `code/analysis_config.R`. The input file intentionally does
not contain precomputed variability columns. The script then rebuilds all
category labels from the continuous metrics.
