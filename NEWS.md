# inough 0.1.0

Initial public release.

## Pipeline

* `inough_signals()` — fits a probit GLMM on accuracy and computes
  response-bias indicators and per-participant LZ complexity.
* `inough_detect()` — three-stage detection: bail-out screening,
  rolling-window candidate localization (`screen()`), and formal
  Welch t-tests with BH-FDR / LZ-informed Bayesian local FDR
  (`test_chunks()`).
* `flags()` — tidy data frame of flagged trials, suitable for
  anti-joining against the input data.
* `report()` — self-contained HTML report with embedded JSON.
* `lz_complexity()` — entropy-corrected LZ76 complexity (exported).

## Configuration

* `inough_control()` and `inough_heuristics()` expose all tuning
  parameters with documented defaults.
