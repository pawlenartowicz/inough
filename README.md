# inough

[![R-CMD-check](https://github.com/pawlenartowicz/inough/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pawlenartowicz/inough/actions/workflows/R-CMD-check.yaml)

**Inattention detection for long psychophysical tasks.**

`inough` is an R package that detects episodes of participant inattention in trial-based experiments (200+ trials). It combines response-pattern heuristics for locating candidate regions with formal Welch t-tests (BH-FDR corrected) for statistical testing, and produces a self-contained HTML report for visual inspection.

## Installation

From GitHub (recommended):

```r
# install.packages("remotes")
remotes::install_github("pawlenartowicz/inough")
```

From a local source clone:

```r
install.packages("/path/to/inough", repos = NULL, type = "source")
# or, equivalently:
# R CMD INSTALL /path/to/inough
```

**Dependencies:** `lme4`, `ggplot2`, `patchwork`, `rlang`, `jsonlite`

## Quick start

```r
library(inough)

# 1. Extract signals — fits a probit GLMM on accuracy,
#    computes response bias indicators and LZ complexity
signals <- inough_signals(
  df       = my_data,
  correct  = correct ~ Stim + Weight + Block,
  response = response ~ answ,
  id       = "ID"
)

# 2. Detect inattention episodes
detected <- inough_detect(signals)

# 3. Inspect
print(detected)
summary(detected)

# 4. Get flagged trials for downstream filtering
flagged <- flags(detected)
clean_data <- my_data[!interaction(my_data$ID, seq_len(nrow(my_data))) %in%
                       interaction(flagged$id, flagged$trial_idx), ]

# 5. Generate interactive HTML report
report(detected, file = "inough_report.html")

# 6. Diagnostic plots for specific participants
plot(detected, id = "10001")
```

### `inough_control()` — detection parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `window_size` | 3 | Half-width of rolling window (total = 2w+1 = 7 trials) |
| `sd_threshold` | 2 | SDs above chance for `\|roll_resp\|` screening threshold |
| `window_weight` | `"uniform"` | Window weighting: `"uniform"` or `"triangular"` |
| `min_chunk` | 6 | Minimum chunk length (trials) to retain after merge |
| `lz_threshold` | 0.2 | LZ complexity below this triggers whole-participant bail-out |
| `comparison` | `"clean"` | t-test reference: `"clean"` (non-flagged trials) or `"rest"` |

The screening threshold is computed analytically: `sd_threshold × sqrt(sum(w²)) / sum(w)` where `w` are the window weights. For uniform weights this simplifies to `sd_threshold / sqrt(2 × window_size + 1)`.

```r
# Example: explicit control (these are the defaults)
detected <- inough_detect(
  signals,
  control = inough_control(
    window_size   = 3,
    sd_threshold  = 2,
    window_weight = "uniform",
    min_chunk     = 6,
    lz_threshold  = 0.2,
    comparison    = "clean"
  )
)
```

### `inough_heuristics()` — post-detection refinements

| Parameter | Default | Description |
|-----------|---------|-------------|
| `boundary_mode` | `"heuristic"` | `"heuristic"` (default): end = fixed extension, start = walk back to first matching stereotype trial. `"fixed"`: extend symmetrically (±w for uniform, ±⅔w for triangular). |
| `spurious` | `TRUE` | Enable spurious-start accuracy trimming |
| `spurious_n` | 6 | Number of trials at chunk start to inspect |
| `spurious_k` | `NULL` | Explicit: flag if >= k correct out of `spurious_n`. When `NULL`, computed from `spurious_n` (see below). |
| `min_left` | 7 | Minimum trials remaining after trimming |

When `spurious_k` is not set, it is computed as `ceiling((0.5 + 1.5 × sqrt(0.25/n)) × n)` (≈ 1.5 SE above chance for `n = spurious_n`).

```r
# Example: override the spurious-trimming threshold
detected <- inough_detect(
  signals,
  heuristics = inough_heuristics(
    spurious   = TRUE,
    spurious_k = 5,
    spurious_n = 6,
    min_left   = 7
  )
)
```

## How it works

The pipeline has three stages: **bail-out** for grossly disengaged participants, **screening** to locate candidate inattention regions using response-pattern heuristics, and **formal testing** with LZ-informed Bayesian local FDR.

```
 Input data (trials x participants)
              |
              v
 ┌──────────────────────────────┐
 │  Signal extraction           │
 │  Probit GLMM -> acc resids   │
 │  Lag-1 autocorrelation       │
 │  LZ complexity (per person)  │
 └──────────────────────────────┘
              |
              v
 ┌──────────────────────────────┐
 │  Bail-out                    │
 │  LZ < threshold -> flag all  │
 │  Acc <= chance+SE -> flag all │
 └──────────────────────────────┘
              |  (passed)
              v
 ┌──────────────────────────────┐
 │  Stage 1: Screening          │
 │  Rolling |resp_lag1| mean    │
 │  Threshold -> candidate runs │
 │  Merge + min-length filter   │
 │  Boundary extension          │
 └──────────────────────────────┘
              |
              v
 ┌──────────────────────────────┐
 │  Heuristics (optional)       │
 │  Spurious accuracy trimming  │
 └──────────────────────────────┘
              |
              v
 ┌──────────────────────────────┐
 │  Stage 2: Formal testing     │
 │  Welch t-test per chunk      │
 │  LZ-informed local FDR       │
 │  -> flagged episodes         │
 └──────────────────────────────┘
              |
              v
 ┌──────────────────────────────┐
 │  HTML report                 │
 │  Per-participant browser     │
 │  Diagnostic plots & export   │
 └──────────────────────────────┘
```

### Signal extraction

`inough_signals()` prepares two independent behavioral signals from raw trial data:

- **Accuracy residuals** — Pearson residuals from a probit GLMM (`lme4::glmer`) with task-specific design predictors and per-participant random effects. These residuals capture trial-level deviations from expected accuracy, controlling for stimulus difficulty and learning.
- **Lag-1 response indicator** — +1 if the participant repeated their previous response, -1 if they switched. Summarized per participant as **LZ complexity** (Lempel-Ziv, entropy-corrected and permutation-normalized), which indexes response stereotypy: low LZ indicates repetitive button-pressing, high LZ indicates random-like responding.

### Bail-out

Before running the trial-level pipeline, participants whose data is too degraded for meaningful analysis are flagged wholesale:

- **LZ < threshold** — response sequence too stereotyped (e.g., pressing the same button repeatedly)
- **Accuracy <= chance + SE** — performing at or below chance across the session

These participants skip the screening and testing stages; all their trials are flagged with `flag_type = "bailout"`.

### Stage 1: Heuristic screening

Response stereotypy locates *where* to look. A centered rolling window computes the weighted mean of lag-1 response indicators. Regions where `|rolling mean| > sd_threshold × window_SD` are candidate inattention zones — the participant is repeating or alternating responses beyond what chance predicts.

Candidates are merged if overlapping, filtered by minimum length, then extended via the configured boundary mode:

- **Fixed** (default): extend ±`window_size` (uniform) or ±⅔`window_size` (triangular)
- **Heuristic**: end boundary uses fixed extension; start boundary walks backwards contiguously, following resp_lag1 trials that match the chunk's stereotype direction (sign of mean rolling response)

### Spurious accuracy heuristic

When enabled, inspects the first N trials of each extended chunk. If accuracy is suspiciously high (>= k/n correct), those trials are excluded from the t-test but kept in the chunk boundary for reporting. This prevents "lucky starts" from masking an accuracy drop.

### Stage 2: Formal testing

Each candidate chunk gets a one-sided Welch t-test on accuracy residuals (chunk vs comparison set), asking: is accuracy inside this chunk significantly *lower* than outside?

**Local FDR.** Each chunk receives a posterior probability of being a true null (no real inattention) via Bayesian local FDR:

```
lfdr = pi0 * f0(t) / [ pi0 * f0(t) + pi1 * f1(t) ]
```

- **pi0, pi1** — LZ-informed prior. `pi0 = min(1, LZ)`, `pi1 = 1 - pi0`. Low LZ (stereotyped) gives a permissive prior; high LZ (random-like) gives a skeptical prior.
- **f0** — null density: central t-distribution at observed df
- **f1** — alternative density: non-central t with noncentrality `ncp = -(acc_clean - 0.5) / SE_theory`, encoding the expectation that inattentive responding drops to chance (50%)

A chunk is flagged when `lfdr < fdr_alpha`. BH-adjusted p-values are also computed for reference.

### HTML report

`report()` generates a self-contained HTML file with:

- Participant browser with color-coded status (bail-out / flagged / clean)
- Sortable and filterable participant list
- Per-participant diagnostic plots (accuracy residuals, rolling response signal, flagged regions)
- Inline chunk table with all tested chunks (significant highlighted, non-significant dimmed), lfdr, effect size, spurious trimming info
- Summary statistics and distribution plots
- CSV export of all flagged chunks

## Output structure

`inough_detect()` returns an `inough_detected` object containing:

| Field | Description |
|-------|-------------|
| `$signals` | The input `inough_signals` object (trial + participant data) |
| `$bailout` | Data frame of bailed-out participants (id, reason, lz_value, accuracy) |
| `$candidates` | Raw candidate regions from screening |
| `$chunks` | Merged/filtered/extended chunks passed to testing |
| `$results` | All tested chunks with t_stat, df, p_raw, p_adj, lfdr, effect_size, significant |
| `$control` | The `inough_control` parameters used |
| `$heuristics` | The `inough_heuristics` parameters used |
| `$fdr_alpha` | The FDR threshold used |

`flags()` extracts a tidy data frame of all flagged trials (both bail-out and chunk-based), suitable for `dplyr::anti_join()` or similar filtering.

## License

GPL (>= 3)

## Authors

**Pawe&#x142; Lenartowicz** ([ORCID](https://orcid.org/0000-0002-6906-7217)) & **Maja Willard** ([ORCID](https://orcid.org/0009-0005-3482-1119))
