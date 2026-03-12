# Master Thesis — ProUni Labour Market Effects in Brazil

LaTeX + R econometrics thesis repo. Topic: causal effects of Brazil's ProUni scholarship programme on formal labour-market outcomes, with heterogeneity analysis by gender, ethnicity, and region.

---

## Project Structure

```
Thesis_Main.tex        # Master LaTeX document
Intro.tex              # Introduction
Context.tex            # Institutional context
Overview.tex           # Literature review
Experiment.tex         # Empirical strategy & data
Results.tex            # Main results
Discussion.tex         # Discussion & heterogeneity
Conclusion.tex         # Conclusion
Appendix_1.tex         # Robustness / supplementary
Bibliography.bib       # BibTeX references
untitled4.R            # FROZEN backup — do not modify
untitled5.R            # Active improved script — build from untitled4.R
Figures/               # Output figures (PDF/PNG)
Tables/                # Output tables (TeX)
Documentation/         # Method notes & data docs
```

---

## Tech Stack

- **LaTeX** (Overleaf-compatible) — thesis writing
- **R** — econometric analysis
  - `did` — Callaway & Sant'Anna (2021) staggered DiD
  - `contdid` — continuous-treatment DiD (Callaway, Goodman-Bacon, Sant'Anna)
  - `fixest` — TWFE benchmarks
  - `basedosdados` — BigQuery access to RAIS & ProUni admin data
  - `geobr` — Brazilian geographic shapefiles
  - `tidyverse`, `data.table`, `ggplot2`, `modelsummary`, `HonestDiD`

---

## Econometric Framework

**Primary estimator:** Callaway & Sant'Anna (2021) via `did::att_gt()`
- Treatment: ProUni scholarship roll-out across Brazilian microregions (staggered, 2005+)
- Outcomes: log formal wages (RAIS), formal employment probability, employment count
- Heterogeneity: gender, ethnicity (branco / pardo / preto), macro-region
- Dose bins: Low (≤ median, 167 microregions) vs High (> median, 140 microregions); 112 never/not-yet-treated controls

**Continuous treatment:** `contdid` for dose-response; fallback to `fixest` spline/polynomial specs

**Benchmark:** TWFE via `fixest::feols()` — document sign-reversal bias vs CS estimates

**Data sources:** RAIS employer–employee panel (2005–2019), ProUni admin records, IBGE Census 2010

---

## Critical Identification Assumptions

1. **Parallel trends** (conditional on covariates) — always present pre-trends tests; all current p-values > 0.10 ✓
2. **No anticipation** — ProUni announced before implementation; discuss in Section 4
3. **SUTVA** — buffer-zone exclusion test already passes (Low/Buffer ≈ Low/Wages) ✓

---

## untitled5.R — Task Specification

`untitled5.R` is a **surgical improvement** of `untitled4.R`. It must NOT re-download data from BigQuery (reuse objects from untitled4.R or saved RDS files). It must NOT change the identification strategy or sample. It must produce outputs to `Figures_Final/` and `Tables_Final/` with a `_v5` suffix to avoid overwriting backup results.

### Fix 1 — contdid (CRITICAL)
The current error is: `Assertion on 'target_parameter' failed: Must be element of set {'level','slope'}, but is not atomic scalar.`
- Pass `target_parameter = "level"` (scalar string, not a vector) in one call
- Pass `target_parameter = "slope"` in a second call
- The `contdid` panel has 558 units, dose range [0, 226.6], 7812 obs — already confirmed valid
- If contdid still fails after the scalar fix, wrap in `tryCatch` and fall back gracefully to the fixest continuous results already computed
- Output: dose-response plot (ATT as function of dose D) saved to `Figures_Final/fig_contdid_level_v5.pdf`

### Fix 2 — HonestDiD sensitivity bounds (CRITICAL)
Current state: all HonestDiD outputs are NA. Steps to fix:
```r
library(HonestDiD)
# Extract event-study coefficients and vcov from aggte() output
# Use sensitivity_results <- HonestDiD::createSensitivityResults(
#   betahat = post_period_coefs,
#   sigma   = event_study_vcov,
#   numPrePeriods  = n_pre,
#   numPostPeriods = n_post,
#   Mbar = seq(0, 0.05, by = 0.01)
# )
# Run for: Low/Wages, Low/DR
# Save sensitivity plot to Figures_Final/fig_honestdid_low_v5.pdf
```
Do NOT run HonestDiD on High-Dose — singular matrix risk.

### Fix 3 — Enrich DR covariates (IMPORTANT)
Current DR spec uses 5 covariates. Add from Census 2010 data already in memory:
- `literacy_rate`, `share_in_school`, `avg_hh_income_pc_mw`
- Use `xformla = ~ pop_18_24 + avg_hh_income_pc_mw + literacy_rate + share_in_school` for Low-Dose DR only
- For High-Dose DR (currently failing due to singular matrix): reduce to `xformla = ~ pop_18_24 + avg_hh_income_pc_mw` (2 covariates max)
- Rerun `att_gt(..., est_method = "dr")` and `aggte()` for both
- Report new ATT and SE; goal is to push Low/DR from p=0.073 toward significance

### Fix 4 — Age-cohort triple interaction (NEW FINDING)
Data already downloaded: cohort RAIS with 1,107,684 rows split by age group (young 22–28, old 29–35).
- Build a new outcome variable: `log_wage_young_minus_log_wage_old` (gap) by microregion-year
- Run `att_gt()` separately for: young-only Low, young-only High, as already done
- NEW: add a **triple-interaction** subgroup — run `att_gt()` on young workers in Low-Dose microregions, conditional on `avg_hh_income_pc_mw` below median (poorest microregions only)
- This tests whether ProUni's wage effect on young workers concentrates in the most credit-constrained areas
- Save event-study figure to `Figures_Final/fig_triple_young_poor_v5.pdf`

### Fix 5 — Drop Wald LATE section
The Wald LATE estimate (LATE = 9.598, SE = 16.560) has a near-zero first stage (FS = 0.0024) and is not credible. Remove it from output tables. Add a comment `# Wald LATE dropped: weak first stage (F < 1)` in the code.

### Fix 6 — Activate placebo cohort test (ages 36–50)
Section in untitled4.R says `# PLACEHOLDER — uncomment to run`. Uncomment and execute it. This tests whether the effect exists for workers who were too old to benefit from ProUni. Expected result: ATT ≈ 0 (confirmation of identification).

### Output discipline
- All new figures: `Figures_Final/fig_*_v5.pdf`
- All new tables: `Tables_Final/tab_*_v5.tex`
- Append a `_v5` audit summary block at end of script
- Do NOT overwrite any `_v4` outputs

---

## Writing Conventions

- European academic economics tone — formal, precise, no overclaiming
- Causal language: "the estimates suggest", not "the programme caused"
- Demographic terminology: **ethnicity** (not "race"), **gender** (not "sex" unless specifically biological)
- LaTeX tables via `modelsummary` exported to `Tables/`; figures to `Figures/`
- BibTeX key format: `AuthorYYYYkeyword` (e.g. `Callaway2021staggered`)

---

## Key Commands

```r
# Run backup (never modify)
source("untitled4.R")

# Run improved script
source("untitled5.R")

# Compile thesis
# pdflatex Thesis_Main.tex && bibtex Thesis_Main && pdflatex Thesis_Main.tex
```

---

## What NOT To Do

- Do **not** modify `untitled4.R` — it is the frozen backup
- Do **not** re-download data from BigQuery in untitled5.R — reuse objects or RDS files
- Do **not** change the identification strategy, sample, or dose bins
- Do **not** hard-code BigQuery credentials — use `basedosdados::set_billing_id()`
- Do **not** commit raw RAIS microdata
- Do **not** interpret TWFE coefficients causally without noting the heterogeneous-treatment-effect bias
- Do **not** pass a vector to `target_parameter` in `contdid` — use a scalar string

---

## Progressive Disclosure

- `Documentation/` — data pipeline notes, variable definitions
- [Callaway & Sant'Anna did vignette](https://bcallaway11.github.io/did/articles/TWFE.html)
- [contdid package site](https://bcallaway11.github.io/contdid)
- [HonestDiD package](https://github.com/asheshrambachan/HonestDiD)
- [basedosdados R docs](https://basedosdados.org)
