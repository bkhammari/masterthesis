# Master Thesis — ProUni Labour Market Effects in Brazil

LaTeX + R econometrics thesis repo. Topic: causal effects of Brazil's ProUni
scholarship programme on formal labour-market outcomes, with heterogeneity
analysis by gender, ethnicity, and region.

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
Figures_Final/         # Output figures (PDF/PNG) — _v5 suffix for new outputs
Tables_Final/          # Output tables (TeX) — _v5 suffix for new outputs
Documentation/         # Method notes, data docs, CLAUDE.md
```

---

## Tech Stack

- **LaTeX** (Overleaf-compatible) — thesis writing
- **R** — econometric analysis
  - `did` — Callaway & Sant'Anna (2021) staggered DiD
  - `contdid` — continuous-treatment DiD (Callaway, Goodman-Bacon, Sant'Anna)
  - `fixest` — TWFE benchmarks + continuous dose specifications
  - `basedosdados` — BigQuery access to RAIS & ProUni admin data
  - `geobr` — Brazilian geographic shapefiles
  - `tidyverse`, `data.table`, `ggplot2`, `modelsummary`, `HonestDiD`

---

## Data Sources

### ProUni Administrative Records
- **Table:** `basedosdados.br_mec_prouni.cursos`
- **Key vars:** `ano`, `id_municipio`, `bolsas_integrais`, `bolsas_parciais`, `id_ies`
- Aggregated to microregion via municipality directorios crosswalk

### RAIS (Relação Anual de Informações Sociais)
- **Table:** `basedosdados.br_me_rais.microdados_vinculos`
- **Panel:** 2005–2019, downloaded year-by-year
- **Key slices:**
  - Main panel (wages + employment, aggregated to microregion × year)
  - Age-cohort split young 22–28 vs old 29–35: 1,107,684 rows
  - Education level split: 628,124 rows
  - Birth-cohort split (born 1987–93 vs 1975–82): 1,066,743 rows
- **Key RAIS codes:**
  - `grau_instrucao = 7` → Superior Completo (university degree)
  - `raca_cor`: 2=Branca, 4=Preta, 8=Parda, 6=Amarela, 1=Indígena, 9=Ignorado
  - Non-white defined as raca_cor ∈ {1, 4, 8}
  - `sexo`: 1=Masculino, 3=Feminino

### IBGE Censo Demográfico 2010
- Youth population denominator (`pop_18_24`), literacy rate, household income,
  informality share, school enrolment — all aggregated to microregion

### Directorios: Município (`br_bd_diretorios_brasil.municipio`)
- 5,571 municipalities, 27 columns, 100% coverage
- **Primary join key:** `id_municipio` (7-digit IBGE) → links RAIS, Censo, ProUni
- **Aggregation key:** `id_municipio` → `id_microrregiao` (558 microregions)
- **Heterogeneity dimension:** `nome_regiao` (5 grandes regiões)
- **Useful flags:** `amazonia_legal`, `capital_uf`, `id_regiao_metropolitana`
- **Spatial data:** `centroide` as WGS84 POINT for distance-based instruments
- Cached as: `bquxjob_7c0fafc2_19cdfc494e9.csv`

### Directorios: IES (`br_bd_diretorios_brasil.instituicao_ensino_superior`)
- 6,589 IES (unique), 7 columns, 100% coverage
- **ProUni-eligible universe:** 2,859 active private IES in 723 municipalities
- **Geographic concentration:** SP(647) > MG(390) > PR(198) > BA(155) > SC(133)
- **Warning:** 3,378/6,589 IES are Inactive (51%) — entry/exit may be endogenous
  to local labour markets; flag in robustness section
- **Validation use:** confirm ProUni scholarship counts = 0 for Federal/Estadual IES
- Cached as: `bquxjob_505d8763_19cdfc6c2d3.csv`

---

## Econometric Framework

### Treatment Variable
```
D_rt = (Σ_{m∈r} ProUni_mt) / pop_18_24_r × 1000
```
- Population denominator: fixed at Censo 2010 (time-invariant)
- `g` = first year D_rt > 0 per microregion
- **Dose thresholds (2019 cross-section):** τ_20 = 20.16 | τ_50 = 48.42 | τ_75 = 78.87
- **Bins:** Low = D_rt_2019 ∈ (0, τ_50] → 167 microregions;
  High = D_rt_2019 > τ_50 → 140 microregions;
  Never-treated (comparison) → 112 microregions
- **Treatment year range:** g ∈ [2010, 2019] | Dose range: [0.0, 226.6]
- All sanity checks passed: no duplicates, D_rt = 0 pre-2009, weakly increasing ✓

### Primary Estimator: Callaway & Sant'Anna (2021)
- `did::att_gt()` with `control_group = "notyettreated"`
- Aggregated via `aggte(type = "simple")`
- **Primary specification:** doubly-robust (DR) with covariates:
  `xformla = ~ pop_18_24 + avg_hh_income_pc_mw + literacy_rate + share_in_school`
- Balanced panel confirmed; `allow_unbalanced_panel = FALSE`

### Continuous Treatment
- `contdid` for dose-response (currently broken — see Fix 1 below)
- Fallback: `fixest::feols()` with D_rt × post interactions (linear, log, quadratic)

### Benchmark
- `fixest::feols()` TWFE — documented sign-reversal bias vs CS estimates
- TWFE wages: −0.65% vs CS Low/DR: +2.6% → do not interpret causally

---

## Key Empirical Results

### Headline ATTs (Callaway & Sant'Anna)

| Model | ATT | SE | Pre-trend p |
|---|---|---|---|
| **Low/DR (primary)** | **+0.0259** | **0.0145** | **0.535** |
| Low/Wages (uncond.) | +0.0128 | 0.0128 | 0.359 |
| High/Wages (uncond.) | −0.0070 | 0.0074 | 0.916 |
| Low/Emp | +0.0182 | 0.0332 | 0.120 |
| High/DR | FAILED (singular matrix) | — | — |

**Central finding:** Consistent asymmetry — all Low-dose estimates ≥ 0,
all High-dose estimates ≤ 0, across all 45 models. Consistent with scarcity-rent
hypothesis in low-saturation markets and crowding-out/saturation in high-dose regions.

### Gender Heterogeneity
| Model | ATT | SE |
|---|---|---|
| Male/Low | +0.0167 | 0.0169 |
| Female/Low | −0.0003 | 0.0007 |
| Male/High | −0.0083 | 0.0092 |

Wage effects concentrated among male workers in low-dose regions. Female effect
is essentially zero. Persistent gender gap post-ProUni.

### Ethnicity Heterogeneity
| Model | ATT | SE |
|---|---|---|
| Non-white/Low | +0.0128 | 0.0119 |
| White/Low | −0.0022 | 0.0088 |
| Non-white/High | −0.0070 | 0.0071 |

Non-white workers in low-dose regions benefit (+1.3%); white workers do not.

### Continuous Dose (fixest, N=7,812)
| Spec | Coef | p |
|---|---|---|
| Linear D×post | −0.000129 | 0.054 |
| Log(1+D)×post | −0.002124 | 0.048 |
| Quadratic D×post | −0.000256 | 0.033 |
| Linear+DR D×post | +0.000063 | 0.434 |

Inverted-U dose-response consistent with binned CS results.

### TWFE Benchmark (methodology illustration)
| Model | Coef | SE |
|---|---|---|
| TWFE/Wages | −0.0065 | 0.0037 |
| TWFE/Employment | −0.0467** | 0.0180 |

TWFE employment estimate has **wrong sign** vs CS (+1.8%). Use prominently in
methodology section to motivate CS over TWFE.

### Abstract Numbers
1. Primary wage effect (Low/DR): **+2.6%** (p = 0.073)
2. Continuous dose (log): **−0.0021** per log-unit (p = 0.048, N=7,812)
3. TWFE employment bias: **−4.7%\*** (wrong sign vs CS +1.8%)
4. Pre-trend p-values: all > 0.08 across all 45 models ✓

---

## Critical Identification Assumptions

1. **Parallel trends** (conditional on covariates) — all pre-trend p-values > 0.08
   across all 45 models ✓
2. **No anticipation** — ProUni announced before implementation; discuss in Section 4
3. **SUTVA** — buffer-zone exclusion (13 border microregions removed) passes:
   Low/Buffer ≈ Low/Wages ✓

---

## untitled5.R — Task Specification

`untitled5.R` is a **surgical improvement** of `untitled4.R`. Rules:
- Do NOT re-download from BigQuery — reuse objects or saved RDS files
- Do NOT change identification strategy, sample, or dose bins
- All outputs use `_v5` suffix; do NOT overwrite `_v4` outputs

### Fix 1 — contdid [CRITICAL] 🔴
**Error:** `Assertion on 'target_parameter' failed: Must be element of set
{'level','slope'}, but is not atomic scalar.`
```r
# CORRECT usage — scalar string, NOT a vector:
out_level <- contdid(..., target_parameter = "level")
out_slope <- contdid(..., target_parameter = "slope")
# Wrap in tryCatch; fall back to fixest results if still failing
```
- Panel: 558 units, dose [0, 226.6], 7,812 obs — confirmed valid
- Output: `Figures_Final/fig_contdid_level_v5.pdf`

### Fix 2 — HonestDiD sensitivity bounds [CRITICAL] 🔴
**Current state:** All HonestDiD outputs are NA.
```r
library(HonestDiD)
sensitivity_results <- HonestDiD::createSensitivityResults(
  betahat       = post_period_coefs,   # from aggte() event-study
  sigma         = event_study_vcov,
  numPrePeriods  = n_pre,
  numPostPeriods = n_post,
  Mbar = seq(0, 0.05, by = 0.01)
)
# Run for: Low/Wages and Low/DR only
# Do NOT run on High-Dose (singular matrix risk)
```
- Output: `Figures_Final/fig_honestdid_low_v5.pdf`

### Fix 3 — Enrich DR covariates [IMPORTANT] 🟡
- Current DR spec: 5 covariates
- Enriched spec:
  `xformla = ~ pop_18_24 + avg_hh_income_pc_mw + literacy_rate + share_in_school`
- High-Dose DR (singular matrix): reduce to `xformla = ~ pop_18_24 + avg_hh_income_pc_mw`
- Goal: push Low/DR from p = 0.073 toward conventional significance
- Rerun `att_gt(..., est_method = "dr")` + `aggte()` for both bins

### Fix 4 — Age-cohort triple interaction [NEW] 🟡
```r
# NEW: young workers (22–28) × Low-Dose × below-median income microregions only
# Tests: does ProUni wage effect concentrate in most credit-constrained areas?
young_poor_data <- data_low |>
  filter(avg_hh_income_pc_mw < median(avg_hh_income_pc_mw))
# Run att_gt() on this subsample
```
- Output: `Figures_Final/fig_triple_young_poor_v5.pdf`

### Fix 5 — Drop Wald LATE [COSMETIC] 🟢
```r
# Wald LATE dropped: weak first stage (F < 1)
# LATE = 9.598, SE = 16.560, FS = 0.0024 — not credible
```
Remove from all output tables. Mention in footnote only if reviewer asks.

### Fix 6 — Activate placebo cohort test ages 36–50 [IMPORTANT] 🟡
- Section in `untitled4.R` is `# PLACEHOLDER — uncomment to run`
- Uncomment and execute
- Expected: ATT ≈ 0 (workers too old to benefit from ProUni)
- Positive result would confirm identification; non-zero would indicate
  general equilibrium spillovers to older workers

---

## Known Issues

| Issue | Severity | Status |
|---|---|---|
| contdid `target_parameter` assertion | 🔴 Critical | Fix in untitled5.R |
| High/DR singular matrix | 🔴 Critical | Reduce covariates in untitled5.R |
| HonestDiD all NA | 🔴 Critical | Fix in untitled5.R |
| Central/Wages skipped (0 controls) | 🟡 Important | Remove from paper |
| Placebo cohort 36–50 not run | 🟡 Important | Activate in untitled5.R |
| IES entry/exit endogeneity | 🟡 Important | Add footnote in Data section |
| Figure 2 (regional divergence) skipped | 🟢 Minor | Remove from paper outline |
| `inf`→`NA` coercion on `g` column | 🟢 Minor | Cosmetic, no effect on results |

---

## Writing Conventions

- European academic economics tone — formal, precise, no overclaiming
- Causal language: "the estimates suggest", not "the programme caused"
- Demographic terminology: **ethnicity** (not "race"), **gender** (not "sex"
  unless specifically biological)
- LaTeX tables via `modelsummary` exported to `Tables_Final/`; figures to `Figures_Final/`
- BibTeX key format: `AuthorYYYYkeyword` (e.g. `Callaway2021staggered`)

---

## Key Commands

```r
# Run frozen backup (never modify)
source("untitled4.R")

# Run improved script
source("untitled5.R")

# Set BigQuery credentials (never hard-code)
basedosdados::set_billing_id("YOUR_PROJECT_ID")
```

```bash
# Compile thesis
pdflatex Thesis_Main.tex && bibtex Thesis_Main && pdflatex Thesis_Main.tex
```

---

## What NOT To Do

- Do **not** modify `untitled4.R` — frozen backup
- Do **not** re-download from BigQuery in `untitled5.R` — reuse objects or RDS
- Do **not** change identification strategy, sample, or dose bins
- Do **not** hard-code BigQuery credentials — use `basedosdados::set_billing_id()`
- Do **not** commit raw RAIS microdata
- Do **not** interpret TWFE coefficients causally without noting
  heterogeneous-treatment-effect bias
- Do **not** pass a vector to `target_parameter` in `contdid` — scalar string only
- Do **not** report the Wald LATE (weak first stage, F < 1)

---

## Progressive Disclosure

- `CLAUDE.md` — full empirical context, all 45 ATT results, data pipeline, known issues
- [Callaway & Sant'Anna did vignette](https://bcallaway11.github.io/did/articles/TWFE.html)
- [contdid package site](https://bcallaway11.github.io/contdid)
- [HonestDiD package](https://github.com/asheshrambachan/HonestDiD)
- [basedosdados R docs](https://basedosdados.org)
