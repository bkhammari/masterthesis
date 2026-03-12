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
untitled4.R            # Main R analysis script
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
  - `tidyverse`, `data.table`, `ggplot2`, `modelsummary`

---

## Econometric Framework

**Primary estimator:** Callaway & Sant'Anna (2021) via `did::att_gt()`
- Treatment: ProUni scholarship roll-out across Brazilian microregions (staggered, 2005+)
- Outcomes: log formal wages (RAIS), formal employment probability, employment count
- Heterogeneity: gender, ethnicity (branco / pardo / preto), macro-region (Norte, Nordeste, Sul, Sudeste, Centro-Oeste)

**Continuous treatment:** `contdid` for dose-response (scholarship intensity per microregion)

**Benchmark:** TWFE via `fixest::feols()` — use with caution, document weighting bias

**Data sources:** RAIS (employer–employee panel), ProUni admin records, IBGE Census 2010

---

## Critical Identification Assumptions

1. **Parallel trends** (conditional on covariates) — always present pre-trends tests
2. **No anticipation** — ProUni announced before implementation; discuss carefully
3. **Stable unit treatment value assumption (SUTVA)** — microregion spillovers are a known concern; discuss in Section 6

---

## Writing Conventions

- European academic economics tone — formal, precise, no overclaiming
- Causal language: "the estimates suggest", not "the programme caused"
- Demographic terminology: **ethnicity** (not "race"), **gender** (not "sex" unless specifically biological)
- LaTeX tables via `modelsummary` exported to `Tables/`; figures exported to `Figures/`
- BibTeX key format: `AuthorYYYYkeyword` (e.g. `Callaway2021staggered`)
- Equation numbering in LaTeX: use `\label{eq:att_gt}` convention

---

## Key Commands

```r
# Install / load core packages
pacman::p_load(did, contdid, fixest, basedosdados, geobr,
               tidyverse, data.table, ggplot2, modelsummary)

# Run main DiD
source("untitled4.R")

# Compile thesis (terminal)
# pdflatex Thesis_Main.tex && bibtex Thesis_Main && pdflatex Thesis_Main.tex
```

---

## What NOT To Do

- Do **not** rename `untitled4.R` without updating this file
- Do **not** hard-code BigQuery credentials — use `basedosdados::set_billing_id()`
- Do **not** commit raw RAIS microdata (confidential; use aggregated outputs only)
- Do **not** interpret TWFE coefficients causally without noting the heterogeneous-treatment-effect bias

---

## Progressive Disclosure

For deeper method guidance, see:
- `Documentation/` — data pipeline notes, variable definitions
- [Callaway & Sant'Anna did vignette](https://bcallaway11.github.io/did/articles/TWFE.html)
- [contdid package site](https://bcallaway11.github.io/contdid)
- [basedosdados R docs](https://basedosdados.org)
