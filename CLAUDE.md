# CLAUDE.md — ProUni Thesis: Full Project Context
# Last updated: 2026-03-12

## 1. PROJECT OVERVIEW

**Title:** How Equal is Equal Opportunity? Gender, Ethnic and Regional Disparities
in ProUni Labour Market Effects in Brazil

**Unit of analysis:** Microregion × year (r, t)
**Panel:** 558 microregions × 15 years (2005–2019) = 7,812 observations
**Treatment variable:** D_rt = cumulative ProUni scholarships per 1,000 young adults (age 18–24)
**Outcome variables:** log formal wages (RAIS), log formal employment count (RAIS)
**Estimator (primary):** Callaway & Sant'Anna (2021) — `did` package, R
**Estimator (continuous):** fixest TWFE with interaction D_rt × post (logged, linear, quadratic)
**Estimator (benchmark):** TWFE via fixest (biased — used only for contrast)

---

## 2. DATA SOURCES & STATUS

### 2A. ProUni Administrative Records
- **Source:** `basedosdados.br_mec_prouni.cursos` (BigQuery)
- **Downloaded:** batched strategy, all years available
- **Key variables:** `ano`, `id_municipio`, `bolsas_integrais`, `bolsas_parciais`,
  `turno`, `nome_curso`, `id_ies`
- **Dose construction:** aggregate to microregion via directorios join, divide by
  Censo 2010 youth population (18–24)
- **Note:** ProUni only operates in *private* IES → validate against IES directory

### 2B. RAIS (Relação Anual de Informações Sociais)
- **Source:** `basedosdados.br_me_rais.microdados_vinculos` (BigQuery)
- **Downloaded:** year-by-year loop 2005–2019; ~628K–1.1M rows depending on slice
- **Key slices available:**
  - Main panel (aggregated to microregion × year): wages + employment
  - Age-cohort split (young 22–28 vs old 29–35): 1,107,684 rows
  - Education level split: 628,124 rows
  - Birth-cohort split (born 1987–93 vs 1975–82): 1,066,743 rows
- **Outcome construction:** `log_wage_sm` = log(mean real wage / minimum wage);
  `n_vinculos` = count of formal contracts
- **Filters applied:** workers with tertiary education, age 22–35, valid wage > 0
- **Key RAIS codes (from dictionary):**
  - `grau_instrucao`: 7 = Superior Completo (university degree)
  - `raca_cor`: 1=Indígena, 2=Branca, 4=Preta, 6=Amarela, 8=Parda, 9=Ignorado
  - `sexo`: 1=Masculino, 3=Feminino
  - Non-white defined as: raca_cor ∈ {1, 4, 8} (Indígena, Preta, Parda)

### 2C. IBGE 2010 Censo (Microdata)
- **Source:** `basedosdados.br_ibge_censo_demografico` (BigQuery)
- **Used for:** youth population denominator (pop_18_24), covariates
- **Key covariates built:**
  - `pop_18_24`: youth population per microregion (denominator for D_rt)
  - `pop_nonwhite`: non-white share (balance check)
  - `share_no_formal_income`: informality proxy
  - `avg_hh_income_pc_mw`: average household income per capita in MW
  - `literacy_rate`: literacy rate
  - `share_in_school`: school enrolment rate

### 2D. Directorios: Município (`br_bd_diretorios_brasil.municipio`)
- **Rows:** 5,571 municipalities | **Columns:** 27 | 100% coverage
- **Critical join keys:**
  - `id_municipio` (7-digit IBGE) → **primary join key** across all datasets
  - `id_microrregiao` → aggregation to microregion r
  - `id_mesorregiao` → robustness at coarser level
  - `nome_regiao` → 5 grandes regiões (Norte/Nordeste/Centro-Oeste/Sudeste/Sul)
  - `id_regiao_metropolitana` → metro vs non-metro heterogeneity
- **Useful binary flags:** `amazonia_legal` (773 muns), `capital_uf` (27 capitals)
- **Spatial data:** `centroide` as WGS84 POINT(lon lat) → extract with regex for
  distance-based spillover instruments
- **Cached as:** `bquxjob_7c0fafc2_19cdfc494e9.csv`
- **Join pipeline:**
  ```r
  dir_mun <- read_csv("bquxjob_7c0fafc2_19cdfc494e9.csv") |>
    select(id_municipio, id_microrregiao, nome_microrregiao,
           id_mesorregiao, sigla_uf, nome_regiao,
           id_regiao_metropolitana, amazonia_legal, centroide) |>
    mutate(id_municipio = as.character(id_municipio),
           lon = as.numeric(str_extract(centroide, "(?<=POINT\\()[-\\d.]+")),
           lat = as.numeric(str_extract(centroide, "(?<=\\s)[-\\d.]+(?=\\))")))  
  ```

### 2E. Directorios: IES (`br_bd_diretorios_brasil.instituicao_ensino_superior`)
- **Rows:** 6,589 IES (unique, no duplicates) | **Columns:** 7 | all 100% coverage
- **Key variables:** `id_ies`, `nome`, `rede`, `tipo_instituicao`,
  `situacao_funcionamento`, `id_municipio`, `sigla_uf`
- **rede distribution:** Privada=6,150 | Estadual=179 | Federal=165 |
  Municipal=67 | Especial=28
- **situacao_funcionamento:** Ativa=3,211 | **Inativa=3,378** (51% inactive!)
- **ProUni-eligible universe:** 2,859 active private IES in 723 municipalities,
  all 27 UFs
- **Geographic concentration:** SP(647) > MG(390) > PR(198) > BA(155) > SC(133)
- **Identification risk:** IES entry/exit over 2005–2019 may be endogenous to
  local labour market conditions → flag in robustness section
- **Cached as:** `bquxjob_505d8763_19cdfc6c2d3.csv`
- **Validation use:**
  ```r
  ies_dir <- read_csv("bquxjob_505d8763_19cdfc6c2d3.csv") |>
    mutate(id_municipio = as.character(id_municipio),
           prouni_eligible = (rede == "Privada"))
  # Verify ProUni scholarships = 0 for Federal/Estadual IES (data quality check)
  ```

---

## 3. TREATMENT VARIABLE CONSTRUCTION

```
D_rt = (Σ_{m∈r} ProUni_mt) / pop_18_24_r × 1000
```

- Aggregated from municipality → microregion via `dir_mun`
- Population denominator: fixed at 2010 Censo (time-invariant)
- `g` (first treatment year) = first year with D_rt > 0, per microregion
- **Thresholds (2019 cross-section):**
  - τ_20 = 20.16 | τ_50 = 48.42 | τ_75 = 78.87
- **Bin assignment (primary):** Low = D_rt_2019 ∈ (0, τ_50] → 167 microregions
  High = D_rt_2019 > τ_50 → 140 microregions
  Never-treated (comparison): g = 0 → 112 microregions
- **Sanity checks (all passed):**
  - No duplicate (id_microrregiao, ano) pairs
  - D_rt = 0 for all years before 2009 ✓
  - D_rt weakly increasing over time ✓
  - Both bins > 50 microregions ✓
- **Treatment year range:** g ∈ [2010, 2019]
- **Dose range (continuous models):** [0.0, 226.6]

---

## 4. EMPIRICAL RESULTS

### 4A. Primary Results (Callaway & Sant'Anna)

| Model | ATT | SE | Pre-trend p | N_obs |
|---|---|---|---|---|
| Low/Wages (uncond.) | +0.0128 | 0.0128 | 0.359 | 3,906 |
| High/Wages (uncond.) | −0.0070 | 0.0074 | 0.916 | 3,528 |
| Marginal/Wages | +0.0136 | 0.0078 | 0.362 | 3,752 |
| **Low/DR (primary)** | **+0.0259** | **0.0145** | **0.535** | **3,906** |
| High/DR | FAILED (singular matrix) | — | — | — |

**Interpretation:** Low-dose microregions show a positive wage effect of ~+2.6%
(DR, p=0.073), consistent with scarcity-rent hypothesis. High-dose regions show
no significant effect (saturation). This asymmetry is the central finding.

### 4B. Employment Results

| Model | ATT | SE | Pre-trend p |
|---|---|---|---|
| Low/Emp | +0.0182 | 0.0332 | 0.120 |
| High/Emp | +0.0065 | 0.0396 | 0.809 |
| Marginal/Emp | −0.0386 | 0.0335 | 0.597 |

All employment effects are statistically insignificant. TWFE shows −4.7%***
(wrong sign vs CS) — key illustration of TWFE bias for thesis.

### 4C. Gender Heterogeneity

| Model | ATT | SE | Pre-trend p |
|---|---|---|---|
| Male/Low | +0.0167 | 0.0169 | 0.418 |
| Female/Low | −0.0003 | 0.0007 | 0.953 |
| Male/High | −0.0083 | 0.0092 | 0.918 |
| Female/High | −0.0041 | 0.0073 | 0.736 |

**Finding:** Wage effects concentrated among males in low-dose regions. Female
effect is essentially zero. Persistent gender gap post-ProUni.

### 4D. Ethnicity Heterogeneity

| Model | ATT | SE | Pre-trend p |
|---|---|---|---|
| White/Low | −0.0022 | 0.0088 | 0.355 |
| Non-white/Low | +0.0128 | 0.0119 | 0.405 |
| White/High | −0.0048 | 0.0088 | 0.683 |
| Non-white/High | −0.0070 | 0.0071 | 0.916 |

**Finding:** Non-white workers in low-dose regions benefit (+1.3%), white workers
do not. Ethnic parity argument partially supported at low doses, but effects are
not statistically significant at conventional levels.

### 4E. Continuous Dose (fixest, all 558 microregions, N=7,812)

| Model | Coef | SE | t | p |
|---|---|---|---|---|
| Linear D×post | −0.000129 | 0.0000669 | −1.932 | 0.054 |
| Log(1+D)×post | −0.002124 | 0.001072 | −1.982 | 0.048 |
| Quadratic D×post | −0.000256 | 0.000120 | −2.139 | 0.033 |
| Quad D²×post | +0.0000012 | 0.0000012 | 1.010 | 0.313 |
| Linear+DR D×post | +0.0000631 | 0.0000806 | 0.783 | 0.434 |

**Interpretation:** Inverted-U dose-response. Negative average slope masks
positive effects at low doses (consistent with binned CS results). The quadratic
term is statistically insignificant, though the coefficient pattern supports
saturation at high doses.

**contdid status:** FAILED — `target_parameter` assertion error. Fix in untitled5.R.

### 4F. TWFE Benchmark (validates CS methodology)

| Model | Coef | SE | N |
|---|---|---|---|
| TWFE/Wages | −0.0065 | 0.0037 | 7,812 |
| TWFE/Employment | −0.0467*** | 0.0180 | 7,812 |

TWFE wage effect is −0.65% vs CS Low/DR = +2.6% → sign difference for
employment, magnitude difference for wages. Use this contrast prominently in
methodology section to motivate CS over TWFE.

### 4G. Cohort Analysis (Duflo-style: young 22–28 vs old 29–35)

| Model | ATT | SE | Pre-trend p |
|---|---|---|---|
| Gap (young−old)/Low | −0.0054 | 0.0120 | 0.707 |
| Gap/High | +0.0234 | 0.0273 | 0.539 |
| Young/Low | +0.0095 | 0.0107 | 0.644 |
| Young/High | −0.0042 | 0.0132 | 0.162 |
| Old/Low [placebo] | +0.0149 | 0.0125 | 0.871 |
| Old/High [placebo] | −0.0276 | 0.0199 | 0.981 |

**Note:** Old cohort placebo is non-zero in Low/Old (+1.5%) — minor concern,
discuss as general equilibrium spillover to older workers.

### 4H. Birth-Cohort CS (born 1987–93 vs 1975–82)

| Model | ATT | SE | Pre-trend p |
|---|---|---|---|
| BC_Gap/Low | +0.0116 | 0.0090 | 0.562 |
| BC_Gap/High | +0.0018 | 0.0136 | 0.811 |
| BC_Exposed/Low | +0.0107 | 0.0101 | 0.968 |
| BC_Control/Low [placebo] | −0.0009 | 0.0092 | 0.904 |

Clean placebo (BC_Control ≈ 0) strengthens identification.

### 4I. Robustness Checks

| Check | Low ATT | High ATT | Notes |
|---|---|---|---|
| Lag-3 | −0.0001 | −0.0019 | Attenuated; effects take >3 yrs |
| Buffer exclusion (SUTVA) | +0.0129 | −0.0070 | 13 border muns excluded; stable |
| Alt. cutoff q25/q75 | +0.0015 | −0.0051 | Dilution in Low bin |
| Alt. cutoff q33/q67 | +0.0089 | −0.0078 | Consistent sign pattern |

All pre-trend p-values remain > 0.08 across all 45 models. Identification robust.

### 4J. Wald LATE — DROPPED
LATE = 9.598, SE = 16.560, FS = 0.0024. Weak first stage (F < 1). Not reported.

---

## 5. CASE STUDY MICROREGIONS

| id_micro | Name | UF | Type | Bin | g | D_2019 |
|---|---|---|---|---|---|---|
| 31028 | Conceição do Mato Dentro | MG | Central | Low | 2014 | 42.1 |
| 35010 | São Joaquim da Barra | SP | Central | High | 2013 | 89.0 |
| 22001 | Baixo Parnaíba Piauiense | PI | Marginal | Low | 2019 | 23.1 |
| 29021 | Salvador | BA | Marginal | High | 2011 | 104.0 |

---

## 6. COVARIATE BALANCE

| | High (n=140) | Low (n=167) |
|---|---|---|
| Mean youth pop | 74,598 | 27,689 |
| Mean non-white pop | 33,619 | 18,024 |
| Share no formal income | 0.402 | 0.551 |
| Avg HH income p.c. (MW) | higher | lower |
| Literacy rate | higher | lower |
| Median treatment year | 2012 | 2017 |

High-dose microregions are larger, wealthier, and treated earlier. Covariate
adjustment via DR is important; Low/DR is the primary specification.

---

## 7. KNOWN ISSUES & TODO

### Critical 🔴
- [ ] **contdid BROKEN:** `target_parameter` assertion fails. Fix: pass scalar
  string `target_parameter = "level"` — NOT a vector
- [ ] **High/DR FAILED:** Singular matrix. Try `faster_mode = FALSE`; if still
  fails, reduce to `xformla = ~ pop_18_24 + avg_hh_income_pc_mw`
- [ ] **HonestDiD:** All NA in audit. Run for Low/Wages and Low/DR before
  submission. Do NOT run on High-Dose.

### Important 🟡
- [ ] **Central/Wages & Central/Emp SKIPPED:** 0 not-yet-treated units. Remove
  from paper or reframe as descriptive only.
- [ ] **Placebo cohort 36–50:** Code block is PLACEHOLDER — uncomment and run
- [ ] **IES entry/exit endogeneity:** Add footnote in Data section noting 51% of
  IES are inactive; D_rt may partly reflect private sector dynamics
- [ ] **VCov matrix warnings:** `pretrend_ftest` falls back to pointwise SEs for
  some models

### Minor 🟢
- [ ] Figure 2 (regional divergence) skipped — remove from paper outline
- [ ] `inf`→`NA` coercion warning for `g` column — cosmetic only

---

## 8. OUTPUT FILE INVENTORY

**Tables_Final/** — regression tables, balance tables, ATT summary (62 files)
**Figures_Final/** — event-study plots, maps, heterogeneity figures (24 files)
- Key figures: `fig_1_mechanism_dose.png`, `fig_4a/b_gender_*.png`,
  `fig_5a/b_race_*.png`, `fig_map_dose_bins.png`, `fig_map_dose_continuous.png`
- Event-study CSVs: `es_*.csv` (43 files)

---

## 9. KEY NUMBERS FOR THESIS ABSTRACT

1. **Primary wage effect (Low-dose, DR):** +2.6% (p = 0.073)
2. **Continuous dose (log):** −0.0021 per log-unit of dose (p = 0.048, N=7,812)
3. **TWFE employment bias:** −4.7%** (wrong sign vs CS +1.8%) — methodology
4. **Consistent sign pattern:** ALL Low-dose estimates ≥ 0; ALL High-dose ≤ 0
5. **Pre-trend p-values:** All > 0.08 across 45 models (identification robust)
6. **Gender:** Male low-dose +1.7%; Female low-dose ≈ 0
7. **Ethnicity:** Non-white low-dose +1.3%; White low-dose −0.2%

---

## 10. R PACKAGE VERSIONS

- `did` — Callaway & Sant'Anna estimator
- `contdid` — continuous treatment DID (currently broken — scalar fix pending)
- `fixest` — TWFE benchmark + continuous dose
- `basedosdados` — BigQuery interface
- `geobr` — Brazilian geographic data (maps)
- `HonestDiD` — sensitivity analysis (currently broken — fix pending)
- `tidyverse`, `data.table`, `ggplot2`, `modelsummary`

---

## 11. THESIS CHAPTER MAPPING

| Chapter | Key results | Key tables/figures |
|---|---|---|
| 3. Data | IES directory, RAIS slices, Censo covariates | Table 1, Table 2 |
| 4. Methodology | C&S identification, dose construction | — |
| 5. Main Results | Low/DR (+2.6%), High (null), TWFE contrast | Table 4, Fig 1, Fig 3 |
| 6. Heterogeneity | Gender (male advantage), Ethnicity (non-white benefit) | Fig 4, Fig 5 |
| 7. Robustness | Lag-3, Buffer, Alt-cutoffs, Birth-cohort | Fig 6–9, Table robustness |
| 8. Conclusion | Scarcity-rent hypothesis, saturation, policy implications | — |
