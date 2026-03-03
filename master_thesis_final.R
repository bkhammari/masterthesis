# ==============================================================================
# MASTER THESIS: THE DISTRIBUTIONAL CONSEQUENCES OF HIGHER EDUCATION EXPANSION
# ==============================================================================
# Author:      Baha Eddine Khammari
# Date:        2026-03
# Methodology: Callaway & Sant'Anna (2021) Staggered DiD with continuous
#              treatment intensity (Callaway, Goodman-Bacon & Sant'Anna 2024).
#
# IDENTIFICATION STRATEGY (Chapter 4):
#   - No genuine never-treated group: all 558 microregions eventually receive
#     some ProUni exposure over 2005-2019.
#   - Comparison group is exclusively NOT-YET-TREATED units at each time t:
#     C_t = {r : G_r^d > t}, enforced via control_group = "notyettreated".
#   - FIX: Dose bins are defined from the INITIAL TREATMENT-YEAR distribution
#     of D_{r,g} (dose at first crossing), NOT the 2019 cross-section.
#     Low = (tau_20, tau_50], High = D_at_entry > tau_75.
#   - Treatment timing g = first year D_{r,t} > tau_20.
#
# OUTPUTS:
#   - 18 baseline CS models (wages + employment x dose/region/gender/race)
#   - 2 lag-robustness models (alternative lag = 3 years)
#   - NEW: Conditional (doubly-robust) models with Census covariates
#   - NEW: HonestDiD sensitivity analysis on primary wage models
#   - NEW: Buffer-zone exclusion test (SUTVA robustness)
#   - NEW: Placebo cohort test placeholder (ages 36-50)
#   - NEW: TWFE benchmark for comparison with CS estimator
#   - NEW: Alternative cutoff robustness (tau_25/75, tau_33/67)
#   - Publication figures (Figures_Final/)
#   - ATT summary table, event-study CSVs, balance table (Tables_Final/)
#   - Audit log with pre-trend tests (Documentation_Final/audit_log.rds)
# ==============================================================================

rm(list = ls())
gc()

# ==============================================================================
# SECTION 0: SETUP — PARAMETERS, LIBRARIES, DIRECTORIES, THEME
# ==============================================================================

# Central parameter store — every tuneable constant lives here
params <- list(
  age_min        = 22L,
  age_max        = 35L,
  lag_years      = 4L,
  q_low          = 0.20,
  q_mid          = 0.50,
  q_high         = 0.75,
  set_seed       = 7364599L,
  years          = 2005L:2019L,
  central_states = c("SP", "RJ", "MG", "ES", "RS", "SC", "PR")
)

set.seed(params$set_seed)

library(basedosdados)
library(dplyr)
library(tidyr)
library(ggplot2)
library(did)        # Callaway & Sant'Anna (2021) estimator
library(readr)
library(xtable)
library(scales)
library(broom)

# NEW: Additional dependencies for robustness and sensitivity analyses
# Wrapped in tryCatch so the script does not fail if a package is missing.
for (pkg in c("fixest", "sf", "geobr")) {
  tryCatch(library(pkg, character.only = TRUE),
           error = function(e) message(sprintf("  WARNING: package '%s' not available — related sections will be skipped.", pkg)))
}
# NEW: HonestDiD — install from GitHub if not available
tryCatch(library(HonestDiD),
         error = function(e) {
           message("  WARNING: HonestDiD not installed. Attempting GitHub install...")
           tryCatch({
             if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
             devtools::install_github("asheshrambachan/HonestDiD")
             library(HonestDiD)
           }, error = function(e2) message("  WARNING: HonestDiD install failed — sensitivity section will be skipped."))
         })

basedosdados::set_billing_id("microdadosbrazil")

dirs <- c("Figures_Final", "Tables_Final", "Documentation_Final")
sapply(dirs, function(x) if (!dir.exists(x)) dir.create(x))

theme_set(theme_classic() + theme(
  plot.title         = element_text(face = "bold", size = 14, color = "black"),
  plot.subtitle      = element_text(size = 11, color = "grey40"),
  axis.text          = element_text(size = 10, color = "black"),
  legend.position    = "bottom",
  panel.grid.major.y = element_line(color = "grey95")
))

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS
# ==============================================================================

# Helper 0 — escape curly braces in error messages for cli compatibility
# The cli package intercepts message() and parses {}-expressions. BigQuery
# errors often contain type info like {INT64, STRING} which crashes cli.
safe_err_msg <- function(e) gsub("\\{", "{{", gsub("\\}", "}}", conditionMessage(e)))

# Helper 1 — join nome_microrregiao into any data frame with id_microrregiao
# Defensive: skips join if columns already exist.
add_micro_names <- function(df, geo) {
  cols_needed <- c("nome_microrregiao", "sigla_uf")
  cols_missing <- setdiff(cols_needed, names(df))
  if (length(cols_missing) == 0) return(df)
  lkp <- geo %>%
    distinct(id_microrregiao, .keep_all = TRUE) %>%
    select(id_microrregiao, all_of(cols_missing))
  df %>% left_join(lkp, by = "id_microrregiao")
}

# Helper 2 — Wald chi-squared pre-trend joint test on ATT_ES(e < 0)
# Receives an att_gt object, computes the event study, then tests the joint
# null that all pre-treatment dynamic effects are zero.
pretrend_ftest <- function(att_gt_obj) {
  es       <- did::aggte(att_gt_obj, type = "dynamic", na.rm = TRUE)
  pre_idx  <- which(es$egt < 0)
  if (length(pre_idx) == 0) {
    warning("No pre-treatment periods available for Wald test.")
    return(list(statistic = NA_real_, df = 0L, p.value = NA_real_))
  }
  ests <- es$att.egt[pre_idx]

  # Defensively find the VCov matrix — field name may vary across did versions
  vcv_candidates <- intersect(
    c("V_analytical.aggte", "V.analytical.aggte", "V_analytical"),
    names(es)
  )
  if (length(vcv_candidates) > 0) {
    vcv_full <- es[[vcv_candidates[1]]]
    vcv <- vcv_full[pre_idx, pre_idx, drop = FALSE]
  } else {
    vcv <- NULL
  }

  if (is.null(vcv) || any(is.na(vcv))) {
    warning("VCov matrix not available; falling back to pointwise SEs.")
    ses <- es$se.egt[pre_idx]
    W   <- sum((ests / ses)^2, na.rm = TRUE)
    df  <- sum(!is.na(ests / ses))
  } else {
    W  <- tryCatch(
      as.numeric(t(ests) %*% solve(vcv) %*% ests),
      error = function(e) { warning("VCov inversion failed."); NA_real_ }
    )
    df <- length(ests)
  }
  list(statistic = W, df = df, p.value = pchisq(W, df = df, lower.tail = FALSE))
}

# Helper 3 — run Callaway & Sant'Anna estimation
# ALWAYS uses control_group = "notyettreated" (Chapter 4 identification).
# Returns: att_gt object, dynamic event study, simple ATT aggregate,
#          tidy df, pre-trend test, label, sample sizes.
#
# agg controls how municipality-level cells are aggregated to microregions:
#   "wmean"  — weighted mean using n_vinculos (for log wages, hourly wages, etc.)
#   "logsum" — log of total count (for employment: log(sum(n_vinculos)))
# NEW: xformla — formula for doubly-robust conditioning covariates (default ~1 = unconditional)
# NEW: covars_df — optional data frame with id_num + covariate columns to merge
run_cs <- function(d, outcome_var, label = "", agg = c("wmean", "logsum"),
                   xformla = ~1, covars_df = NULL) {
  agg <- match.arg(agg)

  d_agg <- if (agg == "wmean") {
    # Population-weighted average: larger municipalities contribute proportionally
    d %>%
      group_by(id_num, ano, g) %>%
      summarise(y = weighted.mean(.data[[outcome_var]], w = n_vinculos, na.rm = TRUE),
                .groups = "drop")
  } else {
    # Log of total employment count across all municipalities in the microregion
    d %>%
      group_by(id_num, ano, g) %>%
      summarise(y = log(sum(.data[[outcome_var]], na.rm = TRUE)),
                .groups = "drop")
  }

  # NEW: Merge pre-treatment covariates if provided (for doubly-robust estimation)
  if (!is.null(covars_df)) {
    d_agg <- d_agg %>% left_join(covars_df, by = "id_num")
  }

  n_treated <- d_agg %>% filter(g > 0) %>% pull(id_num) %>% n_distinct()
  n_notyet  <- d_agg %>% filter(g == 0) %>% pull(id_num) %>% n_distinct()
  n_obs     <- nrow(d_agg)

  if (n_treated < 20) stop(sprintf("[%s] Only %d treated units (need >= 20).", label, n_treated))
  if (n_notyet  < 20) stop(sprintf("[%s] Only %d not-yet-treated units (need >= 20).", label, n_notyet))
  message(sprintf("   [%s] %d treated | %d not-yet-treated | %d obs", label, n_treated, n_notyet, n_obs))

  # --- Core estimator: HARDCODED notyettreated, no adaptive logic ---
  att <- did::att_gt(
    yname                  = "y",
    tname                  = "ano",
    idname                 = "id_num",
    gname                  = "g",
    xformla                = xformla,
    data                   = d_agg,
    control_group          = "notyettreated",
    allow_unbalanced_panel = TRUE,
    print_details          = FALSE
  )

  # Dynamic event study
  es      <- did::aggte(att, type = "dynamic", na.rm = TRUE)
  es_tidy <- broom::tidy(es)

  # Simple ATT aggregate — the headline number for the abstract
  att_simple <- did::aggte(att, type = "simple", na.rm = TRUE)

  pretrend <- pretrend_ftest(att)
  message(sprintf("   [%s] ATT = %.4f (SE = %.4f) | Pre-trend Wald p = %.4f",
                  label, att_simple$overall.att, att_simple$overall.se, pretrend$p.value))

  list(
    att_gt   = att,
    es       = es,
    tidy     = es_tidy,
    simple   = att_simple,
    pretrend = pretrend,
    label    = label,
    n_obs    = n_obs,
    n_treated = n_treated,
    n_notyet  = n_notyet
  )
}

# Helper 4 — build a one-row summary tibble from a run_cs result
build_att_row <- function(res) {
  s <- res$simple
  tibble(
    Model      = res$label,
    ATT        = round(s$overall.att, 4),
    SE         = round(s$overall.se, 4),
    CI_low     = round(s$overall.att - 1.96 * s$overall.se, 4),
    CI_high    = round(s$overall.att + 1.96 * s$overall.se, 4),
    Pretrend_p = round(res$pretrend$p.value, 4),
    N_obs      = res$n_obs,
    N_treated  = res$n_treated,
    N_notyet   = res$n_notyet
  )
}

# Helper 5 — publication event-study plot (two overlaid series)
# ylab parameter allows different labels for wages vs employment outcomes.
plot_event_study <- function(df1, df2, l1, l2, c1, c2, title_txt,
                             ylab = "ATT on log(wages in min. wage units)") {
  ggplot() +
    geom_hline(yintercept = 0, color = "black",  linetype = "solid",  linewidth = 0.4) +
    geom_vline(xintercept = -1, color = "grey40", linetype = "dashed", linewidth = 0.5) +
    geom_ribbon(data = df2, aes(x = event.time, ymin = conf.low, ymax = conf.high),
                alpha = 0.12, fill = c2) +
    geom_line(data = df2, aes(x = event.time, y = estimate, color = l2), linewidth = 1.1) +
    geom_point(data = df2, aes(x = event.time, y = estimate, color = l2), size = 1.8) +
    geom_ribbon(data = df1, aes(x = event.time, ymin = conf.low, ymax = conf.high),
                alpha = 0.20, fill = c1) +
    geom_line(data = df1, aes(x = event.time, y = estimate, color = l1), linewidth = 1.1) +
    geom_point(data = df1, aes(x = event.time, y = estimate, color = l1), size = 1.8) +
    scale_color_manual(name = "", values = setNames(c(c1, c2), c(l1, l2))) +
    scale_x_continuous(breaks = seq(-10, 10, 2)) +
    labs(
      title    = title_txt,
      subtitle = "95% confidence bands | CS (2021), not-yet-treated comparison",
      x        = "Event time relative to first treatment year",
      y        = ylab
    )
}

# ==============================================================================
# SECTION 2: DOCUMENTATION & PROVENANCE (APPENDIX)
# ==============================================================================
message(" [1/12] Fetching Official Dictionaries...")

tryCatch({
  dict_prouni <- basedosdados::read_sql(
    "SELECT nome_coluna, chave, valor
     FROM `basedosdados.br_mec_prouni.dicionario`
     WHERE id_tabela = 'microdados'")
  write_csv(dict_prouni, "Documentation_Final/appendix_prouni_codes.csv")
}, error = function(e) message("  Warning: ProUni dict skipped."))

tryCatch({
  dict_rais <- basedosdados::read_sql(
    "SELECT nome_coluna, chave, valor
     FROM `basedosdados.br_me_rais.dicionario`
     WHERE id_tabela = 'microdados_vinculos'")
  write_csv(dict_rais, "Documentation_Final/appendix_rais_codes.csv")
}, error = function(e) message("  Warning: RAIS dict skipped."))

tryCatch({
  dict_census <- basedosdados::read_sql(
    "SELECT nome_coluna, chave, valor
     FROM `basedosdados.br_ibge_censo_demografico.dicionario`
     WHERE id_tabela = 'microdados_pessoa_2010'")
  write_csv(dict_census, "Documentation_Final/appendix_census_codes.csv")
}, error = function(e) message("  Warning: Census dict skipped."))

# ==============================================================================
# SECTION 3: DATA EXTRACTION FROM BIGQUERY
# ==============================================================================
message(" [2/12] Extracting Data from BigQuery...")

# 3A. Geography crosswalk (municipality -> microregion)
df_geo <- basedosdados::read_sql("
  SELECT
    id_municipio,
    id_microrregiao,
    nome_microrregiao,
    sigla_uf,
    nome_regiao AS regiao
  FROM `basedosdados.br_bd_diretorios_brasil.municipio`
")

# 3B. Census 2010 — youth population denominator + covariates
# NOTE: id_microrregiao is NOT selected from Census microdata — that column has
# unreliable one-digit values. The correct municipality->microregion mapping
# comes from df_geo (the official IBGE directory) via id_municipio joins.
df_census <- basedosdados::read_sql("
  SELECT
    id_municipio,
    v0601 AS sexo,
    v0606 AS raca_cor,

    SUM(peso_amostral) AS pop_18_24,

    /* Census IBGE v0606: 1=Branca, 2=Preta, 3=Amarela, 4=Parda, 5=Indigena */
    SUM(CASE WHEN v0606 IN ('2', '4')
             THEN peso_amostral ELSE 0 END) AS pop_nonwhite,

    SUM(CASE WHEN v6400 = '3'
              AND (v0633 IS NULL OR CAST(v0633 AS INT64) < 11)
             THEN peso_amostral ELSE 0 END) AS pop_18_24_eligible,

    AVG(CASE WHEN CAST(v6511 AS FLOAT64) = 0 OR v6511 IS NULL
             THEN 1.0 ELSE 0.0 END) AS share_no_formal_income,
    AVG(v6532) AS avg_hh_income_pc_mw,
    AVG(CASE WHEN v0627 = '1' THEN 1.0 ELSE 0.0 END) AS literacy_rate,
    AVG(CASE WHEN v0628 = '1' THEN 1.0 ELSE 0.0 END) AS share_in_school,
    AVG(CASE WHEN v0616 IN ('3','4') THEN 1.0 ELSE 0.0 END) AS share_mobility_disability

  FROM `basedosdados.br_ibge_censo_demografico.microdados_pessoa_2010`
  WHERE
    v6036 BETWEEN 18 AND 24
    AND v0601 NOT IN ('9')
    AND v0606 NOT IN ('9')
  GROUP BY
    id_municipio,
    v0601,
    v0606
")

# 3C. ProUni scholarships — batched download with demographic stratification
message("   -> Downloading ProUni (Batched Strategy)...")
# FIX: basedosdados uses 'sexo' and 'raca_cor' (not sexo_beneficiario/raca_cor_beneficiario)
# Verified from: github.com/basedosdados/mais/blob/master/bases/br_mec_prouni/microdados/table_config.yaml
q_prouni_base <- "
  SELECT
    ano,
    id_municipio,
    raca_cor,
    sexo,
    tipo_bolsa,
    COUNT(*) AS n_bolsas
  FROM `basedosdados.br_mec_prouni.microdados`
"
q_prouni_group <- "GROUP BY ano, id_municipio, raca_cor, sexo, tipo_bolsa"

p1 <- basedosdados::read_sql(paste(q_prouni_base, "WHERE ano BETWEEN 2005 AND 2009", q_prouni_group))
p2 <- basedosdados::read_sql(paste(q_prouni_base, "WHERE ano BETWEEN 2010 AND 2014", q_prouni_group))
p3 <- basedosdados::read_sql(paste(q_prouni_base, "WHERE ano BETWEEN 2015 AND 2019", q_prouni_group))
df_prouni_raw <- bind_rows(p1, p2, p3)
rm(p1, p2, p3); gc()

# 3D. RAIS outcomes — safe year-by-year loop
message("   -> Downloading RAIS (Safe Loop)...")
rais_list <- list()
for (y in params$years) {
  q <- paste0("
    SELECT
      ano,
      id_municipio,
      raca_cor,
      sexo,
      COUNT(*)                                             AS n_vinculos,
      AVG(LOG(GREATEST(valor_remuneracao_media_sm, 0.1)))  AS log_wage_sm
    FROM `basedosdados.br_me_rais.microdados_vinculos`
    WHERE
      ano = ", y, "
      AND idade BETWEEN ", params$age_min, " AND ", params$age_max, "
      AND valor_remuneracao_media_sm > 0
      AND tipo_vinculo = '10'
      AND vinculo_ativo_3112 = '1'
    GROUP BY ano, id_municipio, raca_cor, sexo
  ")
  tryCatch({
    rais_list[[as.character(y)]] <- basedosdados::read_sql(q)
    message(sprintf("     Year %d OK", y))
  }, error = function(e) {
    message(sprintf("     Error fetching year %d: %s", y, safe_err_msg(e)))
  })
  gc()
}

# CHECK: all years downloaded before we lose the list
missing_years <- setdiff(as.character(params$years), names(rais_list))
if (length(missing_years) > 0) {
  stop(sprintf("RAIS download incomplete — missing years: %s", paste(missing_years, collapse = ", ")))
}

df_rais <- bind_rows(rais_list)
rm(rais_list); gc()

# ==============================================================================
# SECTION 4: TREATMENT CONSTRUCTION
# ==============================================================================
# This section implements the Chapter 4 identification strategy:
#   D_{r,t} = 1000 * sum_{k=2005}^{t-4} scholarships_{r,k} / youth_pop_{r,2010}
#   FIX: Dose bins from INITIAL TREATMENT-YEAR distribution (D_at_entry),
#        NOT the 2019 cross-section (which caused look-ahead bias).
#   Low = (tau_20, tau_50], High > tau_75 of D_at_entry
#   Timing g = first year D_{r,t} > tau_20, else g = 0 (not-yet-treated)
# ==============================================================================
message(" [3/12] Defining Treatment Dose & Timing...")

# 4.1 Aggregate ProUni to microregion x year
df_prouni_agg <- df_prouni_raw %>%
  mutate(ano = as.integer(ano)) %>%
  inner_join(df_geo %>% distinct(id_municipio, id_microrregiao), by = "id_municipio") %>%
  group_by(id_microrregiao, ano) %>%
  summarise(new_schol = sum(n_bolsas, na.rm = TRUE), .groups = "drop")

# 4.2 Youth population denominator at microregion level
df_pop_micro <- df_census %>%
  inner_join(df_geo %>% distinct(id_municipio, id_microrregiao), by = "id_municipio") %>%
  group_by(id_microrregiao) %>%
  summarise(pop = sum(pop_18_24, na.rm = TRUE), .groups = "drop")

# 4.3 Build balanced panel and compute lagged density D_rt
df_panel <- expand_grid(
    id_microrregiao = unique(df_pop_micro$id_microrregiao),
    ano             = params$years
  ) %>%
  left_join(df_prouni_agg, by = c("id_microrregiao", "ano")) %>%
  mutate(new_schol = replace_na(new_schol, 0L)) %>%
  left_join(df_pop_micro, by = "id_microrregiao") %>%
  group_by(id_microrregiao) %>%
  arrange(ano) %>%
  mutate(
    cum_schol = cumsum(new_schol),
    # D_{r,t}: cumulative scholarship density per 1,000 youth, lagged 4 years
    D_rt = (dplyr::lag(cum_schol, params$lag_years, default = 0) / pop) * 1000
  ) %>%
  ungroup() %>%
  filter(pop > 0)

# 4.4 FIX: Compute dose thresholds from INITIAL TREATMENT-YEAR distribution
# Previously, bins were defined from the 2019 cross-section of D_rt. This
# conditions on post-treatment information and creates look-ahead bias.
# Now: (a) use a preliminary low threshold to identify when each microregion
# first enters treatment, (b) record D_at_entry at that moment, (c) compute
# tau cutoffs from THAT distribution.

# Step 1: Preliminary threshold — p10 of all nonzero D_rt values
# Deliberately low: only needs to identify "any meaningful exposure".
nonzero_D   <- df_panel$D_rt[df_panel$D_rt > 0]
tau_prelim  <- quantile(nonzero_D, 0.10, na.rm = TRUE)
message(sprintf("   FIX: Preliminary threshold (p10 of nonzero D_rt) = %.4f", tau_prelim))

# Step 2: For each microregion that eventually crosses tau_prelim,
# record its dose at the year it first crosses.
df_first_cross <- df_panel %>%
  filter(D_rt > tau_prelim) %>%
  group_by(id_microrregiao) %>%
  slice_min(ano, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(id_microrregiao, first_year = ano, D_at_entry = D_rt)

# Step 3: Define actual bin cutoffs from D_at_entry distribution
tau_20 <- quantile(df_first_cross$D_at_entry, params$q_low,  na.rm = TRUE)
tau_50 <- quantile(df_first_cross$D_at_entry, params$q_mid,  na.rm = TRUE)
tau_75 <- quantile(df_first_cross$D_at_entry, params$q_high, na.rm = TRUE)

message(sprintf("   FIX: Entry-based thresholds: tau_20=%.4f | tau_50=%.4f | tau_75=%.4f",
                tau_20, tau_50, tau_75))
message(sprintf("   FIX: %d microregions cross preliminary threshold (of %d total)",
                nrow(df_first_cross), n_distinct(df_panel$id_microrregiao)))

# 4.5 FIX: Classify dose bins from D_at_entry and compute treatment timing g
# Bins are assigned once, at entry, and never change — no look-ahead bias.
# D_rt_2019 is kept for descriptive tables only — NOT used for bin assignment.
df_2019 <- df_panel %>% filter(ano == max(params$years))

df_final <- df_panel %>%
  left_join(df_first_cross, by = "id_microrregiao") %>%
  left_join(
    df_2019 %>% select(id_microrregiao, D_rt_2019 = D_rt),
    by = "id_microrregiao"
  ) %>%
  group_by(id_microrregiao) %>%
  mutate(
    # FIX: Bins from D_at_entry (dose at first crossing), not D_rt_2019
    dose_bin = case_when(
      D_at_entry >  tau_20 & D_at_entry <= tau_50 ~ "Low",
      D_at_entry >  tau_75                         ~ "High",
      TRUE                                         ~ NA_character_
    ),
    # g = first year D_rt crosses tau_20 (the entry-distribution threshold)
    g = {
      crossing <- ano[D_rt > tau_20]
      if (length(crossing) == 0) 0L else as.integer(min(crossing))
    }
  ) %>%
  ungroup() %>%
  mutate(g = ifelse(is.infinite(g), 0L, g)) %>%
  add_micro_names(df_geo)

message(sprintf("   Microregions: %d total | %d with g > 0 | %d with g = 0",
                n_distinct(df_final$id_microrregiao),
                n_distinct(df_final$id_microrregiao[df_final$g > 0]),
                n_distinct(df_final$id_microrregiao[df_final$g == 0])))

# ==============================================================================
# SECTION 5: SANITY CHECKS (fail-fast)
# ==============================================================================
message(" [4/12] Running Sanity Checks...")

# CHECK 1: No duplicate (id_microrregiao, ano) in the treatment panel
n_dup <- df_final %>% count(id_microrregiao, ano) %>% filter(n > 1) %>% nrow()
if (n_dup > 0) stop(sprintf("FAIL: %d duplicate (id_microrregiao, ano) pairs.", n_dup))
message("   [OK] No duplicate (id_microrregiao, ano) pairs.")

# CHECK 2: D_rt = 0 before first meaningful year (2005 + lag_years = 2009)
pre_lag_bad <- df_final %>%
  filter(ano < (min(params$years) + params$lag_years), D_rt != 0) %>%
  nrow()
if (pre_lag_bad > 0) {
  stop(sprintf("FAIL: %d rows with non-zero D_rt before %d.",
               pre_lag_bad, min(params$years) + params$lag_years))
}
message(sprintf("   [OK] D_rt = 0 for all years before %d.", min(params$years) + params$lag_years))

# CHECK 3: Each dose bin has > 50 microregions
bin_sizes <- df_final %>%
  filter(!is.na(dose_bin)) %>%
  distinct(id_microrregiao, dose_bin) %>%
  count(dose_bin, name = "n_regions")
message("   Bin sizes:")
print(bin_sizes)
for (b in c("Low", "High")) {
  n_b <- bin_sizes$n_regions[bin_sizes$dose_bin == b]
  if (length(n_b) == 0 || n_b <= 50) {
    stop(sprintf("FAIL: Bin '%s' has %s microregions (need > 50).", b,
                 ifelse(length(n_b) == 0, "0", as.character(n_b))))
  }
}
message("   [OK] Both bins have > 50 microregions.")

# CHECK 4: D_rt is weakly increasing over time
non_monotone <- df_final %>%
  group_by(id_microrregiao) %>%
  arrange(ano) %>%
  mutate(d_diff = D_rt - dplyr::lag(D_rt)) %>%
  filter(!is.na(d_diff), d_diff < -1e-10) %>%
  nrow()
if (non_monotone > 0) {
  warning(sprintf("NOTE: %d rows where D_rt decreased — check lag logic.", non_monotone))
} else {
  message("   [OK] D_rt is weakly increasing over time.")
}

message(" -> All sanity checks passed.")

# ==============================================================================
# SECTION 6: DESCRIPTIVE STATISTICS
# ==============================================================================
message(" [5/12] Generating Descriptive Tables...")

# Table 1: Bin summary
tab1 <- df_final %>%
  filter(ano == max(params$years), !is.na(dose_bin)) %>%
  group_by(dose_bin) %>%
  summarise(
    N_Regions          = n_distinct(id_microrregiao),
    Avg_Youth_Pop      = round(mean(pop, na.rm = TRUE), 0),
    Avg_Density_2019   = round(mean(D_rt_2019, na.rm = TRUE), 2),
    Avg_D_at_entry     = round(mean(D_at_entry, na.rm = TRUE), 2),
    Median_Treatment_Yr = median(g[g > 0], na.rm = TRUE),
    .groups = "drop"
  )
write_csv(tab1, "Tables_Final/table1_bin_summary.csv")
print(xtable(tab1, caption = "Summary Statistics by Dose Bin (Low / High)"),
      file = "Tables_Final/table1_bin_summary.tex")
message("   Table 1 saved.")

# Table 2: Covariate balance at microregion level
df_covars <- df_census %>%
  inner_join(df_geo %>% distinct(id_municipio, id_microrregiao), by = "id_municipio") %>%
  group_by(id_microrregiao) %>%
  summarise(
    pop_18_24              = sum(pop_18_24, na.rm = TRUE),
    pop_nonwhite           = sum(pop_nonwhite, na.rm = TRUE),
    share_no_formal_income = mean(share_no_formal_income, na.rm = TRUE),
    avg_hh_income_pc_mw    = mean(avg_hh_income_pc_mw, na.rm = TRUE),
    literacy_rate          = mean(literacy_rate, na.rm = TRUE),
    share_in_school        = mean(share_in_school, na.rm = TRUE),
    .groups = "drop"
  )

tab2 <- df_final %>%
  filter(ano == max(params$years), !is.na(dose_bin)) %>%
  distinct(id_microrregiao, dose_bin) %>%
  inner_join(df_covars, by = "id_microrregiao") %>%
  group_by(dose_bin) %>%
  summarise(
    across(
      c(pop_18_24, pop_nonwhite, share_no_formal_income,
        avg_hh_income_pc_mw, literacy_rate, share_in_school),
      list(mean = ~round(mean(.x, na.rm = TRUE), 3),
           sd   = ~round(sd(.x, na.rm = TRUE), 3))
    ),
    .groups = "drop"
  )
write_csv(tab2, "Tables_Final/table2_balance.csv")
print(xtable(tab2, caption = "Pre-Treatment Covariate Balance by Dose Bin"),
      file = "Tables_Final/table2_balance.tex")
message("   Table 2 saved.")

# NEW: Build covariate dataframe for doubly-robust estimation (CHANGE 2)
# These are pre-treatment (Census 2010), time-invariant covariates at the
# microregion level. id_num is the integer key used by run_cs().
df_covars_dr <- df_covars %>%
  mutate(
    id_num             = as.integer(id_microrregiao),
    log_pop_18_24      = log(pmax(pop_18_24, 1)),
    share_nonwhite     = ifelse(pop_18_24 > 0, pop_nonwhite / pop_18_24, 0)
  ) %>%
  select(id_num, log_pop_18_24, share_nonwhite,
         share_no_formal_income, avg_hh_income_pc_mw, literacy_rate)

# NEW: Formula for conditional (DR) specification
dr_xformla <- ~ log_pop_18_24 + share_nonwhite + share_no_formal_income +
                 avg_hh_income_pc_mw + literacy_rate
message(sprintf("   NEW: DR covariate dataframe: %d microregions, %d covariates",
                nrow(df_covars_dr), ncol(df_covars_dr) - 1L))

# ==============================================================================
# SECTION 7: CASE STUDIES (4 auto-selected representative microregions)
# ==============================================================================
message(" [6/12] Selecting Case Study Microregions...")

candidates <- df_final %>%
  filter(ano == max(params$years), !is.na(dose_bin), g > 0) %>%
  mutate(region_type = if_else(sigla_uf %in% params$central_states, "Central", "Marginal"))

case_micro <- bind_rows(
  candidates %>% filter(region_type == "Central",  dose_bin == "Low")  %>% slice_sample(n = 1),
  candidates %>% filter(region_type == "Central",  dose_bin == "High") %>% slice_sample(n = 1),
  candidates %>% filter(region_type == "Marginal", dose_bin == "Low")  %>% slice_sample(n = 1),
  candidates %>% filter(region_type == "Marginal", dose_bin == "High") %>% slice_sample(n = 1)
) %>%
  select(id_microrregiao, nome_microrregiao, sigla_uf, region_type, dose_bin, g, D_at_entry, D_rt_2019)

stopifnot("Need exactly 4 case study units" = nrow(case_micro) == 4)
write_csv(case_micro, "Tables_Final/table3_casestudy_units.csv")
message("   Selected case study microregions:")
print(case_micro)

# ==============================================================================
# SECTION 8: PREPARE REGRESSION DATASETS
# ==============================================================================
message(" [7/12] Building Regression Datasets...")

# 8.1 Build the FULL merged dataset (all races, both sexes)
# No hard-coded race filter here — subsets are created below for heterogeneity.
df_reg_all <- df_rais %>%
  inner_join(
    df_geo %>% select(id_municipio, id_microrregiao),
    by = "id_municipio"
  ) %>%
  mutate(
    ano    = as.integer(ano),
    id_num = as.integer(id_microrregiao)
  ) %>%
  inner_join(
    df_final %>% select(id_microrregiao, ano, g, dose_bin, sigla_uf),
    by = c("id_microrregiao", "ano")
  )

# 8.2 Demographic subsets
# RAIS uses MTE race coding (NOT IBGE):
#   1 = Indígena, 2 = Branca, 4 = Preta, 6 = Amarela, 8 = Parda, 9 = Não informado
# Non-white = codes 1 (Indígena), 4 (Preta), 8 (Parda)
# White     = code 2 (Branca)
df_nonwhite <- df_reg_all %>% filter(raca_cor %in% c("1", "4", "8"))
df_white    <- df_reg_all %>% filter(raca_cor == "2")

# Gender subsets (from non-white, matching main results population)
# RAIS sexo: 1 = Masculino, 2 = Feminino
df_male   <- df_nonwhite %>% filter(sexo == "1")
df_female <- df_nonwhite %>% filter(sexo == "2")

# 8.3 Dose-bin subsets for MAIN RESULTS (non-white, both genders)
# Each includes the dose-bin's treated units + all g=0 (not-yet-treated) units.
data_low  <- df_nonwhite %>% filter(dose_bin == "Low"  | g == 0)
data_high <- df_nonwhite %>% filter(dose_bin == "High" | g == 0)

# 8.4 Regional heterogeneity subsets (non-white, both genders)
data_cen <- df_nonwhite %>% filter(sigla_uf %in% params$central_states,  !is.na(dose_bin) | g == 0)
data_mar <- df_nonwhite %>% filter(!sigla_uf %in% params$central_states, !is.na(dose_bin) | g == 0)

# 8.5 Gender heterogeneity subsets
data_low_male   <- df_male   %>% filter(dose_bin == "Low"  | g == 0)
data_low_female <- df_female %>% filter(dose_bin == "Low"  | g == 0)
data_high_male  <- df_male   %>% filter(dose_bin == "High" | g == 0)
data_high_female <- df_female %>% filter(dose_bin == "High" | g == 0)

# 8.6 Race heterogeneity subsets
data_low_white    <- df_white    %>% filter(dose_bin == "Low"  | g == 0)
data_low_nonwhite <- df_nonwhite %>% filter(dose_bin == "Low"  | g == 0)
data_high_white   <- df_white    %>% filter(dose_bin == "High" | g == 0)
data_high_nonwhite <- df_nonwhite %>% filter(dose_bin == "High" | g == 0)

message(sprintf("   Main:    data_low=%d | data_high=%d rows", nrow(data_low), nrow(data_high)))
message(sprintf("   Region:  data_cen=%d | data_mar=%d rows",  nrow(data_cen), nrow(data_mar)))
message(sprintf("   Gender:  low_m=%d | low_f=%d | high_m=%d | high_f=%d",
                nrow(data_low_male), nrow(data_low_female),
                nrow(data_high_male), nrow(data_high_female)))
message(sprintf("   Race:    low_w=%d | low_nw=%d | high_w=%d | high_nw=%d",
                nrow(data_low_white), nrow(data_low_nonwhite),
                nrow(data_high_white), nrow(data_high_nonwhite)))

# 8.7 Case study wage trajectory plot
cs_ids <- case_micro$id_microrregiao
df_case_plot <- df_nonwhite %>%
  filter(id_microrregiao %in% cs_ids) %>%
  inner_join(
    case_micro %>% select(id_microrregiao, nome_microrregiao, dose_bin, region_type, g),
    by = "id_microrregiao"
  ) %>%
  mutate(event_time = ano - g) %>%
  group_by(id_microrregiao, nome_microrregiao, dose_bin, region_type, event_time) %>%
  summarise(avg_log_wage = mean(log_wage_sm, na.rm = TRUE), .groups = "drop") %>%
  filter(event_time >= -5, event_time <= 5)

p_cases <- ggplot(df_case_plot, aes(
    x = event_time, y = avg_log_wage,
    color = interaction(dose_bin, region_type, sep = " / "),
    group = id_microrregiao
  )) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0,  linetype = "dotted", color = "black") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  facet_wrap(~nome_microrregiao, scales = "free_y") +
  labs(
    title    = "Case Study: Individual Wage Trajectories",
    subtitle = "Event time relative to first treatment year (e = 0)",
    x = "Event time (e)", y = "Avg. log wage (min. wages)", color = "Bin / Region"
  ) +
  scale_x_continuous(breaks = -5:5)
ggsave("Figures_Final/fig_0_casestudy_trajectories.png", p_cases, width = 10, height = 7, dpi = 300)
message("   Case study figure saved.")

# ==============================================================================
# SECTION 9A: MAIN ESTIMATION — WAGES (non-white, 4 models)
# ==============================================================================
message(" [8/12] Running Main Wage Models (Callaway & Sant'Anna)...")

message("   -> Low Dose / Wages...")
res_low  <- run_cs(data_low,  "log_wage_sm", label = "Low/Wages")

message("   -> High Dose / Wages...")
res_high <- run_cs(data_high, "log_wage_sm", label = "High/Wages")

message("   -> Central / Wages...")
res_cen  <- run_cs(data_cen,  "log_wage_sm", label = "Central/Wages")

message("   -> Marginal / Wages...")
res_mar  <- run_cs(data_mar,  "log_wage_sm", label = "Marginal/Wages")

# ==============================================================================
# SECTION 9B: EMPLOYMENT MARGIN — n_vinculos (non-white, 4 models)
# ==============================================================================
message("   -> Low Dose / Employment...")
res_low_emp  <- run_cs(data_low,  "n_vinculos", label = "Low/Emp",  agg = "logsum")

message("   -> High Dose / Employment...")
res_high_emp <- run_cs(data_high, "n_vinculos", label = "High/Emp", agg = "logsum")

message("   -> Central / Employment...")
res_cen_emp  <- run_cs(data_cen,  "n_vinculos", label = "Central/Emp", agg = "logsum")

message("   -> Marginal / Employment...")
res_mar_emp  <- run_cs(data_mar,  "n_vinculos", label = "Marginal/Emp", agg = "logsum")

# ==============================================================================
# SECTION 9C: GENDER HETEROGENEITY — log wages (4 models)
# ==============================================================================
message(" [9/12] Running Gender Heterogeneity Models...")

message("   -> Male / Low Dose...")
res_male_low   <- run_cs(data_low_male,   "log_wage_sm", label = "Male/Low")

message("   -> Female / Low Dose...")
res_female_low <- run_cs(data_low_female, "log_wage_sm", label = "Female/Low")

message("   -> Male / High Dose...")
res_male_high  <- run_cs(data_high_male,  "log_wage_sm", label = "Male/High")

message("   -> Female / High Dose...")
res_female_high <- run_cs(data_high_female, "log_wage_sm", label = "Female/High")

# ==============================================================================
# SECTION 9D: RACE HETEROGENEITY — log wages (4 models)
# ==============================================================================
message("   Running Race Heterogeneity Models...")

message("   -> White / Low Dose...")
res_white_low    <- run_cs(data_low_white,    "log_wage_sm", label = "White/Low")

message("   -> Non-white / Low Dose...")
res_nonwhite_low <- run_cs(data_low_nonwhite, "log_wage_sm", label = "Nonwhite/Low")

message("   -> White / High Dose...")
res_white_high   <- run_cs(data_high_white,   "log_wage_sm", label = "White/High")

message("   -> Non-white / High Dose...")
res_nonwhite_high <- run_cs(data_high_nonwhite, "log_wage_sm", label = "Nonwhite/High")

# ==============================================================================
# SECTION 9E: ROBUSTNESS — ALTERNATIVE LAG (lag = 3 years, 2 models)
# ==============================================================================
message(" [10/12] Running Robustness: Alternative Lag (3 years)...")

# Recompute D_rt with lag = 3 on the same panel infrastructure
lag_alt <- 3L
df_panel_lag3 <- df_panel %>%
  group_by(id_microrregiao) %>%
  arrange(ano) %>%
  mutate(
    D_rt_alt = (dplyr::lag(cum_schol, lag_alt, default = 0) / pop) * 1000
  ) %>%
  ungroup()

# FIX: Apply same entry-based bin definition for lag=3 (no look-ahead bias)
nonzero_D_lag3   <- df_panel_lag3$D_rt_alt[df_panel_lag3$D_rt_alt > 0]
tau_prelim_lag3  <- quantile(nonzero_D_lag3, 0.10, na.rm = TRUE)

df_first_cross_lag3 <- df_panel_lag3 %>%
  filter(D_rt_alt > tau_prelim_lag3) %>%
  group_by(id_microrregiao) %>%
  slice_min(ano, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(id_microrregiao, first_year_alt = ano, D_at_entry_alt = D_rt_alt)

tau_20_alt <- quantile(df_first_cross_lag3$D_at_entry_alt, params$q_low,  na.rm = TRUE)
tau_50_alt <- quantile(df_first_cross_lag3$D_at_entry_alt, params$q_mid,  na.rm = TRUE)
tau_75_alt <- quantile(df_first_cross_lag3$D_at_entry_alt, params$q_high, na.rm = TRUE)

message(sprintf("   FIX: Lag-3 entry-based thresholds: tau_20=%.4f | tau_50=%.4f | tau_75=%.4f",
                tau_20_alt, tau_50_alt, tau_75_alt))

df_final_lag3 <- df_panel_lag3 %>%
  left_join(df_first_cross_lag3, by = "id_microrregiao") %>%
  group_by(id_microrregiao) %>%
  mutate(
    # FIX: Bins from D_at_entry_alt, not D_rt_2019
    dose_bin_alt = case_when(
      D_at_entry_alt >  tau_20_alt & D_at_entry_alt <= tau_50_alt ~ "Low",
      D_at_entry_alt >  tau_75_alt                                ~ "High",
      TRUE                                                         ~ NA_character_
    ),
    g_alt = {
      crossing <- ano[D_rt_alt > tau_20_alt]
      if (length(crossing) == 0) 0L else as.integer(min(crossing))
    }
  ) %>%
  ungroup() %>%
  mutate(g_alt = ifelse(is.infinite(g_alt), 0L, g_alt))

# Build robustness regression data (non-white only)
df_rob <- df_nonwhite %>%
  select(-g, -dose_bin) %>%
  inner_join(
    df_final_lag3 %>% select(id_microrregiao, ano, g_alt, dose_bin_alt),
    by = c("id_microrregiao", "ano")
  ) %>%
  rename(g = g_alt, dose_bin = dose_bin_alt)

data_low_lag3  <- df_rob %>% filter(dose_bin == "Low"  | g == 0)
data_high_lag3 <- df_rob %>% filter(dose_bin == "High" | g == 0)

message("   -> Low Dose / Lag=3...")
res_low_lag3  <- run_cs(data_low_lag3,  "log_wage_sm", label = "Low/Lag3")

message("   -> High Dose / Lag=3...")
res_high_lag3 <- run_cs(data_high_lag3, "log_wage_sm", label = "High/Lag3")

message(" -> All 18 baseline + 2 lag-robustness models completed.")

# ==============================================================================
# SECTION 9F: NEW — CONDITIONAL (DOUBLY-ROBUST) MODELS (CHANGE 2)
# ==============================================================================
# Run the primary wage models with pre-treatment Census covariates in xformla.
# This implements the doubly-robust estimator of Sant'Anna & Zhao (2020).
# Both unconditional (already estimated above) and conditional are stored;
# if results diverge substantially, it is flagged in the audit log.
# ==============================================================================
message(" [NEW] Running Conditional (Doubly-Robust) Wage Models...")

res_low_dr  <- tryCatch({
  message("   -> Low Dose / DR Wages...")
  run_cs(data_low, "log_wage_sm", label = "Low/DR",
         xformla = dr_xformla, covars_df = df_covars_dr)
}, error = function(e) { message(sprintf("   WARNING: Low/DR failed: %s", safe_err_msg(e))); NULL })

res_high_dr <- tryCatch({
  message("   -> High Dose / DR Wages...")
  run_cs(data_high, "log_wage_sm", label = "High/DR",
         xformla = dr_xformla, covars_df = df_covars_dr)
}, error = function(e) { message(sprintf("   WARNING: High/DR failed: %s", safe_err_msg(e))); NULL })

# NEW: Flag divergence between unconditional and conditional estimates
dr_divergence_flag <- FALSE
if (!is.null(res_low_dr) && !is.null(res_high_dr)) {
  for (pair in list(
    list(unc = res_low,  cond = res_low_dr,  lbl = "Low"),
    list(unc = res_high, cond = res_high_dr, lbl = "High")
  )) {
    diff_att <- abs(pair$unc$simple$overall.att - pair$cond$simple$overall.att)
    se_pool  <- sqrt(pair$unc$simple$overall.se^2 + pair$cond$simple$overall.se^2)
    if (se_pool > 0 && diff_att / se_pool > 2) {
      message(sprintf("   FLAG: %s unconditional vs DR ATTs diverge by %.1f pooled SEs!",
                      pair$lbl, diff_att / se_pool))
      dr_divergence_flag <- TRUE
    }
  }
  if (!dr_divergence_flag) message("   OK: Unconditional and DR estimates are broadly consistent.")
} else {
  message("   WARNING: DR comparison skipped (one or both models failed).")
}

# ==============================================================================
# SECTION 9G: NEW — HonestDiD SENSITIVITY ANALYSIS (CHANGE 3)
# ==============================================================================
# Runs HonestDiD relative-magnitudes sensitivity on the primary wage models.
# Computes breakdown values: largest Mbar for which the CI still excludes zero.
# ==============================================================================
message(" [NEW] Running HonestDiD Sensitivity Analysis...")

# NEW: Helper function for HonestDiD
run_honest <- function(es_obj, label) {
  betahat     <- es_obj$att.egt
  event_times <- es_obj$egt
  post_start  <- min(which(event_times >= 0))

  # Find VCov matrix — field name varies across did versions
  vcv_candidates <- intersect(
    c("V_analytical.aggte", "V.analytical.aggte", "V_analytical"),
    names(es_obj)
  )
  sigma <- if (length(vcv_candidates) > 0) es_obj[[vcv_candidates[1]]] else NULL

  if (is.null(sigma)) {
    warning(sprintf("HonestDiD: no VCov matrix found for %s", label))
    return(list(label = label, results = NULL, breakdown = NA_real_))
  }

  tryCatch({
    results <- HonestDiD::createSensitivityResults_relativeMagnitudes(
      betahat        = betahat,
      sigma          = sigma,
      numPrePeriods  = post_start - 1,
      numPostPeriods = length(betahat) - post_start + 1,
      Mbarvec        = c(0, 0.5, 1, 2)
    )
    bd <- max(c(results$Mbar[results$CI_lower > 0], 0), na.rm = TRUE)
    message(sprintf("   [%s] HonestDiD breakdown Mbar = %.2f", label, bd))
    list(label = label, results = results, breakdown = bd)
  }, error = function(e) {
    warning(sprintf("HonestDiD failed for %s: %s", label, safe_err_msg(e)))
    list(label = label, results = NULL, breakdown = NA_real_)
  })
}

honest_results <- list()
if (requireNamespace("HonestDiD", quietly = TRUE)) {
  for (spec in list(
    list(es = res_low$es,  lbl = "Low/Wages"),
    list(es = res_high$es, lbl = "High/Wages")
  )) {
    honest_results[[spec$lbl]] <- run_honest(spec$es, spec$lbl)
  }
  # NEW: Also run on DR models if available
  if (!is.null(res_low_dr))  honest_results[["Low/DR"]]  <- run_honest(res_low_dr$es,  "Low/DR")
  if (!is.null(res_high_dr)) honest_results[["High/DR"]] <- run_honest(res_high_dr$es, "High/DR")

  # NEW: Generate sensitivity plots
  for (nm in names(honest_results)) {
    hr <- honest_results[[nm]]
    if (!is.null(hr$results)) {
      tryCatch({
        p_honest <- ggplot(hr$results, aes(x = Mbar, y = estimate)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, fill = "#0072B2") +
          geom_line(color = "#0072B2", linewidth = 1) +
          geom_point(color = "#0072B2", size = 2) +
          labs(
            title    = sprintf("HonestDiD Sensitivity: %s", hr$label),
            subtitle = sprintf("Breakdown Mbar = %.2f", hr$breakdown),
            x = expression(bar(M) ~ "(relative magnitudes)"),
            y = "ATT (robust CI)"
          )
        fname <- sprintf("Figures_Final/fig_honest_%s.png",
                         gsub("[/ ]", "_", tolower(hr$label)))
        ggsave(fname, p_honest, width = 7, height = 5, dpi = 300)
        message(sprintf("   Saved: %s", fname))
      }, error = function(e) message(sprintf("   WARNING: HonestDiD plot failed for %s: %s", nm, safe_err_msg(e))))
    }
  }
} else {
  message("   SKIPPED: HonestDiD package not available.")
}

# ==============================================================================
# SECTION 9H: NEW — BUFFER-ZONE EXCLUSION TEST FOR SUTVA (CHANGE 4)
# ==============================================================================
# Exclude comparison microregions that border High-Dose treated units.
# If these untreated neighbours are affected by spillovers, dropping them
# should leave baseline results roughly unchanged (SUTVA holds).
# ==============================================================================
message(" [NEW] Running Buffer-Zone Exclusion Test (SUTVA)...")

res_low_buffer  <- NULL
res_high_buffer <- NULL
border_ids      <- integer(0)

if (requireNamespace("sf", quietly = TRUE) && requireNamespace("geobr", quietly = TRUE)) {
  tryCatch({
    # 1. Load microregion geometries
    micro_sf <- geobr::read_micro_region(year = 2010, simplified = TRUE)

    # 2. Build adjacency
    neighbors <- sf::st_touches(micro_sf)

    # 3. Find comparison units bordering High-Dose microregions
    high_dose_ids <- unique(df_final$id_microrregiao[df_final$dose_bin == "High"])
    for (i in seq_along(micro_sf$code_micro)) {
      code_i <- as.character(micro_sf$code_micro[i])
      if (code_i %in% high_dose_ids) next
      neighbor_codes <- as.character(micro_sf$code_micro[neighbors[[i]]])
      if (any(neighbor_codes %in% high_dose_ids)) {
        border_ids <- c(border_ids, code_i)
      }
    }
    border_ids <- unique(border_ids)
    message(sprintf("   %d comparison microregions border High-Dose units (excluded).",
                    length(border_ids)))

    # 4. Re-run Low and High on restricted sample
    data_low_nobuffer  <- data_low  %>% filter(!(id_microrregiao %in% border_ids))
    data_high_nobuffer <- data_high %>% filter(!(id_microrregiao %in% border_ids))

    message("   -> Low Dose / Buffer exclusion...")
    res_low_buffer  <- run_cs(data_low_nobuffer,  "log_wage_sm", label = "Low/Buffer")

    message("   -> High Dose / Buffer exclusion...")
    res_high_buffer <- run_cs(data_high_nobuffer, "log_wage_sm", label = "High/Buffer")

  }, error = function(e) message(sprintf("   WARNING: Buffer-zone test failed: %s", safe_err_msg(e))))
} else {
  message("   SKIPPED: sf/geobr packages not available.")
}

# ==============================================================================
# SECTION 9I: NEW — PLACEBO COHORT TEST (AGES 36-50) (CHANGE 5)
# ==============================================================================
# Workers aged 36-50 in RAIS are too old to have been ProUni beneficiaries.
# Non-zero effects on this cohort would indicate the results capture local
# economic shocks rather than the ProUni human capital channel.
#
# NOTE: Re-querying BigQuery for the full 2005-2019 panel is expensive.
# The exact query and pipeline are shown below as a placeholder.
# To run: uncomment the block and execute. Results feed into the same
# run_cs() + build_att_row() pipeline as all other models.
# ==============================================================================
message(" [NEW] Placebo Cohort Test (ages 36-50): PLACEHOLDER")
message("   To run this test, uncomment the block below and re-execute.")

res_placebo_low  <- NULL
res_placebo_high <- NULL

# --- UNCOMMENT TO RUN PLACEBO COHORT ---
# message("   -> Downloading RAIS placebo cohort (ages 36-50)...")
# placebo_list <- list()
# for (y in params$years) {
#   q_placebo <- paste0("
#     SELECT
#       ano,
#       id_municipio,
#       raca_cor,
#       sexo,
#       COUNT(*)                                             AS n_vinculos,
#       AVG(LOG(GREATEST(valor_remuneracao_media_sm, 0.1)))  AS log_wage_sm
#     FROM `basedosdados.br_me_rais.microdados_vinculos`
#     WHERE
#       ano = ", y, "
#       AND idade BETWEEN 36 AND 50
#       AND valor_remuneracao_media_sm > 0
#       AND tipo_vinculo = '10'
#       AND vinculo_ativo_3112 = '1'
#     GROUP BY ano, id_municipio, raca_cor, sexo
#   ")
#   tryCatch({
#     placebo_list[[as.character(y)]] <- basedosdados::read_sql(q_placebo)
#     message(sprintf("     Placebo year %d OK", y))
#   }, error = function(e) message(sprintf("     Error fetching placebo year %d: %s", y, safe_err_msg(e))))
#   gc()
# }
#
# df_rais_placebo <- bind_rows(placebo_list)
# rm(placebo_list); gc()
#
# # Build regression data for placebo cohort (non-white, same dose bins as main)
# df_placebo_reg <- df_rais_placebo %>%
#   filter(raca_cor %in% c("1", "4", "8")) %>%
#   inner_join(df_geo %>% select(id_municipio, id_microrregiao), by = "id_municipio") %>%
#   mutate(ano = as.integer(ano), id_num = as.integer(id_microrregiao)) %>%
#   inner_join(df_final %>% select(id_microrregiao, ano, g, dose_bin), by = c("id_microrregiao", "ano"))
#
# data_placebo_low  <- df_placebo_reg %>% filter(dose_bin == "Low"  | g == 0)
# data_placebo_high <- df_placebo_reg %>% filter(dose_bin == "High" | g == 0)
#
# message("   -> Placebo Low Dose / Wages...")
# res_placebo_low  <- run_cs(data_placebo_low,  "log_wage_sm", label = "Placebo/Low")
# message("   -> Placebo High Dose / Wages...")
# res_placebo_high <- run_cs(data_placebo_high, "log_wage_sm", label = "Placebo/High")
# --- END PLACEBO BLOCK ---

# ==============================================================================
# SECTION 9J: NEW — TWFE BENCHMARK (CHANGE 6)
# ==============================================================================
# Estimate the naive two-way fixed effects (TWFE) specification for comparison
# with the heterogeneity-robust CS estimator. TWFE is biased under staggered
# treatment timing, so we expect divergence.
# ==============================================================================
message(" [NEW] Running TWFE Benchmark...")

twfe_fit <- NULL
if (requireNamespace("fixest", quietly = TRUE)) {
  tryCatch({
    twfe_data <- df_nonwhite %>%
      group_by(id_num, ano, g) %>%
      summarise(y = weighted.mean(log_wage_sm, w = n_vinculos, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(treated = as.integer(g > 0 & ano >= g))

    twfe_fit <- fixest::feols(y ~ treated | id_num + ano, data = twfe_data, cluster = ~id_num)
    twfe_coef <- fixest::coeftable(twfe_fit)
    message(sprintf("   TWFE coeff = %.4f (SE = %.4f, p = %.4f)",
                    twfe_coef[1, "Estimate"], twfe_coef[1, "Std. Error"], twfe_coef[1, "Pr(>|t|)"]))
  }, error = function(e) message(sprintf("   WARNING: TWFE failed: %s", safe_err_msg(e))))
} else {
  message("   SKIPPED: fixest package not available.")
}

# ==============================================================================
# SECTION 9K: NEW — ALTERNATIVE DISCRETIZATION CUTOFFS (CHANGE 7)
# ==============================================================================
# Re-run the primary wage models with alternative bin boundaries to verify
# that results are not an artifact of the specific percentile choice.
# ==============================================================================
message(" [NEW] Running Alternative Cutoff Robustness...")

# NEW: Helper to recompute bins and re-run for a given pair of quantiles
run_alt_cutoffs <- function(q_lo, q_hi, suffix) {
  tau_lo <- quantile(df_first_cross$D_at_entry, q_lo, na.rm = TRUE)
  tau_hi <- quantile(df_first_cross$D_at_entry, q_hi, na.rm = TRUE)
  message(sprintf("   Cutoffs %s: tau_lo=%.4f | tau_hi=%.4f", suffix, tau_lo, tau_hi))

  df_alt <- df_final %>%
    mutate(
      dose_bin_alt = case_when(
        D_at_entry >  tau_20 & D_at_entry <= tau_lo ~ "Low",
        D_at_entry >  tau_hi                         ~ "High",
        TRUE                                         ~ NA_character_
      )
    )

  data_low_alt  <- df_nonwhite %>%
    select(-dose_bin) %>%
    inner_join(df_alt %>% select(id_microrregiao, ano, dose_bin_alt), by = c("id_microrregiao", "ano")) %>%
    rename(dose_bin = dose_bin_alt) %>%
    filter(dose_bin == "Low" | g == 0)

  data_high_alt <- df_nonwhite %>%
    select(-dose_bin) %>%
    inner_join(df_alt %>% select(id_microrregiao, ano, dose_bin_alt), by = c("id_microrregiao", "ano")) %>%
    rename(dose_bin = dose_bin_alt) %>%
    filter(dose_bin == "High" | g == 0)

  res_lo <- tryCatch(run_cs(data_low_alt,  "log_wage_sm", label = paste0("Low/", suffix)),
                     error = function(e) { message(sprintf("   WARNING: %s Low failed: %s", suffix, safe_err_msg(e))); NULL })
  res_hi <- tryCatch(run_cs(data_high_alt, "log_wage_sm", label = paste0("High/", suffix)),
                     error = function(e) { message(sprintf("   WARNING: %s High failed: %s", suffix, safe_err_msg(e))); NULL })
  list(low = res_lo, high = res_hi)
}

message("   -> tau_25 / tau_75 cutoffs...")
alt_25_75 <- run_alt_cutoffs(0.25, 0.75, "q25_75")
res_low_25_75  <- alt_25_75$low
res_high_25_75 <- alt_25_75$high

message("   -> tau_33 / tau_67 cutoffs...")
alt_33_67 <- run_alt_cutoffs(1/3, 2/3, "q33_67")
res_low_33_67  <- alt_33_67$low
res_high_33_67 <- alt_33_67$high

message(" -> All new robustness and sensitivity sections completed.")

# ==============================================================================
# SECTION 9L: ATT SUMMARY TABLE (UPDATED)
# ==============================================================================
message("   Building ATT summary table...")

# Collect all results — NULLs are silently dropped by Filter()
all_results <- Filter(Negate(is.null), list(
  res_low, res_high, res_cen, res_mar,
  res_low_emp, res_high_emp, res_cen_emp, res_mar_emp,
  res_male_low, res_female_low, res_male_high, res_female_high,
  res_white_low, res_nonwhite_low, res_white_high, res_nonwhite_high,
  res_low_lag3, res_high_lag3,
  res_low_dr, res_high_dr,
  res_low_buffer, res_high_buffer,
  res_placebo_low, res_placebo_high,
  res_low_25_75, res_high_25_75,
  res_low_33_67, res_high_33_67
))

tab_att <- bind_rows(lapply(all_results, build_att_row))

# NEW: Add TWFE row if available
if (!is.null(twfe_fit)) {
  tc <- fixest::coeftable(twfe_fit)
  tab_att <- bind_rows(tab_att, tibble(
    Model      = "TWFE (naive)",
    ATT        = round(tc[1, "Estimate"], 4),
    SE         = round(tc[1, "Std. Error"], 4),
    CI_low     = round(tc[1, "Estimate"] - 1.96 * tc[1, "Std. Error"], 4),
    CI_high    = round(tc[1, "Estimate"] + 1.96 * tc[1, "Std. Error"], 4),
    Pretrend_p = NA_real_,
    N_obs      = as.integer(fixest::nobs(twfe_fit)),
    N_treated  = NA_integer_,
    N_notyet   = NA_integer_
  ))
}

write_csv(tab_att, "Tables_Final/table4_summary_att.csv")
print(xtable(tab_att, caption = "Summary of ATT Estimates Across All Models",
             digits = c(0, 0, 4, 4, 4, 4, 4, 0, 0, 0)),
      file = "Tables_Final/table4_summary_att.tex")
message("   Table 4 (ATT summary) saved.")
print(tab_att)

# ==============================================================================
# SECTION 10: VISUALIZATION
# ==============================================================================
message(" [11/12] Generating Publication Figures...")

# Figure 1: Dose-Response — wages (High vs Low)
p1 <- plot_event_study(
  res_low$tidy, res_high$tidy,
  "Low Dose", "High Dose",
  "#E69F00", "#0072B2",
  "Dose-Response: High vs. Low Exposure Intensity"
)
ggsave("Figures_Final/fig_1_mechanism_dose.png", p1, width = 9, height = 5.5, dpi = 300)

# Figure 2: Regional Divergence — wages (Central vs Marginal)
p2 <- plot_event_study(
  res_cen$tidy, res_mar$tidy,
  "Central (South/Southeast)", "Marginal (North/Northeast)",
  "#2980b9", "#e74c3c",
  "Regional Divergence: Center vs. Periphery"
)
ggsave("Figures_Final/fig_2_result_region.png", p2, width = 9, height = 5.5, dpi = 300)

# Figure 3: Dose-Response — employment (High vs Low)
p3 <- plot_event_study(
  res_low_emp$tidy, res_high_emp$tidy,
  "Low Dose", "High Dose",
  "#E69F00", "#0072B2",
  "Employment Margin: High vs. Low Exposure",
  ylab = "ATT on log(formal employment)"
)
ggsave("Figures_Final/fig_3_employment_dose.png", p3, width = 9, height = 5.5, dpi = 300)

# Figure 4: Gender heterogeneity — wages (Male vs Female, Low and High)
p4_low <- plot_event_study(
  res_male_low$tidy, res_female_low$tidy,
  "Male", "Female",
  "#1b9e77", "#d95f02",
  "Gender Gap: Low Dose"
)
p4_high <- plot_event_study(
  res_male_high$tidy, res_female_high$tidy,
  "Male", "Female",
  "#1b9e77", "#d95f02",
  "Gender Gap: High Dose"
)
ggsave("Figures_Final/fig_4a_gender_low.png",  p4_low,  width = 9, height = 5.5, dpi = 300)
ggsave("Figures_Final/fig_4b_gender_high.png", p4_high, width = 9, height = 5.5, dpi = 300)

# Figure 5: Race heterogeneity — wages (White vs Non-white, Low and High)
p5_low <- plot_event_study(
  res_white_low$tidy, res_nonwhite_low$tidy,
  "White", "Non-white",
  "#7570b3", "#e7298a",
  "Race Gap: Low Dose"
)
p5_high <- plot_event_study(
  res_white_high$tidy, res_nonwhite_high$tidy,
  "White", "Non-white",
  "#7570b3", "#e7298a",
  "Race Gap: High Dose"
)
ggsave("Figures_Final/fig_5a_race_low.png",  p5_low,  width = 9, height = 5.5, dpi = 300)
ggsave("Figures_Final/fig_5b_race_high.png", p5_high, width = 9, height = 5.5, dpi = 300)

# Figure 6: Robustness — main spec (lag=4) vs alternative (lag=3)
p6 <- plot_event_study(
  res_low$tidy, res_low_lag3$tidy,
  "Baseline (lag = 4)", "Robustness (lag = 3)",
  "#0072B2", "#D55E00",
  "Robustness: Sensitivity to Lag Length (Low Dose)"
)
ggsave("Figures_Final/fig_6_robustness_lag.png", p6, width = 9, height = 5.5, dpi = 300)

# NEW: Figure 7 — Unconditional vs Doubly-Robust (if DR models succeeded)
if (!is.null(res_low_dr)) {
  p7 <- plot_event_study(
    res_low$tidy, res_low_dr$tidy,
    "Unconditional", "Doubly-Robust (DR)",
    "#E69F00", "#009E73",
    "Sensitivity: Unconditional vs. DR (Low Dose)"
  )
  ggsave("Figures_Final/fig_7_dr_low.png", p7, width = 9, height = 5.5, dpi = 300)
}
if (!is.null(res_high_dr)) {
  p7b <- plot_event_study(
    res_high$tidy, res_high_dr$tidy,
    "Unconditional", "Doubly-Robust (DR)",
    "#0072B2", "#009E73",
    "Sensitivity: Unconditional vs. DR (High Dose)"
  )
  ggsave("Figures_Final/fig_7b_dr_high.png", p7b, width = 9, height = 5.5, dpi = 300)
}

# NEW: Figure 8 — Buffer-zone exclusion comparison (if buffer models succeeded)
if (!is.null(res_low_buffer)) {
  p8 <- plot_event_study(
    res_low$tidy, res_low_buffer$tidy,
    "Baseline", "Buffer Excluded",
    "#E69F00", "#CC79A7",
    "SUTVA Robustness: Buffer-Zone Exclusion (Low Dose)"
  )
  ggsave("Figures_Final/fig_8_buffer_low.png", p8, width = 9, height = 5.5, dpi = 300)
}
if (!is.null(res_high_buffer)) {
  p8b <- plot_event_study(
    res_high$tidy, res_high_buffer$tidy,
    "Baseline", "Buffer Excluded",
    "#0072B2", "#CC79A7",
    "SUTVA Robustness: Buffer-Zone Exclusion (High Dose)"
  )
  ggsave("Figures_Final/fig_8b_buffer_high.png", p8b, width = 9, height = 5.5, dpi = 300)
}

# NEW: Figure 9 — Alternative cutoffs comparison
if (!is.null(res_low_25_75)) {
  p9 <- plot_event_study(
    res_low$tidy, res_low_25_75$tidy,
    "Baseline (p20/p50)", "Alt (p25/p75)",
    "#E69F00", "#56B4E9",
    "Cutoff Robustness: Low Dose"
  )
  ggsave("Figures_Final/fig_9_altcut_low.png", p9, width = 9, height = 5.5, dpi = 300)
}

# Save all event-study estimates as CSVs for appendix
es_exports <- list(
  es_low_wage = res_low$tidy, es_high_wage = res_high$tidy,
  es_cen_wage = res_cen$tidy, es_mar_wage = res_mar$tidy,
  es_low_emp = res_low_emp$tidy, es_high_emp = res_high_emp$tidy,
  es_cen_emp = res_cen_emp$tidy, es_mar_emp = res_mar_emp$tidy,
  es_male_low = res_male_low$tidy, es_female_low = res_female_low$tidy,
  es_male_high = res_male_high$tidy, es_female_high = res_female_high$tidy,
  es_white_low = res_white_low$tidy, es_nonwhite_low = res_nonwhite_low$tidy,
  es_white_high = res_white_high$tidy, es_nonwhite_high = res_nonwhite_high$tidy,
  es_low_lag3 = res_low_lag3$tidy, es_high_lag3 = res_high_lag3$tidy
)
# NEW: Add conditional and robustness event-study CSVs
if (!is.null(res_low_dr))      es_exports$es_low_dr      <- res_low_dr$tidy
if (!is.null(res_high_dr))     es_exports$es_high_dr     <- res_high_dr$tidy
if (!is.null(res_low_buffer))  es_exports$es_low_buffer  <- res_low_buffer$tidy
if (!is.null(res_high_buffer)) es_exports$es_high_buffer <- res_high_buffer$tidy
if (!is.null(res_low_25_75))   es_exports$es_low_25_75   <- res_low_25_75$tidy
if (!is.null(res_high_25_75))  es_exports$es_high_25_75  <- res_high_25_75$tidy
if (!is.null(res_low_33_67))   es_exports$es_low_33_67   <- res_low_33_67$tidy
if (!is.null(res_high_33_67))  es_exports$es_high_33_67  <- res_high_33_67$tidy

for (nm in names(es_exports)) {
  write_csv(es_exports[[nm]], sprintf("Tables_Final/%s.csv", nm))
}
message("   All figures and event-study CSVs saved.")

# ==============================================================================
# SECTION 11: AUDIT LOG
# ==============================================================================
message(" [12/12] Saving Audit Log...")

audit_log <- list(
  timestamp  = Sys.time(),
  params     = params,
  n_units    = n_distinct(df_final$id_microrregiao),
  n_years    = length(params$years),
  # FIX: Entry-based thresholds (not 2019 cross-section)
  bin_method     = "entry-based (D_at_entry)",
  tau_prelim     = tau_prelim,
  thresholds     = list(tau_20 = tau_20, tau_50 = tau_50, tau_75 = tau_75),
  thresholds_lag3 = list(tau_20 = tau_20_alt, tau_50 = tau_50_alt, tau_75 = tau_75_alt),
  bin_sizes  = list(
    Low  = bin_sizes$n_regions[bin_sizes$dose_bin == "Low"],
    High = bin_sizes$n_regions[bin_sizes$dose_bin == "High"]
  ),
  g_range = list(
    min = min(df_final$g[df_final$g > 0], na.rm = TRUE),
    max = max(df_final$g[df_final$g > 0], na.rm = TRUE)
  ),
  att_summary = tab_att,
  pretrend = list(
    Low_Wages      = res_low$pretrend,
    High_Wages     = res_high$pretrend,
    Central_Wages  = res_cen$pretrend,
    Marginal_Wages = res_mar$pretrend,
    Low_Emp        = res_low_emp$pretrend,
    High_Emp       = res_high_emp$pretrend,
    Central_Emp    = res_cen_emp$pretrend,
    Marginal_Emp   = res_mar_emp$pretrend,
    Male_Low       = res_male_low$pretrend,
    Female_Low     = res_female_low$pretrend,
    Male_High      = res_male_high$pretrend,
    Female_High    = res_female_high$pretrend,
    White_Low      = res_white_low$pretrend,
    Nonwhite_Low   = res_nonwhite_low$pretrend,
    White_High     = res_white_high$pretrend,
    Nonwhite_High  = res_nonwhite_high$pretrend,
    Low_Lag3       = res_low_lag3$pretrend,
    High_Lag3      = res_high_lag3$pretrend
  ),
  # NEW: Doubly-robust results
  dr_divergence_flag = dr_divergence_flag,
  dr_pretrend = list(
    Low_DR  = if (!is.null(res_low_dr))  res_low_dr$pretrend  else NA,
    High_DR = if (!is.null(res_high_dr)) res_high_dr$pretrend else NA
  ),
  # NEW: HonestDiD breakdown values
  honest_did = lapply(honest_results, function(x) x$breakdown),
  # NEW: Buffer-zone info
  buffer_zone = list(
    n_excluded = length(border_ids),
    excluded_ids = border_ids,
    pretrend_low  = if (!is.null(res_low_buffer))  res_low_buffer$pretrend  else NA,
    pretrend_high = if (!is.null(res_high_buffer)) res_high_buffer$pretrend else NA
  ),
  # NEW: TWFE benchmark
  twfe = if (!is.null(twfe_fit)) {
    tc <- fixest::coeftable(twfe_fit)
    list(coef = tc[1, "Estimate"], se = tc[1, "Std. Error"], p = tc[1, "Pr(>|t|)"])
  } else NA,
  # NEW: Alternative cutoff pretrends
  alt_cutoff_pretrend = list(
    Low_25_75  = if (!is.null(res_low_25_75))  res_low_25_75$pretrend  else NA,
    High_25_75 = if (!is.null(res_high_25_75)) res_high_25_75$pretrend else NA,
    Low_33_67  = if (!is.null(res_low_33_67))  res_low_33_67$pretrend  else NA,
    High_33_67 = if (!is.null(res_high_33_67)) res_high_33_67$pretrend else NA
  ),
  case_studies = case_micro
)
saveRDS(audit_log, "Documentation_Final/audit_log.rds")

# Print summary
n_total_models <- nrow(tab_att)
message("\n=== AUDIT SUMMARY ===")
message(sprintf("  Microregions: %d | Years: %d", audit_log$n_units, audit_log$n_years))
message(sprintf("  Bin method: %s (tau_prelim = %.4f)", audit_log$bin_method, audit_log$tau_prelim))
message(sprintf("  Bins: Low = %d | High = %d", audit_log$bin_sizes$Low, audit_log$bin_sizes$High))
message(sprintf("  Treatment year range: %d - %d", audit_log$g_range$min, audit_log$g_range$max))
message("  --- Headline ATTs (non-white, wages) ---")
message(sprintf("  Low:  %.4f (SE %.4f) | High: %.4f (SE %.4f)",
                res_low$simple$overall.att, res_low$simple$overall.se,
                res_high$simple$overall.att, res_high$simple$overall.se))
message(sprintf("  Central: %.4f | Marginal: %.4f",
                res_cen$simple$overall.att, res_mar$simple$overall.att))
message("  --- Pre-trend p-values (wages, main) ---")
message(sprintf("  Low=%.4f | High=%.4f | Central=%.4f | Marginal=%.4f",
                res_low$pretrend$p.value, res_high$pretrend$p.value,
                res_cen$pretrend$p.value, res_mar$pretrend$p.value))
if (!is.null(res_low_dr)) {
  message("  --- DR estimates ---")
  message(sprintf("  Low/DR: %.4f (SE %.4f) | High/DR: %.4f (SE %.4f)",
                  res_low_dr$simple$overall.att, res_low_dr$simple$overall.se,
                  if (!is.null(res_high_dr)) res_high_dr$simple$overall.att else NA,
                  if (!is.null(res_high_dr)) res_high_dr$simple$overall.se  else NA))
}
if (length(honest_results) > 0) {
  message("  --- HonestDiD breakdown Mbar ---")
  for (nm in names(honest_results)) {
    message(sprintf("  %s: %.2f", nm, honest_results[[nm]]$breakdown))
  }
}
if (!is.null(twfe_fit)) {
  tc <- fixest::coeftable(twfe_fit)
  message(sprintf("  --- TWFE benchmark: %.4f (SE %.4f) ---", tc[1, "Estimate"], tc[1, "Std. Error"]))
}
message(sprintf("  --- Models estimated: %d total ---", n_total_models))
message("=====================\n")

message("DONE. All outputs saved to Figures_Final/, Tables_Final/, Documentation_Final/.")
