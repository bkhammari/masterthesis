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
#   - Dose bins are defined from the 2019 cross-section of the lagged density
#     D_{r,t}: Low = (tau_20, tau_50], High = D > tau_75.
#   - Treatment timing g = first year D_{r,t} > tau_20.
#
# OUTPUTS:
#   - 18 Callaway & Sant'Anna models (wages + employment x dose/region/gender/race)
#   - 2 robustness models (alternative lag = 3 years)
#   - 7 publication figures (Figures_Final/)
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
run_cs <- function(d, outcome_var, label = "") {
  d_agg <- d %>%
    group_by(id_num, ano, g) %>%
    summarise(y = mean(.data[[outcome_var]], na.rm = TRUE), .groups = "drop")

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
df_census <- basedosdados::read_sql("
  SELECT
    id_municipio,
    id_microrregiao,
    v0601 AS sexo,
    v0606 AS raca_cor,

    SUM(peso_amostral) AS pop_18_24,

    SUM(CASE WHEN v0606 IN ('2', '3')
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
    id_microrregiao,
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
      AND vinculo_ativo_3112 = 1
    GROUP BY ano, id_municipio, raca_cor, sexo
  ")
  tryCatch({
    rais_list[[as.character(y)]] <- basedosdados::read_sql(q)
    message(sprintf("     Year %d OK", y))
  }, error = function(e) message(sprintf("     Error fetching year %d: %s", y, e$message)))
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
#   Dose bins from 2019 cross-section: Low = (tau_20, tau_50], High > tau_75
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

# 4.4 Compute dose thresholds from 2019 cross-section
df_2019 <- df_panel %>% filter(ano == max(params$years))

tau_20 <- quantile(df_2019$D_rt, params$q_low,  na.rm = TRUE)
tau_50 <- quantile(df_2019$D_rt, params$q_mid,  na.rm = TRUE)
tau_75 <- quantile(df_2019$D_rt, params$q_high, na.rm = TRUE)

message(sprintf("   Thresholds: tau_20=%.4f | tau_50=%.4f | tau_75=%.4f", tau_20, tau_50, tau_75))

# 4.5 Classify dose bins and compute treatment timing g
# Only Low and High bins. g = first year D_rt crosses tau_20.
# Regions that never cross get g = 0 (not-yet-treated comparisons in did::att_gt).
df_final <- df_panel %>%
  left_join(
    df_2019 %>% select(id_microrregiao, D_rt_2019 = D_rt),
    by = "id_microrregiao"
  ) %>%
  group_by(id_microrregiao) %>%
  mutate(
    dose_bin = case_when(
      D_rt_2019 >  tau_20 & D_rt_2019 <= tau_50 ~ "Low",
      D_rt_2019 >  tau_75                        ~ "High",
      TRUE                                        ~ NA_character_
    ),
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
  select(id_microrregiao, nome_microrregiao, sigla_uf, region_type, dose_bin, g, D_rt_2019)

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
# Non-white = codes 2 (Preta), 4 (Parda), 8 (Indigena)
# White     = code 1 (Branca)
df_nonwhite <- df_reg_all %>% filter(raca_cor %in% c("2", "4", "8"))
df_white    <- df_reg_all %>% filter(raca_cor == "1")

# Gender subsets (from non-white, matching main results population)
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
res_low_emp  <- run_cs(data_low,  "n_vinculos", label = "Low/Emp")

message("   -> High Dose / Employment...")
res_high_emp <- run_cs(data_high, "n_vinculos", label = "High/Emp")

message("   -> Central / Employment...")
res_cen_emp  <- run_cs(data_cen,  "n_vinculos", label = "Central/Emp")

message("   -> Marginal / Employment...")
res_mar_emp  <- run_cs(data_mar,  "n_vinculos", label = "Marginal/Emp")

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

df_2019_lag3 <- df_panel_lag3 %>% filter(ano == max(params$years))
tau_20_alt <- quantile(df_2019_lag3$D_rt_alt, params$q_low,  na.rm = TRUE)
tau_50_alt <- quantile(df_2019_lag3$D_rt_alt, params$q_mid,  na.rm = TRUE)
tau_75_alt <- quantile(df_2019_lag3$D_rt_alt, params$q_high, na.rm = TRUE)

df_final_lag3 <- df_panel_lag3 %>%
  left_join(
    df_2019_lag3 %>% select(id_microrregiao, D_rt_2019_alt = D_rt_alt),
    by = "id_microrregiao"
  ) %>%
  group_by(id_microrregiao) %>%
  mutate(
    dose_bin_alt = case_when(
      D_rt_2019_alt > tau_20_alt & D_rt_2019_alt <= tau_50_alt ~ "Low",
      D_rt_2019_alt > tau_75_alt                               ~ "High",
      TRUE                                                      ~ NA_character_
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

message(" -> All 18 + 2 robustness models completed.")

# ==============================================================================
# SECTION 9F: ATT SUMMARY TABLE
# ==============================================================================
message("   Building ATT summary table...")

all_results <- list(
  res_low, res_high, res_cen, res_mar,
  res_low_emp, res_high_emp, res_cen_emp, res_mar_emp,
  res_male_low, res_female_low, res_male_high, res_female_high,
  res_white_low, res_nonwhite_low, res_white_high, res_nonwhite_high,
  res_low_lag3, res_high_lag3
)

tab_att <- bind_rows(lapply(all_results, build_att_row))
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
  ylab = "ATT on formal employment (n_vinculos)"
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
  thresholds = list(tau_20 = tau_20, tau_50 = tau_50, tau_75 = tau_75),
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
  case_studies = case_micro
)
saveRDS(audit_log, "Documentation_Final/audit_log.rds")

# Print summary
message("\n=== AUDIT SUMMARY ===")
message(sprintf("  Microregions: %d | Years: %d", audit_log$n_units, audit_log$n_years))
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
message("  --- Models estimated: 18 main + 2 robustness = 20 total ---")
message("=====================\n")

message("DONE. All outputs saved to Figures_Final/, Tables_Final/, Documentation_Final/.")
