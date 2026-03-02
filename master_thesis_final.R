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
add_micro_names <- function(df, geo) {
  lkp <- geo %>%
    distinct(id_microrregiao, nome_microrregiao, sigla_uf, .keep_all = FALSE)
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
  # V.analytical contains the full variance-covariance matrix of ATT_ES
  vcv  <- es$V_analytical.aggte[pre_idx, pre_idx, drop = FALSE]
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
# Returns a list: att_gt object, aggte event study, tidy df, pre-trend test.
run_cs <- function(d, outcome_var, label = "") {
  d_agg <- d %>%
    group_by(id_num, ano, g) %>%
    summarise(y = mean(.data[[outcome_var]], na.rm = TRUE), .groups = "drop")

  n_treated <- d_agg %>% filter(g > 0) %>% pull(id_num) %>% n_distinct()
  n_notyet  <- d_agg %>% filter(g == 0) %>% pull(id_num) %>% n_distinct()

  if (n_treated < 20) stop(sprintf("[%s] Only %d treated units (need >= 20).", label, n_treated))
  if (n_notyet  < 20) stop(sprintf("[%s] Only %d not-yet-treated units (need >= 20).", label, n_notyet))
  message(sprintf("   [%s] %d treated | %d not-yet-treated units", label, n_treated, n_notyet))

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

  es      <- did::aggte(att, type = "dynamic", na.rm = TRUE)
  es_tidy <- broom::tidy(es)
  # NOTE: NAs in estimates are preserved — NOT replaced with 0

  pretrend <- pretrend_ftest(att)
  message(sprintf("   [%s] Pre-trend Wald p = %.4f", label, pretrend$p.value))

  list(att_gt = att, es = es, tidy = es_tidy, pretrend = pretrend, label = label)
}

# Helper 4 — publication event-study plot (two overlaid series)
plot_event_study <- function(df1, df2, l1, l2, c1, c2, title_txt) {
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
      subtitle = "95% confidence bands | Callaway & Sant'Anna (2021), not-yet-treated comparison",
      x        = "Event time relative to first treatment year",
      y        = "ATT (log wages, minimum wages)"
    )
}

# ==============================================================================
# SECTION 2: DOCUMENTATION & PROVENANCE (APPENDIX)
# ==============================================================================
message(" [1/10] Fetching Official Dictionaries...")

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
  # BUG FIX: old code wrote dict_rais here instead of dict_census
  write_csv(dict_census, "Documentation_Final/appendix_census_codes.csv")
}, error = function(e) message("  Warning: Census dict skipped."))

# ==============================================================================
# SECTION 3: DATA EXTRACTION FROM BIGQUERY
# ==============================================================================
message(" [2/10] Extracting Data from BigQuery...")

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
# FIX: alias is pop_18_24 (was pop_18_24_total, causing downstream mismatch)
# NEW: pop_nonwhite defined as raca_cor IN ('2','3') = Preta + Parda
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
q_prouni_base <- "
  SELECT
    ano,
    id_municipio,
    raca_cor_beneficiario,
    sexo_beneficiario,
    tipo_bolsa,
    COUNT(*) AS n_bolsas
  FROM `basedosdados.br_mec_prouni.microdados`
"
q_prouni_group <- "GROUP BY ano, id_municipio, raca_cor_beneficiario, sexo_beneficiario, tipo_bolsa"

p1 <- basedosdados::read_sql(paste(q_prouni_base, "WHERE ano BETWEEN 2005 AND 2009", q_prouni_group))
p2 <- basedosdados::read_sql(paste(q_prouni_base, "WHERE ano BETWEEN 2010 AND 2014", q_prouni_group))
p3 <- basedosdados::read_sql(paste(q_prouni_base, "WHERE ano BETWEEN 2015 AND 2019", q_prouni_group))
df_prouni_raw <- bind_rows(p1, p2, p3)
rm(p1, p2, p3); gc()

# 3D. RAIS outcomes — safe year-by-year loop
# NEW: adds sexo, renames to n_vinculos, adds tipo_vinculo + status_vinculo filters
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
      AND tipo_vinculo  = 10
      AND status_vinculo = 10
    GROUP BY ano, id_municipio, raca_cor, sexo
  ")
  tryCatch({
    rais_list[[as.character(y)]] <- basedosdados::read_sql(q)
    message(sprintf("     Year %d OK", y))
  }, error = function(e) message(sprintf("     Error fetching year %d: %s", y, e$message)))
  gc()
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
message(" [3/10] Defining Treatment Dose & Timing...")

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
# REMOVED: "Control" group, "Medium Dose" bin — only Low and High remain.
# g = first year D_rt crosses tau_20; g = 0 for regions that never cross
#     (these serve as not-yet-treated comparisons in did::att_gt, NOT as a
#      permanent "control" group).
df_final <- df_panel %>%
  left_join(
    df_2019 %>% select(id_microrregiao, D_rt_2019 = D_rt),
    by = "id_microrregiao"
  ) %>%
  group_by(id_microrregiao) %>%
  mutate(
    # 2-bin classification from 2019 cross-section of lagged density
    dose_bin = case_when(
      D_rt_2019 >  tau_20 & D_rt_2019 <= tau_50 ~ "Low",
      D_rt_2019 >  tau_75                        ~ "High",
      TRUE                                        ~ NA_character_
    ),
    # Treatment timing: first year D_rt crosses tau_20
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
message(" [4/10] Running Sanity Checks...")

# CHECK 1: No duplicate (id_microrregiao, ano) in the treatment panel
n_dup <- df_final %>% count(id_microrregiao, ano) %>% filter(n > 1) %>% nrow()
if (n_dup > 0) stop(sprintf("FAIL: %d duplicate (id_microrregiao, ano) pairs.", n_dup))
message("   [OK] No duplicate (id_microrregiao, ano) pairs.")

# CHECK 2: D_rt = 0 before first meaningful year (2005 + lag_years - 1 = 2008)
# With a 4-year lag, D_rt for 2005-2008 should be 0 (no prior cumulative to lag).
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

# CHECK 4: D_rt is weakly increasing over time for each microregion
# (cumulative density should not decrease)
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
message(" [5/10] Generating Descriptive Tables...")

# Table 1: Bin summary — replaces old 4-group table
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
message(" [6/10] Selecting Case Study Microregions...")

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

write_csv(case_micro, "Tables_Final/table3_casestudy_units.csv")
message("   Selected case study microregions:")
print(case_micro)

# Case study trajectory plot is generated AFTER df_reg_all is built (Section 8)

# ==============================================================================
# SECTION 8: PREPARE REGRESSION DATASETS
# ==============================================================================
message(" [7/10] Building Regression Datasets...")

# Merge RAIS outcomes with geography and treatment variables
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
  ) %>%
  # FIX: non-white = codes 2 (Preta), 4 (Parda), 8 (Indigena)
  # Old code incorrectly included code "1" (Branca/White)
  filter(raca_cor %in% c("2", "4", "8"))

# --- Dose-bin subsets ---
# Each subset includes the dose-bin's treated units PLUS all g=0 units
# (not-yet-treated comparisons). The did package with control_group="notyettreated"
# uses g=0 units as the comparison pool.
data_low  <- df_reg_all %>% filter(dose_bin == "Low"  | g == 0)
data_high <- df_reg_all %>% filter(dose_bin == "High" | g == 0)

# --- Regional heterogeneity subsets ---
# Include all binned microregions (Low + High) within each geographic area.
# Regions not in either bin are excluded.
data_cen <- df_reg_all %>% filter(sigla_uf %in% params$central_states,  !is.na(dose_bin) | g == 0)
data_mar <- df_reg_all %>% filter(!sigla_uf %in% params$central_states, !is.na(dose_bin) | g == 0)

message(sprintf("   data_low:  %d rows | data_high: %d rows", nrow(data_low), nrow(data_high)))
message(sprintf("   data_cen:  %d rows | data_mar:  %d rows", nrow(data_cen), nrow(data_mar)))

# --- Case study wage trajectory plot ---
cs_ids <- case_micro$id_microrregiao
df_case_plot <- df_reg_all %>%
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
# SECTION 9: CALLAWAY & SANT'ANNA ESTIMATION
# ==============================================================================
message(" [8/10] Running Callaway & Sant'Anna Models...")

message("   -> Low Dose Model...")
res_low  <- run_cs(data_low,  "log_wage_sm", label = "Low")

message("   -> High Dose Model...")
res_high <- run_cs(data_high, "log_wage_sm", label = "High")

message("   -> Central Region Model...")
res_cen  <- run_cs(data_cen,  "log_wage_sm", label = "Central")

message("   -> Marginal Region Model...")
res_mar  <- run_cs(data_mar,  "log_wage_sm", label = "Marginal")

message(" -> All models completed.")

# ==============================================================================
# SECTION 10: VISUALIZATION & AUDIT LOG
# ==============================================================================
message(" [9/10] Generating Publication Figures...")

# Figure 1: Dose-Response (High vs Low)
p1 <- plot_event_study(
  res_low$tidy, res_high$tidy,
  "Low Dose", "High Dose",
  "#E69F00", "#0072B2",
  "Dose-Response: High vs. Low Exposure Intensity"
)
ggsave("Figures_Final/fig_1_mechanism_dose.png", p1, width = 9, height = 5.5, dpi = 300)

# Figure 2: Regional Divergence (Central vs Marginal)
p2 <- plot_event_study(
  res_cen$tidy, res_mar$tidy,
  "Central (South/Southeast)", "Marginal (North/Northeast)",
  "#2980b9", "#e74c3c",
  "Regional Divergence: Center vs. Periphery"
)
ggsave("Figures_Final/fig_2_result_region.png", p2, width = 9, height = 5.5, dpi = 300)

# Save event-study estimates as CSVs for appendix
write_csv(res_low$tidy,  "Tables_Final/es_low_dose.csv")
write_csv(res_high$tidy, "Tables_Final/es_high_dose.csv")
write_csv(res_cen$tidy,  "Tables_Final/es_central.csv")
write_csv(res_mar$tidy,  "Tables_Final/es_marginal.csv")

message(" [10/10] Saving Audit Log...")

# ---------- AUDIT LOG ----------
audit_log <- list(
  timestamp  = Sys.time(),
  params     = params,
  n_units    = n_distinct(df_final$id_microrregiao),
  n_years    = length(params$years),
  thresholds = list(tau_20 = tau_20, tau_50 = tau_50, tau_75 = tau_75),
  bin_sizes  = list(
    Low  = bin_sizes$n_regions[bin_sizes$dose_bin == "Low"],
    High = bin_sizes$n_regions[bin_sizes$dose_bin == "High"]
  ),
  g_range = list(
    min = min(df_final$g[df_final$g > 0], na.rm = TRUE),
    max = max(df_final$g[df_final$g > 0], na.rm = TRUE)
  ),
  pretrend = list(
    Low      = res_low$pretrend,
    High     = res_high$pretrend,
    Central  = res_cen$pretrend,
    Marginal = res_mar$pretrend
  ),
  case_studies = case_micro
)
saveRDS(audit_log, "Documentation_Final/audit_log.rds")

# Print summary for console
message("\n=== AUDIT SUMMARY ===")
message(sprintf("  Microregions: %d | Years: %d", audit_log$n_units, audit_log$n_years))
message(sprintf("  Bins: Low = %d | High = %d", audit_log$bin_sizes$Low, audit_log$bin_sizes$High))
message(sprintf("  Treatment year range: %d - %d", audit_log$g_range$min, audit_log$g_range$max))
message(sprintf("  Pre-trend p-values: Low=%.4f | High=%.4f | Central=%.4f | Marginal=%.4f",
                res_low$pretrend$p.value, res_high$pretrend$p.value,
                res_cen$pretrend$p.value, res_mar$pretrend$p.value))
message("=====================\n")

message("DONE. All outputs saved to Figures_Final/, Tables_Final/, Documentation_Final/.")
