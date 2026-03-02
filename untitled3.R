# ==============================================================================
# MASTER THESIS: THE DISTRIBUTIONAL CONSEQUENCES OF HIGHER EDUCATION EXPANSION
# AUTHOR: Baha Eddine Khammari
# DATE: 2026-01-27
# METHODOLOGY: Callaway, Goodman-Bacon & Sant'Anna (2024) - Continuous DiD
# ==============================================================================

rm(list = ls())
gc()

# 1. SETUP & LIBRARIES
# ------------------------------------------------------------------------------
library(basedosdados)
library(dplyr)
library(tidyr)
library(ggplot2)
library(did)        # Callaway & Sant'Anna Estimator
library(contdid)
library(readr)      # Data I/O
library(xtable)     # LaTeX Table Generation
library(scales)     # Plot Formatting
library(broom)      # Tidy Model Outputs
library(geobr)

# Set Project ID
basedosdados::set_billing_id("microdadosbrazil")

# Create Directory Structure
dirs <- c("Figures_Final", "Tables_Final", "Documentation_Final")
sapply(dirs, function(x) if(!dir.exists(x)) dir.create(x))

# Set Professional Theme
theme_set(theme_classic() + theme(
  plot.title = element_text(face = "bold", size = 14, color = "black"),
  plot.subtitle = element_text(size = 11, color = "grey40"),
  axis.text = element_text(size = 10, color = "black"),
  legend.position = "bottom",
  panel.grid.major.y = element_line(color = "grey95")
))

# ==============================================================================
# PART 1: DOCUMENTATION & PROVENANCE (APPENDIX)
# ==============================================================================
message(" [1/7] Fetching Official Dictionaries...")

# A. ProUni DICTIONARY
tryCatch({
  dict_prouni <- basedosdados::read_sql("SELECT nome_coluna, chave, valor FROM `basedosdados.br_mec_prouni.dicionario` WHERE id_tabela = 'microdados'")
  write_csv(dict_prouni, "Documentation_Final/appendix_prouni_codes.csv")
}, error = function(e) message("Warning: ProUni dict skipped."))

# B. RAIS DICTIONARY
tryCatch({
  dict_rais <- basedosdados::read_sql("SELECT nome_coluna, chave, valor FROM `basedosdados.br_me_rais.dicionario` WHERE id_tabela = 'microdados_vinculos'")
  write_csv(dict_rais, "Documentation_Final/appendix_rais_codes.csv")
}, error = function(e) message("Warning: RAIS dict skipped."))

# C. Census DICTIONARY

tryCatch({
  dict_census <- basedosdados::read_sql("SELECT nome_coluna, chave, valor FROM `basedosdados.br_ibge_censo_demografico.dicionario` WHERE id_tabela = 'microdados_pessoa_2010'")
  write_csv(dict_rais, "Documentation_Final/appendix_census_codes.csv")
}, error = function(e) message("Warning: Census dict skipped."))

# ==============================================================================
# PART 2: ROBUST DATA EXTRACTION (SAFE MODE)
# ==============================================================================
message(" [2/7] Extracting Data from BigQuery...")

# A. Geography
df_geo <- basedosdados::read_sql("
  SELECT 
    id_municipio, 
    id_microrregiao, 
    nome_microrregiao, 
    sigla_uf, 
    nome_regiao 
  FROM `basedosdados.br_bd_diretorios_brasil.municipio`
")

# B. Census (Denominator - 2010 Youth Population, Stratified)
df_census <- basedosdados::read_sql("
  SELECT 
    id_municipio,
    id_microrregiao,
    v0601  AS sexo,
    v0606  AS raca_cor,

    -- PRIMARY DENOMINATOR: all 18-24 youth (matches thesis Eq. 1 exactly)
    SUM(peso_amostral) AS pop_18_24_total,

    -- ROBUSTNESS DENOMINATOR: HS complete, not currently enrolled
    SUM(CASE WHEN v6400 = '3' 
              AND (v0633 IS NULL OR CAST(v0633 AS INT64) < 11)
             THEN peso_amostral ELSE 0 END) AS pop_18_24_eligible,

    -- PRE-TREATMENT COVARIATE 1: Share with zero formal labour income
    -- Captures baseline labour market exclusion (key for your North/NE heterogeneity)
    AVG(CASE WHEN CAST(v6511 AS FLOAT64) = 0 OR v6511 IS NULL 
             THEN 1.0 ELSE 0.0 END) AS share_no_formal_income,

    -- PRE-TREATMENT COVARIATE 2: Average household income per capita (in MW)
    -- Use as balance-table control, NOT as denominator filter
    AVG(v6532) AS avg_hh_income_pc_mw,

    -- PRE-TREATMENT COVARIATE 3: Literacy rate (proxy for baseline human capital)
    AVG(CASE WHEN v0627 = '1' THEN 1.0 ELSE 0.0 END) AS literacy_rate,

    -- PRE-TREATMENT COVARIATE 4: Share currently attending school
    -- Captures pre-existing educational momentum in region
    AVG(CASE WHEN v0628 = '1' THEN 1.0 ELSE 0.0 END) AS share_in_school,

    -- PRE-TREATMENT COVARIATE 5: Share with disability (marginalisation proxy)
    -- v0616 = difficulty walking/climbing stairs (most common Census disability var)
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



# C. ProUni (Numerator - 3 Chunks)
message("   -> Downloading ProUni (Batched Strategy)...")
q_prouni <- "SELECT ano, id_municipio FROM `basedosdados.br_mec_prouni.microdados`"
p1 <- basedosdados::read_sql(paste(q_prouni, "WHERE ano BETWEEN 2005 AND 2009"))
p2 <- basedosdados::read_sql(paste(q_prouni, "WHERE ano BETWEEN 2010 AND 2014"))
p3 <- basedosdados::read_sql(paste(q_prouni, "WHERE ano BETWEEN 2015 AND 2019"))
df_prouni_raw <- rbind(p1, p2, p3); rm(p1, p2, p3); gc()

# D. RAIS (Outcomes - Safe Loop)
message("   -> Downloading RAIS (Safe Loop)...")
rais_list <- list()
years <- 2005:2019

for(y in years) {
  q <- paste0("
    SELECT ano, id_municipio, raca_cor,
    COUNT(*) as emp_count,
    AVG(LOG(GREATEST(valor_remuneracao_media_sm, 0.1))) as log_wage_sm
    FROM `basedosdados.br_me_rais.microdados_vinculos`
    WHERE ano = ", y, " AND idade BETWEEN 22 AND 35 AND valor_remuneracao_media_sm > 0
    GROUP BY ano, id_municipio, raca_cor")
  
  tryCatch({
    rais_list[[paste0(y)]] <- basedosdados::read_sql(q)
  }, error = function(e) message(paste("Error fetching year:", y)))
  gc() 
}
df_rais <- bind_rows(rais_list)

# ==============================================================================
# PART 3: TREATMENT DEFINITION (FIXED FOR FULL COVERAGE)
# ==============================================================================
message(" [3/7] Defining Treatment Dose & Timing...")

# 3.1 Aggregate ProUni
df_prouni_agg <- df_prouni_raw %>%
  mutate(ano = as.numeric(ano)) %>%
  inner_join(df_geo, by="id_municipio") %>%
  group_by(id_microrregiao, ano) %>%
  summarise(new_schol = n(), .groups="drop")

# 3.2 Build Panel & Continuous Density
df_pop_micro <- df_census %>% inner_join(df_geo, by="id_municipio") %>% 
  group_by(id_microrregiao) %>% summarise(pop=sum(pop_18_24, na.rm=T)) %>% ungroup()

df_panel <- expand_grid(id_microrregiao = unique(df_pop_micro$id_microrregiao), ano = 2005:2019) %>%
  left_join(df_prouni_agg, by=c("id_microrregiao","ano")) %>%
  mutate(new_schol = replace_na(new_schol, 0)) %>%
  left_join(df_pop_micro, by="id_microrregiao") %>%
  group_by(id_microrregiao) %>% arrange(ano) %>%
  mutate(
    cum_schol = cumsum(new_schol),
    # CONTINUOUS VARIABLE: Density of Graduates (Lagged 4 Years)
    dens = (dplyr::lag(cum_schol, 4, default=0) / pop) * 1000
  ) %>% ungroup() %>% filter(pop > 0)

# 3.3 Define Dose Bins (THE FIX: Force Bottom 20% to be Control)
df_2019 <- df_panel %>% filter(ano == 2019)

# We use the 20th Percentile as the "Baseline" cutoff
tau_20 <- quantile(df_2019$dens, 0.20, na.rm=T) 
tau_50 <- quantile(df_2019$dens, 0.50, na.rm=T) 
tau_75 <- quantile(df_2019$dens, 0.75, na.rm=T) 

df_final <- df_panel %>% group_by(id_microrregiao) %>%
  mutate(
    final_dens = max(dens),
    
    # NEW DOSE DEFINITION:
    dose_group = case_when(
      final_dens <= tau_20 ~ "Control", # Bottom 20% is the new Baseline
      final_dens > tau_20 & final_dens < tau_50 ~ "Low Dose",
      final_dens >= tau_50 & final_dens < tau_75 ~ "Medium Dose",
      final_dens >= tau_75 ~ "High Dose"
    ),
    
    # NEW TIMING DEFINITION:
    # If you are in the "Control" group, we force g = 0 (Never Treated)
    # Otherwise, g is the year you crossed the Baseline (tau_20)
    g = if_else(dose_group == "Control", 0, min(ano[dens > tau_20], na.rm=T))
    
  ) %>% ungroup() %>%
  # Handle Infinites if regions never crossed threshold despite being non-control
  mutate(g = ifelse(is.infinite(g), 0, g)) %>%
  inner_join(df_geo %>% distinct(id_microrregiao, sigla_uf, nome_microrregiao), 
             by="id_microrregiao")

message("   -> Control Group Redefined as Bottom 20% (Baseline Intensity).")
# ==============================================================================
# PART 4: DESCRIPTIVE STATISTICS & CASE STUDY
# ==============================================================================
message(" [4/7] Generating Tables...")

# Table 1: Treated vs Control
tab1 <- df_final %>%
  filter(ano == 2010) %>%
  group_by(dose_group) %>%
  summarise(
    N_Regions = n(),
    Avg_Youth_Pop = round(mean(pop), 0),
    Avg_Scholarship_Density = round(mean(final_dens), 2)
  )
write_csv(tab1, "Tables_Final/table1_descriptives.csv")
print(xtable(tab1, caption="Summary Statistics by Intensity Group"), file="Tables_Final/table1_descriptives.tex")

# Table 2: Case Study
# The column 'nome_microrregiao' now exists in df_final, so filter will work.
targets <- c("São Paulo", "Fortaleza")
tab2 <- df_final %>%
  filter(ano == 2019) %>%
  filter(grepl(paste(targets, collapse="|"), nome_microrregiao, ignore.case=TRUE)) %>%
  select(Name=nome_microrregiao, State=sigla_uf, Treated_Year=g, Intensity=dose_group, Density=dens)

write_csv(tab2, "Tables_Final/table2_casestudy.csv")
print(tab2) # Verify output in console

# ==============================================================================
# PART 5: PREPARE REGRESSION DATASETS
# ==============================================================================
message(" [5/7] Splitting Data Subsets...")

# Merge Outcomes with Geography FIRST
df_reg_all <- df_rais %>%
  # 1. Join Geo (Only IDs to prevent dupes)
  inner_join(df_geo %>% select(id_municipio, id_microrregiao), by = "id_municipio") %>% 
  # 2. Numeric IDs
  mutate(ano = as.numeric(ano), id_num = as.numeric(id_microrregiao)) %>%
  # 3. Join Treatment (Brings 'g', 'dose_group', 'sigla_uf')
  inner_join(df_final %>% select(id_microrregiao, ano, g, dose_group, sigla_uf), 
             by = c("id_microrregiao", "ano")) %>%
  # 4. Filter Mechanism (Non-White)
  filter(raca_cor %in% c("1", "4", "8")) 

# Subsets for Analysis
data_low  <- df_reg_all %>% filter(dose_group %in% c("Control", "Low Dose"))
data_high <- df_reg_all %>% filter(dose_group %in% c("Control", "High Dose"))

# Define Regions: INCLUDES 'ES' (Espirito Santo) in Central
central_states <- c("SP", "RJ", "MG", "ES", "RS", "SC", "PR")

data_cen  <- df_reg_all %>% filter(sigla_uf %in% central_states)
data_mar  <- df_reg_all %>% filter(!sigla_uf %in% central_states)

# ==============================================================================
# PART 6: ESTIMATION (CS REGRESSIONS - FIXED)
# ==============================================================================
message(" [6/7] Running Callaway & Sant'Anna Models...")

# Wrapper Function
run_cs <- function(d, outcome_var) {
  # 1. Check if we have a valid Control Group (g=0)
  has_control <- any(d$g == 0)
  
  if(!has_control) {
    message("   ! Warning: No pure control group in this subset. Switching to 'notyettreated'.")
    ctrl_opt <- "notyettreated"
  } else {
    ctrl_opt <- "nevertreated"
  }
  
  d_agg <- d %>% group_by(id_num, ano, g) %>% summarise(y = mean(.data[[outcome_var]], na.rm=T), .groups='drop')
  
  # 2. Run CS Estimator
  att <- did::att_gt(yname="y", tname="ano", idname="id_num", gname="g", 
                     data=d_agg, 
                     control_group=ctrl_opt, # Adaptive Control
                     allow_unbalanced_panel=TRUE, 
                     print_details=FALSE)
  
  # 3. Event Study Aggregation
  es <- did::aggte(att, type="dynamic", na.rm=TRUE)
  
  broom::tidy(es) %>% mutate(
    estimate = ifelse(is.na(estimate), 0, estimate),
    conf.low = ifelse(is.na(conf.low), 0, conf.low),
    conf.high = ifelse(is.na(conf.high), 0, conf.high)
  )
}

# Run 4 Core Models
message("   -> Running Low Dose Model...")
es_low  <- run_cs(data_low, "log_wage_sm")

message("   -> Running High Dose Model...")
es_high <- run_cs(data_high, "log_wage_sm")

message("   -> Running Central Regional Model...")
es_cen  <- run_cs(data_cen, "log_wage_sm")

message("   -> Running Marginal Regional Model...")
es_mar  <- run_cs(data_mar, "log_wage_sm")

message("   -> Models Completed Successfully.")

# ==============================================================================
# PART 7: VISUALIZATION (FINAL PLOTS)
# ==============================================================================
message(" [7/7] Generating Publication Figures...")

plot_div <- function(df1, df2, l1, l2, c1, c2, title_txt) {
  ggplot() +
    geom_vline(xintercept = -1, color = 'grey', linetype = "dotted") + 
    geom_hline(yintercept = 0, color="black", linetype = "dotted") +
    geom_vline(xintercept = 4, linetype = "dashed", alpha=0.5) +
    
    geom_ribbon(data = df2, aes(x=event.time, ymin=conf.low, ymax=conf.high), alpha=0.1, fill=c2) +
    geom_line(data = df2, aes(x=event.time, y=estimate, color=l2), linewidth=1.2) +
    
    geom_ribbon(data = df1, aes(x=event.time, ymin=conf.low, ymax=conf.high), alpha=0.2, fill=c1) +
    geom_line(data = df1, aes(x=event.time, y=estimate, color=l1), linewidth=1.2) +
    
    scale_color_manual(name="", values=setNames(c(c1, c2), c(l1, l2))) +
    labs(title=title_txt, x="Years relative to First Cohort", y="Effect on Log Wages") +
    scale_x_continuous(breaks=seq(-10, 8, 2))
}

# Figure 1: Dose-Response (Mechanism)
p1 <- plot_div(es_low, es_high, "Low Dose", "High Dose", "#E69F00", "darkblue", "Dose-Response: High vs Low Intensity")
ggsave("Figures_Final/fig_1_mechanism_dose.png", p1, width=9, height=5)

# Figure 2: Regional Divergence (Main Result)
p2 <- plot_div(es_cen, es_mar, "Central (South/SE)", "Marginal (North/NE)", "#2980b9", "#e74c3c", "Regional Divergence: Center vs Periphery")
ggsave("Figures_Final/fig_2_result_region.png", p2, width=9, height=5)

message("DONE. All files saved to 'Figures_Final' and 'Tables_Final'.")