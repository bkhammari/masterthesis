# How Equal is Equal Opportunity? 🎓🇧🇷
## Gender, Ethnic and Regional Disparities in Labor Market Outcomes of Higher Education Policies in Brazil

Welcome to the repository for my Master's Thesis in Economics. This project empirically investigates the distributional and causal impacts of Brazil's **ProUni** (Programa Universidade para Todos) higher education policy on early-career labor market outcomes.

### 📌 Overview

This thesis utilizes modern Difference-in-Differences (DiD) econometrics—specifically the **Callaway & Sant'Anna (2021) estimator with continuous treatment**—to analyze the heterogeneous effects of scholarship intensity across 558 Brazilian microregions. 

By merging comprehensive administrative microdata from the Ministry of Education (ProUni) and the Ministry of Labor (RAIS), this research uncovers critical insights into how higher education expansion intersects with existing structural disparities in Brazil's formal labor market.

### 🎯 Key Objectives & Contributions
- **Causal Inference**: Estimates the true causal effect of ProUni on early-career formal earnings (log wages) and employment margins.
- **Heterogeneity Analysis**: Explores disparities across:
  - **Gender** 
  - **Demographic Background** (Ethnicity)
  - **Geography** (North/Northeast vs. South/Southeast)
- **Methodological Rigor**: Addresses staggered treatment timing and heterogeneous treatment effects using robust continuous-treatment DiD designs (`did` and `contdid` packages in R).

### 📊 Data & Environment

This repository contains the data pipelines and analysis scripts required to reproduce the study.

**Primary Data Sources:**
- **RAIS** (Relação Anual de Informações Sociais): Employer-employee administrative panel (2005-2019).
- **ProUni**: Administrative scholarship allocation records.
- **IBGE Census 2010**: Demographic baseline data.

**Data Access:**
Data is efficiently extracted and processed directly from Google BigQuery via the Brazilian public data initiative [`basedosdados`](https://basedosdados.org/). Geographic variables are managed using `geobr`.

### 🛠️ Tech Stack & Methods
- **Language**: R
- **Econometric Packages**: `did`, `contdid`, `fixest`, `modelsummary`
- **Data Manipulation**: `tidyverse`, `data.table`, `basedosdados`
- **Visualization**: `ggplot2`
- **Document Rendering**: LaTeX (Overleaf template structure included)

### 📂 Repository Structure

unknown
├── Thesis_Main/             # LaTeX thesis template and text files
└── README.md
