# preprocess_hypertension_data.R
# ----------------------------------------------
# Preprocessing UK Biobank data for hypertensive participants on losartan
# This script filters participants, applies quality control, and builds a losartan user cohort.

# Required packages
required_packages <- c(
  "data.table", "tidyverse", "qqman", "gridExtra", "grid", "rms", "genpwr",
  "mice", "rvest", "bigsnpr", "devtools", "dplyr"
)

# Install missing packages
installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load libraries
library(data.table)
library(tidyverse)
library(qqman)
library(gridExtra)
library(grid)
library(rms)
library(genpwr)
library(mice)
library(rvest)
library(bigsnpr)
library(devtools)
library(dplyr)

# Load helper functions
source("relevant_functions.R")  # <-- Replace with actual path or GitHub source if public

# Load data
bd <- as_tibble(fread("PLACEHOLDER_PATH/bd_enc.csv"))  # Replace with path to encoded UKB data
withdrawn_consent <- as_tibble(fread("PLACEHOLDER_PATH/withdrawn_consent2.csv"))

# Filter out participants who withdrew consent
bd <- bd %>% filter(!f.eid %in% withdrawn_consent$V1)

# Load data dictionary
dict <- as_tibble(fread("PLACEHOLDER_PATH/Data_Dictionary_Showcase.csv"))

# Identify participants with hypertension
bd_htn <- select(bd, contains(c("f.eid", "f.53.", "f.20002."))) %>%
  ukb_reshape_long(dict) %>%
  filter(`Non-cancer illness code, self-reported` %in% c(1065, 1072)) %>%
  filter(I == ".0.") %>%
  select(f.eid) %>%
  distinct()

# Touchscreen data (field 6150)
bd_htn_ts <- select(bd, contains(c("f.eid", "f.53.", "f.6150."))) %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0." & `Vascular/heart problems diagnosed by doctor` == "High blood pressure") %>%
  select(f.eid) %>%
  distinct()

# Combine HTN cohorts
bd_htn <- bind_rows(bd_htn, bd_htn_ts) %>% distinct()

# QC: Ethnicity, sex discrepancies, aneuploidy, outliers
bd_check <- select(bd, contains(c("f.eid", "f.53.", "f.31.", "f.22001.", "f.22006."))) %>%
  filter(f.eid %in% bd_htn$f.eid) %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0.")

bd_caucasian <- bd_check %>%
  filter(`Genetic ethnic grouping` == "Caucasian") %>%
  select(-`Genetic ethnic grouping`)

bd_sex_checked <- bd_caucasian %>%
  filter(Sex == `Genetic sex`) %>%
  select(f.eid)

x_aneup <- select(bd, contains(c("f.eid", "f.22019.")))
x_aneuploidy <- x_aneup %>%
  filter(!is.na(`f.22019.0.0`)) %>%
  select(f.eid)

bd_sex_checked_ane <- bd_sex_checked %>%
  filter(!f.eid %in% x_aneuploidy$f.eid)

# Remove heterozygosity/missing rate outliers
outliers <- as_tibble(fread("PLACEHOLDER_PATH/outliers.csv")) %>%
  filter(!is.na(f.22027.0.0))

bd_not_outliers <- bd_sex_checked_ane %>%
  filter(!f.eid %in% outliers$f.eid)

# Relatedness filtering
related <- as_tibble(fread("PLACEHOLDER_PATH/UKBB_relatedness.csv")) %>%
  filter(Kinship > 0.0884) %>%
  select(ID1, ID2)

missingness <- as_tibble(fread("PLACEHOLDER_PATH/missingness.csv")) %>%
  rename(ID1 = eid, miss_ID1 = `22005-0.0`)

related_ID1 <- left_join(related, missingness, by = "ID1")

missingness <- missingness %>%
  rename(ID2 = ID1, miss_ID2 = miss_ID1)

related_ID2 <- left_join(related_ID1, missingness, by = "ID2")

# Keep participant with lower missingness
new_related <- ifelse(related_ID2$miss_ID1 < related_ID2$miss_ID2,
                      related_ID2$ID1, related_ID2$ID2)
new_related_df <- tibble(f.eid = unique(new_related))

# Participants to exclude due to relatedness
related_ids <- unique(c(related$ID1, related$ID2))
for_exclusion <- tibble(f.eid = related_ids) %>%
  filter(!f.eid %in% new_related_df$f.eid)

bd_unrelated <- bd_not_outliers %>%
  filter(!f.eid %in% for_exclusion$f.eid)

# Medication info (field 20003)
bd_drugs <- select(bd, contains(c("f.eid", "f.53.", "f.20003."))) %>%
  filter(f.eid %in% bd_unrelated$f.eid) %>%
  ukb_reshape_long(dict) %>%
  distinct()

# Medication coding
dict_meds <- as_tibble(fread("PLACEHOLDER_PATH/Codings.csv")) %>%
  filter(Coding == 4) %>%
  select(Value, Meaning) %>%
  rename(`Treatment/medication code` = Value,
         `Treatment/medication` = Meaning) %>%
  mutate(`Treatment/medication code` = parse_integer(`Treatment/medication code`))

# Join with medication names
bd_drugs <- bd_drugs %>%
  left_join(dict_meds, by = "Treatment/medication code") %>%
  filter(I == ".0.") %>%
  select(-`Treatment/medication code`)

# Identify losartan users
losartan_drug_pattern <- "losartan|cozaar"
bd_los <- bd_drugs %>%
  filter(str_detect(tolower(`Treatment/medication`), losartan_drug_pattern)) %>%
  select(f.eid) %>%
  distinct()

bd_non_los <- bd_drugs %>%
  filter(!f.eid %in% bd_los$f.eid) %>%
  select(f.eid) %>%
  distinct()

# Final losartan-labeled HTN dataset
bd_htn3 <- bd_drugs %>%
  select(-`Treatment/medication`) %>%
  mutate(losartan = "other")

bd_htn3$losartan[bd_htn3$f.eid %in% bd_non_los$f.eid] <- "no losartan"
bd_htn3$losartan[bd_htn3$f.eid %in% bd_los$f.eid] <- "losartan"
bd_htn3$losartan <- factor(bd_htn3$losartan)

# Save intermediate output if desired
# write_csv(bd_htn3, "outputs/processed_losartan_htn.csv")

# Covariate selection
# Placeholder: physical activity, genotyping batch, etc. to be added in second chunk
# ----------------------------------------------------------
# CONTINUED IN SECOND SCRIPT CHUNK

# ===============================================
# Missing Data Summary and Covariate Preparation
# ===============================================

# --- Missing Data Summary for bd_HTN_3 ---

# Calculate percentage of missing values per column
missing_counts <- colSums(is.na(bd_HTN_3))
total_rows <- nrow(bd_HTN_3)
percentage_missing <- (missing_counts / total_rows) * 100

print(percentage_missing)

# ==========================
# Covariates Dataframe Setup
# ==========================

# Select relevant fields for covariate analysis
# Fields: sex (31), date of assessment (53), assessment center (54), BMI (21001),
# age (21003), smoking (20116, 20117), alcohol (1558), physical activity (6164),
# medications (20003), education, genotype array (22000), Townsend deprivation (189), etc.

bd_long <- bd %>%
  select(contains(c(
    "f.eid", "f.31.", "f.53.", "f.189.", "f.6164.", "f.21001.", "f.21003.",
    "f.20116.", "f.20117.", "f.30520.", "f.30530.", "f.30740.", "f.30880.",
    "f.22000.", "f.93.", "f.4080."
  ))) %>%
  filter(f.eid %in% bd_HTN_3$f.eid) %>%
  mutate_at(vars(contains("21003")), as.integer) %>%
  ukb_reshape_long(dict) %>%
  filter(I == ".0.") %>%
  distinct()

# --- Unique Participant Count ---
bd_long %>%
  select(f.eid) %>%
  distinct() %>%
  nrow()

# --- Missing Data Summary for bd_long ---
missing_counts <- colSums(is.na(bd_long))
total_rows <- nrow(bd_long)
percentage_missing <- (missing_counts / total_rows) * 100

print(percentage_missing)

# Note: Use `nrow()` for dataframes, `length()` for vectors
# Use `distinct()` or `unique()` to remove duplicated rows when one participant has multiple entries

# ==================================
# Physical Activity (PA) Processing
# ==================================

# Define PA levels (ordered by intensity)
pa_levels <- c(
  "Strenuous sports",
  "Other exercises (eg: swimming, cycling, keep fit, bowling)",
  "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)",
  "Light DIY (eg: pruning, watering the lawn)",
  "Walking for pleasure (not as a means of transport)",
  "None of the above",
  "Prefer not to answer"
)

# Assign PA levels and summarize most intense per participant
bd_pa <- bd_long %>%
  mutate(pa = factor(`Types of physical activity in last 4 weeks`, levels = pa_levels)) %>%
  select(f.eid, I, `Date of attending assessment centre`, pa) %>%
  group_by(f.eid) %>%
  summarise(pa = pa_levels[min(as.numeric(pa), na.rm = TRUE)], .groups = "drop")

# --- Missing Data Summary for bd_pa ---
missing_counts <- colSums(is.na(bd_pa))
total_rows <- nrow(bd_pa)
percentage_missing <- (missing_counts / total_rows) * 100

print(percentage_missing)

# =============================
# Genotype Array (Batch) Setup
# =============================

# Load and prepare genotype batch dictionary
dict_array <- fread("Codings.csv") %>%
  as_tibble() %>%
  filter(Coding == 22000) %>%
  select(Value, Meaning) %>%
  rename(array = Meaning) %>%
  mutate(Value = as.integer(Value))

# Join genotype array info to participants
bd_array <- bd_long %>%
  select(f.eid, I, `Date of attending assessment centre`, `Genotype measurement batch`) %>%
  distinct() %>%
  rename(
    date = `Date of attending assessment centre`,
    Value = `Genotype measurement batch`
  ) %>%
  left_join(dict_array, by = "Value") %>%
  select(f.eid, array) %>%
  mutate(array = case_when(
    startsWith(array, "UKBiLEVE") ~ "bileve",
    startsWith(array, "Batch") ~ "axiom",
    TRUE ~ NA_character_
  ))

# --- Unique Participant Count for bd_array ---
bd_array %>%
  select(f.eid) %>%
  distinct() %>%
  nrow()  # ~104,418 participants

# 03_preprocessing_covariates_outcomes.R
# ---------------------------------------
# Prepare covariates, medications, and outcomes for analysis
# Author: [Your Name]
# Date: [YYYY-MM-DD]

library(tidyverse)
library(data.table)
library(mice)
library(forcats)

# Extract and rename relevant covariates
bd_cov <- bd_long %>%
  select(f.eid, I, `Date of attending assessment centre`,
         `Age when attended assessment centre`, Sex,
         `Body mass index (BMI)`, `Smoking status`,
         `Alcohol drinker status`, `Townsend deprivation index at recruitment`) %>%
  unique() %>%
  rename(date = `Date of attending assessment centre`,
         age = `Age when attended assessment centre`,
         sex = Sex,
         bmi = `Body mass index (BMI)`,
         smoking = `Smoking status`,
         alcohol = `Alcohol drinker status`,
         townsend = `Townsend deprivation index at recruitment`)

# Add physical activity and genotyping array
bd_cov <- bd_cov %>%
  left_join(bd_pa, by = "f.eid") %>%
  left_join(bd_array, by = "f.eid")

# Prepare principal components (PCs 1–40)
pcs <- bd %>%
  select(contains(c("f.eid", "f.22009."))) %>%
  filter(f.eid %in% bd_HTN_3$f.eid) %>%
  rename_with(~ c("f.eid", paste0("C", 1:40)), everything())

# Merge PCs with covariates
bd_cov <- left_join(bd_cov, pcs, by = "f.eid")

# -----------------------------------------------
# Process medication usage
# -----------------------------------------------

# Prepare antihypertensive medication status
bd_HTN_3 <- bd_HTN_3 %>%
  select(f.eid, losartan)

bd_meds <- bd %>%
  select(contains(c("f.eid", "f.53.", "f.20003."))) %>%
  filter(f.eid %in% bd_HTN_3$f.eid) %>%
  ukb_reshape_long(dict) %>%
  rename(Date = `Date of attending assessment centre`,
         Coding = `Treatment/medication code`) %>%
  filter(I == ".0.")

dict_meds <- fread("medication_GWAS_dict.csv") %>%
  as_tibble()

bd_meds <- bd_meds %>%
  left_join(dict_meds, by = "Coding") %>%
  select(-Date, -Coding)

# Define medication classes by ATC codes
get_med_class <- function(atc_pattern) {
  bd_meds %>% filter(str_detect(`Medication ATC code`, atc_pattern)) %>%
    pull(f.eid) %>% unique()
}

# Medication groups
losartan_ids <- get_med_class("^C09CA01")
raas_ids <- setdiff(get_med_class("^C09"), losartan_ids)
diuretic_ids <- get_med_class("^C03")
bb_ids <- get_med_class("^C07")
ccb_ids <- get_med_class("^C08")

# Annotate medication use in bd_HTN_3
bd_HTN_3 <- bd_HTN_3 %>%
  mutate(
    losartan = factor(ifelse(losartan == "losartan", "yes", "no")),
    raas = if_else(f.eid %in% raas_ids, "yes", "no"),
    diuretics = if_else(f.eid %in% diuretic_ids, "yes", "no"),
    bb = if_else(f.eid %in% bb_ids, "yes", "no"),
    ccb = if_else(f.eid %in% ccb_ids, "yes", "no")
  )

# Combine medication data with covariates
bd_cov <- left_join(bd_cov, bd_HTN_3, by = "f.eid")

# -----------------------------------------------
# Recode factor variables
# -----------------------------------------------

# Collapse physical activity categories
bd_cov <- bd_cov %>%
  mutate(pa = fct_collapse(pa,
                           strenuous = "Strenuous sports",
                           moderate = c("Other exercises (eg: swimming, cycling, keep fit, bowling)",
                                        "Heavy DIY (eg: weeding, lawn mowing, carpentry, digging)"),
                           mild = c("Light DIY (eg: pruning, watering the lawn)",
                                    "Walking for pleasure (not as a means of transport)"),
                           none = "None of the above"))

# Define factor levels
factors <- c("sex", "smoking", "alcohol", "pa", "losartan", "array", "raas", "diuretics", "bb", "ccb")
levels <- list(
  c("Female", "Male"),
  c("Never", "Previous", "Current", "Prefer not to answer"),
  c("Never", "Previous", "Current", "Prefer not to answer"),
  c("none", "mild", "moderate", "strenuous", "Prefer not to answer"),
  c("no", "yes"),
  c("axiom", "bileve"),
  rep(list(c("no", "yes")), 4)
)

bd_cov[factors] <- bd_cov %>%
  select(all_of(factors)) %>%
  map2(levels, factor) %>%
  ukb_recode_factor_na()

# -----------------------------------------------
# Prepare outcomes: Systolic Blood Pressure (SBP)
# -----------------------------------------------

bd_outcomes <- bd_long %>%
  select(f.eid, I, `Date of attending assessment centre`,
         `Systolic blood pressure, manual reading`,
         `Systolic blood pressure, automated reading`) %>%
  unique() %>%
  rename(date = `Date of attending assessment centre`,
         manual = `Systolic blood pressure, manual reading`,
         automated = `Systolic blood pressure, automated reading`) %>%
  group_by(f.eid) %>%
  mutate(SBP = if_else(is.na(manual), automated, manual)) %>%
  summarise(SBP = mean(SBP, na.rm = TRUE))

# Combine covariates and outcome
bd_all_df <- left_join(bd_cov, bd_outcomes, by = "f.eid")

# Check missingness
missing_counts <- colSums(is.na(bd_all_df))
print(missing_counts)

# -----------------------------------------------
# Linear regression
# -----------------------------------------------

# Crude model
lm_crude <- lm(SBP ~ losartan, data = bd_all_df)
summary(lm_crude)

# Adjusted models
lm_adj1 <- lm(SBP ~ losartan + age + sex, data = bd_all_df)
summary(lm_adj1)

lm_adj2 <- lm(SBP ~ losartan + age + sex + bmi + townsend + smoking + pa + alcohol +
                array + raas + diuretics + bb + ccb, data = bd_all_df)
summary(lm_adj2)

# -----------------------------------------------
# Multiple Imputation
# -----------------------------------------------

tempData <- bd_all_df %>%
  select(-f.eid, -I, -date) %>%
  mice(m = 1, maxit = 10, seed = 7)

bd_imp <- complete(tempData, 1) %>%
  as_tibble() %>%
  mutate(f.eid = bd_all_df$f.eid)

# Save imputed dataset
write_rds(bd_imp, "bd_all_df_imp.rds")

# Check missingness post-imputation
print(colSums(is.na(bd_imp)))

# Repeat regression with imputed data
lm_imp <- lm(SBP ~ losartan + age + sex + bmi + townsend + smoking + pa + alcohol +
               array + raas + diuretics + bb + ccb, data = bd_imp)
summary(lm_imp)

# Export coefficients
coeff_df <- summary(lm_imp)$coefficients[, 1:4] %>%
  as.data.frame() %>%
  rename(Estimate = 1, Std.Error = 2, t.value = 3, `Pr(>|t|)` = 4)

write.csv(coeff_df, file = "coefficients_df_adjusted.csv")

# Confidence intervals
conf_intervals <- confint(lm_imp)

# Export full sample file
write.table(bd_imp, file = "losartan_adj.sample", sep = " ", row.names = FALSE, quote = FALSE)


######BASH code
#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Example pipeline for genome-wide QTL, vQTL and GEI analyses
# This script processes chromosomes 1-22 using placeholder directories.
# Replace placeholders with actual paths when deploying.

# === Configuration ===

cohort="losartan"
outcome="sbp"

# Base directories (placeholders)
data_dir="data/genotypes"
sample_dir="data/samples"
keep_dir="data/keep_files"
output_dir="results"

# OSCA executable placeholder (add your version and install instructions)
osca_exec="tools/osca-0.46.1-linux-x86_64/osca-0.46.1"
# TODO: Add README note or setup script to install OSCA and plink binaries

# PLINK2 executable placeholder
plink2_exec="tools/plink2"

# Loop over chromosomes 1 to 22
for chr in {1..22}; do

  echo "Processing chromosome $chr ..."

  # Step 1: Convert BGEN to PLINK binary format with filters
  $plink2_exec \
    --bgen "${data_dir}/ukb22828_c${chr}_b0_v3.bgen" ref-first \
    --sample "${sample_dir}/ukb22828_c12_b0_v3_s487159.sample" \
    --keep "${keep_dir}/keep_${cohort}_${outcome}.txt" \
    --make-bed \
    --out "${output_dir}/${cohort}_chr_${chr}" \
    --geno 0.05 \
    --maf 0.05 \
    --hwe 0.00001 \
    --minimac3-r2-filter 0.3 \
    --hard-call-threshold 0.1 \
    --rm-dup force-first

  # Step 2: vQTL analysis with OSCA
  $osca_exec \
    --vqtl \
    --bfile "${output_dir}/${cohort}_chr_${chr}" \
    --pheno "${keep_dir}/${outcome}_${cohort}_osca_sens.txt" \
    --vqtl-mtd 2 \
    --out "${output_dir}/osca_${cohort}_chr_${chr}_${outcome}"

  # Step 3: QTL analysis with PLINK2 (GWAS)
  $plink2_exec \
    --bfile "${output_dir}/${cohort}_chr_${chr}" \
    --pheno "${keep_dir}/${outcome}_${cohort}_osca.txt" \
    --glm allow-no-covars \
    --out "${output_dir}/plink_${cohort}_chr_${chr}_${outcome}"

done

# Step 4: Merge OSCA vQTL results for Manhattan plot input
sed -n 1p "${output_dir}/osca_${cohort}_chr_1_${outcome}.vqtl" > "${output_dir}/osca_${cohort}_${outcome}.txt"
for chr in {1..22}; do
  sed '1d' "${output_dir}/osca_${cohort}_chr_${chr}_${outcome}.vqtl" >> "${output_dir}/osca_${cohort}_${outcome}.txt"
done
bgzip --threads 32 "${output_dir}/osca_${cohort}_${outcome}.txt"

# Step 5: Merge PLINK QTL GWAS results for Manhattan plot input
head -n 1 "${output_dir}/plink_${cohort}_chr_1_${outcome}.SBP.glm.linear" > "${output_dir}/plink_${cohort}_${outcome}_Manhattan.txt"
for chr in {1..22}; do
  grep -w ADD "${output_dir}/plink_${cohort}_chr_${chr}_${outcome}.SBP.glm.linear" >> "${output_dir}/plink_${cohort}_${outcome}_Manhattan.txt"
done
bgzip --threads 32 "${output_dir}/plink_${cohort}_${outcome}_Manhattan.txt"

echo "Pipeline completed."


######post processing of genetic association files:
# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(readr)
# NOTE: Add your custom plotting functions in the main repo or here:
# source('path_to/ukb_manhattan_plots.R')
# source('path_to/ukb_qq_gwas_plots.R')
# source('path_to/ukb_lambda.R')

# === Configuration ===
cohort <- "losartan"
outcome <- "sbp"
p_threshold <- 1e-4

# === Read and process vQTL results ===
vqtl_file <- paste0("results/osca_", cohort, "_", outcome, ".txt.gz")
vqtl_results <- fread(vqtl_file) %>%
  rename(CHR = Chr, MAF = freq, BP = bp)

# Filter significant vQTL SNPs
sig_vqtl <- filter(vqtl_results, P < p_threshold)
write.csv(sig_vqtl, file = paste0("results/osca_", cohort, "_", outcome, "_significant.csv"), row.names = FALSE)

cat("Total SNPs analyzed (vQTL):", nrow(vqtl_results), "\n")
cat("Significant vQTL SNPs:", nrow(sig_vqtl), "\n")
cat("Sample size summary (vQTL):\n")
print(summary(vqtl_results$NMISS))

# Plotting (requires your custom functions)
ukb_manhattan_plots(vqtl_results, cohort, outcome, "vqtl", genomewideP = p_threshold, suggestive = FALSE)
ukb_qq_gwas_plots(vqtl_results, cohort, outcome, "vqtl")
cat("Genomic inflation factor (lambda):", ukb_lambda(vqtl_results), "\n")

# === Read and process QTL results ===
qtl_file <- paste0("results/plink_", cohort, "_", outcome, "_Manhattan.txt.gz")
qtl_results <- fread(qtl_file)
colnames(qtl_results) <- c("CHR", "BP", "SNP", "alleleA", "alleleB", "PROVISIONAL_REF", "A1", "OMITTED", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")

sig_qtl <- filter(qtl_results, P < p_threshold)
write.csv(sig_qtl, file = paste0("results/plink_", cohort, "_", outcome, "_significant.csv"), row.names = FALSE)

cat("Total SNPs analyzed (QTL):", nrow(qtl_results), "\n")
cat("Significant QTL SNPs:", nrow(sig_qtl), "\n")
cat("Sample size summary (QTL):\n")
print(summary(qtl_results$OBS_CT))

ukb_manhattan_plots(qtl_results, cohort, outcome, "qtl", genomewideP = p_threshold, suggestive = FALSE)
ukb_qq_gwas_plots(qtl_results, cohort, outcome, "qtl")
cat("Genomic inflation factor (lambda):", ukb_lambda(qtl_results), "\n")

# === Note: SNP handling error to be added later ===
# TODO: Document SNP error exclusions when merging files.

# === Downstream analysis for vQTL (top SNP extraction) ===
top_snps <- sig_vqtl %>% select(SNP)
write.table(top_snps, "results/top_snps_losartan_osca.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Placeholder for PLINK extraction and merge commands (see README for details)

# === GEI analysis for top independent SNPs ===

# Load phenotype and covariate data
phenotype <- fread("data/phenotypes/sbp_losartan_osca_sens.txt")
covariates <- fread("data/covariates/covar_losartan_sbp.txt")

pheno <- left_join(phenotype, covariates, by = c("FID", "IID"))

# Load genotype data using bigsnpr (install if needed)
library(bigsnpr)
obj_bigSNP <- snp_attach("results/losartan_sbp.rds")  # placeholder path
genotype <- obj_bigSNP$genotypes
fam_order <- as.data.table(obj_bigSNP$fam) %>%
  rename(FID = family.ID, IID = sample.ID)

# Merge phenotype with genotype fam order to ensure matching
y <- left_join(fam_order, pheno, by = c("FID", "IID"))

# Load SNP list and filter out known problematic SNP (e.g. rs1609685)
snps_list <- fread("results/osca_losartan_sbp_ld.txt") %>%
  filter(SNP != "rs1609685")

# Load list of near-independent SNPs
snps_ind <- fread("results/osca_losartan_sbp_vqtl.snps", header = FALSE)
snps_list <- filter(snps_list, SNP %in% snps_ind$V1)

# Add columns for GEI results
snps_list <- snps_list %>%
  mutate(
    int_p = NA_real_,
    n_losartan = NA_integer_,
    beta_losartan = NA_real_,
    se_losartan = NA_real_,
    n_not_losartan = NA_integer_,
    beta_not_losartan = NA_real_,
    se_not_losartan = NA_real_
  )

# Loop over SNPs to fit interaction models
for (i in seq_len(nrow(snps_list))) {
  snp_pos <- which(obj_bigSNP$map$marker.ID == snps_list$SNP[i])
  y2 <- mutate(y, snp = genotype[, snp_pos])

  model <- lm(as.formula(paste0(outcome, " ~ snp + losartan + snp:losartan")), data = y2)
  coef_summary <- summary(model)$coefficients

  # Interaction p-value (snp:losartan)
  snps_list$int_p[i] <- coef_summary["snp:losartan", "Pr(>|t|)"]

  # Stratified effects
  snps_list$n_losartan[i] <- sum(y2$losartan == 1)
  snps_list$beta_losartan[i] <- coef_summary["snp", "Estimate"] + coef_summary["snp:losartan", "Estimate"]
  snps_list$se_losartan[i] <- coef_summary["snp", "Std. Error"]  # Note: This is an approximation

  snps_list$n_not_losartan[i] <- sum(y2$losartan == 0)
  snps_list$beta_not_losartan[i] <- coef_summary["snp", "Estimate"]
  snps_list$se_not_losartan[i] <- coef_summary["snp", "Std. Error"]
}

write.csv(snps_list, "results/osca_losartan_sbp_interaction_results.csv", row.names = FALSE)

# === End of analysis ===
