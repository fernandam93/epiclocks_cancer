#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 12 October 2021

# This script runs two-sample MR for epigenetic clock acceleration
# measures and cancer in consortiums (additional subtypes)

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("your_working_directory") 

# Install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot", "openxlsx", "readxl")

#---------------------------------------------------------------------#
#                              Read exposure                          #----
#---------------------------------------------------------------------#

exp_dat <- read.table(file = "exp_data.txt", header = T)

#---------------------------------------------------------------------#
#                           Outcomes                              #----
#---------------------------------------------------------------------#

# Read BCAC data on breast cancer subtypes
out_func_bc <- function(name, beta, se)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read.xlsx("subtype_BC_raw.xlsx", 
                           colNames = T) %>% format_data(snps = exp_dat$SNP, 
                                                        type = "outcome", 
                                                        snp_col = "SNP.iCOGs", 
                                                        beta_col = beta, 
                                                        se_col = se, 
                                                        effect_allele_col = "Effect.Meta", 
                                                        other_allele_col = "Baseline.Meta",
                                                        chr_col = "chr.iCOGs",
                                                        pos_col = "Position.iCOGs",
                                                        eaf_col = "EAFcontrols.iCOGs", 
                                                        pval_col = "MTOP_p_value_meta")
  outcome_var$outcome <- name
  return(outcome_var)
}

luminal_A <- out_func_bc("Luminal A", "Luminal_A_log_or_meta", "Luminal_A_se_meta")
luminal_B <- out_func_bc("Luminal B", "Luminal_B_log_or_meta", "Luminal_B_se_meta")
luminal_B_HER2Neg <- out_func_bc("Luminal B HER2 Negative", "Luminal_B_HER2Neg_log_or_meta", "Luminal_B_HER2Neg_se_meta")
HER2_Enriched <- out_func_bc("HER2 Enriched", "HER2_Enriched_log_or_meta", "HER2_Enriched_se_meta")
Triple_Neg <- out_func_bc("Triple Negative", "Triple_Neg_log_or_meta", "Triple_Neg_se_meta")


#Read GECCO data on colorectal cancer subtypes
out_func_gecco_sub <- function(sheet, name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read_excel("GECCO_all.xlsx", sheet = sheet, col_names = T)  
  outcome_var$SNP <- exp_dat$SNP[match(outcome_var$Position, exp_dat$pos.exposure)]
  outcome_var <- outcome_var %>%  format_data(snps = exp_dat$SNP,
                                                type = "outcome",
                                                snp_col = "SNP",
                                                beta_col = "Effect",
                                                se_col = "StdErr",
                                                effect_allele_col = "Allele1",
                                                other_allele_col = "Allele2",
                                                chr_col = "Chr",
                                                pos_col = "Position",
                                                eaf_col = "Freq1",
                                                pval_col = "P.value")
  outcome_var$outcome <- name
  return(outcome_var)
}


colon <- out_func_gecco_sub("colon", "Colon cancer")
distal <- out_func_gecco_sub("distal", "Distal colon cancer")
female <- out_func_gecco_sub("female", "Female colorectal cancer")
male <- out_func_gecco_sub("male", "Male colorectal cancer")
proximal <- out_func_gecco_sub("proximal", "Proximal colon cancer")
rectal <- out_func_gecco_sub("rectal", "Rectal cancer")

# Upload outcome datasets including proxies (these are already in two-sample MR format)
# CIMBA subtypes
BRCA1_BC_new_proxies <- read.table("CIMBA_BRCA1_BC_replaced_proxies.txt", header = T)
BRCA2_BC_new_proxies <- read.table("CIMBA_BRCA2_BC_replaced_proxies.txt", header = T)
BRCA1_OC_new_proxies <- read.table("CIMBA_BRCA1_OC_replaced_proxies.txt", header = T)
BRCA2_OC_new_proxies <- read.table("CIMBA_BRCA2_OC_replaced_proxies.txt", header = T)
# PRACTICAL subtypes
PRACT_adv_proxies <- read.table("PRACTICAL_adv_replaced_proxies.txt", header = T)
PRACT_age55_proxies <- read.table("PRACTICAL_age55_replaced_proxies.txt", header = T)
PRACT_caseonly_proxies <- read.table("PRACTICAL_caseonly_replaced_proxies.txt", header = T)
PRACT_gleason_proxies <- read.table("PRACTICAL_gleason_replaced_proxies.txt", header = T)
PRACT_highvslow_proxies <- read.table("PRACTICAL_highvslow_replaced_proxies.txt", header = T)
PRACT_highvslowint_proxies <- read.table("PRACTICAL_highvslowint_replaced_proxies.txt", header = T)

# Combine outcomes
subtype_out_dat <- rbind(luminal_A, luminal_B, luminal_B_HER2Neg, HER2_Enriched, 
                         Triple_Neg, BRCA1_BC_new_proxies, BRCA2_BC_new_proxies, BRCA1_OC_new_proxies, BRCA2_OC_new_proxies, 
                         PRACT_adv_proxies, PRACT_age55_proxies, PRACT_caseonly_proxies, PRACT_gleason_proxies, 
                         PRACT_highvslow_proxies, PRACT_highvslowint_proxies, colon, rectal, male, female, proximal, distal)

#---------------------------------------------------------------------#
#                            Harmonisation                             #----
#---------------------------------------------------------------------#

# Harmonise exposure and outcome datasets
data_subtype <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = subtype_out_dat
)

#---------------------------------------------------------------------#
#                            Results                                   #----
#---------------------------------------------------------------------#
#MR
results_subtype <- mr(data_subtype, method_list = c("mr_ivw", "mr_egger_regression", 
                                            "mr_weighted_median", "mr_weighted_mode"))

#---------------------------------------------------------------------#
#                            Save results                             #----
#---------------------------------------------------------------------#

#write.table(results_subtype, "results_subtypes.txt", row.names = F)
