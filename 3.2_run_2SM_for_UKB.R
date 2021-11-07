#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 17 August 2021

# This script runs two-sample MR for epigenetic clock acceleration
# measures and cancer in UK Biobank

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("your_working_directory") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot", "openxlsx")

#---------------------------------------------------------------------#
#                            Read exposure                             #----
#---------------------------------------------------------------------#
#read exposure
exp_dat <- read.table("data/exp_data.txt", header = T)

#---------------------------------------------------------------------#
#                            Outcomes                                  #----
#---------------------------------------------------------------------#

# Function for UKB outcome datasets
out_func_UKB <- function(file, name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read.table(paste("folder_containing_data", file,".txt", sep = ""), 
                            header = T) %>% format_data(snps = exp_dat$SNP, 
                                                        type = "outcome", 
                                                        snp_col = "SNP", 
                                                        beta_col = "beta", 
                                                        se_col = "se", 
                                                        effect_allele_col = "ALLELE1", 
                                                        other_allele_col = "ALLELE0",
                                                        chr_col = "CHR",
                                                        pos_col = "GENPOS",
                                                        eaf_col = "A1FREQ", 
                                                        pval_col = "P_BOLT_LMM_INF",
                                                        ncase_col = "ncase", 
                                                        ncontrol_col = "ncontrol")
  outcome_var$outcome <- name
  return(outcome_var)
}

lung_cancer <- out_func_UKB("lung_cancer_data", "lung cancer")
lung_cancer_unadj <- out_func_UKB("lung_cancer_unadj_data", "lung cancer unadjusted")
ovarian_cancer <- out_func_UKB("ovarian_cancer_data", "ovarian cancer")
breast_cancer <- out_func_UKB("breast_cancer_data", "breast cancer")
prostate_cancer <- out_func_UKB("prostate_cancer_data", "prostate cancer")
colorectal_cancer <- out_func_UKB("colorectal_cancer_data", "colorectal cancer")
pan_cancer <- out_func_UKB("pan_cancer_data", "pan cancer")
pan_cancer_incC44 <- out_func_UKB("pan_inclc44_cancer_data", "pan cancer inc C44")

# Combine outcome datasets
UKB_out_dat <- rbind(lung_cancer, lung_cancer_unadj, ovarian_cancer, breast_cancer, 
                     prostate_cancer, colorectal_cancer, pan_cancer, pan_cancer_incC44)

#---------------------------------------------------------------------#
#                            Harmonisation                             #----
#---------------------------------------------------------------------#
# Harmonise exposure and outcome datasets
data_UKB <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = UKB_out_dat
)
#---------------------------------------------------------------------#
#                            Results                                   #----
#---------------------------------------------------------------------#
# Run two-sample MR
results_UKB <- mr(data_UKB, method_list = c("mr_ivw", "mr_egger_regression", 
                                        "mr_weighted_median", "mr_weighted_mode"))

# NOTE: BOLT-LMM transformation was not required here, as the data had already been transformed when we obtained it

#---------------------------------------------------------------------#
#                            Save results                             #----
#---------------------------------------------------------------------#

#write.table(results_UKB, "results/results/results_UKB_new.txt", row.names = F)

