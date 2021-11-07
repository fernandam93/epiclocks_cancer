#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 17 August 2021

# This script creates exposure datasets for epigenetic clock acceleration
# measures, which can later be used in two-sample MR analyses
# This script also estimates R2 and F-statistics for each of the exposures

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("your_working_directory") 

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools")

#---------------------------------------------------------------------#
#                           Exposures                                  #----
#---------------------------------------------------------------------#
# Create function to extract GWAS significant SNPs and perform LD-clumping 
# Output in two-sample MR format
# NOTE: this takes really long if run in RStudio
exposure_func <- function(data_name, file, SNP) {
  x <- read.table(file = file, header = T)
  x$pheno <- "DNA methylation ageing"
  x <- format_data(x, type = "exposure", 
                   phenotype_col = "pheno",
                   snp_col = SNP,
                   beta_col = "Effect", 
                   se_col = "SE", 
                   eaf_col = "Freq1", 
                   effect_allele_col = "A1", 
                   other_allele_col = "A2", 
                   pval_col = "P", 
                   samplesize_col = "N", 
                   chr_col = "chr", 
                   pos_col = "bp")
  x$id.exposure <- data_name
  x <- clump_data(x, clump_p1 = 5e-08, clump_p2 = 5e-08)
}

# Apply function to raw epigenetic age acceleration datasets
GrimAge_exp_dat <- exposure_func("GrimAge","GrimAge_EUR_summary_statistics.txt", "rsID")
Hannum_exp_dat <- exposure_func("Hannum","Hannum_EUR_summary_statistics.txt", "roblrsID")
IEAA_exp_dat <- exposure_func("IEAA","IEAA_EUR_summary_statistics.txt", "rsID")
PhenoAge_exp_dat <- exposure_func("PhenoAge","PhenoAge_EUR_summary_statistics.txt", "rsID")

#combine exposures
exp_dat <- rbind(GrimAge_exp_dat, Hannum_exp_dat, IEAA_exp_dat, PhenoAge_exp_dat)

#save unique list of SNPs
#write.table(unique(exp_dat$SNP), "SNP_list.txt", row.names = F, col.names = F)

#---------------------------------------------------------------------#
#                          R2 and F-statistic                         #----
#---------------------------------------------------------------------#

# Calculate R2 and F statistics for each exposure dataset
#method 1
exp_dat$r2 <- (2 * (exp_dat$beta.exposure^2) * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure)) /
  (2 * (exp_dat$beta.exposure^2) * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure) +
     2 * exp_dat$samplesize.exposure * exp_dat$eaf.exposure * 
     (1 - exp_dat$eaf.exposure) * exp_dat$se.exposure^2)
exp_dat$F <- exp_dat$r2 * (exp_dat$samplesize.exposure - 2) / (1 - exp_dat$r2)
#method 2
# exp_dat$F_stat <- exp_dat$beta.exposure^2 / exp_dat$se.exposure^2
# exp_dat$R2_stat <- exp_dat$F_stat/(exp_dat$samplesize.exposure-2+exp_dat$F_stat)

# Calculate total R2 for each exposure dataset 
r2_func <- function(id)
{
  x <- exp_dat[which(exp_dat$id.exposure==id),]
  sum(x$r2, na.rm = T)
}

variance_GrimAge <- r2_func("GrimAge") # 0.47%
variance_Hannum <- r2_func("Hannum") # 1.48%
variance_IEAA <- r2_func("IEAA") # 4.41%
variance_PhenoAge <- r2_func("PhenoAge") #1.86%

# Calculate minimum F-statistic for each exposure dataset 
Fmin_func <- function(id)
{
  x <- exp_dat[which(exp_dat$id.exposure==id),]
  min(x$F, na.rm = T)
}

Fmin_GrimAge <- Fmin_func("GrimAge") # 31
Fmin_Hannum <- Fmin_func("Hannum") # 31
Fmin_IEAA <- Fmin_func("IEAA") # 31
Fmin_PhenoAge <- Fmin_func("PhenoAge") # 32 

# Calculate maximum F-statistic for each exposure dataset 
Fmax_func <- function(id)
{
  x <- exp_dat[which(exp_dat$id.exposure==id),]
  max(x$F, na.rm = T)
}

Fmax_GrimAge <- Fmax_func("GrimAge") # 45
Fmax_Hannum <- Fmax_func("Hannum") # 99
Fmax_IEAA <- Fmax_func("IEAA") # 240
Fmax_PhenoAge <- Fmax_func("PhenoAge") # 89 

# Calculate median F-statistic for each exposure dataset 
Fmedian_func <- function(id)
{
  x <- exp_dat[which(exp_dat$id.exposure==id),]
  median(x$F, na.rm = T)
}

Fmedian_GrimAge <- Fmedian_func("GrimAge") # 36
Fmedian_Hannum <- Fmedian_func("Hannum") # 38
Fmedian_IEAA <- Fmedian_func("IEAA") # 47
Fmedian_PhenoAge <- Fmedian_func("PhenoAge") # 45 

#save
#write.table(exp_dat, "exp_data.txt", row.names = F)
