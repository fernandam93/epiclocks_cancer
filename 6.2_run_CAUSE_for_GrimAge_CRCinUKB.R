
#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 18 August 2021

# This script runs CAUSE for GrimAge acceleration and colorectal cancer
# in UK Biobank 
# This needs to be run using the terminal

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Install packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "readr", "ieugwasr")
pacman::p_load_gh("jean997/cause@v1.2.0", "explodecomputer/genetics.binaRies")

# Available outcomes
ao <- available_outcomes()
#---------------------------------------------------------------------#
#                            Read exposure                             #----
#---------------------------------------------------------------------#

# Read exposure
GrimAge <- read.table("GrimAge_EUR_summary_statistics.txt", header = T)

#---------------------------------------------------------------------#
#                            Outcomes                                  #----
#---------------------------------------------------------------------#

# Read outcome
# Extract outcome SNPs matching the SNPs in the exposure dataset (output is in the two sample MR package format)
out_func_UKB <- function(file, name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read.table(paste("your_folder/", file,".txt", sep = ""), 
                            header = T) %>% format_data(snps = GrimAge$rsID, 
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

colorectal_ca <- out_func_UKB("maf_overall_colorectal_cancer_imputed", "colorectal cancer")


#---------------------------------------------------------------------#
#                            Merge GWAS data                          #----
#---------------------------------------------------------------------#

X <- gwas_merge(GrimAge, colorectal_ca, 
                snp_name_cols = c("rsID", "SNP"), 
                beta_hat_cols = c("Effect", "beta.outcome"), 
                se_cols = c("SE", "se.outcome"), 
                A1_cols = c("A1", "effect_allele.outcome"), 
                A2_cols = c("A2", "other_allele.outcome"))

#---------------------------------------------------------------------#
#                    Calculate nuisance parameters                    #----
#---------------------------------------------------------------------#

# only > 100,000 variants
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)
head(params$mix_grid)

#---------------------------------------------------------------------#
#                                Clump data                           #----
#---------------------------------------------------------------------#
X$p_value <- 2*pnorm(abs(X$beta_hat_1/X$seb1), lower.tail=FALSE)
X_clump <- X %>% rename(rsid = snp,
                        pval = p_value) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = 0.01,
                     clump_p = 1e-03, #should use larger p value, 1e-03 as used for posteriors, not 5e-08
                    # plink_bin = genetics.binaRies::get_plink_binary(),
                     #bfile = "~/EUR"
                    )
keep_snps <- X_clump$rsid
#---------------------------------------------------------------------#
#                    MR-CAUSE analysis                                #----
#---------------------------------------------------------------------#

# X is unclumped data and variants clumped data
res <- cause(X=X, variants = keep_snps, param_ests = params)
plot(res$sharing)
plot(res$causal)
summary(res, ci_size=0.95)
plot(res)
plot(res, type="data")

png('CAUSE/CAUSE_GrimAge_CRC1.png', res=300, height=2000, width=3500)
plot(res)
dev.off()

png('CAUSE/CAUSE_GrimAge_CRC2.png', res=300, height=2000, width=3500)
plot(res, type="data")
dev.off()
