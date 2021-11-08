#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 18 August 2021

# This script runs CAUSE for GrimAge acceleration and prostate cancer
# in PRACTICAL (data for prostate cancer can be obtained from MR-Base)
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
#read exposure
GrimAge <- read.table("GrimAge_EUR_summary_statistics.txt", header = T)

#---------------------------------------------------------------------#
#                            Outcomes                                  #----
#---------------------------------------------------------------------#

# Option 1: obtain complete summary statistics from MR-Base
# Extract outcome SNPs matching the SNPs in the exposure dataset (output is in the two sample MR package format)
prostate_ca <- extract_outcome_data(
  snps = GrimAge$rsID,
  outcomes = c("ieu-b-85"), proxies = F
)

# Save
#write.table(prostate_ca, "prostate_ca_complete_summary_stats.txt", row.names = F, col.names = T)

# Option 2: Upload complete summary statistics (preferred option)
prostate_ca <- read.table("prostate_ca_complete_summary_stats.txt", header = T)

#---------------------------------------------------------------------#
#                            Merge GWAS data                          #----
#---------------------------------------------------------------------#

X <- gwas_merge(GrimAge, prostate_ca, 
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

png('CAUSE/CAUSE_GrimAge_PrC1.png', res=300, height=2000, width=3500)
plot(res)
dev.off()

png('CAUSE/CAUSE_GrimAge_PrC2.png', res=300, height=2000, width=3500)
plot(res, type="data")
dev.off()
