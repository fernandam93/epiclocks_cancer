
#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 12 November 2021

# This script runs CAUSE for GrimAge acceleration and colorectal cancer in GECCO
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
colorectal_ca <- read.table("GECCO_summary_statistics.txt", header = T) #already in two-sample MR format

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
