####################################################################################
#                         EPIGENETIC CLOCKS AND MULTIPLE CANCERS
####################################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 17 August 2021

# This script finds LD-proxies for cancer subtypes in PRACTICAL, which can then be
# used in two-sample MR analyses
# Output in two-sample MR format 

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("your_working_directory") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools", "LDlinkR")

#---------------------------------------------------------------------#
#             Read exposure and outcome datasets                      #----
#---------------------------------------------------------------------#

# Read exposure
exp_dat <- read.table("exp_data.txt", header = T)

# Function to read outcome
out_func_practical <- function(var, name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read.table(paste("PRACTICAL_", var, ".txt", sep = ""), 
                            header = T) %>% format_data(snps = exp_dat$SNP, 
                                                        type = "outcome", 
                                                        snp_col = "SNP", 
                                                        beta_col = "Effect", 
                                                        se_col = "StdErr", 
                                                        effect_allele_col = "Allele1", 
                                                        other_allele_col = "Allele2",
                                                        chr_col = "Chr",
                                                        pos_col = "position",
                                                        eaf_col = "Freq1", 
                                                        pval_col = "Pvalue")
  outcome_var$outcome <- name
  return(outcome_var)
}


PRACT_adv <- out_func_practical("adv", "Advanced prostate cancer")
PRACT_age55 <- out_func_practical("age55", "Early onset prostate cancer")
PRACT_caseonly <- out_func_practical("caseonly", "Advanced prostate cancer (cases only)")
PRACT_gleason <- out_func_practical("gleason", "Gleason Score")
PRACT_highvslow <- out_func_practical("highvslow", "High risk prostate cancer (vs low risk)")
PRACT_highvslowint <- out_func_practical("highvslowint", "High risk prostate cancer (vs low and intermediate risk)")

#---------------------------------------------------------------------#
#              Identify SNPs that need proxies                        #----
#---------------------------------------------------------------------#

# Function to find list of snps in exposure dataset that are missing from the outcome dataset
find_missing_SNP <- function(out_dat) {
                                        snps_need_proxy <- subset(exp_dat, !(exp_dat$SNP %in% out_dat$SNP))
}

s <- find_missing_SNP(PRACT_adv) 

count(s) #check how many snps are in list
s$SNP #see list of snps
s[1,1] #see snp missing n#1
s[2,1] #see snp missing n#2
s[3,1] #see snp missing n#3
s[4,1] #see snp missing n#4
s[5,1] #see snp missing n#4

#---------------------------------------------------------------------#
#                    Find proxies for these snps                      #----
#---------------------------------------------------------------------#

# Function to find LD proxy using LDLINK
find_LD_proxy <- function(snps_need_proxy) {
                              proxy <- (LDproxy(snps_need_proxy[1,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[c(1,4),] 
                              proxy$original <- snps_need_proxy[1,1]
                              proxy2 <- (LDproxy(snps_need_proxy[2,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[c(1,3),] 
                              proxy2$original <- snps_need_proxy[2,1]
                              proxy3 <- (LDproxy(snps_need_proxy[3,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[1:2,] 
                              proxy3$original <- snps_need_proxy[3,1]
                              proxy4 <- (LDproxy(snps_need_proxy[4,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[c(1,9),] 
                              proxy4$original <- snps_need_proxy[4,1]
                              proxy5 <- (LDproxy(snps_need_proxy[5,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[1:2,] 
                              proxy5$original <- snps_need_proxy[5,1]
                              proxies <- rbind(proxy, proxy2, proxy3, proxy4, proxy5)
                              proxies
                              # we could change number of proxies we want to find
}

a <- find_LD_proxy(s)
a[2,1] #see proxy snp (for snp missing n#1)
a[4,1] #see proxy snp (for snp missing n#2)
a[6,1] #see proxy snp (for snp missing n#3)
a[8,1] #see proxy snp (for snp missing n#4)
a[10,1] #see proxy snp (for snp missing n#5)

# We need to make sure the identified proxy SNPs are available in outcome dataset before continuing 
# Here, we used the terminal to do this (e.g., zcat finngen_R5_C3_BREAST_EXALLC.gz | grep rs290794)
# If SNPs aren't available, you need to find the next best proxy for the missing SNP

# List all SNPs included in the outcome dataset, including proxies for those that are missing
# This can then be used to extract data related to these SNPs using grep in the terminal
list_all_snps<- function(out_dat, proxy) {
                                          exp_snps <- out_dat$SNP
                                          proxy_snp <- proxy[c(2,4,6,8,10),1]
                                          all_snps <- c(exp_snps, proxy_snp)
}
all <- list_all_snps(PRACT_adv, a)

write.table(all, "PRACTICAL_SNP_list_inc_proxies.txt", quote = F, sep = " ", col.names = F, row.names = F)

#---------------------------------------------------------------------#
#            NOW USE THE TERMINAL TO EXTRACT DATA FOR ALL SNPS        #----
#---------------------------------------------------------------------#
#Example:
#zcat finngen_R5_C3_BRONCHUS_LUNG_EXALLC.gz | grep -w -F -f FINNGEN_SNP_list_inc_proxies.txt -e rsids > finngen_lung_exallc_inc_proxies.txt 

#---------------------------------------------------------------------#
#                       Back to R -> Format data first                #----
#---------------------------------------------------------------------#
# Function to format data (two-sample MR format)
# Here, it is important to use the format_data function using the SNP list including proxies, not the original SNP list
out_func_practical_proxies <- function(var, name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read.table(paste("PRACTICAL_", var, "_inc_proxies.txt", sep = ""), 
                            header = T) %>% format_data(snps = all, 
                                                        type = "outcome", 
                                                        snp_col = "SNP", 
                                                        beta_col = "Effect", 
                                                        se_col = "StdErr", 
                                                        effect_allele_col = "Allele1", 
                                                        other_allele_col = "Allele2",
                                                        chr_col = "Chr",
                                                        pos_col = "position",
                                                        eaf_col = "Freq1", 
                                                        pval_col = "Pvalue")
  outcome_var$outcome <- name
  return(outcome_var)
}

PRACT_adv_proxies <- out_func_practical_proxies("adv", "Advanced prostate cancer")
PRACT_age55_proxies <- out_func_practical_proxies("age55", "Early onset prostate cancer")
PRACT_caseonly_proxies <- out_func_practical_proxies("caseonly", "Advanced prostate cancer (cases only)")
PRACT_gleason_proxies <- out_func_practical_proxies("gleason", "Gleason Score")
PRACT_highvslow_proxies <- out_func_practical_proxies("highvslow", "High risk prostate cancer (vs low risk)")
PRACT_highvslowint_proxies <- out_func_practical_proxies("highvslowint", "High risk prostate cancer (vs low and intermediate risk)")


# Find list of SNPs in the proxy dataset that are missing from the outcome dataset, 
# just to make sure that we aren't missing any SNPs before moving on to the next step
find_missing_SNP_proxy <- function(out_dat) {
  snps_need_proxy <- subset(as.data.frame(all), !(all %in% out_dat$SNP))
}

s2 <- find_missing_SNP_proxy(PRACT_adv_proxies) 

#---------------------------------------------------------------------#
#  Modify outcome dataset so that proxy SNP "rsid" is replaced by original SNP "rsid" #----
#---------------------------------------------------------------------#
# Function to add columns for the replacement data to the proxy dataset 
replacement_cols_proxy <- function(proxy) {
  proxy <- proxy %>% separate(Coord, c("chr","pos"), sep = "([:])") 
  proxy$chr <- gsub("chr", "", proxy$chr)
  proxy <- proxy %>% separate(Correlated_Alleles, c("A1_original","A1_proxy", "A2_original", "A2_proxy"), sep = "([,=])")
  proxy$original_chr <- c(proxy[1,2], proxy[1,2], proxy[3,2], proxy[3,2], proxy[5,2], proxy[5,2], proxy[7,2], proxy[7,2], proxy[9,2], proxy[9,2])
  proxy$original_pos <- c(proxy[1,3], proxy[1,3], proxy[3,3], proxy[3,3], proxy[5,3], proxy[5,3], proxy[7,3], proxy[7,3], proxy[9,3], proxy[9,3])
  proxy
  
}

b <- replacement_cols_proxy(a)

# Function to replace proxy SNP details by original SNP details
replace_proxy_by_original <- function(proxy, out_dat_proxies) {
  #proxy 1
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[2,1] & out_dat_proxies$effect_allele.outcome == proxy[2,10]] <- proxy[2,9] #A1 effect allele
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[2,1] & out_dat_proxies$effect_allele.outcome == proxy[2,12]] <- proxy[2,11] #A1 effect allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[2,1] & out_dat_proxies$other_allele.outcome == proxy[2,10]] <- proxy[2,9] #A2 other allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[2,1] & out_dat_proxies$other_allele.outcome == proxy[2,12]] <- proxy[2,11] #A2 other allele
  out_dat_proxies$chr.outcome[out_dat_proxies$SNP==proxy[2,1]] <- proxy[2,16] #change chr
  out_dat_proxies$pos.outcome[out_dat_proxies$SNP==proxy[2,1]] <- proxy[2,17] #change pos
  out_dat_proxies$SNP[out_dat_proxies$SNP==proxy[2,1]] <- proxy[2,15] #change rsid
  #proxy 2
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[4,1] & out_dat_proxies$effect_allele.outcome == proxy[4,10]] <- proxy[4,9] #A1 effect allele
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[4,1] & out_dat_proxies$effect_allele.outcome == proxy[4,12]] <- proxy[4,11] #A1 effect allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[4,1] & out_dat_proxies$other_allele.outcome == proxy[4,10]] <- proxy[4,9] #A2 other allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[4,1] & out_dat_proxies$other_allele.outcome == proxy[4,12]] <- proxy[4,11] #A2 other allele
  out_dat_proxies$chr.outcome[out_dat_proxies$SNP==proxy[4,1]] <- proxy[4,16] #change chr
  out_dat_proxies$pos.outcome[out_dat_proxies$SNP==proxy[4,1]] <- proxy[4,17] #change pos
  out_dat_proxies$SNP[out_dat_proxies$SNP==proxy[4,1]] <- proxy[4,15] #change rsid
  #proxy 3
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[6,1] & out_dat_proxies$effect_allele.outcome == proxy[6,10]] <- proxy[6,9] #A1 effect allele
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[6,1] & out_dat_proxies$effect_allele.outcome == proxy[6,12]] <- proxy[6,11] #A1 effect allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[6,1] & out_dat_proxies$other_allele.outcome == proxy[6,10]] <- proxy[6,9] #A2 other allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[6,1] & out_dat_proxies$other_allele.outcome == proxy[6,12]] <- proxy[6,11] #A2 other allele
  out_dat_proxies$chr.outcome[out_dat_proxies$SNP==proxy[6,1]] <- proxy[6,16] #change chr
  out_dat_proxies$pos.outcome[out_dat_proxies$SNP==proxy[6,1]] <- proxy[6,17] #change pos
  out_dat_proxies$SNP[out_dat_proxies$SNP==proxy[6,1]] <- proxy[6,15] #change rsid
  #proxy 4
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[8,1] & out_dat_proxies$effect_allele.outcome == proxy[8,10]] <- proxy[8,9] #A1 effect allele
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[8,1] & out_dat_proxies$effect_allele.outcome == proxy[8,12]] <- proxy[8,11] #A1 effect allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[8,1] & out_dat_proxies$other_allele.outcome == proxy[8,10]] <- proxy[8,9] #A2 other allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[8,1] & out_dat_proxies$other_allele.outcome == proxy[8,12]] <- proxy[8,11] #A2 other allele
  out_dat_proxies$chr.outcome[out_dat_proxies$SNP==proxy[8,1]] <- proxy[8,16] #change chr
  out_dat_proxies$pos.outcome[out_dat_proxies$SNP==proxy[8,1]] <- proxy[8,17] #change pos
  out_dat_proxies$SNP[out_dat_proxies$SNP==proxy[8,1]] <- proxy[8,15] #change rsid
  out_dat_proxies
  #proxy 5
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[10,1] & out_dat_proxies$effect_allele.outcome == proxy[10,10]] <- proxy[10,9] #A1 effect allele
  out_dat_proxies$effect_allele.outcome[out_dat_proxies$SNP==proxy[10,1] & out_dat_proxies$effect_allele.outcome == proxy[10,12]] <- proxy[10,11] #A1 effect allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[10,1] & out_dat_proxies$other_allele.outcome == proxy[10,10]] <- proxy[10,9] #A2 other allele
  out_dat_proxies$other_allele.outcome[out_dat_proxies$SNP==proxy[10,1] & out_dat_proxies$other_allele.outcome == proxy[10,12]] <- proxy[10,11] #A2 other allele
  out_dat_proxies$chr.outcome[out_dat_proxies$SNP==proxy[10,1]] <- proxy[10,16] #change chr
  out_dat_proxies$pos.outcome[out_dat_proxies$SNP==proxy[10,1]] <- proxy[10,17] #change pos
  out_dat_proxies$SNP[out_dat_proxies$SNP==proxy[10,1]] <- proxy[10,15] #change rsid
  out_dat_proxies
}

PRACT_adv_new <- replace_proxy_by_original(b, PRACT_adv_proxies)
PRACT_age55_new <- replace_proxy_by_original(b, PRACT_age55_proxies)
PRACT_caseonly_new <- replace_proxy_by_original(b, PRACT_caseonly_proxies)
PRACT_gleason_new <- replace_proxy_by_original(b, PRACT_gleason_proxies)
PRACT_highvslow_new <- replace_proxy_by_original(b, PRACT_highvslow_proxies)
PRACT_highvslowint_new <- replace_proxy_by_original(b, PRACT_highvslowint_proxies)



write.table(PRACT_adv_new, "PRACTICAL_adv_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(PRACT_age55_new, "PRACTICAL_age55_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(PRACT_caseonly_new, "PRACTICAL_caseonly_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(PRACT_highvslow_new, "PRACTICAL_highvslow_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(PRACT_highvslowint_new, "PRACTICAL_highvslowint_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(PRACT_gleason_new, "PRACTICAL_gleason_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)

