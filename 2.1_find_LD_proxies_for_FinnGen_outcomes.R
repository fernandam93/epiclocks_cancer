####################################################################################
#                         EPIGENETIC CLOCKS AND MULTIPLE CANCERS
####################################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 17 August 2021

# This script finds LD-proxies for cancer outcomes in FinnGen, which can then be
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
#                 Read exposure and outcome datasets                  #----
#---------------------------------------------------------------------#

# Read exposure
exp_dat <- read.table("exp_data.txt", header = T)

# Function to read outcome
out_func_FINNGEN <- function(file, name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read.table(paste("finngen_", file,".txt", sep = ""), 
                            header = T) %>% format_data(snps = exp_dat$SNP, 
                                                        type = "outcome", 
                                                        snp_col = "rsids", 
                                                        beta_col = "beta", 
                                                        se_col = "sebeta", 
                                                        effect_allele_col = "alt", 
                                                        other_allele_col = "ref",
                                                        chr_col = "chrom",
                                                        pos_col = "pos",
                                                        eaf_col = "maf", 
                                                        pval_col = "pval")
  outcome_var$outcome <- name
  return(outcome_var)
}

# Read FinnGen outcome (all finngen datasets are missing the same SNPs, 
# so we can find proxies for only one of the outcomes)
lung_cancer_exallc <- out_func_FINNGEN("lung_exallc", "lung cancer (excluding cancer in controls)")

#---------------------------------------------------------------------#
#              Identify SNPs that need proxies                        #----
#---------------------------------------------------------------------#

# Function to find list of snps in exposure dataset that are missing from the outcome dataset
find_missing_SNP <- function(out_dat) {
                                        snps_need_proxy <- subset(exp_dat, !(exp_dat$SNP %in% out_dat$SNP))
}

s <- find_missing_SNP(lung_cancer_exallc) 

count(s) #check how many snps are in list
s$SNP #see list of snps
s[1,1] #see snp missing n#1
s[2,1] #see snp missing n#2
s[3,1] #see snp missing n#3
s[4,1] #see snp missing n#4

#---------------------------------------------------------------------#
#                    Find proxies for these snps                      #----
#---------------------------------------------------------------------#

# Function to find LD proxy using LDLINK
find_LD_proxy <- function(snps_need_proxy) {
                              proxy <- (LDproxy(snps_need_proxy[1,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[c(1,4),] 
                              proxy$original <- snps_need_proxy[1,1]
                              proxy2 <- (LDproxy(snps_need_proxy[2,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[1:2,] 
                              proxy2$original <- snps_need_proxy[2,1]
                              proxy3 <- (LDproxy(snps_need_proxy[3,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[1:2,] 
                              proxy3$original <- snps_need_proxy[3,1]
                              proxy4 <- (LDproxy(snps_need_proxy[4,1], "EUR", "r2", token = Sys.getenv("LDLINK_TOKEN"), file = F))[1:2,] 
                              proxy4$original <- snps_need_proxy[4,1]
                              proxies <- rbind(proxy, proxy2, proxy3, proxy4)
                              proxies
                              # we could change number of proxies we want to find
}

a <- find_LD_proxy(s)
a[2,1] #see proxy snp (for snp missing n#1)
a[4,1] #see proxy snp (for snp missing n#2)
a[6,1] #see proxy snp (for snp missing n#3)
a[8,1] #see proxy snp (for snp missing n#4)

# We need to make sure the identified proxy SNPs are available in outcome dataset before continuing 
# Here, we used the terminal to do this (e.g., zcat finngen_R5_C3_BREAST_EXALLC.gz | grep rs290794)
# If SNPs aren't available, you need to find the next best proxy for the missing SNP


# List all SNPs included in the outcome dataset, including proxies for those that are missing
# This can then be used to extract data related to these SNPs using grep in the terminal
list_all_snps<- function(out_dat, proxy) {
                                          exp_snps <- out_dat$SNP
                                          proxy_snp <- proxy[c(2,4,6,8),1]
                                          all_snps <- c(exp_snps, proxy_snp)
}
all <- list_all_snps(lung_cancer_exallc, a)

write.table(all, "FINNGEN_SNP_list_inc_proxies.txt", quote = F, sep = " ", col.names = F, row.names = F)
all <- read.table("FINNGEN_SNP_list_inc_proxies.txt", header = F)
#---------------------------------------------------------------------#
#            NOW USE THE TERMINAL TO EXTRACT DATA FOR ALL SNPS        #----
#---------------------------------------------------------------------#
#Example:
#zcat finngen_R5_C3_BRONCHUS_LUNG_EXALLC.gz | grep -w -F -f FINNGEN_SNP_list_inc_proxies.txt -e rsids > finngen_lung_exallc_inc_proxies.txt 

        # Use nano finngen_lung_exallc_inc_proxies.txt to remove hashtag in column names before using the grep command

#---------------------------------------------------------------------#
#                       Back to R -> Format data first                #----
#---------------------------------------------------------------------#
# Function to format data (two-sample MR format)
# Here, it is important to use the format_data function using the SNP list including proxies, not the original SNP list
out_func_FINNGEN_proxies <- function(file, name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read.table(paste("finngen_", file,"_inc_proxies.txt", sep = ""), 
                            header = T) %>% format_data(snps = all, 
                                                        type = "outcome", 
                                                        snp_col = "rsids", 
                                                        beta_col = "beta", 
                                                        se_col = "sebeta", 
                                                        effect_allele_col = "alt", 
                                                        other_allele_col = "ref",
                                                        chr_col = "chrom",
                                                        pos_col = "pos",
                                                        eaf_col = "maf", 
                                                        pval_col = "pval")
  outcome_var$outcome <- name
  return(outcome_var)
}

lung_cancer_exallc_proxies <- out_func_FINNGEN_proxies("lung_exallc", "lung cancer (excluding cancer in controls)")
breast_cancer_exallc_proxies <- out_func_FINNGEN_proxies("breast_exallc", "breast cancer (excluding cancer in controls)")
colorectal_cancer_exallc_proxies <- out_func_FINNGEN_proxies("colorectal_exallc", "colorectal cancer (excluding cancer in controls)")
ovarian_cancer_exallc_proxies <- out_func_FINNGEN_proxies("ovary_exallc", "ovarian cancer (excluding cancer in controls)")
prostate_cancer_exallc_proxies <- out_func_FINNGEN_proxies("prostate_exallc", "prostate cancer (excluding cancer in controls)")


# Find list of SNPs in the proxy dataset that are missing from the outcome dataset, 
# just to make sure that we aren't missing any SNPs before moving on to the next step
find_missing_SNP_proxy <- function(out_dat) {
  snps_need_proxy <- subset(as.data.frame(all), !(all %in% out_dat$SNP))
}

s2 <- find_missing_SNP_proxy(lung_cancer_exallc_proxies) 

#---------------------------------------------------------------------#
#  Modify outcome dataset so that the proxy SNP "rsid" is replaced by the original SNP "rsid" #----
#---------------------------------------------------------------------#
# Function to add columns for the replacement data to the proxy dataset 
replacement_cols_proxy <- function(proxy) {
  proxy <- proxy %>% separate(Coord, c("chr","pos"), sep = "([:])") 
  proxy$chr <- gsub("chr", "", proxy$chr)
  proxy <- proxy %>% separate(Correlated_Alleles, c("A1_original","A1_proxy", "A2_original", "A2_proxy"), sep = "([,=])")
  proxy$original_chr <- c(proxy[1,2], proxy[1,2], proxy[3,2], proxy[3,2], proxy[5,2], proxy[5,2], proxy[7,2], proxy[7,2])
  proxy$original_pos <- c(proxy[1,3], proxy[1,3], proxy[3,3], proxy[3,3], proxy[5,3], proxy[5,3], proxy[7,3], proxy[7,3])
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
}

lung_cancer_exallc_new <- replace_proxy_by_original(b, lung_cancer_exallc_proxies)
breast_cancer_exallc_new <- replace_proxy_by_original(b, breast_cancer_exallc_proxies)
colorectal_cancer_exallc_new <- replace_proxy_by_original(b, colorectal_cancer_exallc_proxies)
ovarian_cancer_exallc_new <- replace_proxy_by_original(b, ovarian_cancer_exallc_proxies)
prostate_cancer_exallc_new <- replace_proxy_by_original(b, prostate_cancer_exallc_proxies)


write.table(breast_cancer_exallc_new, "finngen_breast_exallc_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(ovarian_cancer_exallc_new, "finngen_ovarian_exallc_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(prostate_cancer_exallc_new, "finngen_prostate_exallc_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(lung_cancer_exallc_new, "finngen_lung_exallc_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)
write.table(colorectal_cancer_exallc_new, "finngen_colorectal_exallc_replaced_proxies.txt", sep = " ", row.names = F, col.names = T)



