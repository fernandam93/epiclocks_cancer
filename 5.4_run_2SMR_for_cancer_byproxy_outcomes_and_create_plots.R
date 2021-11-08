#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 12 October 2021

# This script creates plots for two-sample MR results 
# for epigenetic age acceleration and parental cancer in UK Biobank
# This script uses the forestplot function to create colourful plots 

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("your_working_directory") 

# Install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot", "openxlsx")

#---------------------------------------------------------------------#
#                             Read Exposures                          #----
#---------------------------------------------------------------------#

exp_dat <- read.table(file = "exp_data.txt", header = T)

#---------------------------------------------------------------------#
#                           Outcomes                              #----
#---------------------------------------------------------------------#

# Function to read and format outcome data
out_func_proxies <- function(name, file)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- read.table(paste("folder_containing_data/" , file, sep = ""), 
                           header = T) %>% format_data(snps = exp_dat$SNP, 
                                                         type = "outcome", 
                                                         snp_col = "SNP", 
                                                         beta_col = "BETA", 
                                                         se_col = "SE", 
                                                         effect_allele_col = "ALLELE1", 
                                                         other_allele_col = "ALLELE0",
                                                         chr_col = "CHR",
                                                         pos_col = "BP",
                                                         eaf_col = "A1FREQ", 
                                                         pval_col = "P_BOLT_LMM_INF", 
                                                         info_col = "INFO")
  outcome_var$outcome <- name
  return(outcome_var)
}

proxy_breast <- out_func_proxies("Breast cancer", "proxy_breast.txt")
proxy_prostate <- out_func_proxies("Prostate cancer", "proxy_prostate.txt")
proxy_lung <- out_func_proxies("Lung cancer", "proxy_lung.txt")
proxy_colorectal <- out_func_proxies("Colorectal cancer", "proxy_bowel.txt")

# Combine outcomes
proxies_out_dat <- rbind(proxy_colorectal, proxy_breast, proxy_lung, proxy_prostate)

#---------------------------------------------------------------------#
#                            Harmonisation                             #----
#---------------------------------------------------------------------#

# Harmonise exposure and outcome datasets
data_proxies <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = proxies_out_dat
)

#---------------------------------------------------------------------#
#                            Results                                   #----
#---------------------------------------------------------------------#

# Run MR analysis
results_proxies <- mr(data_proxies, method_list = c("mr_ivw", "mr_egger_regression", 
                                                    "mr_weighted_median", "mr_weighted_mode"))

#---------------------------------------------------------------------#
#                           Transform UKB outcomes                      #----
#---------------------------------------------------------------------#

# Create dataframe with data on number of cases and controls for each outcome
cancer <- c("Colorectal cancer", "Lung cancer", "Prostate cancer", "Breast cancer")
ncase <- c(45213, 51073, 31527, 35356)
ncontrol <- c(412429, 404606, 160579, 206992)
UKB_cc <- data.frame(cancer, ncase, ncontrol)

# Function to transform betas from risk difference (BOLT-LMM output) to log odds scale before exponentiation in UK Biobank
UKB_func <- function(outcome) {
  ncase <- UKB_cc$ncase[UKB_cc$cancer==outcome]
  ncontrol <- UKB_cc$ncontrol[UKB_cc$cancer==outcome]
  u <- ncase/(ncase+ncontrol)
}

# Create u column based on number of cases and controls for each outcome
results_proxies$u[results_proxies$outcome=="Colorectal cancer"] <- UKB_func("Colorectal cancer")
results_proxies$u[results_proxies$outcome=="Lung cancer"] <- UKB_func("Lung cancer")
results_proxies$u[results_proxies$outcome=="Prostate cancer"] <- UKB_func("Prostate cancer")
results_proxies$u[results_proxies$outcome=="Breast cancer"] <- UKB_func("Breast cancer")

# Correct UK Biobank betas and SEs using the newly created u column
results_proxies$b <- results_proxies$b/(results_proxies$u*(1-results_proxies$u))
results_proxies$se <- results_proxies$se/(results_proxies$u*(1-results_proxies$u))

# Remove u column from UK Biobank results
results_proxies <- subset(results_proxies, select = -c(u) )

#write.table(results_proxies, "results_proxies.txt", row.names = F)
#results_proxies <- read.table("results_proxies.txt", header = T)
#---------------------------------------------------------------------#
#                       Prepare for visualization                     #----
#---------------------------------------------------------------------#

# Reorder methods as factors
as.factor(results_proxies$method)
results_proxies$method <- factor(results_proxies$method, levels = c("Weighted mode", "Weighted median",    
                                                                    "MR Egger", "Inverse variance weighted"))

# Reporder outcomes as factors
as.factor(results_proxies$outcome)
results_proxies$outcome <- factor(results_proxies$outcome, levels = c("Breast cancer",
                                                                      "Prostate cancer",
                                                                      "Lung cancer",
                                                                      "Colorectal cancer"
                                                                       ))

results_proxies <- arrange(results_proxies, outcome)

#---------------------------------------------------------------------#
#                            Forest function                           #----
#---------------------------------------------------------------------#

myforestplot_3 <- function(df, exp_dataset, exp_dataset2, xlab)
{
  x <- forestplot(
    df = df[which(df$id.exposure==exp_dataset & df$method!="MR Egger"),],
    estimate = b,
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = T,
    colour = method,
    title = exp_dataset2,
    xlab = xlab,
    xlim= c(0.9,1.15)
  ) 
  colours_BP <- c("#6DC6BD", "#2E7EBB", "#08306B")
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.3, 1, 3, 5, 8, 10))
  print(x)
}

#---------------------------------------------------------------------#
#                            Forest plots                             #----
#---------------------------------------------------------------------#

# Create plots using forest plot functions above
p_grim <- myforestplot_3(results_proxies, "GrimAge", "GrimAge", "Odds ratio (95% CI) per year increase in GrimAge acceleration")
#p_grim <- ggarrange(p_grim, legend = "bottom")
#p_grim
p_pheno <- myforestplot_3(results_proxies, "PhenoAge", "PhenoAge", "Odds ratio (95% CI) per year increase in PhenoAge acceleration")
#p_pheno <- ggarrange(p_pheno, legend = "bottom")
#p_pheno
p_hannum <- myforestplot_3(results_proxies, "Hannum", "HannumAge", "Odds ratio (95% CI) per year increase in HannumAge acceleration")
#p_hannum <- ggarrange(p_hannum, legend = "bottom")
#p_hannum
p_ieaa <- myforestplot_3(results_proxies, "IEAA", "Intrinsic HorvathAge", "Odds ratio (95% CI) per year increase in Intrinsic HorvathAge acceleration")
#p_ieaa <- ggarrange(p_ieaa, legend = "bottom")
#p_ieaa

# Combine plots
combine_func <- function(A, B, C, D){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B, C, D, labels=c('A', 'B', 'C', 'D'),
                 ncol = 2, nrow = 2, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}

proxies_all <- combine_func(p_grim, p_pheno, p_hannum, p_ieaa)

# Save
save_func22 <- function(file_name, plot_name)
{
  png(file_name, res=330, height=3000, width=5100)
  print(plot_name)
  dev.off()
}
save_func22('proxies_grim_pheno_hannum_ieaa.png', proxies_all)
