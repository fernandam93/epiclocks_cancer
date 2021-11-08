#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 12 October 2021

# This script creates plots for the two-sample MR results 
# for epigenetic age acceleration and cancer subtypes in international consortiums
# This script uses the forestplot function to create colourful plots 

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("your_working_directory") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot")


#---------------------------------------------------------------------#
#                            Read results                             #----
#---------------------------------------------------------------------#

results_IEU <- read.table("results/results/results_IEU.txt", header = T)
results_subtype <- read.table("results/results/results_subtypes.txt", header = T)
results_subtype <- results_subtype[which(results_subtype$outcome!="Gleason Score"),]
results_IEU <- rbind(results_IEU, results_subtype)
#---------------------------------------------------------------------#
#                      Prepare for analyses                           #----
#---------------------------------------------------------------------#

#create groups for IEU outcomes
results_IEU$group <- "Cancers"
results_IEU$group[results_IEU$id.outcome=="ieu-a-1121"| results_IEU$id.outcome=="ieu-a-1122"| results_IEU$id.outcome=="ieu-a-1123"| results_IEU$id.outcome=="ieu-a-1124"| results_IEU$id.outcome=="ieu-a-1125"] <- "Ovarian cancer subtypes" 
results_IEU$group[results_IEU$id.outcome=="ieu-a-1127"| results_IEU$id.outcome=="ieu-a-1128" | results_IEU$outcome=="BRCA1 breast cancer"| results_IEU$outcome=="BRCA2 breast cancer" | results_IEU$outcome=="Luminal B HER2 Negative" | results_IEU$outcome=="Triple Negative" | results_IEU$outcome=="Luminal A" | results_IEU$outcome=="Luminal B" | results_IEU$outcome=="HER2 Enriched"] <- "Breast cancer subtypes" 
results_IEU$group[results_IEU$id.outcome=="ieu-a-967"| results_IEU$id.outcome=="ieu-a-965"] <- "Lung cancer subtypes" 
results_IEU$group[results_IEU$outcome=="BRCA1 ovarian cancer"| results_IEU$outcome=="BRCA2 ovarian cancer"] <- "Ovarian cancer subtypes"
results_IEU$group[results_IEU$outcome=="Early onset prostate cancer"|results_IEU$outcome=="High risk prostate cancer (vs low and intermediate risk)" | results_IEU$outcome=="Advanced prostate cancer" | results_IEU$outcome=="High risk prostate cancer (vs low risk)" | results_IEU$outcome=="Advanced prostate cancer (cases only)"] <- "Prostate cancer subtypes"
results_IEU$group[results_IEU$outcome=="Colon cancer"|results_IEU$outcome=="Rectal cancer" | results_IEU$outcome=="Male colorectal cancer" | results_IEU$outcome=="Female colorectal cancer" | results_IEU$outcome=="Distal colon cancer" | results_IEU$outcome=="Proximal colon cancer"] <- "Colorectal cancer subtypes"

#split outcome so that id doesnt appear in plot
results_IEU <- split_outcome(results_IEU)

results_IEU$outcome[results_IEU$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"] <- "Breast cancer"
results_IEU$outcome[results_IEU$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"] <- "ER+"
results_IEU$outcome[results_IEU$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"] <- "ER-"

#rename subtypes so that they fit in plot
results_IEU$outcome[results_IEU$outcome=="High grade serous ovarian cancer"] <- "High grade serous"
results_IEU$outcome[results_IEU$outcome=="Low grade serous ovarian cancer"] <- "Low grade serous"
results_IEU$outcome[results_IEU$outcome=="Invasive mucinous ovarian cancer"] <- "Invasive mucinous"
results_IEU$outcome[results_IEU$outcome=="Clear cell ovarian cancer"] <- "Clear cell"
results_IEU$outcome[results_IEU$outcome=="Endometrioid ovarian cancer"] <- "Endometrioid"
results_IEU$outcome[results_IEU$outcome=="Lung adenocarcinoma"] <- "Adenocarcinoma"
results_IEU$outcome[results_IEU$outcome=="Squamous cell lung cancer"] <- "Squamous cell"
results_IEU$outcome[results_IEU$outcome=="High risk prostate cancer (vs low and intermediate risk)"] <- "High risk (vs low and intermediate risk)"
results_IEU$outcome[results_IEU$outcome=="High risk prostate cancer (vs low risk)"] <- "High risk (vs low risk)"
results_IEU$outcome[results_IEU$outcome=="Advanced prostate cancer (cases only)"] <- "Advanced (cases only)"
results_IEU$outcome[results_IEU$outcome=="Advanced prostate cancer"] <- "Advanced"
results_IEU$outcome[results_IEU$outcome=="Early onset prostate cancer"] <- "Early onset"

# Reorder groups as factors
as.factor(results_IEU$group)
results_IEU$group <- factor(results_IEU$group, levels = c("Cancers", "Breast cancer subtypes", "Ovarian cancer subtypes", "Prostate cancer subtypes", "Lung cancer subtypes", "Colorectal cancer subtypes"))

# Reorder outcomes as factors
as.factor(results_IEU$outcome)
results_IEU$outcome <- factor(results_IEU$outcome, levels = c("Breast cancer" ,"ER+" , "ER-"  ,
                                                              "Triple Negative", "Luminal B HER2 Negative", 
                                                               "HER2 Enriched", "Luminal A", "Luminal B", "BRCA1 breast cancer", 
                                                              "BRCA2 breast cancer", 
                                                              "Ovarian cancer", "High grade serous", 
                                                              "Low grade serous" ,"Invasive mucinous", 
                                                              "Clear cell", "Endometrioid",
                                                              "BRCA1 ovarian cancer", "BRCA2 ovarian cancer",
                                                              "Prostate cancer",
                                                              "Advanced",
                                                              "Advanced (cases only)",
                                                              "Early onset",
                                                              "High risk (vs low risk)",
                                                              "High risk (vs low and intermediate risk)",
                                                              "Adenocarcinoma" , "Lung cancer" , 
                                                              "Squamous cell",
                                                                      "Colorectal cancer",
                                                              "Colon cancer",
                                                              "Proximal colon cancer",
                                                              "Distal colon cancer",
                                                              "Rectal cancer",
                                                              "Male colorectal cancer",
                                                              "Female colorectal cancer"))

results_IEU <- arrange(results_IEU, outcome)

# Reorder methods as factors
as.factor(results_IEU$method)
results_IEU$method <- factor(results_IEU$method, levels = c("Weighted mode", "Weighted median", "MR Egger",  "Inverse variance weighted"))



#---------------------------------------------------------------------#
#                            Forest function                           #----
#---------------------------------------------------------------------#

# Function for IVW results only
myforestplot <- function(df, exp_dataset, clock_label, xlab)
{
  x <- forestplot(
    df = df[which(df$id.exposure==exp_dataset & df$method=="Inverse variance weighted"),],
    estimate = b,
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = T,
    #colour = method,
    title = clock_label,
    xlab = xlab,
    #xlim= c(0.5,1.5)
  ) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  #colours_BP <- c("#FE7A25", "#FFCE41", "#00bfff", "#98C409")
  #x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.3, 1, 3, 5, 8, 10))
  print(x)
}

# Function for sensitivity analyses and IVW MR
myforestplot_3 <- function(df, exp_dataset, clock_label, xlab)
{
  x <- forestplot(
    df = df[which(df$id.exposure==exp_dataset & df$group!="Cancers" & df$method!="MR Egger" ),],
    estimate = b,
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = T,
    colour = method,
    title = clock_label,
    xlab = xlab,
    xlim= c(0.5,1.5)
  ) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  colours_BP <- c("#6DC6BD", "#2E7EBB", "#08306B")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.50, 0.75, 1.00, 1.25, 1.5, 1.75, 2))
  print(x)
}

# Function for prostate and colorectal cancer subtypes only
myforestplot_4 <- function(df, exp_dataset, clock_label, xlab)
{
  x <- forestplot(
    df = df[which(df$id.exposure==exp_dataset & df$group!="Cancers" & df$group!="Breast cancer subtypes" & df$group!="Lung cancer subtypes" & df$group!="Ovarian cancer subtypes" & df$method!="MR Egger" ),],
    estimate = b,
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = T,
    colour = method,
    title = clock_label,
    xlab = xlab,
    xlim= c(0.5,1.5)
  ) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  colours_BP <- c("#6DC6BD", "#2E7EBB", "#08306B")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.50, 0.75, 1.00, 1.25, 1.5, 1.75, 2))
  print(x)
}

#---------------------------------------------------------------------#
#                            Forest plots                             #----
#---------------------------------------------------------------------#

# Create plots using forest plot functions above
p1 <- myforestplot(results_IEU, "GrimAge", "GrimAge", "Odds ratio (95% CI)")
p2 <- myforestplot(results_IEU, "IEAA", "Intrinsic HorvathAge", "Odds ratio (95% CI)")
p3 <- myforestplot(results_IEU, "PhenoAge", "PhenoAge", "Odds ratio (95% CI)")
p4 <- myforestplot(results_IEU, "Hannum", "HannumAge", "Odds ratio (95% CI)")

p_grim <- myforestplot_3(results_IEU, "GrimAge", "GrimAge", "Odds ratio (95% CI) per year increase in GrimAge acceleration")
p_pheno <- myforestplot_3(results_IEU, "PhenoAge", "PhenoAge", "Odds ratio (95% CI) per year increase in PhenoAge acceleration")
p_hannum <- myforestplot_3(results_IEU, "Hannum", "HannumAge", "Odds ratio (95% CI) per year increase in HannumAge acceleration")
p_ieaa <- myforestplot_3(results_IEU, "IEAA", "Intrinsic HorvathAge", "Odds ratio (95% CI) per year increase in Intrinsic HorvathAge acceleration")

p_grim2 <- myforestplot_4(results_IEU, "GrimAge", "GrimAge", "Odds ratio (95% CI) per year increase in GrimAge acceleration")
p_grim2 <- ggarrange(p_grim2, legend = "bottom")
p_grim2

# Combine plots
combine_func2 <- function(A, B){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B, labels=c('A', 'B'),
                 ncol = 2, nrow = 1, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}

p_grim_pheno <- combine_func2(p_grim, p_pheno)
p_hannum_ieaa <- combine_func2(p_hannum, p_ieaa)

# Save plots
save_func20 <- function(file_name, plot_name)
{
  png(file_name, res=300, height=5000, width=3000)
  print(plot_name)
  dev.off()
}
save_func20('results/plots/subtype_consortia.png', p_grim)

save_func2000 <- function(file_name, plot_name)
{
  png(file_name, res=330, height=3000, width=3000)
  print(plot_name)
  dev.off()
}
save_func2000('results/plots/subtype_grim_prc_crc.png', p_grim2)

save_func21 <- function(file_name, plot_name)
{
  png(file_name, res=300, height=5000, width=5000)
  print(plot_name)
  dev.off()
}
save_func21('results/plots/subtype_grim_pheno.png', p_grim_pheno)
save_func21('results/plots/subtype_hannum_ieaa.png', p_hannum_ieaa)
