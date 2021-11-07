#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 11 October 2021

# This script creates a plot for the meta-analysis of two-sample MR results 
# across UK Biobank, FinnGen and international consortiums
# This script uses the forestplot function to create colourful plots 

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("your_working_directory") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("meta", "metafor", "openxlsx", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggforestplot", "cowplot")

#---------------------------------------------------------------------#
#                   Load meta-analysis results                        #----
#---------------------------------------------------------------------#

results <- read.xlsx("MA_clocks_incFDR.xlsx")
GrimAge <- results[which(results$exposure=="GrimAge"),]
PhenoAge <- results[which(results$exposure=="PhenoAge"),]
Hannum <- results[which(results$exposure=="Hannum"),]
IEAA <- results[which(results$exposure=="IEAA"),]

#---------------------------------------------------------------------#
#                 Format meta-analysis results                        #----
#---------------------------------------------------------------------#

cancers <- c('Breast cancer', 'Ovarian cancer', 'Prostate cancer', 'Lung cancer', 'Colorectal cancer')

as.factor(GrimAge$outcome)
GrimAge$outcome <- factor(GrimAge$outcome, levels = cancers)
GrimAge <- arrange(GrimAge, outcome)

as.factor(PhenoAge$outcome)
PhenoAge$outcome <- factor(PhenoAge$outcome, levels = cancers)
PhenoAge <- arrange(PhenoAge, outcome)

as.factor(Hannum$outcome)
Hannum$outcome <- factor(Hannum$outcome, levels = cancers)
Hannum <- arrange(Hannum, outcome)

as.factor(IEAA$outcome)
IEAA$outcome <- factor(IEAA$outcome, levels = cancers)
IEAA <- arrange(IEAA, outcome)

#reorder methods
as.factor(GrimAge$method)
GrimAge$method <- factor(GrimAge$method, levels = c("Weighted mode", "Weighted median",
                                                            "MR Egger", "Inverse variance weighted"))
#reorder methods
as.factor(PhenoAge$method)
PhenoAge$method <- factor(PhenoAge$method, levels = c("Weighted mode", "Weighted median",
                                                        "MR Egger", "Inverse variance weighted"))
#reorder methods
as.factor(Hannum$method)
Hannum$method <- factor(Hannum$method, levels = c("Weighted mode", "Weighted median",
                                                        "MR Egger", "Inverse variance weighted"))
#reorder methods
as.factor(IEAA$method)
IEAA$method <- factor(IEAA$method, levels = c("Weighted mode", "Weighted median",
                                                                "MR Egger", "Inverse variance weighted"))

#---------------------------------------------------------------------#
#                            Plot Function                            #----
#---------------------------------------------------------------------#

myforestplot <- function(df, title, xlab)
{
  x <- forestplot(
    df = df[which(df$method!="MR Egger"),],
    estimate = b,    # b and se have already been multiplied by 10 in the MA stage
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = TRUE,
    colour = method,
    title = title,
    xlab = xlab,
    xlim= c(0.88,1.22)
  ) 
  colours <- c("#6DC6BD", "#2E7EBB", "#08306B")
  #colours <- c("#3FBC73FF", "#238A8DFF", "#3A548CFF", "black")
  x <- x + scale_color_manual(values=colours)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  print(x)
}

#---------------------------------------------------------------------#
#                            Forest plots                             #----
#---------------------------------------------------------------------#

#main analyses 
GrimAge_MA <- myforestplot(GrimAge, "GrimAge", "Odds ratio (95% CI) per 1 year increase in GrimAge acceleration")
PhenoAge_MA <- myforestplot(PhenoAge, "PhenoAge", "Odds ratio (95% CI) per 1 year increase in PhenoAge acceleration")
Hannum_MA <- myforestplot(Hannum, "HannumAge", "Odds ratio (95% CI) per 1 year increase in HannumAge acceleration")
IEAA_MA <- myforestplot(IEAA, "Intrinsic HorvathAge", "Odds ratio (95% CI) per 1 year increase in Intrinsic HorvathAge acceleration")


#in one plot
comb_func <- function(A, B, C, D){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B, C, D, labels=c('A', 'B', 'C', 'D'),
                 ncol = 2, nrow = 2, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}

all <- comb_func(GrimAge_MA, PhenoAge_MA, Hannum_MA, IEAA_MA)

#---------------------------------------------------------------------#
#                              Save plot                              #----
#---------------------------------------------------------------------#

save_func2 <- function(file_name, plot_name)
{
  png(file_name, res=330, height=3000, width=5100)
  print(plot_name)
  dev.off()
}

#save 
save_func2('MA_all.png', all)

