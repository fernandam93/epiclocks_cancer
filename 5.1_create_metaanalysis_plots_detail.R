#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 11 October 2021

# This script creates plots for the meta-analysis of two-sample MR results 
# across UK Biobank, FinnGen and international consortiums
# This script uses the forest.meta function to create plots with details on 
# individual study ORs and CIs, between study heterogeneity, etc

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
#                   Load two-sample MR results                        #----
#---------------------------------------------------------------------#

# Function to load MR results 
load_data_func <- function(file_name)
{
  results <- read.table(file = file_name, header = T)
  results <- split_outcome(results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
}

results_IEU <- load_data_func("results_IEU.txt")
results_UKB <- load_data_func("results_UKB_new.txt")
results_FINNGEN <- load_data_func("results_FINNGEN.txt")

#---------------------------------------------------------------------#
#                           Format MR results                         #----
#---------------------------------------------------------------------#

# Add new column with study name
results_IEU$study[results_IEU$outcome=="Colorectal cancer"] <- "GECCO"
results_IEU$study[results_IEU$outcome=="Lung cancer"] <- "ILCCO"
results_IEU$study[results_IEU$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"] <- "BCAC"
results_IEU$study[results_IEU$outcome=="Ovarian cancer"] <- "OCAC"
results_IEU$study[results_IEU$outcome=="Prostate cancer"] <- "PRACTICAL"
results_UKB$study <- "UK Biobank"
results_FINNGEN$study <- "FinnGen"

# Rename IEU results and filter to only include some of them in the meta-analysis
results_IEU$outcome[results_IEU$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"] <- "Breast cancer"
results_IEU <- results_IEU[which(results_IEU$outcome=='Breast cancer' | results_IEU$outcome=='Ovarian cancer' | results_IEU$outcome=='Prostate cancer' | results_IEU$outcome=='Lung cancer' | results_IEU$outcome=='Colorectal cancer'),]

# Rename UKB results and filter to only include some of them in the meta-analysis
results_UKB$outcome[results_UKB$outcome=="breast cancer"] <- 'Breast cancer'
results_UKB$outcome[results_UKB$outcome=="ovarian cancer"] <- 'Ovarian cancer'
results_UKB$outcome[results_UKB$outcome=="prostate cancer"] <- 'Prostate cancer'
results_UKB$outcome[results_UKB$outcome=="lung cancer unadjusted"] <- 'Lung cancer'
results_UKB$outcome[results_UKB$outcome=="colorectal cancer"] <- 'Colorectal cancer'
results_UKB <- results_UKB[which(results_UKB$outcome=='Breast cancer' | results_UKB$outcome=='Ovarian cancer' | results_UKB$outcome=='Prostate cancer' | results_UKB$outcome=='Lung cancer' | results_UKB$outcome=='Colorectal cancer'),]

# Rename FinnGen results and filter to only include some of them in the meta-analysis
results_FINNGEN$outcome[results_FINNGEN$outcome=="breast cancer (excluding cancer in controls)"] <- 'Breast cancer'
results_FINNGEN$outcome[results_FINNGEN$outcome=="ovarian cancer (excluding cancer in controls)"] <- 'Ovarian cancer'
results_FINNGEN$outcome[results_FINNGEN$outcome=="prostate cancer (excluding cancer in controls)"] <- 'Prostate cancer'
results_FINNGEN$outcome[results_FINNGEN$outcome=="lung cancer (excluding cancer in controls)"] <- 'Lung cancer'
results_FINNGEN$outcome[results_FINNGEN$outcome=="colorectal cancer (excluding cancer in controls)"] <- 'Colorectal cancer'
results_FINNGEN <- results_FINNGEN[which(results_FINNGEN$outcome=='Breast cancer' | results_FINNGEN$outcome=='Ovarian cancer' | results_FINNGEN$outcome=='Prostate cancer' | results_FINNGEN$outcome=='Lung cancer' | results_FINNGEN$outcome=='Colorectal cancer'),]

# Combine results
results <- rbind(results_IEU, results_UKB, results_FINNGEN)

# Sort results by outcome
cancers <- c("Breast cancer", "Ovarian cancer", "Prostate cancer", "Lung cancer", "Colorectal cancer")
results<-results[order(match(results$outcome, cancers)),]

#---------------------------------------------------------------------#
#                           Run meta-analysis                         #----
#---------------------------------------------------------------------#

#Function to run a fixed effect meta-analysis
meta_func2 <- function(method_varname, exp_varname, out1, out2="", out3="", out4="", 
                      out5="", out6="", out7="", out8="", out9="", out10="", 
                      out11="", out12="", out13="", out14="", out15="", 
                      out16="", out17="", out18="", out19="", out20="", out21="")
{
  input <- results[which(results$method==method_varname),]
  input <- input[which(input$id.exposure==exp_varname),]
  input <- input[which(input$outcome==out1 | input$outcome==out2 | input$outcome==out3 | 
                         input$outcome==out4 | input$outcome==out5 | input$outcome==out6 | 
                         input$outcome==out7 | input$outcome==out8 | input$outcome==out9 | 
                         input$outcome==out10 | input$outcome==out11 | input$outcome==out12 |
                         input$outcome==out13 | input$outcome==out14 | input$outcome==out15 | 
                         input$outcome==out16 | input$outcome==out17 | input$outcome==out18 | 
                         input$outcome==out19 | input$outcome==out20 | input$outcome==out21),]
  #meta-analysis
  a <- metagen(TE = b, seTE = se, data = input, 
               studlab = paste(study), sm = "OR",
               hakn = FALSE, byvar = c(outcome),
               method.tau="DL", comb.fixed = T, comb.random = F, exclude = id.outcome %in% c("J2M7ze")) #excluding colorectal cancer in UKB from MA
  print(a)
  
  #forest plot
  forest.meta(a, studlab = TRUE, 
              comb.fixed = a$comb.fixed,
              type.study="square",
              squaresize=0.5,
              lty.fixed = 2,
              type.fixed="diamond",
              bylab = "",
              text.fixed = "Total", # write anything 
              text.fixed.w = "Total",
              col.study="black", 
              col.square="black", 
              col.diamond="white", 
              col.diamond.lines="black",
              col.label.right="black",
              col.label.left="black", 
              colgap.right = "0.5cm",
              colgap.forest.left ="2.2cm",
              col.by = "black",
              smlab="Odds ratio [95% CI]", 
              leftcols=c("studlab"),# To remove "logHR" and "seHR" from plot
              leftlabs = c("Study"),
              rightcols=c("w.fixed", "effect", "ci"),
              rightlabs=c("Weight", "OR","[95% CI]"),
              test.overall = F,
              lwd=0.9,
              print.I2 = a$comb.fixed,
              plotwidth="4.5cm",
              print.I2.ci = a$comb.fixed, 
              print.tau2 = F, 
              print.Q = FALSE,
              digits = 2, 
              fontsize = 9,
              overall = FALSE,
              overall.hetstat = FALSE)
  m <- recordPlot()
}

#IVW
IVW_GrimAge <- meta_func2("Inverse variance weighted", "GrimAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                         "Lung cancer", "Colorectal cancer")
IVW_PhenoAge <- meta_func2("Inverse variance weighted", "PhenoAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                          "Lung cancer", "Colorectal cancer")
IVW_Hannum <- meta_func2("Inverse variance weighted", "Hannum", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                        "Lung cancer", "Colorectal cancer")
IVW_IEAA <- meta_func2("Inverse variance weighted", "IEAA", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                      "Lung cancer", "Colorectal cancer")


#Egger
egger_GrimAge <- meta_func2("MR Egger", "GrimAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                           "Lung cancer", "Colorectal cancer")
egger_PhenoAge <- meta_func2("MR Egger", "PhenoAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                            "Lung cancer", "Colorectal cancer")
egger_Hannum <- meta_func2("MR Egger", "Hannum", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                          "Lung cancer", "Colorectal cancer")
egger_IEAA <- meta_func2("MR Egger", "IEAA", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                        "Lung cancer", "Colorectal cancer")


#Mode based
mode_GrimAge <- meta_func2("Weighted mode", "GrimAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                          "Lung cancer", "Colorectal cancer")
mode_PhenoAge <- meta_func2("Weighted mode", "PhenoAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                           "Lung cancer", "Colorectal cancer")
mode_Hannum <- meta_func2("Weighted mode", "Hannum", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                         "Lung cancer", "Colorectal cancer")
mode_IEAA <- meta_func2("Weighted mode", "IEAA", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                       "Lung cancer", "Colorectal cancer")

#median based
median_GrimAge <- meta_func2("Weighted median", "GrimAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                            "Lung cancer", "Colorectal cancer")
median_PhenoAge <- meta_func2("Weighted median", "PhenoAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                             "Lung cancer", "Colorectal cancer")
median_Hannum <- meta_func2("Weighted median", "Hannum", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                           "Lung cancer", "Colorectal cancer")
median_IEAA <- meta_func2("Weighted median", "IEAA", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                         "Lung cancer", "Colorectal cancer")

#---------------------------------------------------------------------#
#                     Save meta-analysis plots                        #----
#---------------------------------------------------------------------#

save_func_MA <- function(file_name, plot_name)
{
  png(file_name, res=330, height=3000, width=2000)
  print(plot_name)
  dev.off()
}

#save 
save_func_MA('MA_IVW_GrimAge.png', IVW_GrimAge)
save_func_MA('MA_IVW_PhenoAge.png', IVW_PhenoAge)
save_func_MA('MA_IVW_Hannum.png', IVW_Hannum)
save_func_MA('MA_IVW_IEAA.png', IVW_IEAA)

save_func_MA('MA_egger_GrimAge.png', egger_GrimAge)
save_func_MA('MA_egger_PhenoAge.png', egger_PhenoAge)
save_func_MA('MA_egger_Hannum.png', egger_Hannum)
save_func_MA('MA_egger_IEAA.png', egger_IEAA)

save_func_MA('MA_median_GrimAge.png', median_GrimAge)
save_func_MA('MA_median_PhenoAge.png', median_PhenoAge)
save_func_MA('MA_median_Hannum.png', median_Hannum)
save_func_MA('MA_median_IEAA.png', median_IEAA)

save_func_MA('MA_mode_GrimAge.png', mode_GrimAge)
save_func_MA('MA_mode_PhenoAge.png', mode_PhenoAge)
save_func_MA('MA_mode_Hannum.png', mode_Hannum)
save_func_MA('MA_mode_IEAA.png', mode_IEAA)

