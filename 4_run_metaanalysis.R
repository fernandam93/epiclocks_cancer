#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 11 October 2021

# This script runs a meta-analysis of two-sample MR results 
# across UK Biobank, FinnGen and international consortiums
# This script also applies an FDR correction to main meta-analysed IVW results

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

#---------------------------------------------------------------------#
#                           Run meta-analysis                         #----
#---------------------------------------------------------------------#

# Function to run a fixed effect meta-analysis 
meta_func <- function(method_varname, exp_varname, out1, out2="", out3="", out4="", 
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

  #extract values from meta output
  TE.tibble <- as_tibble(a$TE.fixed.w)
  se.tibble <- as_tibble(a$seTE.fixed.w)
  p.tibble <- as_tibble(a$pval.fixed.w)
  bylevs.tibble <- as_tibble(a$bylevs)
  #combine tibbles and change column names
  tibble <- cbind(TE.tibble, se.tibble, p.tibble, bylevs.tibble)
  colnames(tibble) <- c("b", "se", "pval", "outcome")
  #add columns for exposure and method
  tibble$exposure <- exp_varname
  tibble$method <- method_varname
  tibble
  
}

#IVW
IVW_GrimAge <- meta_func("Inverse variance weighted", "GrimAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                         "Lung cancer", "Colorectal cancer")
IVW_PhenoAge <- meta_func("Inverse variance weighted", "PhenoAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                          "Lung cancer", "Colorectal cancer")
IVW_Hannum <- meta_func("Inverse variance weighted", "Hannum", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                        "Lung cancer", "Colorectal cancer")
IVW_IEAA <- meta_func("Inverse variance weighted", "IEAA", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                      "Lung cancer", "Colorectal cancer")


#Egger
egger_GrimAge <- meta_func("MR Egger", "GrimAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                         "Lung cancer", "Colorectal cancer")
egger_PhenoAge <- meta_func("MR Egger", "PhenoAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                          "Lung cancer", "Colorectal cancer")
egger_Hannum <- meta_func("MR Egger", "Hannum", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                        "Lung cancer", "Colorectal cancer")
egger_IEAA <- meta_func("MR Egger", "IEAA", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                      "Lung cancer", "Colorectal cancer")


#Mode based
mode_GrimAge <- meta_func("Weighted mode", "GrimAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                           "Lung cancer", "Colorectal cancer")
mode_PhenoAge <- meta_func("Weighted mode", "PhenoAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                            "Lung cancer", "Colorectal cancer")
mode_Hannum <- meta_func("Weighted mode", "Hannum", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                          "Lung cancer", "Colorectal cancer")
mode_IEAA <- meta_func("Weighted mode", "IEAA", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                        "Lung cancer", "Colorectal cancer")

#median based
median_GrimAge <- meta_func("Weighted median", "GrimAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                          "Lung cancer", "Colorectal cancer")
median_PhenoAge <- meta_func("Weighted median", "PhenoAge", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                           "Lung cancer", "Colorectal cancer")
median_Hannum <- meta_func("Weighted median", "Hannum", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                         "Lung cancer", "Colorectal cancer")
median_IEAA <- meta_func("Weighted median", "IEAA", "Breast cancer", "Ovarian cancer", "Prostate cancer", 
                       "Lung cancer", "Colorectal cancer")

# COMBINE METHODS
GrimAge <- rbind(IVW_GrimAge, egger_GrimAge, median_GrimAge, mode_GrimAge)
PhenoAge <- rbind(IVW_PhenoAge, egger_PhenoAge, median_PhenoAge, mode_PhenoAge)
Hannum <- rbind(IVW_Hannum, egger_Hannum, median_Hannum, mode_Hannum)
IEAA <- rbind(IVW_IEAA, egger_IEAA, median_IEAA, mode_IEAA)

##################################################################
#               Multiple testing correction
##################################################################

# Function to apply multiple testing correction to main meta-analysed IVW results
multiple_testing_func <- function(results_IVW){
  results_IVW$p.fdr<- p.adjust(results_IVW$pval, method = "fdr")
  results_IVW$p.bon<- p.adjust(results_IVW$pval, method = "bonferroni")
  results_IVW$p.hoch <- p.adjust(results_IVW$pval, method = "hochberg")
  results_IVW$p.sig <- ifelse(results_IVW$pval < .05, "*", "")
  results_IVW$p.fdr.sig <- ifelse(results_IVW$p.fdr < .05, "*", "")
  results_IVW$p.bon.sig <- ifelse(results_IVW$p.bon < .05, "*", "")
  results_IVW$p.hoch.sig <- ifelse(results_IVW$p.hoch < .05, "*", "")
  return(results_IVW)
}

IVW_GrimAge <- multiple_testing_func(IVW_GrimAge) 
IVW_PhenoAge <- multiple_testing_func(IVW_PhenoAge)
IVW_Hannum <- multiple_testing_func(IVW_Hannum)
IVW_IEAA <- multiple_testing_func(IVW_IEAA)


# Add FDR p-values
GrimAge <- merge(GrimAge, IVW_GrimAge[, c("p.fdr", "method", "outcome")], by=c("method", "outcome"), all.x = T)
PhenoAge <- merge(PhenoAge, IVW_PhenoAge[, c("p.fdr", "method", "outcome")], by=c("method", "outcome"), all.x = T)
Hannum <- merge(Hannum, IVW_Hannum[, c("p.fdr", "method", "outcome")], by=c("method", "outcome"), all.x = T)
IEAA <- merge(IEAA, IVW_IEAA[, c("p.fdr", "method", "outcome")], by=c("method", "outcome"), all.x = T)

# Combine datasets and save meta-analysis results
MA_fdr_corrected <- rbind(GrimAge, PhenoAge, Hannum, IEAA)
write.xlsx(MA_fdr_corrected, "MA_clocks_incFDR.xlsx")

