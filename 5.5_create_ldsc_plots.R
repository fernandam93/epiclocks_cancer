#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 14 October 2021

# This script creates a figure for LD score regression results 
# for epigenetic age acceleration and cancer in international consortiums,
# UK Biobank and FinnGen 

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("your_working_directory") 

library("openxlsx")

#---------------------------------------------------------------------#
#                            LD Score plots                           #----
#---------------------------------------------------------------------#

# Read results stored in an xlsx file (sheet "rg" contains genetic correlations and sheet "p" contains p-values)
rg <- as.matrix(read.xlsx("ldsc_results_clocks_cancers.xlsx", sheet = "rg", colNames = T, rowNames = T))
p<- as.matrix(read.xlsx("ldsc_results_clocks_cancers.xlsx", sheet = "p", colNames = T, rowNames = T))

# Set colours
col2 <- colorRampPalette(c('firebrick3', 'white', 'deepskyblue3'))

# Create plot
corrplot::corrplot(rg, method = "color", p.mat = p, is.corr = F, insig = "label_sig", 
                   na.label = "-", 
                   tl.col = "black", tl.cex = 0.8, tl.srt = 45, 
                   #addCoef.col = "black", 
                   number.cex = 0.8, sig.level = c(0.001, 0.01, 0.05), cl.ratio = 0.2, col = col2(37), 
                   pch.col = 'black', pch.cex = 1.5, #cl.pos = "n", 
                   outline = "black", 
                   mar = c(1,0,0,4), 
                   cl.align.text = "l") 

m <- recordPlot()

m

# Save
save_func_MA <- function(file_name, plot_name)
{
  png(file_name, res=330, height=2500, width=2100)
  print(plot_name)
  dev.off()
}
save_func_MA('ldsc.png', m)
