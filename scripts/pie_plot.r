# this script builds plots to visualize the data from NGS sequencing

#set the working directory
setwd("/Volumes/@Home/R/Nature\ com/")

#load the libraries
#install.packages("devtools")
#install.packages("Cairo")
library(Cairo) #package to save the data in different file format (.pdf)
library(ggplot2)

##--------------------------------------------------
## DATA PREPROCESSING
##--------------------------------------------------

#reading the dataframe in the path
setwd("/Volumes/@Home/R/Nature\ com/4_cultures/")
df1 <- read.table("final_results_4cultures.txt", header =  T, stringsAsFactors=FALSE)
setwd("/Volumes/@Home/R/Nature\ com/6_subcultures/")
#df2 <- read.table("6subcultures_filtrated.txt", header =  T, stringsAsFactors=FALSE)
df2 <- read.table("6_subcultures_filtrated.txt", header =  T, stringsAsFactors=FALSE)

###########################################
#pie plot plotting
##########################################

cairo_pdf(file = "wt1_1.txt.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_0$VarFreq, labels = NA, col = culture_0$colors, init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "wt1_2.txt", cex.main = 1.6)
legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()