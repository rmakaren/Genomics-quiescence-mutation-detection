#load neccecary libraries
library(ggplot2) # plotting the data
library(grid)
library(gridExtra)
library(plyr)
library(dplyr)
require(scales)
library(data.table)
library(stringr)
library(tidyverse)
library(gsubfn)
library(Cairo) #package to save the data in different file format (.pdf)
require(graphics)

######################################
#functions for data preprocessing ####
######################################

#transforming the percentage (factor) from VarScan output into numeric in Varfreq column
percentage_transform <- function(x) (as.numeric(sub("%", "", x))/100)

#function in order to sort the dataframe by genes infollowing order: sgf73, win1, wis1, sty1, mkh1, pek1, pmk1, pmc1, tif452
gene_order <- function(x) {
  x$gene_order <- 
    ifelse((x$gene =="sgf73"), "1",
           ifelse((x$gene =="win1"), "2", 
                  ifelse((x$gene =="wis1"), "3",
                         ifelse((x$gene =="sty1"), "4",
                                ifelse((x$gene =="mkh1"), "5",
                                       ifelse((x$gene =="pek1"), "6", 
                                              ifelse((x$gene =="pmk1"), "7",
                                                     ifelse((x$gene =="pmc1"), "8",
                                                            ifelse((x$gene =="tif452"), "9", "no")))))))))
  return(x)
}

#diverse 2 months and 3 months 
time_in_quiescence <- function(x) {
  x$time_in_quiescence <- 
    ifelse(grepl("2m.txt", x$Sample), "2_months",
           ifelse(grepl("3m.txt", x$Sample), "3_months", "2_months"))
  return(x)
}

#############################################################################
#functions for filtering variants from NGS and other data analysis
#############################################################################

#subsetting based on the base quality of the alternative variant in comparison to the reference. 
# The maximum difference is 7%
base_quality <- function(x){
  x <- subset(x, ((abs(x$Qual1 - x$Qual2) < 0.07*((x$Qual1 + x$Qual2/2))) == TRUE))
  return(x)
}

#this function filtrate the SNPs base on the proportion of the SNPs on both read directions with 20% range
filter_reads_proportion <- function(x) {
  tmp <- subset(x, mutation_type == "SNP")
  tmp1 <- subset(x, mutation_type != "SNP")
  tmp2<- subset(tmp, (tmp$Reads1Plus/tmp$Reads2Plus)/(tmp$Reads1Minus/tmp$Reads2Minus) <1.2 & (tmp$Reads1Plus/tmp$Reads2Plus)/(tmp$Reads1Minus/tmp$Reads2Minus)>0.8)
  x <- rbind(tmp1, tmp2)
  return(x)
}

# #this function filtrate the indels base on the proportion of the indels on both read directions with 20% range
filter_reads_proportion_indels <- function(x) {
  tmp <- subset(x, mutation_type != "SNP")
  tmp1 <- subset(x, mutation_type == "SNP")
  tmp2<- subset(tmp, (tmp$Reads1Plus/tmp$Reads2Plus)/(tmp$Reads1Minus/tmp$Reads2Minus) <1.2 & (tmp$Reads1Plus/tmp$Reads2Plus)/(tmp$Reads1Minus/tmp$Reads2Minus)>0.8)
  x <- rbind(tmp1, tmp2)
  return(x)
}

#dropping the indels of size +/- below frequency 0.5% base pair bacause there is no way to filter them from the sequencing errors 
drop_low_indels <- function(x) {
  x <- subset(x, ((mutation_type == "Insertion" | mutation_type == "Deletion") & nchar(x$VarAllele) != 2) | mutation_type == "SNP" | ((mutation_type == "Insertion" | mutation_type == "Deletion") & VarFreq < 0.005 & nchar(x$VarAllele) != 2) | ((mutation_type == "Insertion" | mutation_type == "Deletion") & VarFreq > 0.005 & nchar(x$VarAllele) == 2))
  return(x)
}

#this function subsets variants in common between and calculates range and mean within mutations between 2 datasets (subcultures 2 months)
range_and_mean <- function(x){
  tmp <- duplo_5(x) #subsetting the duplicated variables within 2 datasets
  tmp2 <- with(tmp, tapply(VarFreq, Position, mean)) # calculating the mean within 2 datasets
  tmp2 <- round(tmp2, 3)
  tmp3 <- duplo(tmp) #subsetting unique variants within the duplicates
  tmp3$V1.x <- NULL # removing unneccecary columns
  tmp3$V1.y <- NULL # removing unneccecary columns
  tmp3$VarFreq <- tmp2 # adding the new frequency value to the dataframe
  
  ranging <- function(x){
    tmp <- subset(x, V1 = 2)
    y <- tmp$VarFreq[seq(1, length(tmp$VarFreq), 2)]
    z <- tmp$VarFreq[seq(2, length(tmp$VarFreq), 2)]
    tmp1 <- paste0(rep("±"), round(abs(y-z)/2, 3)*100, rep("%"))
    return(tmp1)
  }
  
  tmp3$range <- ranging(tmp)
  return(tmp3)
}

######################################
#functions to annotate the variants
######################################

#this function annotate the 9 genes from targeted resequencing data (sgf73, win1, wis1, sty1, mkh1, pek1, pmk1, sgf73, tif452)
gene_annotation <- function(x){
  x$gene <- 
    ifelse((x$Position >=5090441 & x$Position <=5094751), "win1",
           ifelse((x$Position >=206864 & x$Position <=207913), "sty1",
                  ifelse((x$Position >=613498 & x$Position <=616848), "mkh1",
                         ifelse((x$Position >=4313854 & x$Position <=4312763), "pek1",
                                ifelse((x$Position >=730860 & x$Position <=732583), "pmk1",
                                       ifelse((x$Position >=1143967 & x$Position <=1145784), "wis1", 
                                              ifelse((x$Position >=2122751 & x$Position <=2123856), "sgf73",
                                                     ifelse((x$Position >=1134499 & x$Position <=1135377), "tif452",
                                                            ifelse((x$Position >=2019758 & x$Position <=2020488), "ras1",
                                                                   ifelse((x$Position >=2713053 & x$Position <=2716931), "pmc1", "no"))))))))))
  return(x)
}

#cDNA annotation
cDNA_annotation <- function(x){
  x$cDNA <- 
    ifelse((x$gene =="win1" & x$mutation_type == "Insertion" | x$gene =="win1" & x$mutation_type == "Deletion"), abs(x$Position - (5090441-2)), 
           ifelse((x$gene =="sty1"), abs(207913 - x$Position) + 1, 
                  ifelse((x$gene =="mkh1" & x$mutation_type == "SNP"), abs(x$Position - (616848) - 1),
                         ifelse((x$gene =="pek1"), abs(x$Position - (4313854)),
                                ifelse((x$gene =="pmk1"& x$Position > 731026 & x$Position < 731205), abs(x$Position - (730860)),
                                       ifelse((x$gene =="wis1"), abs(x$Position - (1145784)-1), 
                                              ifelse((x$gene =="sgf73" & x$Position <=  2123202), abs(x$Position - (2123856-71)),
                                                     ifelse((x$gene =="tif452"), abs(x$Position - (1134499-1)),
                                                            ifelse((x$gene =="ras1"), abs(x$Position - (2020488)),
                                                                   ifelse((x$gene =="pmc1"), abs((2716931) - x$Position+1),
                                                                          ifelse((x$gene == "sgf73" & x$Position >=  2123202),  abs(x$Position - (2123856)),
                                                                                 ifelse((x$gene == "pmk1" & x$Position > 731479 & x$Position < 731550 & x$mutation_type != "SNP"), abs(x$Position - (730860+100+273+1)),
                                                                                        ifelse((x$gene == "pmk1" & x$Position > 731479 & x$Position < 731550 & x$mutation_type != "SNP"), abs(x$Position - (730860+100+273+2)),
                                                                                               ifelse((x$gene == "pmk1" & x$Position > 731633 & x$Position < 732583 & x$mutation_type != "SNP"), abs(x$Position - (730860+100+273+82-2)),
                                                                                                      ifelse((x$gene == "pmk1" & x$Position > 731479 & x$Position < 731550 & x$mutation_type == "SNP"), abs(x$Position - (730860+100-0)),
                                                                                                             ifelse((x$gene == "pmk1" & x$Position > 731479 & x$Position < 731550 & x$mutation_type == "SNP"), abs(x$Position - (730860+100+273-0)),
                                                                                                                    ifelse((x$gene == "pmk1" & x$Position > 731633 & x$Position < 732583 & x$mutation_type == "SNP"), abs(x$Position - (730860+100+273+82-1)),
                                                                                                                           ifelse((x$gene =="mkh1" & x$mutation_type == "Insertion" | x$gene =="mkh1" & x$mutation_type == "Deletion"), abs(x$Position - (616848)),
                                                                                                                                  ifelse((x$gene =="win1" & x$mutation_type == "SNP"), abs(x$Position - (5090441))+1, NA)))))))))))))))))))
  return(x)
}



#this function annotate the mutations by mutation type
mutation_type <- function(x) {
  x$mutation_type <- 
    ifelse(grepl("\\+", x$VarAllele), "Insertion", 
         ifelse(grepl("\\-", x$VarAllele), "Deletion", "SNP"))
return(x)
}

#this function marks the strain of DNA where gene is located
strain_direction <- function(x) {
  x$strain_direction <- 
    ifelse((x$gene =="win1"), "+",
           ifelse((x$gene =="sty1"), "-", 
                  ifelse((x$gene =="mkh1"), "-",
                         ifelse((x$gene =="pek1"), "+",
                                ifelse((x$gene =="pmk1"), "+",
                                       ifelse((x$gene =="wis1"), "-", 
                                              ifelse((x$gene =="sgf73"), "-",
                                                     ifelse((x$gene =="tif452"), "+",
                                                            ifelse((x$gene =="pmc1"), "-", "no")))))))))
  return(x)
}


#this function annotate the mutations by mutation type
mutation_annotation <- function(x){
  x$annotation <- 
    ifelse((x$mutation_type == "SNP"), paste0(x$gene, rep ("-"), x$Ref, x$cDNA, x$VarAllele, rep (" "), x$VarFreq*100, rep("%")),
           ifelse((x$mutation_type == "Insertion"  & nchar(x$VarAllele) <= 3), paste0(x$gene, rep ("-"), x$cDNA, x$VarAllele, rep (" "), x$VarFreq*100, rep("%")),
                  ifelse((x$mutation_type == "Deletion" & nchar(x$VarAllele) <= 3), paste0(x$gene, rep ("-"), x$cDNA, x$VarAllele, rep (" "), x$VarFreq*100, rep("%") ),
                         ifelse((x$mutation_type == "Insertion"  & nchar(x$VarAllele) > 3), paste0(x$gene, rep ("-"), x$cDNA, rep ("+"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "),  x$VarFreq*100, rep("%")),
                                ifelse((x$mutation_type == "Deletion"  & nchar(x$VarAllele) > 3 & x$strain_direction == "+"), paste0(x$gene, rep ("-"), x$cDNA, rep ("-"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "),  x$VarFreq*100, rep("%")),
                                       ifelse((x$mutation_type == "Deletion"  & nchar(x$VarAllele) > 3 & x$strain_direction == "-"), paste0(x$gene, rep ("-"), x$cDNA-nchar(x$VarAllele)+2, rep ("-"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "),  x$VarFreq*100, rep("%")),
                                              paste0(x$gene, rep ("-"), x$cDNA, rep (" "), "CPX", rep (" "), x$VarFreq*100, rep("%"))))))))
  x$annotation <- ifelse((x$strain_direction == "-"), gsubfn(".", list("C" = "G", "G" = "C", "A" = "T", "T" = "A"), x$annotation), gsubfn(".", list("C" = "C", "G" = "G", "A" = "A", "T" = "T"), x$annotation))
  return(x)
}

#this function annotate the mutations by mutation type with mean and range for (2nd experiment with 2 datasets)
mutation_annotation_range <- function(x){
  x$annotation <- 
    ifelse((x$mutation_type == "SNP" & x$range != "±0%"), paste0(x$gene, rep ("-"), x$Ref, x$cDNA, x$VarAllele, rep (" "), round(x$VarFreq*100, 2), x$range),
           ifelse((x$mutation_type == "Insertion"  & nchar(x$VarAllele) <= 3 & x$range != "±0%"), paste0(x$gene, rep ("-"), x$cDNA, x$VarAllele, rep (" "), round(x$VarFreq*100, 2), x$range),
                  ifelse((x$mutation_type == "Deletion" & nchar(x$VarAllele) <= 3 & x$range != "±0%"), paste0(x$gene, rep ("-"), x$cDNA, x$VarAllele, rep (" "), round(x$VarFreq*100, 2), x$range ),
                         ifelse((x$mutation_type == "Insertion"  & nchar(x$VarAllele) > 3 & x$range != "±0%"), paste0(x$gene, rep ("-"), x$cDNA, rep ("+"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "),  round(x$VarFreq*100, 2), x$range),
                                ifelse((x$mutation_type == "Deletion"  & nchar(x$VarAllele) > 3 & x$strain_direction == "+" & x$range != "±0%"), paste0(x$gene, rep ("-"), x$cDNA, rep ("-"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "),  round(x$VarFreq*100, 2), x$range),
                                       ifelse((x$mutation_type == "Deletion"  & nchar(x$VarAllele) > 3 & x$strain_direction == "-" & x$range != "±0%"), paste0(x$gene, rep ("-"), x$cDNA-nchar(x$VarAllele)+2, rep ("-"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "), round(x$VarFreq*100, 2), x$range),
                                              ifelse((x$mutation_type == "SNP" & x$range == "±0%"), paste0(x$gene, rep ("-"), x$Ref, x$cDNA, x$VarAllele, rep (" "), x$VarFreq*100, rep("%")),
                                                     ifelse((x$mutation_type == "Insertion"  & nchar(x$VarAllele) <= 3 & x$range == "±0%"), paste0(x$gene, rep ("-"), x$cDNA, x$VarAllele, rep (" "),  x$VarFreq*100, rep("%")),
                                                            ifelse((x$mutation_type == "Deletion" & nchar(x$VarAllele) <= 3 & x$range == "±0%"), paste0(x$gene, rep ("-"), x$cDNA, x$VarAllele, rep (" "),  x$VarFreq*100, rep("%")),
                                                                   ifelse((x$mutation_type == "Insertion"  & nchar(x$VarAllele) > 3 & x$range == "±0%"), paste0(x$gene, rep ("-"), x$cDNA, rep ("+"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "),   x$VarFreq*100, rep("%")),
                                                                          ifelse((x$mutation_type == "Deletion"  & nchar(x$VarAllele) > 3 & x$strain_direction != "+" & x$range == "±0%"), paste0(x$gene, rep ("-"), x$cDNA, rep ("-"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "),   x$VarFreq*100, rep("%")),
                                                                                 ifelse((x$mutation_type == "Deletion"  & nchar(x$VarAllele) > 3 & x$strain_direction != "-" & x$range == "±0%"), paste0(x$gene, rep ("-"), x$cDNA, rep ("-"), as.character(nchar(x$VarAllele)-1), rep ("bp"), rep (" "),   x$VarFreq*100, rep("%")),
                                                                                        paste0(x$gene, rep ("-"), x$cDNA, rep (" "), "CPX", rep (" "),  x$VarFreq*100, rep("%"))))))))))))))
  x$annotation <- ifelse((x$strain_direction == "-"), gsubfn(".", list("C" = "G", "G" = "C", "A" = "T", "T" = "A"), x$annotation), gsubfn(".", list("C" = "C", "G" = "G", "A" = "A", "T" = "T"), x$annotation))
  return(x)
}
##############################
#additional minor filters
##############################

#function to remove the rows by index
removeRows <- function(rowNum, data) {
  newData <- data[-rowNum, , drop = FALSE]
  rownames(newData) <- NULL
  newData
}

#diverse 2 months and 3 months 
time_in_quiescence <- function(x) {
  x$time_in_quiescence <- 
    ifelse(grepl("2m.txt", x$Sample), "2_months",
           ifelse(grepl("3m.txt", x$Sample), "3_months", "2_months"))
  return(x)
}

#excluding the variants if consensus is not identical to the reference
consensus_filter <- function(x){
  x$Ref <- as.character(x$Ref)
  x$Cons <- as.character(x$Cons)
  x<- subset(x, x$Ref == x$Cons)
  return(x)
}

#reverse the letters
reverse_the_letters <- function(x){
  x <- gsubfn(".", list("C" = "G", "G" = "C", "A" = "T", "T" = "A"), x)
  return(x)
}

#function to reverse the string (for mutation on the - strain for example)
revString = function(string, index = 1:nchar(string)){
  paste(rev(unlist(strsplit(string, NULL)))[index], collapse = "")
}

strReverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

strReverse(c("abc", "Statistics"))

function(x){
  gsub('.{1}$', '', x)
}

# sorting the data based on the gene order, sample, time in quiescence and position on the cDNA
sorting <- function(x){
  x <- x[order(x$Sample, x$time_in_quiescence, x$gene_order, x$cDNA),]
  rownames(x) <- seq(length=nrow(x))
  return(x)
}

#------------------------------------
#several functions to index duplicates
#------------------------------------

#sorting duplicated variants
duplo <- function(x){
  tmp <- ddply(x, .(Position), nrow) #indexing duplicates and not duplicates betweeen two sequencing samples by Position
  tmp2 <- merge(tmp, x, by=c("Position")) # adding indexes to dataframes
  tmp3 <- tmp2[!(!duplicated(tmp2[, 1:2], fromLast=T) & !duplicated(tmp2[, 1:2])),] #keeping the duplicates
  x <- subset(tmp3, Sample!=tmp3$Sample[1]) # keeping duplicates from one subsample
  return(x)
}

#sorting duplicated variants by month
duplo_1 <- function(x){
  tmp <- ddply(x, .(Position), nrow) #indexing duplicates and not duplicates betweeen two sequencing samples by Position
  tmp2 <- merge(tmp, x, by=c("Position")) # adding indexes to dataframes
  tmp3 <- tmp2[!(!duplicated(tmp2[, 1:2], fromLast=T) & !duplicated(tmp2[, 1:2])),] #keeping the duplicates
  x <- subset(tmp3, Sample!=tmp3$Sample[1]) # keeping duplicates from one subsample
  return(x)
}

# duplo_2 is identical to duplo_1 but uptakes the second variable from the Sample column
duplo_2 <- function(x){
  tmp <- ddply(x, .(Position), nrow) #indexing duplicates and not duplicates betweeen two sequencing samples by Position
  tmp2 <- merge(tmp, x, by=c("Position")) # adding indexes to dataframes
  tmp3 <- tmp2[!(!duplicated(tmp2[, 1:2], fromLast=T) & !duplicated(tmp2[, 1:2])),] #keeping the duplicates
  x <- subset(tmp3, Sample!=tmp3$Sample[2]) # keeping duplicates from one subsample
  return(x)
}

#indexing duplicates without removing them
duplo_3 <- function(x){
  tmp <- ddply(x, .(Position), nrow) #indexing duplicates and not duplicates betweeen two sequencing samples by Position
  x <- merge(tmp, x, by=c("Position")) # adding indexes to dataframes
  return(x)
}

#marking the duplicates and unique
duplo_4 <- function(x){
  tmp <- ddply(x, .(Position), nrow) #indexing duplicates and not duplicates betweeen two sequencing samples by Position
  x <- merge(tmp, x, by=c("Position")) # adding indexes to dataframes
  return(x)
}

#grabbing the unique variants among the 2 months and 3 months of quiescence
unique_variants_per_month <- function (x){
  tmp <- ddply(x, .(Position, Sample), nrow)
  tmp2 <- merge(tmp, x, by=c("Position", "Sample"))
  tmp3 <- test[!(duplicated(test[, 1:2], fromLast=T) & !duplicated(test[, 1:2])),]
  return(x)
}

#sorting unique  variants
not_duplo <- function(x){
  tmp <- ddply(x, .(Position), nrow) #indexing duplicates and not duplicates betweeen two sequencing samples by Position
  tmp2 <- merge(tmp, x, by=c("Position")) # adding indexes to dataframes
  tmp3 <- tmp2[!(!duplicated(tmp2[, 1:2], fromLast=T) & !duplicated(tmp2[, 1:2])),] #keeping the duplicates
  x <- subset(tmp3, Sample!=tmp3$Sample[1]) # keeping duplicates from one subsample
  return(x)
}

#indexing the datatable without the subsetting
duplo_5 <- function(x){
  tmp <- ddply(x, .(Position), nrow) #indexing duplicates and not duplicates betweeen two sequencing samples by Position
  tmp2 <- merge(tmp, x, by=c("Position")) # adding indexes to dataframes
  x <- tmp2[!(!duplicated(tmp2[, 1:2], fromLast=T) & !duplicated(tmp2[, 1:2])),] #keeping the duplicates
  return(x)
}

#rename the samples
sample_names <- function(x) {
  x$Sample <- 
    ifelse((x$Sample == "wt0_2m.txt"), "culture_0",
           ifelse((x$Sample=="wt0_3m.txt"), "culture_0",
                  ifelse((x$Sample=="wt1_2m.txt"), "culture_1",
                         ifelse((x$Sample=="wt1_3m.txt"), "culture_1",
                                ifelse((x$Sample=="wt2_2m.txt"), "culture_2",
                                       ifelse((x$Sample=="wt2_3m.txt"), "culture_2",
                                              ifelse((x$Sample=="wt3_2m.txt"), "culture_3",
                                                     ifelse((x$Sample=="wt3_3m.txt"), "culture_3",
                                                            ifelse((x$Sample=="wt1_1.txt"), "subculture_1",
                                                                   ifelse((x$Sample=="wt1_2.txt"), "subculture_1",
                                                                          ifelse((x$Sample=="wt2_1.txt"), "subculture_2",
                                                                                 ifelse((x$Sample=="wt2_2.txt"), "subculture_2",
                                                                                        ifelse((x$Sample=="wt3_1.txt"), "subculture_3",
                                                                                              ifelse((x$Sample=="wt3_2.txt"), "subculture_3",
                                                                                                    ifelse((x$Sample=="wt4_1.txt"), "subculture_4",
                                                                                                            ifelse((x$Sample=="wt4_2.txt"), "subculture_4",
                                                                                                                  ifelse((x$Sample=="wt5_1.txt"), "subculture_5",
                                                                                                                          ifelse((x$Sample=="wt5_2.txt"), "subculture_5",
                                                                                                                                ifelse((x$Sample=="wt6_1.txt"), "subculture_6",
                                                                                                                                        ifelse((x$Sample=="wt6_2.txt"), "subculture_6", "none"))))))))))))))))))))
  return(x)
}

#annotion for building plots
mutation <- function(x) {
  x$mutation <- 
    ifelse((x$mutation_type == "SNP" & x$strain_direction == "+"), paste0(x$Ref, ">", x$VarAllele),
           ifelse((x$mutation_type == "Insertion"& x$strain_direction == "+"), paste0(x$VarAllele),
                  ifelse((x$mutation_type == "Deletion"& x$strain_direction == "+"), paste0(x$VarAllele),
                         ifelse((x$mutation_type == "SNP" & x$strain_direction == "-"), paste0(reverse_the_letters(x$Ref), ">", reverse_the_letters(x$VarAllele)),
                                ifelse((x$mutation_type == "Insertion" & x$strain_direction == "-"), paste0( "+", gsub('.{1}$', '', reverse_the_letters(strReverse(x$VarAllele)))),
                                       ifelse((x$mutation_type == "Deletion" & x$strain_direction == "-"), paste0( "-", gsub('.{1}$', '', reverse_the_letters(strReverse(x$VarAllele)))),
                                              ifelse((x$mutation_type == "Complex"), paste0(x$VarAllele), "no")))))))
  return(x)
}

#functions to calculate the coverage of the sequencing data: avarege, minimum and maximum
average <- function(x) {
  sum(x$Reads1 + x$Reads2)/length(x$Reads1)
}

min_coverage <- function(x){
  min(x$Reads1 + x$Reads2)
}

max_coverage <- function(x){
  max(x$Reads1 + x$Reads2)
}

#multiple plot building
p <- lapply(tmp, function(tmp) ggplot(data = tmp, aes(x = VarFreq)) + 
              list(
                geom_histogram(binwidth = 0.001),
                expand_limits(y=c(0, 50)),
                coord_cartesian(xlim = c(0, 1), ylim=c(0,130)),
                annotation_logticks(),
                theme_bw(),
                NULL
              ))

#venn diagrams
venn.plot <- draw.pairwise.venn(length(wt2_1.txt$VarFreq), 
                                length(wt2_2.txt$VarFreq), 
                                length(wt2_test$VarFreq), 
                                category = c("wt4_1", "wt4_2"),
                                fill = c("blue", "red"),
                                lty = "blank",
                                cex = 2,
                                cat.cex = 0,
                                cat.pos = c(285, 105),
                                cat.dist = 0.09,
                                cat.just = list(c(-1, -1), c(1, 1)),
                                ext.pos = 30,
                                ext.dist = -0.05,
                                ext.length = 0.85,
                                ext.line.lwd = 0,
                                ext.line.lty = "dashed"
);
grid.draw(venn.plot);


#multiple file processing and saving the data
lapply(files, function(x) {
  t <- read.table(x, header=T, sep ="\t", na.strings=c("","NA")) # load file
  # apply function
  out <- function(t)
    # write to file
    write.table(out, "path/to/output", sep="\t", quote=F, row.names=F, col.names=T)
})

#---------------------
#small function to repair the list into characters
list_to_vectors <- function(x){
  x$Ref <- as.character(x$Ref)
  x$VarAllele <- as.character(x$VarAllele)
  return(x)
}

#adding empt row to a dataframes
addding_row_to_df <- function(x){
  tmp<- rep(NA, ncol(x))
  x <- rbind(x, tmp)
  return(x)
}

####################
#plotting the distribution of the variant frequencies

plot_histogram <- function(x) ggplot(data = x, aes(VarFreq)) + 
  list(
    geom_histogram(binwidth = 0.001, show.legend = T),
    coord_cartesian(xlim = c(0, 0.1), ylim=c(0,130)),
    scale_x_continuous(labels = percent, breaks = seq(0, 0.1, 0.01), minor_breaks = seq(0, 0.1, 0.001), name = "Variant frequency"),
    scale_y_continuous(breaks = seq(0, 200, 10), minor_breaks = seq(0, 200, 1), name = "Number of variants"),
    theme(axis.title=element_text(family = "Times New Roman", face="bold", size=14,color="black"), axis.text=element_text(family = "Times New Roman", face="bold", size=14,color="black")) ,
    NULL
  )

# saving the output into separate tiff file
output_tiff <- function(x) {
  for (i in 1:length(x)) {
    file_name = paste(names(x[i]), ".tiff", sep="")
    tiff(file_name, width = 1000, height = 1000, units="px", res=200)
    print(x[[i]])
    dev.off()
  }
}

#------------------------------

# rest of the population function
rest_pop <- function (x){
  test_test<- (1 -sum(x$VarFreq, na.rm = TRUE))
  return(test_test)
}

rest_of_population <- function(x){
  tmp <- rep(NA, ncol(x)) #creating a new empty row 
  x <- rbind(x, tmp) #binding two rows
  tmp2 <- sum(x$VarFreq, na.rm = TRUE)
  x$VarFreq[ nrow(x) ] <- tmp2
  x$Sample[ nrow(x) ] <- x$Sample[1]
  x$annotation[ nrow(x) ] <- paste0("total", "=", (tmp2*100), "%")
  x$time_in_quiescence[ nrow(x) ] <- x$time_in_quiescence[1]
  x$gene_order[ nrow(x) ] <-99
  return(x)
}

#################
#Specific function to create the Complex mutations

Complex_subculture_1 <- function(x){
  a <- which(grepl(5090804, x$Position) & grepl("A", x$VarAllele))
  b <- which(grepl(5090806, x$Position) & grepl("A", x$VarAllele))
  c <- which(grepl(5090808, x$Position) & grepl("A", x$VarAllele))
  d <- which(grepl(5090796, x$Position) & grepl("-TCGTGCT$", x$VarAllele))
  x <- removeRows(c(a, b, c), x)
  x[d, which( colnames(x)=="VarAllele" )] = "-TCGTGCT*T>A*T>A*C>A"
  x[d, which( colnames(x)=="mutation_type" )] = "Complex"
  x[d, which( colnames(x)=="annotation" )] = "win1-357-Complex-1.2%"
  return(x)
}

Complex_subculture_3 <- function(x){
  a <- which(grepl(1144503, x$Position) & grepl("A", x$VarAllele))
  d <- which(grepl(1144504, x$Position) & grepl("A", x$VarAllele))
  x <- removeRows(c(a), x)
  x[d, which( colnames(x)=="Ref" )] = "GG"
  x[d, which( colnames(x)=="VarAllele" )] = "AA"
  x[d, which( colnames(x)=="mutation_type" )] = "Complex"
  x[d, which( colnames(x)=="annotation" )] = "wis1-C1281T 8.5%"
  return(x)
}

Complex_subculture_6 <- function(x){
  a <- which(grepl(5094224, x$Position) & grepl("-A", x$VarAllele))
  d <- which(grepl(5094226, x$Position) & grepl("C", x$VarAllele))
  x<- x[-a, , drop = FALSE]
  rownames(x) <- NULL
  x[d, which( colnames(x)=="VarAllele" )] = "T>C*-C"
  x[d, which( colnames(x)=="mutation_type" )] = "Complex"
  x[d, which( colnames(x)=="annotation" )] = "win1-3785-Complex-1%"
  return(x)
}

Complex_mkh1 <- function(x){
  a <- which(grepl(614965, x$Position) & grepl("A", x$VarAllele))
  d <- which(grepl(614964, x$Position) & grepl("C", x$VarAllele))
  x<- x[-a, , drop = FALSE]
  rownames(x) <- NULL
  x[d, which( colnames(x)=="Ref" )] = "GA"
  x[d, which( colnames(x)=="VarAllele" )] = "AC"
  x[d, which( colnames(x)=="mutation_type" )] = "Complex"
  x[d, which( colnames(x)=="annotation" )] = "mkh-CT1884TG 0.2%"
  return(x)
}

Complex_pmk1 <- function(x){
  a <- which(grepl(731642, x$Position) & grepl("T", x$VarAllele))
  x<- x[-a, , drop = FALSE]
  rownames(x) <- NULL
  return(x)
}

remove_hot_spots <- function(x){
  a <- which(grepl(616208, x$Position) & grepl("+T", x$VarAllele))
  b <- which(grepl(2122889, x$Position) & grepl("+A", x$VarAllele))
  c <- which(grepl(5093449, x$Position) & grepl("+A", x$VarAllele))
 d <- which(grepl(5091712, x$Position) & grepl("+TATCCTTCACGTCGTTC", x$VarAllele))
 e <- which(grepl(5090833, x$Position) & grepl("+CTACTACATCCTC", x$VarAllele))
#  x<- x[-c(a, b, c), , drop = FALSE]
  x<- x[-c(a, b, c, d, e), , drop = FALSE]
  rownames(x) <- NULL
  return(x)
}

freq_to_per <- function(x){
  x$persentage <- paste0(x$VarFreq*100, rep("%"))
  return(x)
}
