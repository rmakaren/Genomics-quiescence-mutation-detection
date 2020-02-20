# This script analysis NGS data from targeted resequencing 1000 colonies after Day1, 2 months, 3 months of quiescence
# of PCR products of 9 genes that undergo positive selection during quiescence
# all the function used in this script are saved in a separate file called "functions"

#set working directory
setwd("/home/rostyslav/Desktop/R_analysis/NGS_sequencing/6_subcultures")

#load the libraries
library(ggplot2) # plotting the data
library(grid)
library(gridExtra)
library(plyr)
library(dplyr)
library(gsubfn)

# MULTIPLE SAMPLE ANALYSIS
#--------------------------

##########################
## DATA PREPROCESSING
##########################

#----------------------------------------------------------
##1. uploading the files
#----------------------------------------------------------
filesToProcess <- dir(pattern = ".txt$") # identify files to read in the path
listOfFiles <- lapply(filesToProcess, function(x) read.table(x, header = T, sep ="\t", na.strings=c("","NA"), stringsAsFactors = F))# Iteract over each of those file names with lapply
names(listOfFiles) <- filesToProcess # assign the names of the files to the dataframes
listOfFiles <- Map(cbind, listOfFiles, Sample = names(listOfFiles)) # add the column with dataframe name to the dataframes

#-------------------------------------------
##2. removing unnessesary factors variables
#-------------------------------------------

# apply function percentage transform to the each of the Varfreq column
listOfFiles <- lapply(listOfFiles, function(x) {x$VarFreq <- sapply(x$VarFreq, percentage_transform)
x})
# transform all the dataframes names from last column into character
listOfFiles <- lapply(listOfFiles, function(x) {x[20] <- lapply(x[20], as.character) 
x})
# transform all the dataframes names from last column into character
listOfFiles <- lapply(listOfFiles, function(x) {x$Ref <- lapply(x$Ref, as.character)
x})# transform all the dataframes names from last column into character
listOfFiles <- lapply(listOfFiles, function(x) {x$VarAllele <- lapply(x$VarAllele, as.character)
x})

########################
## DATA ANALYSIS
########################

#---------------------------------
##3. allele candidates filtering
#---------------------------------

# keep the variants that present at frequency > 0.1%
listOfFiles <- lapply(listOfFiles, subset, VarFreq > 0.0008) 

# round the value of frequency to 0.001 (0.1%, 3rd number after dot)
listOfFiles <- lapply(listOfFiles, function(x){
  x$VarFreq <- round(x$VarFreq, 3)
  return(x)
})

# keep the variants that comes from coverage > 5000 reads/nucleotide
listOfFiles <- lapply(listOfFiles, subset, Reads1 + Reads2 > 5000) 
#list2env(listOfFiles ,.GlobalEnv) 

# sort the alleles by mutation type SNPs, Indels, Deletions
listOfFiles <- lapply(listOfFiles, mutation_type) 

# keep the variants that have the proper base quality
listOfFiles <- lapply(listOfFiles, base_quality)

# keep the SNPs that are present in the equvalent proportion on the both pair-end reads
listOfFiles <- lapply(listOfFiles, filter_reads_proportion)

# keep the variants that have passed the Fisher`s test
listOfFiles <- lapply(listOfFiles, subset, Pvalue < 0.05)

#dropping the indels of size +/-1 b.p. below frequency 0.5% base pair bacause there is no way to distunguish them from the sequencing errors 
listOfFiles <- lapply(listOfFiles, drop_low_indels)

# keep the indels that are present in the equvalent proportion on the both pair-end reads
#listOfFiles <- lapply(listOfFiles, filter_reads_proportion_indels)

#---------------------------------
##4. allele candidates annotation
#--------------------------------

# assigning the gene names
listOfFiles <- lapply(listOfFiles, gene_annotation)

#marking the strain direction of the gene on the DNA
listOfFiles <- lapply(listOfFiles, strain_direction)

# dropping the variants outside the cDNA of the 9 genes
listOfFiles <- lapply(listOfFiles, subset, gene != "no")

# calculting the position of the mutation on the cDNA of the gene
listOfFiles <- lapply(listOfFiles, cDNA_annotation)

#
listOfFiles <- lapply(listOfFiles, function(x) {x[24] <- lapply(x[24], as.integer)
x})

#annotating the mutations
listOfFiles <- lapply(listOfFiles, mutation_annotation)

# removing unneccecary sequences
listOfFiles <- lapply(listOfFiles, subset, strain_direction != "no")

# adding the indexes that order genes in following direction: sgf73, win1, wis1, sty1, mkh1, pek1, pmk1, pmc1, tif452
listOfFiles <- lapply(listOfFiles, gene_order)

# adding the column with time in quiescence
listOfFiles <- lapply(listOfFiles, time_in_quiescence)

#assigning the names of the samples
#tmp <- sample_names(tmp)

# sorting the data based on the gene order, sample, time in quiescence and position on the cDNA
listOfFiles <- lapply(listOfFiles, sorting)

#5. here I need to unlist the data frames 
list2env(listOfFiles ,.GlobalEnv) 
#list2env(listOfFiles_2 ,.GlobalEnv)

#-----------------------------------------------------------------------------
#updating data after manual check of bam file - adding the complex mutations
#-----------------------------------------------------------------------------

wt1_1.txt <- Complex_subculture_1(wt1_1.txt)
wt1_2.txt <- Complex_subculture_1(wt1_2.txt)

wt3_1.txt <- Complex_subculture_3(wt3_1.txt)
wt3_2.txt <- Complex_subculture_3(wt3_2.txt)

wt6_1.txt <- Complex_subculture_6(wt6_1.txt)
wt6_2.txt <- Complex_subculture_6(wt6_2.txt)

#---------------------------------

wt1 <- rbind(wt1_1.txt, wt1_2.txt)
wt2 <- rbind(wt2_1.txt, wt2_2.txt)
wt3 <- rbind(wt3_1.txt, wt3_2.txt)
wt4 <- rbind(wt4_1.txt, wt4_2.txt)
wt5 <- rbind(wt5_1.txt, wt5_2.txt)
wt6 <- rbind(wt6_1.txt, wt6_2.txt)

# keeping only duplicate variants
all_together <- do.call("rbind", lapply(ls("pattern" = "wt.$"),get))

#calculating the range and mean between the 2 sequencing samples and saving into single dataframe
tmp_list <- lapply(ls(pattern = "wt.$"),get)
names(tmp_list) <- ls(pattern = "wt.$")
tmp_list <- lapply(tmp_list, range_and_mean)
tmp_list <- lapply(tmp_list, rest_of_population)
list2env(tmp_list ,.GlobalEnv) 

# calculating the rest of the population
listOfFiles_2 <- lapply(listOfFiles, rest_of_population)

#calculating the total number of the mutations
all_together <- do.call("rbind", lapply(ls("pattern" = "wt.$"),get))
all_together <- list_to_vectors(all_together)
all_together <- mutation_annotation_range(all_together)
all_together <- sorting(all_together)

#raw sequences without calculating the range and mean
raw_sequences <- lapply(ls(pattern = "\\>.txt"),get)
names(raw_sequences) <- ls(pattern = "\\>.txt")
raw_sequences <- lapply(raw_sequences, rest_of_population)
list2env(raw_sequences ,.GlobalEnv) 

#removing the controls
rm(D1_1.txt, D1_2.txt)
raw_sequences <- do.call("rbind", lapply(ls("pattern" = "\\>.txt"),get))
raw_sequences <- list_to_vectors(raw_sequences)
raw_sequences <- sorting(raw_sequences)

#saving the data
#write.table(raw_sequences, file = "two_sequences_samples.txt", sep="\t", quote=T, row.names=F, col.names=T)

# merging everythinginto single dataframe
#rm(list=ls(pattern="\\>.txt")) # removing the old dataframes
#all_together <- do.call("rbind", lapply(ls(pattern = "wt.$"),get))

#test <- split.data.frame(tmp, tmp$Sample)

#saving the final results as an output file
all_together <- list_to_vectors(all_together)
write.table(all_together, file = "/home/rostyslav/Desktop/R_analysis/NGS_sequencing/6_subcultures_filtrated.txt", sep="\t", quote=T, row.names=F, col.names=T)
#write.table(raw_sequences, file = "6_subcultures_filtrated_2_datasets.txt", sep="\t", quote=T, row.names=F, col.names=T)

# subsetting SNPs
rm(list=ls(pattern="\\>.txt")) # removing the old dataframes
all_together <- do.call("list", lapply(ls(pattern = "wt.$"),get))
ls(pattern = "wt.$")

######################################
#Plotting the variants
#####################################

test <- lapply(ls(pattern = "wt.$"),get)
names(test) <- ls(pattern = "wt.$")

#subsetting SNPs
list_SNP <- lapply(test, subset, mutation_type =="SNP")
for (i in names(list_SNP)){
  names(list_SNP) <- paste(names(list_SNP), "_SNP", sep = "")
  return(names(list_SNP))
}

#subsetting indels
list_indel <- lapply(test, subset, mutation_type !="SNP" & mutation_type !="Complex")
for (i in names(list_indel)){
  names(list_indel) <- paste(names(list_indel), "_indel", sep = "")
  return(names(list_indel))
}

##6. plotting raw variants

#plotting the SNPs
p_SNP <- lapply(list_SNP, plot_histogram)

#plotting the indels
p_indel <- lapply(list_SNP, plot_histogram)

##--------------------------------------------------
#saving the plots to tiff files
##--------------------------------------------------

output_tiff(p_SNP)
output_tiff(p_indel)

#--------------------------------------------------
#Statistical anaylis of the number of mutations
#--------------------------------------------------

unique_variants_per_month <- function (x){
  tmp <- ddply(x, .(Position, Sample), nrow)
  x <- merge(tmp, x, by=c("Position", "Sample"))
  x <- x[!(duplicated(x[, 1:2], fromLast=T) & !duplicated(x[, 1:2])),]
  return(x)
}

all_together <- do.call("rbind", lapply(ls("pattern" = "wt.*"),get))
all_together <- list_to_vectors(all_together)
all_together <- sorting(all_together)
all_together <- sample_names(all_together)
all_together <- unique_variants_per_month(all_together)

mutation_per_gene <- function(x){
  tmp <- x %>%
    select(gene, mutation_type)
  tmp <- as.data.frame(table(tmp))
  tmp2 <- reshape(tmp, direction="wide", timevar="mutation_type",idvar="gene")
  tmp2 <- gene_order(tmp2)
  tmp2 <- tmp2[order(tmp2$gene_order),]
  tmp2$Indels <-  c(tmp2$Freq.Insertion + tmp2$Freq.Deletion)
  tmp2 <- tmp2 %>%
    select(gene, Indels, Freq.SNP)
  tmp2$SNPs <- tmp2$Freq.SNP
  tmp2$Freq.SNP <- NULL
  rownames(tmp2) <- seq(length=nrow(tmp2))
  return(tmp2)
}

remove_hot_spots_homopolymers <- function(x){
  a <- which(grepl(616208, x$Position) & grepl("+T", x$VarAllele))
  b <- which(grepl(2122889, x$Position) & grepl("+A", x$VarAllele))
  c <- which(grepl(5093449, x$Position) & grepl("+A", x$VarAllele))
  x<- x[-c(a, b, c), , drop = FALSE]
  rownames(x) <- NULL
  return(x)
}

all_together_2 <- remove_hot_spots_homopolymers(all_together)
remove_hot_spots_all <- function(x){
  a <- which(grepl(616208, x$Position) & grepl("+T", x$VarAllele))
  b <- which(grepl(2122889, x$Position) & grepl("+A", x$VarAllele))
  c <- which(grepl(5093449, x$Position) & grepl("+A", x$VarAllele))
  d <- which(grepl(5091712, x$Position) & grepl("+TATCCTTCACGTCGTTC", x$VarAllele))
  e <- which(grepl(5090833, x$Position) & grepl("+CTACTACATCCTC", x$VarAllele))
  x<- x[-c(a, b, c, d, e), , drop = FALSE]
  rownames(x) <- NULL
  return(x)
}
all_together_3 <- remove_hot_spots_all(all_together)

temp1 <- mutation_per_gene(all_together)
temp2 <- mutation_per_gene(all_together_2)
temp3 <- mutation_per_gene(all_together_3)

sum(temp1$Indels)
sum(temp1$SNPs)

sum(temp2$Indels)
sum(temp2$SNPs)

sum(temp3$Indels)
sum(temp3$SNPs)
