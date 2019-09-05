library(dplyr)


#adding manually the variants from WGS
#tmp <- c(207717, "I", "C", 0.002, "A", "culture_0", "SNP", "3_months", "sty1", 197, "0.2%", "-", NA) #sty1 29AP(8W)#
#tmp1 <- c(207840, "I", "G", 0.000, "C", "culture_0", "SNP", "3_months", "sty1", 74, "0%", "-", NA) #sty1 49PP(6R)#
#tmp2 <- c(2714003, "I", "T", 0.000, "-AACCCACAA", "culture_0", "Deletion", "3_months", "pmc1", 2920, "0%", "-", NA) #pmc1 38AP(10W)
#tmp3 <- c(613910, "I", "T", 0.000, "-TATGCTTAGA", "culture_0", "Deletion", "3_months", "mkh1", 2929, "0%", "-", NA) #mkh1 48AA(14W)
#tmp4 <- c(4313731, "II", "T", 0.000, "A", "culture_0", "SNP", "3_months", "pek1", 969, "0%", "-", NA) #pek1 51AA(15W)
#tmp5 <- c(731111, "II", "T", 0.000, "+A", "culture_0", "Insertion", "3_months", "pmk1", 153, "0%", "+", NA) #pmk1 55AA(18W)
#tmp6 <- c(2715686, "I", "C", 0.000, "A", "culture_0", "SNP", "3_months", "pmc1", 1247, "0%", "-", NA) #pmc1 60AP(19W)
#tmp7 <- c(2715686, "II", "A", 0.000, "+T", "culture_0", "Insertion", "3_months", "wis1", NA, "0%", "-", NA) #wis1 B(26W)
#tmp8 <- c(5091712, "I", "A", 0.003, "+TATCCTTCACGTCGTTC", "culture_0", "Insertion", "3_months", "win1", 1273, "0.3%", "+", NA) #win1 E(29W)
#tmp9 <- c(614355, "I", "C", 0.000, "T", "culture_0", "SNP", "3_months", "mkh1", 2494, "0%", "-", NA) #mkh1 F(30W)
#tmp10 <- c(2716606, "I", "A", 0.000, "T", "culture_0", "SNP", "3_months", "pmc1", 326, "0%", "-", NA) #pmc1 7PP(1R)
#tmp11 <- c(2714025, "I", "A", 0.000, "-G", "culture_0", "Deletion", "3_months", "pmc1", 2906, "0%", "-", NA) #pmc1 24PP(2R)
#tmp12 <- c(2716550, "I", "G", 0.000, "A", "culture_0", "SNP", "3_months", "pmc1", 382, "0%", "-", NA) #pmc1 71PP(8R)
#tmp13 <- c(2714491, "I", "G", 0.000, "+AATATTATCACCAGTAAC", "culture_0", "Insertion", "3_months", "pmc1", 2441, "0%", "-", NA) #pmc1 I(33W)
#tmp14 <- c(1134749, "II", "G", 0.002, "T", "culture_0", "SNP", "3_months", "tif452", 251, "0.2%", "+", NA) #tif452 63AP(21W)
#tmp15 <- c(5090833, "I", "A", 0.002, "+CTACTACATCCTC", "culture_0", "Insertion", "3_months", "win1", 394, "0.2%", "+", NA) #win1 76AA(24W)
#tmp16 <- c(731537, "II", "T", 0.004, "G", "wt1_2m", "SNP", "3_months", "pmk1", 305, "0.2%", "+", NA) #pmk1
#tmp15 <- c(5091526, "I", "G", 0.000, "-G", "culture_0", "Insertion", "3_months", "win1", 1086, "0%", "+", NA) #win1 76AA(24W)

####################################################
# real code for all together
####################################################

setwd("/home/rostyslav/Desktop/R project/1st_exp/")
df1 <- read.table("final_results_4cultures.txt", header =  T, stringsAsFactors=FALSE)
setwd("/home/rostyslav/Desktop/R project/2nd_exp/")
#df2 <- read.table("final_6_subcultures_3.txt", header =  T, stringsAsFactors=FALSE)
df2 <- read.table("6_subcultures_filtrated.txt", header =  T, stringsAsFactors=FALSE)
df2$range <- NULL
df2 <- gene_order(df2)
df2 <- sample_names(df2)
freq_to_per <- function(x){
  x$persentage <- paste0(x$VarFreq*100, rep("%"))
  return(x)
}
df1 <- freq_to_per(df1)
df2 <- freq_to_per(df2)
df <- rbind(df1, df2)

#fixing
df_tmp <- split(df2, f = list(df2$Sample, df2$time_in_quiescence))
df_tmp <- lapply(df_tmp, function(x) x[-nrow(x),] )
df2 <- do.call("rbind", df_tmp) #transform the list into the dataframe
df2 <- sorting(df2)

#sort everything again
df <- sorting(df)

#######################################################
#subsetting and manipulationg with win1 
#tempo <- df[!(df$gene=="win1" & df$VarFreq < 0.01),]

minor_win1_allels <- function(x){
  #calculating the total number of the win1 alleles <1%
  tmp <- subset(x, gene == "win1" & VarFreq < 0.01)
  tmp1 <- tmp %>%
    select(Sample, time_in_quiescence, gene)
  tmp2 <- as.data.frame(table(tmp1))
  
  #counting the number of the alleles in 
  tmp3 <- tmp %>%
    group_by(Sample, time_in_quiescence) %>%
    summarise(VarFreq = sum(VarFreq))
  
  #merging 2 dataframes
  results <- merge(tmp2, tmp3, by = c("Sample", "time_in_quiescence"))
  tmp3$annotation <- paste0(results$Freq, rep(" win1 alleles<1% ="), results$VarFreq*100, rep("%"))
  tmp3$gene <- "win1"
  
  #removing the minor win1 alleles from the dataframe
  tempo <- x[!(x$gene=="win1" & x$VarFreq < 0.01),]
  x <- rbind.fill(tempo, tmp3)
  return(x)
}

df <- minor_win1_allels(df)
df <- gene_order(df)
df <- sorting(df)
freq_to_per <- function(x){
  x$persentage <- paste0(x$VarFreq*100, rep("%"))
  return(x)
}
df <- freq_to_per(df)
dt_test <- rest_of_population(df)
tmp <- rest_of_population(df1)
df$gene_order <- as.integer(df$gene_order)

########################################
# assigning the colors, based on the gene and cDNA position
# x - is a dataframe, y is a gene, 
#######################################

df_tmp <- split(df, f = list(df$Sample, df$time_in_quiescence))
rest_of_population <- function(x){
  tmp <- rep(NA, ncol(x)) #creating a newempty row 
  x <- rbind(x, tmp) #binding two rows
  tmp2 <- sum(x$VarFreq, na.rm = TRUE)
  x$VarFreq[ nrow(x) ] <- 1 - tmp2
  x$Sample[ nrow(x) ] <- x$Sample[1]
  x$annotation[ nrow(x) ] <- paste0("total", "=", (tmp2*100), "%")
  x$gene[ nrow(x) ] <- paste0("rest of the population")
  x$time_in_quiescence[ nrow(x) ] <- x$time_in_quiescence[1]
  x$gene_order[ nrow(x) ] <-99
  return(x)
}
df_tmp <- df_tmp[sapply(df_tmp, function(x) dim(x)[1]) > 1] #removing empty lists
df_tmp <- lapply(df_tmp, rest_of_population)
df_tmp <- do.call("rbind", df_tmp) #transform the list into the dataframe
pie_plot_colors <- function(x){
  tmp_list <- split(x, f = x$gene) #splitting dataframe into list of dataframes by genes
  tmp_list <- lapply(tmp_list, duplo_4) #counting duplicates by position of cDNA
  tmp_list <- lapply(tmp_list, function(x) {x[order(x$cDNA),]}) # sorting the data by cDNA postion
  tmp_list_2 <- lapply(tmp_list, function(tmp_list) tmp_list[!duplicated(tmp_list$cDNA),]) #subsetting unique variants from the mutations
  coloring <- function(x) {
    x$colors <- 
      ifelse((x$gene == "sgf73"), hsv(0.85, .8, seq(.75,.9,length.out = sum(x$gene == "sgf73"))), #pink
             ifelse((x$gene  == "win1"), hsv(0.65, .7, seq(.2,1,length.out = sum(x$gene == "win1"))), # darkblue
                    ifelse((x$gene == "wis1"), hsv(0.75, .7, seq(.4,1,length.out = sum(x$gene == "wis1"))), # violet
                           ifelse((x$gene == "sty1"), hsv(0.55, .7, seq(.35,1,length.out = sum(x$gene == "sty1"))), # blue
                                  ifelse((x$gene == "mkh1"), hsv(0.3, .7, seq(.4,1,length.out = sum(x$gene == "mkh1"))), # green
                                         ifelse((x$gene == "pek"), hsv(0.08, 1, seq(.8,1,length.out = sum(x$gene == "pek"))), # orange
                                                ifelse((x$gene== "pmk1"), hsv(0.1, .8, seq(.35,.8,length.out = sum(x$gene == "pmk1"))), # brown
                                                       ifelse((x$gene == "pmc1"), hsv(0.15, .7, seq(.9,1,length.out = sum(x$gene == "pmc1"))), # yellow
                                                              ifelse((x$gene == "tif452"), hsv(0.98, .8, seq(.7,1,length.out = sum(x$gene == "tif452"))), # red
                                                                     "#F0FFFF")))))))))
    return(x)
  } #adding the gamma colors based on the gene type and position of the mutation on cDNA
  tmp_list_3 <- map(tmp_list_2, coloring) #adding the gamma colors based on the gene type and position of the mutation on cDNA
  temporary <- function(x, y) {
    x$colors <- rep(y$colors, times = y$V1)
    return(x)
  } #function to add the same color within the variant duplication
  tmp_list_4 <- map2(tmp_list, tmp_list_3, temporary) #adding the colors to the original dataframes
  tmp_list_4 <- do.call("rbind", tmp_list_4) #transform the list into the dataframe
  tmp_list_4 <- sorting(tmp_list_4) #sorting the finl dataframe
  return(tmp_list_4)
}
df_tmp <- pie_plot_colors(df_tmp)

##########################
#pieplots to plot the data
##########################

df_tmp <- split(df_tmp, f = list(df_tmp$Sample, df_tmp$time_in_quiescence))
list2env(df_tmp ,.GlobalEnv) 


#subsetting SNPs

cairo_pdf(file = "subculture_1.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = subculture_1.2_months$VarFreq, labels = NA, col = subculture_1.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "subculture_1.2_months", cex.main = 1.6)
legend(x = .05, y = .45, subculture_1.2_months$annotation, col = "white",  cex = 1.6, fill = subculture_1.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "subculture_2.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = subculture_2.2_months$VarFreq, labels = NA, col = subculture_2.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "subculture_2.2_months", cex.main = 1.6)
legend(x = .05, y = .45, subculture_2.2_months$annotation, col = "white",  cex = 1.6, fill = subculture_2.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "subculture_3.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = subculture_3.2_months$VarFreq, labels = NA, col = subculture_3.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "subculture_3.2_months", cex.main = 1.6)
legend(x = .05, y = .45, subculture_3.2_months$annotation, col = "white",  cex = 1.6, fill = subculture_3.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "subculture_4.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = subculture_4.2_months$VarFreq, labels = NA, col = subculture_4.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "subculture_4.2_months", cex.main = 1.6)
legend(x = .05, y = .45, subculture_4.2_months$annotation, col = "white",  cex = 1.6, fill = subculture_4.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "subculture_5.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = subculture_5.2_months$VarFreq, labels = NA, col = subculture_5.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "subculture_6.2_months", cex.main = 1.6)
legend(x = .05, y = .45, subculture_5.2_months$annotation, col = "white",  cex = 1.6, fill = subculture_5.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "subculture_6.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = subculture_6.2_months$VarFreq, labels = NA, col = subculture_6.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "subculture_6.2_months", cex.main = 1.6)
legend(x = .05, y = .45, subculture_6.2_months$annotation, col = "white",  cex = 1.6, fill = subculture_6.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "culture_0.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_0.2_months$VarFreq, labels = NA, col = culture_0.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "culture_0.2_months", cex.main = 1.6)
legend(x = .05, y = .45, culture_0.2_months$annotation, col = "white",  cex = 1.6, fill = culture_0.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "culture_0.3_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_0.3_months$VarFreq, labels = NA, col = culture_0.3_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "culture_0.3_months", cex.main = 1.6)
legend(x = .05, y = .45, culture_0.3_months$annotation, col = "white",  cex = 1.6, fill = culture_0.3_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "culture_1.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_1.2_months$VarFreq, labels = NA, col = culture_1.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "culture_1.2_months", cex.main = 1.6)
legend(x = .05, y = .45, culture_1.2_months$annotation, col = "white",  cex = 1.6, fill = culture_1.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "culture_1.3_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_1.3_months$VarFreq, labels = NA, col = culture_1.3_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "culture_1.3_months", cex.main = 1.6)
legend(x = .05, y = .45, culture_1.3_months$annotation, col = "white",  cex = 1.6, fill = culture_1.3_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "culture_2.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_2.2_months$VarFreq, labels = NA, col = culture_2.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "culture_2.2_months", cex.main = 1.6)
legend(x = .05, y = .45, culture_2.2_months$annotation, col = "white",  cex = 1.6, fill = culture_2.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "culture_2.3_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_2.3_months$VarFreq, labels = NA, col = culture_2.3_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "culture_2.3_months", cex.main = 1.6)
legend(x = .05, y = .45, culture_2.3_months$annotation, col = "white",  cex = 1.6, fill = culture_2.3_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "culture_3.2_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_3.2_months$VarFreq, labels = NA, col = culture_3.2_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "culture_3.2_months", cex.main = 1.6)
legend(x = .05, y = .45, culture_3.2_months$annotation, col = "white",  cex = 1.6, fill = culture_3.2_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

cairo_pdf(file = "culture_3.3_months.pdf", width = 11, height =13, family = 'Times New Roman')
pie(x = culture_3.3_months$VarFreq, labels = NA, col = culture_3.3_months$colors,  init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
title(line = -16, main = "culture_3.3_months", cex.main = 1.6)
legend(x = .05, y = .45, culture_3.3_months$annotation, col = "white",  cex = 1.6, fill = culture_3.3_months$colors,  trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
dev.off()

#I give up, I am tired to write such big functions


# saving the output into separate tiff file
output_tiff <- function(x) {
  for (i in 1:length(x)) {
    file_name = paste(names(x[i]), ".pdf", sep="")
    cairo_pdf(file = file_name, width = 11, height =13, family = 'Times New Roman')
    pie(x$VarFreq, labels = NA, col = x$colors, init.angle =  90, radius = .4, border = "darkgrey", lty = "blank", cex = 1, angle = 100)
    title(line = -16, main = "wt1_2.txt", cex.main = 1.6)
    legend(x = .05, y = .45, x$annotation, col = "white",  cex = 1.6, fill = x$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = NA)
    dev.off()
  }
}

output_tiff <- function(x) {
  for (i in 1:length(x)) {
    file_name = paste(names(x[i]), ".tiff", sep="")
    tiff(file_name, width = 1000, height = 1000, units="px", res=200)
    print(x[[i]])
    dev.off()
  }
}


#testing_tmp <- unlist(test_tmp, use.names=FALSE) #transforming the lists into character vectors 
#testing_tmp <- lapply(test_tmp, unlist)
#testing_tmp_2 <- lapply(testing_tmp, unname)
#test_tmp_2 <- map2(df_test, testing_tmp_2, cbind)


#subsetting unique variants within 2 months and 3 months per culture
unique_variants_per_month <- function (x){
  tmp <- ddply(x, .(Position, Sample), nrow)
  x <- merge(tmp, x, by=c("Position", "Sample"))
  x <- x[!(duplicated(x[, 1:2], fromLast=T) & !duplicated(x[, 1:2])),]
  return(x)
}
df1 <- unique_variants_per_month(df1)
df1$V1 <- NULL
df <- rbind(df1, df2)
df <- sorting(df)

write.table(total_2, file = "final_variants_and_table_total.txt", sep="\t", quote=T, row.names=F, col.names=T)


#tmp code
###########
#colors for genes
pink = hsv(0.85, .8, seq(.75,.9,length.out = length(x$cDNA))) # sgf73
darkblue = hsv(0.65, .7, seq(.2,1,length.out = length(x$cDNA))) # win1
violet = hsv(0.75, .7, seq(.4,1,length.out = length(x$cDNA))) # wis1
blue = hsv(0.55, .7, seq(.35,1,length.out = length(x$cDNA))) #sty1
green = hsv(0.3, .7, seq(.4,1,length.out = length(x$cDNA))) # mkh1
orange = hsv(0.08, 1, seq(.8,1,length.out = length(x$cDNA))) #pek#
brown = hsv(0.1, .8, seq(.35,.8,length.out = length(x$cDNA))) # pmk1
yellow = hsv(0.15, .7, seq(.9,1,length.out = length(x$cDNA))) #pmc1
red =  hsv(0.98, .8, seq(.7,1,length.out = length(x$cDNA))) # tif452

tmp <- subset(df, gene == "mkh1") #subsetting the gene to assign the color with 
tmp <- duplo_4(tmp) # marking the number of duplicates the variant is present in the dataset
tmp <- tmp[order(tmp$gene, tmp$cDNA),]# sort just in case
tmp1 <- subset(tmp, !duplicated(x = tmp$Position)) # subsetting the uniques occasions of mutations
tmp2 <- rep(tmp1$colors, times = tmp1$V1)
tmp$colors <- tmp2
tmp1$colors <- hsv(0.3, .7, seq(.4,1,length.out = length(tmp1$cDNA))) # assigning the colors to the gene 
barplot(rep(1,length(tmp1$colors)), col = tmp1$colors)