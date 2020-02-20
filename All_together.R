##############################
# merging the two dataframes #
##############################

setwd("/home/rostyslav/Desktop/R_analysis/NGS_sequencing")
df1 <- read.table("final_results_4cultures.txt", header =  T, stringsAsFactors=FALSE)
#setwd("/home/rostyslav/Desktop/R_analysis/NGS_sequencing")
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
df <- rbind(df1, df2)
df <- sorting(df) #sort everything again
#write.table(df, file = "final_variants_and_table_total.txt", sep="\t", quote=T, row.names=F, col.names=T)

#########################################################
# subsetting and manipulationg with win1               ##
# tempo <- df[!(df$gene=="win1" & df$VarFreq < 0.01),] ##
#########################################################

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
#dt_test <- rest_of_population(df)
#tmp <- rest_of_population(df1)
df$gene_order <- as.integer(df$gene_order)

#############################################################
# assigning the colors, based on the gene and cDNA position #
# x - is a dataframe, y is a gene,                          #
#############################################################

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

#########################
# building the pieplots #
#########################

df_tmp <- split(df_tmp, f = list(df_tmp$Sample, df_tmp$time_in_quiescence))
list2env(df_tmp ,.GlobalEnv) 

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

write.table(df, file = "final_variants_and_table_total.txt", sep="\t", quote=T, row.names=F, col.names=T)