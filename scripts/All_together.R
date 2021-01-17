##############################
# merging the two dataframes #
##############################

setwd("home/rostyslav/Desktop/R_analysis/PhD_project/NGS_sequencing/4_cultures")
df1 <- read.table("final_results_4cultures.txt", header =  T, stringsAsFactors=FALSE)
setwd("/home/rostyslav/Desktop/R_analysis/PhD_project/NGS_sequencing/6_subcultures")
df2 <- read.table("6_subcultures_filtrated.txt", header =  T, stringsAsFactors=FALSE)
df2$range <- NULL
df2 <- gene_order(df2)
df1 <- sample_names(df1)
df2 <- sample_names(df2)
#freq_to_per <- function(x){
#  x$persentage <- paste0(x$VarFreq*100, rep("%"))
#  return(x)
#}
#df1 <- freq_to_per(df1)
#df2 <- freq_to_per(df2)
#df <- rbind(df1, df2)

#fixing
df_tmp <- split(df2, f = list(df2$Sample, df2$time_in_quiescence))
df_tmp <- lapply(df_tmp, function(x) x[-nrow(x),] )
df2 <- do.call("rbind", df_tmp) #transform the list into the dataframe
df2 <- sorting(df2)
df <- rbind(df1, df2)
df <- sorting(df) #sort everything again


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