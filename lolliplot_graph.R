#----------------------------------------------
#R setup
#----------------------------------------------

#set working directory
setwd("/Users/rostyslav/Desktop/4_cultures")

#prepare the packages
source("https://bioconductor.org/biocLite.R")
biocLite("trackViewer")

library(trackViewer)
require(reshape2)
library(gsubfn)
library(tidyr)
library(dplyr)

###############################################
#Data input
###############################################

df <- read.table(file = "for_lolliplot_annotation.txt", sep="\t", header = T, stringsAsFactors = F)
df$Ref <- as.character(df$Ref)
df$VarAllele <- as.character(df$VarAllele)
df$gene <- as.character(df$gene)
df$mutation_type <- as.character(df$mutation_type)
df$annotation <- as.character(df$annotation)

df <- df[order(df$Sample, df$time_in_quiescence, df$gene_order, df$cDNA),]

###############################################
#Data preproceccing
###############################################

#annotion for building plots
df$annotation_lolliplot <- 
  ifelse((df$mutation_type == "SNP"), paste0(df$Ref, df$cDNA, df$VarAllele),
         ifelse((df$mutation_type == "Insertion"  & nchar(df$VarAllele) <= 3), paste0(df$cDNA, df$VarAllele),
                ifelse((df$mutation_type == "Deletion" & nchar(df$VarAllele) <= 3), paste0(df$cDNA, df$VarAllele),
                       ifelse((df$mutation_type == "Insertion"  & nchar(df$VarAllele) > 3), paste0(df$cDNA, rep ("+"), as.character(nchar(df$VarAllele)-1), rep ("bp")),
                              ifelse((df$mutation_type == "Deletion"  & nchar(df$VarAllele) > 3), paste0(df$cDNA, rep ("-"), as.character(nchar(df$VarAllele)-1), rep ("bp")), paste0("complex"))))))

#annotion for colors
df$color <- 
  ifelse((df$mutation_type == "SNP" & df$gene == "sty1"), paste0("yellow"),
         ifelse((df$mutation_type == "SNP" & df$gene == "wis1"), paste0("yellow"),
                ifelse((df$mutation_type == "SNP" & df$gene == "tif452"), paste0("yellow"),
                      ifelse((df$mutation_type == "Insertion"), paste0("lightcoral"),
                            ifelse((df$mutation_type == "SNP" & df$gene == "win1"), paste0("black"),
                                   ifelse((df$mutation_type == "Deletion"), paste0("darkviolet"),
                                          ifelse((df$mutation_type == "Complex"), paste0("white"),
                                                 ifelse((df$cDNA == 250), paste0("black"),
                                                        ifelse((df$cDNA == 959), paste0("yellow"),
                                                               ifelse((df$cDNA == 2517), paste0("green"),
                                                                      ifelse((df$cDNA == 3040), paste0("yellow"),
                                                                             ifelse((df$cDNA == 326), paste0("black"),
                                                                                    ifelse((df$cDNA == 382), paste0("black"),
                                                                                           ifelse((df$cDNA == 1247), paste0("yellow"),
                                                                                                  ifelse((df$cDNA == 4001), paste0("yellow"),
                                                                                                         ifelse((df$cDNA == 593), paste0("yellow"), "green"))))))))))))))))

#updating SNPs for mkh1
tmp <- subset(df, gene == "mkh1" & mutation_type == "SNP")
tmp <- tmp[order(tmp$cDNA),]

df[71,"color"] <- "yellow" # mkh1-CT1884TG mkh1-F629V 
df[6,"color"] <- "yellow" # mkh1-T1884G mkh1-F629V 
df[125,"color"] <- "yellow" # mkh1-T1884G mkh1-F629V
df[38,"color"] <- "black" # mkh1-C2378G mkh1-793*
df[21,"color"] <- "yellow" # mkh1-T2494C mkh1-G974R
df[7,"color"] <- "yellow" # mkh1-T2914C mkh1-S974P
df[73,"color"] <- "yellow" # mkh1-G3124T mkh1-A1042S
df[74,"color"] <- "black" # mkh1-T3137G mkh1-1045*
df[39,"color"] <- "yellow" # mkh1-A3138C mkh1-L1045F
df[53,"color"] <- "black" # mkh1-G3145T mkh1-1049*
df[141,"color"] <- "yellow" # mkh1-A3152G mkh1-K1057R
df[8,"color"] <- "yellow" # mkh1-A3351T mkh1-*1117Y

#updating SNPs for win1
tmp <- subset(df, gene == "win1" & mutation_type == "SNP")
tmp <- tmp[order(tmp$cDNA),]

df[176,"color"] <- "yellow" # win1-C1751T win1-S1117F
df[15,"color"] <- "yellow" # win1-C2350A win1-P1117T
df[48,"color"] <- "yellow" # win1-T2381G win1-L1117R

#updating SNPs for pmk1
tmp <- subset(df, gene == "pmk1" & mutation_type == "SNP")
tmp <- tmp[order(tmp$cDNA),]

df[41,"color"] <- "yellow" # pmk1-G497T pmk1-C166F
df[43,"color"] <- "yellow" # pmk1-G577T pmk1-W193G
df[55,"color"] <- "yellow" # pmk1-G592A pmk1-E198K

#updating SNPs for pmc1
tmp <- subset(df, gene == "pmc1" & mutation_type == "SNP")
tmp <- tmp[order(tmp$cDNA),]

df[79,"color"] <- "black" # pmc1-G280T 

tmp <- subset(df, gene == "pmc1" & mutation_type == "Insertion")
df[29,"cDNA"] <- 2434 # pmc1-G280T 
df[29,"annotation"] <- "pmc1-2434+18bp 0%"
df[29,"annotation_lolliplot"] <- "2434+18bp"

#updating SNPs for pek1
tmp <- subset(df, gene == "pek1" & mutation_type == "SNP")
tmp <- tmp[order(tmp$cDNA),]

df[23,"color"] <- "black" #

df[23,"color"] <- "black" #
df[71,"annotation_lolliplot"] <- "CT1884TG"
df[118,"annotation_lolliplot"] <- "CC1281TT"

df$annotation_lolliplot <- ifelse((df$strain_direction == "-"), gsubfn(".", list("C" = "G", "G" = "C", "A" = "T", "T" = "A"), df$annotation_lolliplot), gsubfn(".", list("C" = "C", "G" = "G", "A" = "A", "T" = "T"), df$annotation_lolliplot))

df <- rbind(total, df2)


#----------------------------------------------
#plotting mutations on sty1 gene
#----------------------------------------------

#domain coordinates on cDNA
prokin <- c(60, 899) #protein kinase domain
ATP <- c(78, 105) #ATP-binding region 
ATP_b <- c(147, 149) #ATP binding site
proton <- c(423, 425) #Proton acceptor - active site 
TXY_1 <- c(513, 521) #short sequence motif - TXY
TXY_2 <- c(528, 536) #short sequence motif - TXY

# domain names
block_names <- c("cDNA", "protein kinase domain", "ATP-binding region", "", "", "short sequence motif - TXY", "short sequence motif - TXY")
#block_names <- c("", "", "", "", "", "", "")

tmp <- subset(df, gene == "sty1")
sty <- as.vector(tmp$cDNA)

#sty <- c(74, 175, 188, 197, 418, 481, 535, 545, 553, 558, 659) #mutation positions on cDNA in sty1 gene
seq_names <- tmp$annotation_lolliplot

sty.gr <- GRanges("chr1", IRanges(sty, width=1, names=paste0(seq_names))) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
sty_features <- GRanges("chr1", IRanges(c(1, 60, 78, 147, 423, 513, 528), # start coordinates for features (domains)
                                        width=c(1050, 893, 27, 3, 3, 8, 8), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

## change the height of features
sty_features$height <- c(0.02, 0.025, 0.05, 0.0, 0.0, 0.05, 0.05)

xaxis <- c(1, 200, 400, 600, 800, 1000, 1051)

# domain colours
sty_features$fill <- c("yellow2", "green3", "yellow3", "red1", "red3", "blue", "blue")

## mutation occurance
sty.gr$score <- rep(1:1, each = length(sty))

## Set up the color of lolliplot depending on the type of mutation
#sty.gr$color <- 
#  ifelse((tmp$mutation_type == "SNP"), "yellow2",
#         ifelse((tmp$mutation_type == "Insertion"), "blue1",
#                ifelse((tmp$mutation_type == "Deletion"), "black1", "white")))

sty.gr$color <- tmp$color

## add the legend
#legend <- c("yellow") ## legend fill color
#legend_names <- c("missense mutation")
#names(legend) <- paste0(legend_names) ## legend labels

#Control the labels
sty.gr.rot <- sty.gr
sty.gr.rot$label.parameter.rot <- 60
par(family = "Times New Roman")
lolliplot(sty.gr.rot, sty_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.8, family = "Times New Roman") #plot the graph

cairo_pdf(file = "lolliplot_sty1.pdf", width = max(xaxis)/100/2.54, height =15/2.54, pointsize = 9, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(sty.gr.rot, sty_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

cairo_pdf(file = "lolliplot_sgf73.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 9.5, family = 'Times New Roman', bg = "white", fallback_resolution=70)


#saving the data
library(Cairo)
cairo_pdf(file = "lolliplot_sty1.pdf", units = "cm", width = max(xaxis)/100/1.5, height =8.936, pointsize = 12, family = 'Times New Roman')
lolliplot(sty.gr.rot, sty_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=.8) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")

Cairo(file="lolliplot_sty1.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=11, 
      dpi=50,
      family = 'Times New Roman')
lolliplot(sty.gr.rot, sty_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8) #plot the graph
dev.off()

#----------------------------------------------
#plotting mutations on sgf73 gene
#----------------------------------------------

# domain names
block_names <- c("cDNA", "SCA7")
tmp <- subset(df, gene == "sgf73")

seq_names <- tmp$annotation_lolliplot

# mutations coordinates
#sgf73 <- c(890)
sgf73 <- tmp$cDNA

sgf.gr <- GRanges("chr1", IRanges(sgf73, width=1, names=seq_names)) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
sgf_features <- GRanges("chr1", IRanges(c(1, 570), # start coordinates for features (domains)
                                        width=c(1106, 185), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

xaxis <- c(1, 200, 400, 600, 800, 1000, 1106)

## change the height of features
sgf_features$height <- c(0.02, 0.025)

# domain colours
sgf_features$fill <- c("yellow2", "azure3")

## mutation occurance
#sgf_freq <- c(2)
#sgf.gr$score <- sgf_freq

## Change the color of lollipop.
sgf.gr$color <- c("lightcoral")

## add the legend
#legend <- c("lightcoral") ## legend fill color
#legend_names <- c("Insertion")
#names(legend) <- paste0(legend_names) ## legend labels

#Control the labels
sgf.gr.rot <- sgf.gr
sgf.gr.rot$label.parameter.rot <- 60
lolliplot(sgf.gr.rot, sgf_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.8) #plot the graph

cairo_pdf(file = "lolliplot_sgf73.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 9.5, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(sgf.gr.rot, sgf_features, legend=legend, xaxis=xaxis,ylab = NULL, yaxis=FALSE, cex=.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

#
cairo_pdf(file = "lolliplot_sgf73.pdf", width = max(xaxis)/100/1.5, height =8.936, pointsize = 15, family = 'Times New Roman')
lolliplot(sgf.gr.rot, sgf_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

Cairo(file="lolliplot_sgf73.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=11, 
      dpi=55,
      family = 'Times New Roman')
lolliplot(sgf.gr.rot, sgf_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8) #plot the graph
dev.off()


#----------------------------------------------
#plotting mutations on mkh1 gene
#----------------------------------------------

# domain names
block_names <- c("", "", "", "", "")

tmp <- subset(df, gene == "mkh1")
mkh1 <- tmp$c.DNA.position

mkh <- as.numeric(tmp$cDNA)
mkh_names <- tmp$annotation_lolliplot

mkh.gr <- GRanges("chr1", IRanges(mkh, width=1, names=paste0(mkh_names))) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
mkh_features <- GRanges("chr1", IRanges(c(1, 2475, 2493, 2519, 2865), # start coordinates for features (domains)
                                        width=c(3351, 791, 26, 3, 3), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

xaxis <- c(1, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3351)

## change the height of features
mkh_features$height <- c(0.02, 0.025, 0.05, 0.05, 0.0)

# domain colours
mkh_features$fill <- c("yellow2", "green3", "yellow3", "red1", "red3")

## mutation occurance

#tmp <- subset(df, gene == "sty1")
#sty <- as.vector(tmp$cDNA)

tmp <- subset(tmp, gene == "mkh1")
unique(mkh1)

mkh_freq <- rep.int(1, length(tmp$gene))

#mkh_freq <- c(1, 9, 1, 1, 1, 1, 1, 1, 1)
mkh.gr$score <- mkh_freq

## Change the color of lollipop.
mkh.gr$color <- tmp$color

## add the legend
legend <- c("lightcoral", "darkviolet", "black", "yellow2" , "green") ## legend fill color
legend_names <- c("Insertion", "Deletion", "Nonsence mutation", "Missense mutation", "Synonymous mutation")
names(legend) <- paste0(legend_names) ## legend labels

mkh.gr.rot <- mkh.gr
mkh.gr.rot$label.parameter.rot <- 60
lolliplot(mkh.gr.rot, mkh_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.8) #plot the graph

#
cairo_pdf(file = "lolliplot_mkh1.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 12, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(mkh.gr.rot, mkh_features, legend=legend, xaxis=xaxis,ylab = NULL, yaxis=FALSE, cex=.8) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

?cairo_pdf

pdfFonts("Times")
names(pdfFonts())
par(family= "Times")
Cairo(file="lolliplot_mkh1.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=10, 
      dpi=70,
      family= "serif",
      fonts = "serif")
lolliplot(mkh.gr.rot, mkh_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8, family= "Times") #plot the graph
dev.off()

cairo_pdf(file = "lolliplot_mkh1.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 12, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(mkh.gr.rot, mkh_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

#----------------------------------------------
#plotting mutations on pek1 gene
#----------------------------------------------

# domain names
block_names <- c("", "", "", "", "", "", "")

tmp <- subset(df, gene == "pek1")
pek <- as.vector(tmp$cDNA)

pek.gr <- GRanges("chr1", IRanges(pek, width=1, names=paste0(tmp$annotation_lolliplot))) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
pek_features <- GRanges("chr1", IRanges(c(1, 237, 255, 324, 618, 702, 714), # start coordinates for features (domains)
                                        width=c(1092, 794, 26, 3, 3, 3, 3), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

xaxis <- c(1, 200, 400, 600, 800, 1000, 1092)

## change the height of features
pek_features$height <- c(0.02, 0.025, 0.05, 0.0, 0.0, 0.0, 0.0)

## change the colour of features
pek_features$fill <- c("yellow2", "green3", "yellow3", "red1", "red3", "white", "white")

## Change the color of lollipop.
pek.gr$color <- tmp$color

## add the legend
#legend <- c("lightcoral", "darkviolet", "yellow2") ## legend fill color
#legend_names <- c("frameshift insertion", "frameshift deletion", "missense mutation")
#names(legend) <- paste0(legend_names) ## legend labels

pek.gr.rot <- pek.gr
pek.gr.rot$label.parameter.rot <- 60

cairo_pdf(file = "lolliplot_pek1.pdf", width = max(xaxis)/100/1.5, height =8.936, pointsize = 15, family = 'Times New Roman')
lolliplot(pek.gr.rot, pek_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

Cairo(file="lolliplot_pek1.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=11, 
      dpi=55,
      family = 'Times New Roman')
lolliplot(pek.gr.rot, pek_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8) #plot the graph
dev.off()

cairo_pdf(file = "lolliplot_pek1.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 9, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(pek.gr.rot, pek_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

#----------------------------------------------
#plotting mutations on pmk1 gene
#----------------------------------------------

# domain names
block_names <- c("", "", "", "", "", "", "TXY short sequence motif", "", "Poly-Ser")

tmp_pmk <- subset(df, gene == "pmk1")
pmk <- tmp_pmk$cDNA

pmk.gr <- GRanges("chr1", IRanges(pmk, width=1, names=paste0(tmp_pmk$annotation_lolliplot))) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
pmk_features <- GRanges("chr1", IRanges(c(1, 63, 81, 156, 447, 558, 558, 564, 1167), # start coordinates for features (domains)
                                        width=c(1269, 944, 26, 3, 3, 3, 6, 3, 17), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

xaxis <- c(1, 200, 400, 600, 800, 1000, 1200, 1269)

## change the height of features
pmk_features$height <- c(0.02, 0.025, 0.05, 0.00, 0.00, 0.00, 0.05, 0.00, 0.025)

## change the colour of features
pmk_features$fill <- c("yellow2", "green3", "yellow3", "red1", "red3", "white", "blue", "white", "darkgray")

## Change the color of lollipop.
pmk.gr$color <- tmp_pmk$color

## add the legend
#legend <- c("lightcoral", "darkviolet", "yellow2") ## legend fill color
#legend_names <- c("frameshift insertion", "frameshift deletion", "missense mutation")
#names(legend) <- paste0(legend_names) ## legend labels

pmk.gr.rot <- pmk.gr
pmk.gr.rot$label.parameter.rot <- 60
lolliplot(pmk.gr.rot, pmk_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.8) #plot the graph

cairo_pdf(file = "lolliplot_pmk1.pdf", width = max(xaxis)/100/1.5, height =8.936, pointsize = 15, family = 'Times New Roman')
lolliplot(pmk.gr.rot, pmk_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

Cairo(file="lolliplot_pmk1.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=11, 
      dpi=55,
      family = 'Times New Roman')
lolliplot(pmk.gr.rot, pmk_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8) #plot the graph
dev.off()

cairo_pdf(file = "lolliplot_pmk1.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 10, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(pmk.gr.rot, pmk_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()


#----------------------------------------------
#plotting mutations on pmc1 gene
#----------------------------------------------

block_names <- c("cDNA", "Poly-Asp", "Ser-rich", "Helical-transmembrane region", "Helical-transmembrane region", "Helical-transmembrane region", "Helical-transmembrane region", "",
                 "Helical-transmembrane region", "Helical-transmembrane region", "Helical-transmembrane region", "Helical-transmembrane region", "Helical-transmembrane region")

#block_names <- c("", "", "", "", "", "", "",
#                 "", "", "", "", "")


tmp_pmc1 <- subset(df, gene == "pmc1")
pmc <- tmp_pmc1$cDNA

pmc.gr <- GRanges("chr1", IRanges(pmc, width=1, names=paste0(tmp_pmc1$annotation_lolliplot))) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
pmc_features <- GRanges("chr1", IRanges(c(1, 18, 30, 711, 822, 1347, 1467, 1635, 2817, 2901, 3051, 3255, 3348), # start coordinates for features (domains)
                                        width=c(3879, 11, 287, 60, 62, 62, 62, 3, 62, 62, 62, 62, 62), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

test <- seq(0, 3879, by=200)
test[1] <- test[1] +1
test <- append(test, 3879)
xaxis <- test

## change the height of features
pmc_features$height <- c(0.02, 0.0, 0.025, 0.025, 0.025, 0.025, 0.025, 0.0, 0.025, 0.025, 0.025, 0.025, 0.025)

## change the colour of features
pmc_features$fill <- c("yellow2", "darkgreen", "lightblue", "red1", "red1", "red1", "red1", "black", "red1", "red1", "red1", "red1", "red1")

## mutation occurance
#pmc_freq <-length(tmp$gene)
pmc.gr$score <- 1

pmc.gr$color <- tmp_pmc1$color

## add the legend
#legend <- c("lightcoral", "darkviolet", "black", "yellow2") ## legend fill color
#legend_names <- c("frameshift insertion", "frameshift deletion", "nonsence mutation", "missense mutation")
#names(legend) <- paste0(legend_names) ## legend labels

pmc.gr.rot <- pmc.gr
pmc.gr.rot$label.parameter.rot <- 60
lolliplot(pmc.gr.rot, pmc_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.8) #plot the graph

cairo_pdf(file = "lolliplot_pmc1.pdf", width = max(xaxis)/100/1.5, height =12, pointsize = 15, family = 'Times New Roman')
lolliplot(pmc.gr.rot, pmc_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

Cairo(file="lolliplot_pmc1.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=10, 
      dpi=70)
lolliplot(pmc.gr.rot, pmc_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8) #plot the graph
dev.off()

cairo_pdf(file = "lolliplot_pmc1.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 12, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(pmc.gr.rot, pmc_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=.8) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()


#----------------------------------------------
#plotting mutations on wis1 gene
#----------------------------------------------

#block_names <- c("cDNA", "Phosphoserine - modified residue", "Protein kinase domain", "ATP binding region", "ATP - binding site",
#                 "Proton acceptor - active site", "Phosphoserine - modified residue", "Phosphoserine - modified residue")

block_names <- c("", "", "", "", "",
                 "", "", "")


tmp_wis1 <- subset(df, gene == "wis1" & cDNA != is.na(cDNA))
wis <- tmp_wis1$cDNA

wis.gr <- GRanges("chr1", IRanges(wis, width=1, names=paste0(tmp_wis1$annotation_lolliplot))) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
wis_features <- GRanges("chr1", IRanges(c(1, 759, 960, 978, 1047, 1323, 1407, 1419), # start coordinates for features (domains)
                                        width=c(1818, 3, 779, 44, 3, 3, 3, 3), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

test <- seq(0, 1818, by=200)
test[1] <- test[1] +1
test <- append(test, 1818)
xaxis <- test

## change the height of features
wis_features$height <- c(0.02, 0.0, 0.025, 0.05, 0.0, 0.0, 0.0, 0.0)

## change the colour of features
wis_features$fill <- c("yellow2", "black", "green3", "yellow3", "white", "blue", "white", "lightgray")

## mutation occurance
wis_freq <- rep.int(1, length(tmp_wis1$gene))
wis.gr$score <- wis_freq

## Change the color of lollipop.
wis.gr$color <- tmp_wis1$color

wis.gr.rot <- wis.gr
wis.gr.rot$label.parameter.rot <- 60

## add the legend
legend <- c("lightcoral", "yellow2") ## legend fill color
#legend_names <- c("frameshift insertion", "missense mutation")
#names(legend) <- paste0(legend_names) ## legend labels

cairo_pdf(file = "lolliplot_wis1.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 11, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(wis.gr.rot, wis_features, legend=legend, xaxis=xaxis,ylab = NULL, yaxis=FALSE, cex=.8) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

cairo_pdf(file = "lolliplot_wis1.pdf", width = max(xaxis)/100/1.5, height =8.936, pointsize = 15, family = 'Times New Roman')
lolliplot(wis.gr.rot, wis_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

Cairo(file="lolliplot_wis1.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=11, 
      dpi=60,
      family = 'Times New Roman')
lolliplot(wis.gr.rot, wis_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8) #plot the graph
dev.off()

#----------------------------------------------
#plotting mutations on tif452 gene
#----------------------------------------------

block_names <- c("cDNA")

tmp_tif452 <- subset(df, gene == "tif452")
tif <- tmp_tif452$cDNA

tif.gr <- GRanges("chr1", IRanges(tif, width=1, names=paste0(tmp_tif452$annotation_lolliplot))) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
tif_features <- GRanges("chr1", IRanges(c(1), # start coordinates for features (domains)
                                        width=c(732), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

test <- seq(0, 732, by=200)
test[1] <- test[1] +1
test <- append(test, 732)
xaxis <- test

## change the height of features
tif_features$height <- c(0.02)

## change the colour of features
tif_features$fill <- c("yellow2")

## Change the color of lollipop.
tif.gr$color <- tmp_tif452$color

## add the legend
legend <- c("yellow2") ## legend fill color
legend_names <- c("missense mutation")
#names(legend) <- paste0(legend_names) ## legend labels

tif.gr.rot <- tif.gr
tif.gr.rot$label.parameter.rot <- 60
lolliplot(tif.gr.rot, tif_features, legend=legend, xaxis=xaxis, yaxis=FALSE, cex=0.8) #plot the graph

cairo_pdf(file = "lolliplot_tif452.pdf", width = (max(xaxis)/100/2.7)+3, height =8, pointsize = 12, family = 'Times New Roman')
lolliplot(tif.gr.rot, tif_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()

Cairo(file="lolliplot_tif452.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=11, 
      dpi=50,
      family = 'Times New Roman')
lolliplot(tif.gr.rot, tif_features, legend=legend, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8) #plot the graph
dev.off()

cairo_pdf(file = "lolliplot_tif452.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 8, family = 'Times New Roman', bg = "white", fallback_resolution=55)
lolliplot(tif.gr.rot, tif_features, legend=legend, xaxis=xaxis,ylab = NULL, yaxis=FALSE, cex=.9) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()


#----------------------------------------------
#plotting mutations on win1 gene
#----------------------------------------------

block_names <- c("cDNA", "", "", "interaction with tea4", "protein kinase domain", "ATP binding region", "", "" )

#mutation positions in win1 gene regarding the genome

tmp_win1 <- subset(df, gene ==  "win1")
win <- as.numeric(tmp_win1$cDNA)

#uploading mutations

win.gr <- GRanges("chr1", IRanges(win, width=1, names = paste0(tmp_win1$annotation_lolliplot))) #align SNPs from the vector to the coordinates(?)
#you cannot plot sample without the features!
win_features <- GRanges("chr1", IRanges(c(1, 672, 678, 846, 3360, 3378, 3447, 3732), # start coordinates for features (domains)
                                        width=c(4311, 3, 3, 2525, 860, 44, 3, 3), # lenght for features (domains)
                                        names=paste0(block_names))) # create place for domain names

test <- seq(0, 4311, by=200)
test[1] <- test[1] +1
test <- append(test, 4311)
xaxis <- test

## change the height of features
win_features$height <- c(0.02, 0.0, 0.0, 0.025, 0.025, 0.05, 0.0, 0.0)

## change the colour of features
win_features$fill <- c("yellow2", "green3", "yellow3", "cyan", "green3", "yellow3", "blue", "white")

## Change the color of lollipop.
win.gr$color <- tmp_win1$color

## add the legend
#legend <- c("lightcoral", "darkviolet", "black", "yellow2") ## legend fill color
#legend_names <- c("frameshift insertion", "frameshift deletion", "nonsence mutation", "missense mutation")

win.gr.rot <- win.gr
win.gr.rot$label.parameter.rot <- 60

par(family = 'Times New Roman')
lolliplot(SNP.gr = win.gr.rot, features = win_features, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.9) #plot the graph


cairo_pdf(file = "lolliplot_win11.pdf", width = max(xaxis)/100/2.54, height =15/2.54,  pointsize = 12, family = 'Times New Roman', bg = "white", fallback_resolution=70)
lolliplot(win.gr.rot, win_features, legend=legend, xaxis=xaxis,ylab = NULL, yaxis=FALSE, cex=.8) #plot the graph
#legend(x = .05, y = .45, culture_0$annotation, col = "white",  cex = 1.6, fill = culture_0$colors, trace = T, text.col = "black", bty = "o", box.lwd = 1, box.col = "white")
dev.off()


Cairo(file="lolliplot_win1.pdf", 
      bg="white",
      type="pdf",
      units="cm", 
      width=max(xaxis)/100, 
      height=15, 
      pointsize=10, 
      dpi=70)
lolliplot(SNP.gr = win.gr.rot, features = win_features, xaxis=xaxis, ylab = NULL, yaxis=FALSE, cex=0.8) #plot the graph
dev.off()



#######
#tmp_data
#######

setwd("/Users/rostyslav/Desktop")
tmp_win1 <- read.table(file = "tmp_win1.csv", sep=";", header = T) # read the final exel file 
tmp_win1 <- sapply(tmp_win1, as.character) # preporocessing
tmp_win1 <- data.frame(tmp_win1, stringsAsFactors = F) # preporocessing
tmp <- subset(df, gene == "win1") # subsetting the win1 mutations from output
tmp$match <- match(tmp$cDNA, tmp_win1$c.DNA.position) # check which mutations are matched by cDNA position
tmp <- tmp %>% drop_na()

tmp$match <- match(tmp$cDNA, tmp_win1$c.DNA.position) 
tmp_win1$match <- match(tmp_win1$c.DNA.position, tmp$cDNA) 

tmp[59,"annotation_lolliplot"] <- "complex"
tmp[59,"color"] <- "white"
tmp[1,"color"] <- "yellow"
tmp[20,"color"] <- "yellow"
tmp[106,"annotation"]
tmp<- removeRows(68, tmp)
tmp <- rbind(tmp, c(NA, NA, NA, 0.005, "+ATACTAGTAAT", "wt3_2m", "Insertion", "1st_4cultures", "win1", 383, "0,5%", "+", "win1-383+11bp", "383+11bp", "lightcoral", NA))
tmp <- rbind(tmp, c(NA, NA, NA, 0.002, "+CTACTACATCCT", "wt0_2m", "Insertion", "1st_4cultures", "win1", 394, "0,2%", "+", "win1-394+13bp", "394+13bp", "lightcoral", NA))
tmp <- rbind(tmp, c(NA, NA, NA, 0.002, "+CTACTACATCCT", "wt3_2m", "Insertion", "1st_4cultures", "win1", 394, "0,2%", "+", "win1-394+13bp", "394+13bp", "lightcoral", NA))



tmp$cDNA <- as.numeric(tmp$cDNA)
tmp_win1$cDNA <- as.numeric(tmp_win1$cDNA)
tmp_win1$c.DNA.position <- as.numeric(tmp_win1$c.DNA.position)

tmp <- tibble::rowid_to_column(tmp, "ID")
tmp$ID <- NULL

tmp$match <- NULL
tmp_win1$match <- NULL

tmp_win1 <- tmp

tmp$persentage<- as.character(tmp$persentage)
tmp$annotation<- as.character(tmp$annotation)
tmp$annotation_lolliplot<- as.character(tmp$annotation_lolliplot)
arrange(tmp)
lolliplot_annotation(tmp_win1)
tmp2<-tmp[!(is.an(tmp$match == NA),]

tmp2 <- tmp %>% drop_na()
