# Rare mutation detection from S.pombe quiescence population 

# Background
Precise detection of mutations in Genomics data is a challenging process. Up-to-date sequencing platforms generate the data
with various level of technical errors (Edward J Fox et al 2014). Moreover, sequencing errors are superimposed on the errors obtained during sample preparation or due to imperfect experimental design. Exact and accurate somatic mutations calling is necessary for diagnostic and preventing tumorigenesis at an early stage; the identification of the rare virus population is important to counter infection outbreak. Therefore the development of the sufficient algorithm to reduce the noise in the sequencing data is crucial. 

# Summary
This project contains custom-developed R scripts to detect rare mutations (up to 0.1% of the population) from S.pombe quiescence populations. Briefly, the data represents a targeted resequencing experiment (Illumina) of the 9 genes of S/MAPK pathways. For more details, please follow to Makarenko et al. (2020). If you are using this code, please cite this paper. Scripts also include some minor statistical analysis as a number of the mutations per gene and mutations annotation and visualisation. 

# Usage

The scripts have to be run in the following order:

functions.R

Variants_distribution_independent_cultures.R

Variants_distribution_subcultures.R

All_together.R

mutation_vis_trackViewer.R

# functions.R 

This script uploads the libraries and the custom-developed functions for mutations  preprocessing, filtering, annotation and visualisation technics used in this analysis

# Variants_distribution_independent_cultures.R 

This script takes as an input multiple .txt files generated by VarScan software  in a table-like format. The output .txt file contains all filtered mutations from the independent cultures

The algorithm for detecting mutations is following:

1) keep the variants that origin from position with coverage > 1000 reads/position
2) keep the variants that have the difference of a base quality with the reference sequenced allele no more that 7% 
3) keep the SNPs that are present in the equvalent proportion on the both pair-end reads with 20% difference
4) keep the variants that have passed the Fisher_s p-value test
5) discarb the indels of size +/- 1 b.p. below frequency 0.5% in the population because it is not possible to distinguish them from the sequencing errors 

Then the script add some annotation data for detected mutations: postion on cDNA, strain direction (+/-), gene name, time in quiescence and the name of the culture name, sorting etc.

Finally the script plots the distribution of the frequencies of the final variants  for the SNPs and the indels and saves the resulting plots into .tiff files

# Variants_distribution_subcultures.R 

This script is almost as same as Variants_distribution_independent_cultures.R, but in addition it maintains only the variants that were found in common between sequencing dublicates. This script also  calculates the individual mutation frequency mean and range within the dublicates. The output .txt file contains all filtered mutations from the subcultures

# All_together.R #

This script merges the data from independent cultures and subcultures. Then it builds and annotates the pieplots. It assign the colors based on the gene and mutations position in the gene's cDNA withthe resulting plots saved into separate .pdf files

# mutation_vis_trackViewer.R #

This script visualises the mutations (SNPs and indels) on the gene cDNA sequence.
The black color corresponds to the nonsence SNP mutation, 
the green one - synonymus SNP mutation,
the yellow one - missense SNP mutation,
the lightcoral - insertion,
the violet - deletion,
The resulting plots are saved in separate .pdf files
