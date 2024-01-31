# Evolutionary Genetics and Genomics (HGEN 6092)
# R Session 2 - Detecting Natural Selection I
# Casey Meili

# Some exercises were adapted from: https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html

# load packages
library("rehh")

# load dataset for chromosome 12 
data(haplohh_cgu_bta12)

# determine the number of chromosomes and SNPs
dset <- haplo(haplohh_cgu_bta12)
dset[1:5,1:5]
nrow(dset) # number of chromosomes sequenced
ncol(dset) # number of SNPs on chromosome 12
dset

# calculate and plot unfolded site frequency spectrum
daf <- colSums(dset)/nchr
hist(daf)

# analyzing SNP F1215500
# extract SNP data column 
snp_column <- dset[, "F1215500", drop = FALSE]
print(snp_column)

# calculate allele frequency 
allele_frequency <- mean(snp_column)
allele_frequency

# caluculate EHH for focal SNP
myres <- calc_ehh(haplohh_cgu_bta12,mrk="F1215500",include_nhaplo=T)
print(myres)

# plot EHH values
plot(myres)

# calculate IHS for all SNPs in chromosome 12 dataset
myscan <- scan_hh(haplohh_cgu_bta12)
myscan.ihs <- ihh2ihs(myscan) 

# plot the IHS values along the chromosome and p-values
manhattanplot(myscan.ihs)
manhattanplot(myscan.ihs,pval=TRUE,threshold=4)

# IHS value for SNP of interest
myscan.ihs$ihs["F1215500",]

# find the row with the lowest p-value (highest -log10(p value))
min_p_value_row <- myscan.ihs$ihs[which.min(myscan.ihs$ihs$LOGPVALUE), ]

# extract relevant information
marker_name <- rownames(min_p_value_row)
p_value <- 10^(-min_p_value_row$LOGPVALUE)  # convert -log10(p value) back to p value
ihs_value <- min_p_value_row$IHS
