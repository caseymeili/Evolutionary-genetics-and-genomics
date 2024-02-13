# Evolutionary Genetics and Genomics (HGEN 6092)
# R Session 4 - Population demography and structure (Part 1)
# Casey Meili

# set working directory
setwd("/Users/caseymeili/Desktop/Utah/Courses/Evol Genetics:Genomics")

# import dataset and check dimensions
dset<-read.table("drosophila_minor_allele_counts.txt",head=T)
head(dset) # sum of major and minor alleles for each location = sample size for each population
dim(dset) # each row is a SNP 17219,  

# site frequency spectrum for individual populations 
# count up how many sites occur with each allele count
sfs_zim<-table(dset[,"ZIM_MINOR"]) # Zimbabwe
sfs_zim 
sfs_bei<-table(dset[,"BEI_MINOR"]) # Beijing
sfs_bei 
sfs_ith<-table(dset[,"ITH_MINOR"]) # Ithaca
sfs_ith 
sfs_net<-table(dset[,"NET_MINOR"]) # Netherlands
sfs_net 
sfs_tas<-table(dset[,"TAS_MINOR"]) # Tasmania
sfs_tas 

# print all bar plots 
par(mfrow=c(2,3)) # to plot all five together
barplot(sfs_zim[-1],main="Zimbabwe")
barplot(sfs_bei[-1],main="Beijing")
barplot(sfs_ith[-1],main="Ithica")
barplot(sfs_net[-1],main="Netherlands")
barplot(sfs_tas[-1],main="Tasmania")

# function to calculate proportion of singletons for each population
calculate_singleton_proportion <- function(df, population) {
  # extract the column corresponding to the minor allele counts for the specified population
  minor_allele_counts <- df[[paste0(population, "_MINOR")]]
  
  # calculate total number of variants with at least one copy (sum of minor allele counts >= 1)
  total_variants_with_copy <- sum(minor_allele_counts >= 1, na.rm = TRUE)
  
  # calculate number of singletons (variants with exactly one copy)
  singletons <- sum(minor_allele_counts == 1, na.rm = TRUE)
  
  # calculate proportion of singletons
  proportion_singletons <- singletons / total_variants_with_copy
  
  return(proportion_singletons)
}

# apply the function to the data frame for each population
proportion_singletons_zim <- calculate_singleton_proportion(dset, "ZIM")
proportion_singletons_bei <- calculate_singleton_proportion(dset, "BEI")
proportion_singletons_ith <- calculate_singleton_proportion(dset, "ITH")
proportion_singletons_net <- calculate_singleton_proportion(dset, "NET")
proportion_singletons_tas <- calculate_singleton_proportion(dset, "TAS")

# print out the proportions
proportion_singletons_zim
proportion_singletons_bei
proportion_singletons_ith
proportion_singletons_net
proportion_singletons_tas


# joint frequency spectrums 
jsfs_zim_bei<-table(dset[,"ZIM_MINOR"],dset[,"BEI_MINOR"]) # Zimbabwe + Beijing 
jsfs_zim_bei 
jsfs_zim_bei<-table(dset[,"NET_MINOR"],dset[,"BEI_MINOR"]) # Beijing + Netherlands
jsfs_zim_bei 

# plot heat maps
heatmap(jsfs_zim_bei,Colv=NA,Rowv=NA)
heatmap(jsfs_zim_bei,Colv=NA,Rowv=NA)

