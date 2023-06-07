install.packages("circlize")

install.packages("BiocManager")
library("BiocManager")
BiocManager::install("VariantAnnotation", dependencies = TRUE)
BiocManager::install("GenomicRanges", force = TRUE)

install.packages("dplyr")

library(circlize)
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)

vcfr <- read.vcfR("C:/Users/alexl/OneDrive/Documents/takala_harrison/Pf3D7_all_v3.snp.hardfilt_Q50_MQ40_QD2_FS60_MQRS12.5_RPRS8.vcf")

vcf <- readVcf("C:/Users/alexl/OneDrive/Documents/takala_harrison/Pf3D7_all_v3.snp.hardfilt_Q50_MQ40_QD2_FS60_MQRS12.5_RPRS8.vcf")

chr <- vcfr@fix
chrom <- chr[, "CHROM"]
pos <- chr[, "POS"]

chrom <- vcf$CHROM
pos <- vcf$POS
depth <- vcf$DP

df <- data.frame(chrom, pos, depth)
