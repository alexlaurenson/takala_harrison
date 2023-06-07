install.packages("circlize")

install.packages("BiocManager")
library("BiocManager")
BiocManager::install("VariantAnnotation", dependencies = TRUE)
BiocManager::install("GenomicRanges", force = TRUE)
install.packages("vcfR")
install.packages("dplyr")

library(circlize)
library(VariantAnnotation)
library(GenomicRanges)
library(vcfR)
library(dplyr)

vcfr <- read.vcfR("C:/Users/alexl/OneDrive/Documents/takala_harrison/Pf3D7_01_v3.snp.hardfilt_Q50_MQ40_QD2_FS60_MQRS12.5_RPRS8.vcf")

vcf <- readVcf("C:/Users/alexl/OneDrive/Documents/takala_harrison/Pf3D7_01_v3.snp.hardfilt_Q50_MQ40_QD2_FS60_MQRS12.5_RPRS8.vcf")

chr <- vcfr@fix
chrom <- chr[, "CHROM"]

pos <- start(vcf)[1:4817]
dp <- assays(vcf)$DP[1:4817]

length(chrom)

df <- data.frame(pos, dp)

# Create a circus plot with genome positions and depths
circos.initialize(factors = unique(df$position), xlim = c(0, 1))
circos.trackPlotRegion(ylim = c(0, max(df$depth)), bg.border = NA)

# Add lines representing the depth at each genome position
circos.lines(df$position, df$depth, col = "blue")

# Add labels for genome positions
circos.text(factors = df$position, labels = df$position, 
            sector.index = 1, cex = 0.6, facing = "inside", niceFacing = TRUE)
