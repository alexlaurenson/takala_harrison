#qc visualization post-qc filtering using bcftools
#load tools
library(vcfR)
library(SNPfiltR)
library(pinfsc50)
library(data.table)
library(ape)
#read files into play area
vcfr <- read.vcfR("/local/scratch/matta/icemr_p2/W10026114_W100291417/called_snps/result/hardfilt/Pf3D7_01_v3.snp.hardfilt_Q50_MQ40_QD2_FS60_MQRS12.5_RPRS8.vcf", verbose = FALSE)
dna <- ape::read.dna("/local/scratch/matta/icemr_p2/W10026114_W100291417/called_snps/Pf3D7_01_v3.fasta", format = "fasta")
gff <- fread("/local/scratch/matta/icemr_p2/W10026114_W100291417/called_snps/Pf3D7_01_v3.gff", sep="\t", quote="", fill=TRUE)
#create chrom file using loaded files
chrom <- create.chromR(name='Supercontig', vcf=vcfr, seq=dna, ann=gff)
plot(chrom)
chrom_proc <- proc.chromR(chrom, verbose=TRUE)
plot(chrom_proc)
chromoqc(chrom_proc, boxp=TRUE, dp.alpha=255)
#read csv file for popmap
umbplate <- read.csv("/local/scratch/matta/icemr_p2/W10026114_W100291417/called_snps/key_WGSnames_18sqPCR_popmap.csv")
#missingness by sample
vcfr_miss <-missing_by_sample(vcfR=vcfr, cutoff = .7)
#missingness by snp
vcfr_misssnp <- missing_by_snp(vcfr)
assess_missing_data_pca(vcfR = vcfr, popmap=umbplate, thresholds = c(.7), clustering = FALSE)
# Read the data from the file, skipping the first line
data <- read.csv("/local/scratch/matta/icemr_p2/W10026114_W100291417/called_snps/coverage_all.txt", header = TRUE, sep = ",")
# Rename the first column to "id"
colnames(data)[1] <- "id"
# Read the data from the file, skipping the first line
key <- read.csv("/home/matt.adams/projects/icemr_p2/input/key_WGSnames_18sqPCR.csv", header = TRUE, sep = ",")
# Replace - / \ with _ in the 'id' column
key$id <- gsub("[-/\\\\]", "_", key$umb_rna_id)
data1 <- data
data <- right_join(key, data1, by="id") 
data[1:3,1:20]
data <- data[ , c("id", "umb_acd_pcd", "umb_pf_mean_dna_ngul", "avg_depth","per_5xcov",	"per_10xcov","per_25xcov","per_50xcov","per_100xcov")]
head(data)
dim(data)
data <- arrange(data, desc(avg_depth))
# Install tidyverse if it's not already installed
# install.packages("tidyverse")

# Load the required libraries
library(tidyverse)

# Pivot the data to longer format
long_data <- data %>%
  pivot_longer(
    cols = starts_with("per_5"),
    names_to = "coverage",
    values_to = "value"
  )

# Create a line plot
ggplot(long_data, aes(x = id, y = value, group = coverage, color = coverage)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Coverage values per sample",
       x = "Sample ID",
       y = "Coverage Value",
       color = "Coverage")
# Load the required library
library(ggplot2)

# Pivot the data to longer format including 'column_3'
long_data <- data %>%
  pivot_longer(
    cols = c(starts_with("per_5x"), "avg_depth"),
    names_to = "coverage",
    values_to = "value"
  )

# Create a line plot
p <- ggplot(long_data, aes(x = id)) +
  geom_line(aes(y = value, group = coverage, color = coverage), data = subset(long_data, coverage != "avg_depth")) +
  geom_line(aes(y = value, group = coverage, color = coverage), data = subset(long_data, coverage == "avg_depth")) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Average Depth")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Average Percent 5x Coverage and Average Depth Per Sample",
       x = "Sample",
       y = "Coverage Value",
       color = "Coverage")

# Print the plot
print(p)
# Load the required library
library(ggplot2)

# Pivot the data to longer format including 'column_3'
long_data <- data %>%
  pivot_longer(
    cols = c(starts_with("per_50x"), "avg_depth"),
    names_to = "coverage",
    values_to = "value"
  )

# Create a line plot
p <- ggplot(long_data, aes(x = id)) +
  geom_line(aes(y = value, group = coverage, color = coverage), data = subset(long_data, coverage != "column_3")) +
  geom_line(aes(y = value, group = coverage, color = coverage), data = subset(long_data, coverage == "column_3")) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Average Depth")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Average 50x Coverage and Average Depth Per Sample",
       x = "Sample",
       y = "Coverage Value",
       color = "Coverage")

# Print the plot
print(p)
# Load the required libraries
library(ggplot2)
library(dplyr)

# Pivot the data to longer format including 'avg_depth'
long_data <- data %>%
  pivot_longer(
    cols = c(starts_with("per_5x"), "avg_depth"),
    names_to = "coverage",
    values_to = "value"
  )

# Sort 'id' in the order of 'avg_depth'
sorted_ids <- data %>% 
  arrange(desc(avg_depth)) %>% 
  pull(id)

# Make 'id' a factor with levels in the order of 'sorted_ids'
long_data$id <- factor(long_data$id, levels = sorted_ids)

# Create a line plot
p <- ggplot(long_data, aes(x = id)) +
  geom_line(aes(y = value, group = coverage, color = coverage), data = subset(long_data, coverage != "avg_depth")) +
  geom_line(aes(y = value, group = coverage, color = coverage), data = subset(long_data, coverage == "avg_depth")) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Average Depth")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Average Percent 5x Coverage and Average Depth Per Sample",
       x = "Sample",
       y = "Coverage Value",
       color = "Coverage")

# Print the plot
print(p)
# Load the required libraries
library(ggplot2)
library(dplyr)

# Pivot the data to longer format including 'avg_depth'
long_data <- data %>%
  pivot_longer(
    cols = c(starts_with("per_5x"), "avg_depth"),
    names_to = "coverage",
    values_to = "value"
  )

# Make sure "per_5xcov" is in your data frame
# Adjust the column name as necessary to match your actual data

# Sort 'id' in the order of 'per_5xcov'
sorted_ids <- data %>% 
  arrange(desc(per_5xcov)) %>% 
  pull(id)

# Make 'id' a factor with levels in the order of 'sorted_ids'
long_data$id <- factor(long_data$id, levels = sorted_ids)

# Create a line plot
p <- ggplot(long_data, aes(x = id)) +
  geom_line(aes(y = value, group = coverage, color = coverage), data = subset(long_data, coverage != "avg_depth")) +
  geom_line(aes(y = value, group = coverage, color = coverage), data = subset(long_data, coverage == "avg_depth")) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Average Depth")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Average Percent 5x Coverage and Average Depth Per Sample",
       x = "Sample",
       y = "Coverage Value",
       color = "Coverage")

# Print the plot
print(p)