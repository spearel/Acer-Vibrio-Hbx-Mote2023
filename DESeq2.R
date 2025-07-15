##### STEP 1 - PREPARE ENVIRONMENT AND LOAD DATA
#set working directory
setwd("/Users/lspeare3/Desktop/FLK_2022_16S_raw_data/LS01-421351153/raw-fastq-files/DESeq2_APR16-2025")

#Install necessary stuff
install.packages("devtools")
library("devtools")
BiocManager::install(version = "3.19")
BiocManager::install("dada2")
library(dada2); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())

#Import OTU Table
otu_table = read.table ("ASV_table.csv", header = TRUE, row.names = 1, sep = ",")
otu_table_matrix <- as.matrix(otu_table)
otu_table_phyloseq <- otu_table(otu_table_matrix, taxa_are_rows = FALSE)
physeq_object <- phyloseq(otu_table_phyloseq)
physeq_object
#MetaData
sample_data = read.table ("metadata.csv", header = TRUE, row.names = 1, sep = ",")
sample_data_phyloseq <- (sample_data(sample_data))
#Taxa data
tax_table = read.table ("Taxa.csv", header = TRUE, row.names = 1, sep = ",")
tax_table_matrix <- as.matrix(tax_table)
tax_table_phyloseq <- tax_table(tax_table_matrix)

#Combine all into one phyloseq object
physeq_object <- phyloseq(otu_table_phyloseq, sample_data_phyloseq, tax_table_phyloseq)


#### STEP 2 - Differential abundance with DESeq2

#load DESeq2
BiocManager::install("DESeq2")
library(DESeq2)


#Aggregate data to Order Level
physeq_object_order <- tax_glom(physeq_object, taxrank = "Order")
tax_table(physeq_object_order)

#subset samples by time point
physeq_time0  <- subset_samples(physeq_object_order, Time =="0")
physeq_time1  <- subset_samples(physeq_object_order, Time =="1")
physeq_time4  <- subset_samples(physeq_object_order, Time =="4")


#now do a step-by-step pairwise comparison for each time point.

#Compare all groups to control at Time 1
# Relevel the treatment factor in your sample_data
sample_data(physeq_time1)$Treatment <- relevel(factor(sample_data(physeq_time1)$Treatment), ref = "N")
dds_time1 <- phyloseq_to_deseq2(physeq_time1, ~ Treatment)
dds_time1 <- DESeq(dds_time1)
#extract results for each treatment vs control
res_H <- results(dds_time1, contrast = c("Treatment", "H", "N"))
res_V <- results(dds_time1, contrast = c("Treatment", "V", "N"))
res_B <- results(dds_time1, contrast = c("Treatment", "B", "N"))
# Convert to data frame and add comparison label
res_H_df <- as.data.frame(res_H)
res_H_df$Comparison <- "H_vs_N"
res_H_df$ASV <- rownames(res_H_df)
res_V_df <- as.data.frame(res_V)
res_V_df$Comparison <- "V_vs_N"
res_V_df$ASV <- rownames(res_V_df)
res_B_df <- as.data.frame(res_B)
res_B_df$Comparison <- "B_vs_N"
res_B_df$ASV <- rownames(res_B_df)
#combine results and add taxonomy
combined_results <- rbind(res_H_df, res_V_df, res_B_df)
tax_df <- as.data.frame(tax_table(physeq_time1))
tax_df$ASV <- rownames(tax_df)
combined_results_tax <- merge(combined_results, tax_df, by = "ASV", all.x = TRUE)
#export just sig results
sig_combined <- subset(combined_results_tax, padj < 0.05)
write_xlsx(sig_combined, path = "Order_DESeq2_sig_results_time1.xlsx")

#Compare all groups to control at Time 4
# Relevel the treatment factor in your sample_data
sample_data(physeq_time4)$Treatment <- relevel(factor(sample_data(physeq_time4)$Treatment), ref = "N")
dds_time4 <- phyloseq_to_deseq2(physeq_time4, ~ Treatment)
dds_time4 <- DESeq(dds_time4)
#extract results for each treatment vs control
res_H <- results(dds_time4, contrast = c("Treatment", "H", "N"))
res_V <- results(dds_time4, contrast = c("Treatment", "V", "N"))
res_B <- results(dds_time4, contrast = c("Treatment", "B", "N"))
# Convert to data frame and add comparison label
res_H_df <- as.data.frame(res_H)
res_H_df$Comparison <- "H_vs_N"
res_H_df$ASV <- rownames(res_H_df)
res_V_df <- as.data.frame(res_V)
res_V_df$Comparison <- "V_vs_N"
res_V_df$ASV <- rownames(res_V_df)
res_B_df <- as.data.frame(res_B)
res_B_df$Comparison <- "B_vs_N"
res_B_df$ASV <- rownames(res_B_df)
#combine results and add taxonomy
combined_results <- rbind(res_H_df, res_V_df, res_B_df)
tax_df <- as.data.frame(tax_table(physeq_time4))
tax_df$ASV <- rownames(tax_df)
combined_results_tax <- merge(combined_results, tax_df, by = "ASV", all.x = TRUE)
#export just sig results
sig_combined <- subset(combined_results_tax, padj < 0.05)
write_xlsx(sig_combined, path = "Order_DESeq2_sig_results_time4.xlsx")

#Compare all groups to control at Time 0
# Relevel the treatment factor in your sample_data
sample_data(physeq_time0)$Treatment <- relevel(factor(sample_data(physeq_time0)$Treatment), ref = "N")
dds_time0 <- phyloseq_to_deseq2(physeq_time0, ~ Treatment)
dds_time0 <- DESeq(dds_time0)
#extract results for each treatment vs control
res_H <- results(dds_time0, contrast = c("Treatment", "H", "N"))
res_V <- results(dds_time0, contrast = c("Treatment", "V", "N"))
res_B <- results(dds_time0, contrast = c("Treatment", "B", "N"))
# Convert to data frame and add comparison label
res_H_df <- as.data.frame(res_H)
res_H_df$Comparison <- "H_vs_N"
res_H_df$ASV <- rownames(res_H_df)
res_V_df <- as.data.frame(res_V)
res_V_df$Comparison <- "V_vs_N"
res_V_df$ASV <- rownames(res_V_df)
res_B_df <- as.data.frame(res_B)
res_B_df$Comparison <- "B_vs_N"
res_B_df$ASV <- rownames(res_B_df)
#combine results and add taxonomy
combined_results <- rbind(res_H_df, res_V_df, res_B_df)
tax_df <- as.data.frame(tax_table(physeq_time0))
tax_df$ASV <- rownames(tax_df)
combined_results_tax <- merge(combined_results, tax_df, by = "ASV", all.x = TRUE)
#export just sig results
sig_combined <- subset(combined_results_tax, padj < 0.05)
write_xlsx(sig_combined, path = "Order_DESeq2_sig_results_time0.xlsx")


###now do an individual comparison between V and B at Time points 1 and 4

#Time 1
# Relevel the treatment factor in your sample_data
sample_data(physeq_time1)$Treatment <- relevel(factor(sample_data(physeq_time1)$Treatment), ref = "V")
dds_time1 <- phyloseq_to_deseq2(physeq_time1, ~ Treatment)
dds_time1 <- DESeq(dds_time1)
#extract results for each treatment vs control
res_B <- results(dds_time1, contrast = c("Treatment", "B", "V"))
# Convert to data frame and add comparison label
res_B_df <- as.data.frame(res_B)
res_B_df$Comparison <- "B_vs_V"
res_B_df$ASV <- rownames(res_B_df)
#combine results and add taxonomy
#combined_results <- rbind(res_H_df, res_V_df, res_B_df)
tax_df <- as.data.frame(tax_table(physeq_time1))
tax_df$ASV <- rownames(tax_df)
combined_results_tax <- merge(res_B_df, tax_df, by = "ASV", all.x = TRUE)
#export just sig results
sig_combined <- subset(combined_results_tax, padj < 0.05)
write_xlsx(sig_combined, path = "Order_DESeq2_sig_results_time1_VB.xlsx")

#Time 4
# Relevel the treatment factor in your sample_data
sample_data(physeq_time4)$Treatment <- relevel(factor(sample_data(physeq_time4)$Treatment), ref = "V")
dds_time4 <- phyloseq_to_deseq2(physeq_time4, ~ Treatment)
dds_time4 <- DESeq(dds_time4)
#extract results for each treatment vs control
res_B <- results(dds_time4, contrast = c("Treatment", "B", "V"))
# Convert to data frame and add comparison label
res_B_df <- as.data.frame(res_B)
res_B_df$Comparison <- "B_vs_V"
res_B_df$ASV <- rownames(res_B_df)
#combine results and add taxonomy
#combined_results <- rbind(res_H_df, res_V_df, res_B_df)
tax_df <- as.data.frame(tax_table(physeq_time4))
tax_df$ASV <- rownames(tax_df)
combined_results_tax <- merge(res_B_df, tax_df, by = "ASV", all.x = TRUE)
#export just sig results
sig_combined <- subset(combined_results_tax, padj < 0.05)
write_xlsx(sig_combined, path = "Order_DESeq2_sig_results_time4_VB.xlsx")


