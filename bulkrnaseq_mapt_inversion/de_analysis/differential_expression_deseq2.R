# Differential-expression analysis (DESeq2 + CQN normalization) for the brain
# bulk RNA-seq, run per ancestry/haplotype contrast. Produces the
# `all_genes_result*.csv` tables consumed by the figure scripts.
#
# Requires R with: biomaRt, DESeq2, cqn, limma, edgeR, variancePartition,
# corrplot (and their dependencies).

# ===============================================================
# CONFIG — edit these paths for your environment. Throughout the script,
# replace the "path/to/BulkRNASeq" prefix with the root of your data tree:
#   - setwd() target below: where outputs / .rda checkpoints are written
#   - the per-sample RSEM ".genes.results" directories (DLPFC / MFG / cau)
#   - the QC-metrics CSV and the sample-metadata CSV
# ===============================================================

library(biomaRt)

# Connect to Ensembl version 111 (provides GC content used for CQN).
message("Connecting to Ensembl version 111...")
mart <- useEnsembl(biomart = "ensembl",
                   dataset = "hsapiens_gene_ensembl",
                   version = 111)

# If you are behind a proxy, set it here (these were cluster-internal values;
# leave commented out on a normal network):
# Sys.setenv(http_proxy = "http://YOUR_PROXY:PORT")
# Sys.setenv(https_proxy = "http://YOUR_PROXY:PORT")

# Output / working directory (holds the .rda checkpoints and result CSVs).
setwd("path/to/BulkRNASeq/Brain_data/Routputs_align_partial_grch38_new_rsem/output")



##### Pipeline #####

### (0) Load RawData
### (1) Gene Filtering
### (2) Normalization
### (3) Outlier Removal

##### (0) Load RawData #####


#Read in datSeq_numeric
datSeq_numeric <- read.csv("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/file_metrics_v38_no_wasp_all_with_2_extra_MFG.csv")

# Install and load necessary library
library(tidyverse)


DLPFC <- c("2177-DLPFC", "22390-DLPFC", "HBCC-1111-DLPFC", "HBCC-1193-DLPFC", "HBCC-1314-DLPFC", "HBCC-1633-DLPFC", "HBCC-1753-DLPFC", "HBCC-1887-DLPFC", "HBCC-1987-DLPFC", "HBCC-2266-DLPFC", "HBCC-2542-DLPFC", "HBCC-2543-DLPFC", "HBCC-2791-DLPFC", "HBCC-2871-DLPFC", "HBCC-846-DLPFC", "HBCC-859-DLPFC", "HBCC-984-DLPFC")
DLPFC_1 <- c("NDBB010-2-1-DLPFC", "NDBB021-4-6-DLPFC", "NDBB023-2-2-DLPFC", "NDBB024-3-4-DLPFC", "NDBB040-5-3-DLPFC", "NDBB057-2-1-DLPFC", "NDBB096-1-4-DLPFC", "NDBB186-5-2-DLPFC", "NDBB265-2-5-DLPFC", "NDBB331-5-2-DLPFC", "NPBB028-DLPFC", "P2927z5L-DLPFC")
MFG <- c("105492-MFG", "112976-MFG", "2020-068-MFG", "2021-051-MFG", "2021-119-MFG", "271-MFG", "274-MFG", "331-MFG", "352-MFG", "451-MFG", "578-MFG")
MFG_1 <- c("63-MFG", "6936-MFG", "7126-MFG", "713-MFG", "7182-MFG", "7375-MFG", "HBCC-3134-MFG", "May-15-MFG", "Mayo-10-MFG", "Mayo-11-MFG", "Mayo-12-MFG")
MFG_1_mayo13 <- c("Mayo-13-MFG")
MFG_2 <- c("Mayo-14-MFG", "Mayo-16-MFG", "Mayo-2-MFG", "Mayo-3-MFG", "Mayo-4-MFG", "Mayo-5-MFG", "Mayo-6-MFG", "Mayo-7-MFG", "Mayo-8-MFG", "Mayo-9-MFG", "P2581-MFG", "P2595-MFG", "P2848-MFG", "T2213-MFG")
cau <- c("112976-02-cau", "112976-cau", "1992-030-cau", "1994-075-cau", "1995-078-cau", "1996-113-cau", "2004-052-cau", "2004-072-cau", "2009-025-cau", "2011-102-cau", "2015-105-cau", "2019-050-cau", "2020-054-1-cau")
cau_1 <- c("2020-054-2-cau-put", "2020-068-cau", "2021-051-cau", "2021-119-cau", "2177-cau", "22390-cau", "271-cau", "274-cau", "331-cau", "352-cau", "451-cau", "6210-cau", "63-cau", "6734-cau", "6904-cau", "6936-cau")
cau_2 <- c("7126-cau", "713-cau", "7182-cau", "7375-cau", "752-cau", "HBCC-1111-cau", "HBCC-1193-cau", "HBCC-1314-cau", "HBCC-1633-cau", "HBCC-1887-cau", "HBCC-1987-cau", "HBCC-2266-cau", "HBCC-2542-cau", "HBCC-2543-cau", "HBCC-2791-cau", "HBCC-2871-cau", "HBCC-3134-cau", "HBCC-846-cau", "HBCC-859-cau")
cau_3 <- c("Mayo-10-cau", "Mayo-11-cau", "Mayo-12-cau", "Mayo-13-cau", "Mayo-14-cau", "Mayo-15-cau", "Mayo-16-cau", "Mayo-2-cau", "Mayo-3-cau", "Mayo-4-cau", "Mayo-5-cau", "Mayo-6-cau", "Mayo-7-cau", "Mayo-8-cau", "Mayo-9-cau")
cau_4 <- c("NDBB010-7-5-cau", "NDBB021-4-5-cau", "NDBB023-7-8-cau", "NDBB024-5-6-cau", "NDBB028-5-8-cau", "NDBB040-cau", "NDBB057-6-8-cau", "NDBB096-6-2-cau", "NDBB186-6-4-cau", "NDBB265-4-4-cau", "NDBB331-7-4-cau", "P2581-cau", "P2595-cau", "P2848-cau", "P2927x6L-cau", "T2213-cau")
MFG_2_extra <- c("USP-114-MFG", "USP-2133-MFG")

# Define the new set of sample names
sample_names <- c(cau, cau_1, cau_2, cau_3, cau_4, MFG, MFG_1, MFG_2, MFG_1_mayo13, DLPFC, DLPFC_1, MFG_2_extra)


# Create an empty list for genes df
rsem_gene_dfs <- list()

# Loop through sample names to read each file
for (sample in sample_names) {
  # Determine the correct directory based on the sample name
  if (sample %in% DLPFC) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/DLPFC/rsem/", sample, ".genes.results")
  } else if (sample %in% DLPFC_1) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/DLPFC/rsem/", sample, ".genes.results")
  } else if (sample %in% MFG) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/MFG/rsem/", sample, ".genes.results")
  } else if (sample %in% MFG_1) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/MFG/rsem/", sample, ".genes.results")
  } else if (sample %in% MFG_2) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/MFG/rsem/", sample, ".genes.results")
  } else if (sample %in% MFG_1_mayo13) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/MFG/rsem/", sample, ".genes.results")
  } else if (sample %in% cau) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/cau/rsem/", sample, ".genes.results")
  } else if (sample %in% cau_1) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/cau/rsem/", sample, ".genes.results")
  } else if (sample %in% cau_2) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/cau/rsem/", sample, ".genes.results")
  } else if (sample %in% cau_3) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/cau/rsem/", sample, ".genes.results")
  } else if (sample %in% cau_4) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/cau/rsem/", sample, ".genes.results")
  } else if (sample %in% MFG_2_extra) {
    file_path <- paste0("path/to/BulkRNASeq/Brain_data/no_wasp_outputs/MFG/2_added_samples/rsem/", sample, ".genes.results")
  } else {
    next  # Skip this iteration if the sample doesn't match (optional, based on your logic)
  }
  rsem_gene_dfs[[sample]] <- read_tsv(file_path)
}

# Combine the individual data frames into one
combined_rsem_gene <- bind_rows(rsem_gene_dfs, .id = "sample_name")

## Only keep transcripts on chr1-Y
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype","percentage_gene_gc_content")

#ensembl_db = listDatasets(useEnsembl(biomart="ENSEMBL_MART_ENSEMBL")) # search for "hsapiens_gene_ensembl" got GRCh38
#mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
#mart <- useMart(host = "https://useast.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
#mart <- useMart(host = "https://www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

#combined_rsem_gene = combined_rsem_gene[idx,]

#expression_data <- combined_rsem_gene$expected_count
#cpm <- (expression_data / sum(expression_data)) * 1000000

# Selecting the columns
rsem_gene <- combined_rsem_gene %>%
  dplyr::select(sample_name, gene_id, expected_count) %>%
  tidyr::spread(key = sample_name, value = expected_count) %>%
  tibble::column_to_rownames(var = "gene_id") %>%
  as.data.frame()

# Selecting the columns
rsem_gene_effLen <- combined_rsem_gene %>%
  dplyr::select(sample_name, gene_id, effective_length) %>%
  tidyr::spread(key = sample_name, value = effective_length) %>%
  tibble::column_to_rownames(var = "gene_id")

# Convert to a base R data frame if you want to ensure it's not a tibble
rsem_gene_effLen <- as.data.frame(rsem_gene_effLen)
#only keep the average!
rsem_gene_effLen <- data.frame(average = rowMeans(rsem_gene_effLen, na.rm = TRUE), row.names = rownames(rsem_gene_effLen))
rsem_gene_effLen <- setNames(rsem_gene_effLen$average, rownames(rsem_gene_effLen))

#get gene id's in all samples
#gene_id_reduced <- data.frame(gene_id = unique(rownames(rsem_gene)))
#geneAnno1 <- getBM(attributes = getinfo, filters = "ensembl_gene_id", values = gene_id_reduced, mart = mart)

geneAnno1 <- getBM(attributes = getinfo, filters = c("ensembl_gene_id"), values = substr(rownames(rsem_gene),1,15), mart = mart)
geneAnno2 = geneAnno1[geneAnno1$chromosome_name %in% c(1:22, "X", "Y"),]
idx = which(substr(rownames(rsem_gene),1,15) %in% geneAnno2$ensembl_gene_id)
rsem_gene = rsem_gene[idx,]
rsem_gene_effLen = rsem_gene_effLen[idx]


#save(combined_rsem_gene, filtered, rsem_gene, rsem_gene_effLen, datSeq_numeric, geneAnno1, file = "3_01_RawData_chr1toY.rda")
save(combined_rsem_gene, rsem_gene, rsem_gene_effLen, datSeq_numeric, geneAnno1, file = "3_01_RawData_chr1toY_with_MFG_2_added_samples_Ensembl111.rda")

rm(list = ls())
lnames = load("3_01_RawData_chr1toY.rda")

# Apply CPM normalization to each sample (column-wise)
cpm <- apply(rsem_gene, 2, function(x) (x/sum(x))*1000000)

# Assign row names from the original data to the new cpm matrix
#rownames(cpm) <- rownames(rsem_genes)

df_idx = as.data.frame(matrix(ncol = 3, nrow = 0))
for (cpm_th in seq(0, 1, 0.1)) {
  for (x in seq(0, 0.5, 0.05)) {
    keep = rowSums(cpm > cpm_th) # Count genes above threshold per sample
    idx = which(keep > x * ncol(cpm)) # Genes above threshold in x% of samples
    df_idx = rbind(df_idx, c(cpm_th, x, length(idx)))
  }
}
colnames(df_idx) = c("cpm_th", "x_pctSamples", "n_genes_kept")
pdf("3_01_NumGenesKept_inXpctSamples_withDifferentCpmThreshold.pdf")
df_idx$cpm_th = factor(df_idx$cpm_th, levels = seq(0, 1, 0.1))
ggplot(df_idx, aes(x = x_pctSamples, y = n_genes_kept, col = cpm_th)) +
  geom_point() +
  geom_line()
dev.off()

# Observation: 30% of samples seem to be the general kneeing point. Line cpm_th = 0.3 seems to be stable for n_genes_kept. Also 18733 genes with CPM > 0.3 in at least 30% samples, consistent with Google's suggestion that human have 47k genes, 40% of which are expressed in neurons (~18800 genes).

#plot cpm per sample distribution - ADD
#cpm1 <- as.data.frame(cpm)

#ggplot(cpm1, aes(x='051091-Gln')) + 
#  geom_histogram(bins=30, fill="blue", alpha=0.7)

### Keep genes with cpm > 0.3 in 30% of samples
keep = apply(cpm > 0.3, 1, sum)
idx = which(keep > 0.30*dim(cpm)[2]) ## cpm > 0.3 in 30% of samples
cpm = cpm[idx,]

#Remove any remaining genes with eff gene length <=50bp
rsem_gene_effLen = rsem_gene_effLen[match(rownames(cpm), names(rsem_gene_effLen))]
stopifnot(rownames(cpm) == names(rsem_gene_effLen))

geneAnno <- geneAnno1[match(substr(rownames(cpm),1,15),geneAnno1[,1]),]
geneAnno <- geneAnno[!is.na(geneAnno[,1]),] # 18733 genes have ensembl_gene_id

# Make sure to include key antisense genes at the 17q21 locus
geneAnno$ensembl_gene_id[geneAnno$hgnc_symbol == "MAPT-AS1"] # ENSG00000264589
geneAnno$ensembl_gene_id[geneAnno$hgnc_symbol == "KANSL1-AS1"] # ENSG00000214401
idx1 = which(substr(names(rsem_gene_effLen), 1, 15) == "ENSG00000264589")
idx2 = which(substr(names(rsem_gene_effLen), 1, 15) == "ENSG00000214401")
rsem_gene_effLen[idx1] # 2809.77 #Sarah:2386.191  #YL: 1902.199 #combined: 2185.523
rsem_gene_effLen[idx2] # 98.61.  #Sarah:126.3027  #YL: 110.5579 #combined: 203.6872
hist(rsem_gene_effLen, xlim = c(0,8e3), breaks = 4000) # a bar every 200bp with breaks = 1000

idx = which(rsem_gene_effLen <= 50)
short_remove = rownames(cpm)[idx]
idx_remove = match(short_remove,rownames(cpm)) # 54 genes
cpm = cpm[-idx_remove,]
rsem_gene_effLen = rsem_gene_effLen[match(rownames(cpm), names(rsem_gene_effLen))]
rsem_gene = rsem_gene[match(rownames(cpm), rownames(rsem_gene)),]

dim(rsem_gene)

#RUN FROM JIANI TEST:
# Indices of genes containing '_PAR_Y' in their names
idx_remove = which(grepl("_PAR_Y", rownames(rsem_gene)))
# Gene names containing '_PAR_Y'
genes_with_PAR_Y = rownames(rsem_gene)[idx_remove]
# Print the gene names
print(genes_with_PAR_Y)
# Extract the first 15 characters from each gene name
first_15_chars = substr(genes_with_PAR_Y, start = 1, stop = 15)
# Print the shortened gene names
print(first_15_chars)
# Convert your vector of first 15 characters to full gene IDs as needed (if they are already full IDs, skip this step)
genes_to_remove = first_15_chars
# Remove the rows where the row names (gene IDs) are in the genes_to_remove list
rsem_gene_filtered = rsem_gene[!rownames(rsem_gene) %in% genes_to_remove, ]
# Check the first few rows of the updated dataframe
head(rsem_gene_filtered)

genes_to_remove = first_15_chars
rsem_gene = rsem_gene[!rownames(rsem_gene) %in% genes_to_remove, ]
rsem_gene_effLen = rsem_gene_effLen[!names(rsem_gene_effLen) %in% genes_to_remove]

idx_remove = which(grepl("_PAR_Y",rownames(rsem_gene))) # removed 13 genes
rsem_gene = rsem_gene[-idx_remove,]
rsem_gene_effLen = rsem_gene_effLen[-idx_remove]

geneAnno <- geneAnno[match(substr(rownames(rsem_gene),1,15),geneAnno[,1]),]
save(geneAnno, file = "3_01_BioMart_geneAnno.rda") # save for later use

rm(cpm,keep,short_remove,df_idx,idx,idx1,idx2,idx_remove,lnames,x)
save(list = ls(), file = "3_01_GeneFiltered_Data.rda")

##### (2) Normalization #####
rm(list = ls())
lnames = load("3_01_GeneFiltered_Data.rda")

library(cqn)
library(edgeR)

### Get CQN GC Content and Read Length Adjustment
geneAnno$effective_length = rsem_gene_effLen[match(geneAnno[,1],substr(names(rsem_gene_effLen),1,15))]
rownames(geneAnno) = geneAnno[,1]

cqn.dat <- cqn(rsem_gene, lengths = as.numeric(geneAnno$effective_length), x = as.numeric(geneAnno$percentage_gene_gc_content), lengthMethod = c("smooth"), sqn = FALSE) # run cqn with no quantile normalization
save(cqn.dat, file = "3_02_CQN_results.rda")

### limma-trend normalization, using the cqn offset values
dge2 <- DGEList(counts = rsem_gene)
dge2 <- calcNormFactors(dge2)
logCPM_offset <- cpm(dge2, log = TRUE, offset = cqn.dat$offset, prior.count = 0.5) # prior count helps to dampen variance of low count genes

datExpr = logCPM_offset
rsem_gene_effLen = rsem_gene_effLen[match(rownames(datExpr), names(rsem_gene_effLen))]

#Read in metadata
datMeta <- read.csv("path/to/BulkRNASeq/Ancestry_analysis/metadata_ancestry_astrocytes_fixed_removed_sample_050975.csv")

# List of cell lines to include
include_cell_lines <- c(
  "Astrocytes-050594", 
  "Astrocytes-050675", 
  "Astrocytes-050689", 
  "Astrocytes-050698", 
  #  "Astrocytes-050975", 
  "Astrocytes-051048", 
  "Astrocytes-051061", 
  "Astrocytes-051064", 
  "Astrocytes-051069", 
  "Astrocytes-051076", 
  "Astrocytes-051079", 
  "Astrocytes-051118", 
  "Astrocytes-051118-KRB", 
  "Astrocytes-051182", 
  "Astrocytes-12442", 
  "Astrocytes-12461", 
  "Astrocytes-12462-C7", 
  "Astrocytes-12468", 
  "Astrocytes-13505-1-N", 
  "Astrocytes-300-12-NWBB7", 
  "Astrocytes-BJSIPS", 
  "Astrocytes-ID6", 
  "Astrocytes-ID7", 
  "Astrocytes-iPS40", 
  "Astrocytes-MSN04", 
  "Astrocytes-WTC11")

# Filter the DataFrame
datSeq_numeric <- datSeq_numeric[datSeq_numeric$cell_line %in% include_cell_lines, ]

save(datExpr, rsem_gene_effLen, datMeta, datSeq_numeric, file="3_02_NormalizedGeneFilteredRNAseq.rda")

##### (3) Outlier Removal #####
rm(list = ls())
load("3_02_NormalizedGeneFilteredRNAseq.rda")

table(datMeta$Batch) 
table(datMeta$haplotype) 
table(datMeta[,c("Batch", "haplotype")])

### PCA-based Outlier Removal
norm <- t(scale(t(datExpr),scale=F))
PCA_res <- prcomp(norm,center=FALSE)
varexp <- (PCA_res$sdev)^2 / sum(PCA_res$sdev^2)
sum(varexp[1:12]) # First 12 PCs explain > 80% of total variance
nPC = 12
PCs = as.data.frame(PCA_res$rotation[,1:nPC])

PCs <- PCs[match(datMeta$File.Name_1, rownames(PCs)), ]

stopifnot(rownames(PCs) == datMeta$File.Name_1)
PCs$haplotype = datMeta$haplotype
PCs$Batch = datMeta$Batch
#PCs$diffbatch = datMeta$differentiationbatch
PCs$sex = datMeta$sex
PCs$apoe = datMeta$apoe
PCs$rin = datMeta$rin
PCs$Batch.Info = datMeta$Batch.Info
PCs$Ancestry = datMeta$Ancestry
PCs$read_depth = datSeq_numeric$TOTAL_READS/1e6

# Define axis text size
axis_text_size <- 10

# Create individual plots
p1 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = haplotype)) +
  geom_point(size = 1) + theme(axis.text = element_text(size = axis_text_size))

p2 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = Batch)) +
  geom_point(size = 1) + theme(axis.text = element_text(size = axis_text_size))

p3 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = sex)) +
  geom_point(size = 1) + theme(axis.text = element_text(size = axis_text_size))

p4 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = apoe)) +
  geom_point(size = 1) + theme(axis.text = element_text(size = axis_text_size))

p5 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = rin)) +
  geom_point(size = 1) +
  scale_color_distiller(palette = "RdYlBu") + theme(axis.text = element_text(size = axis_text_size))

p6 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = read_depth)) +
  geom_point(size = 1) +
  scale_color_distiller(palette = "RdYlBu") +
  labs(col = "read_depth\n(million)") + theme(axis.text = element_text(size = axis_text_size))

p8 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = Batch.Info)) +
  geom_point(size = 1) + theme(axis.text = element_text(size = axis_text_size))

p9 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = Ancestry)) +
  geom_point(size = 1) + theme(axis.text = element_text(size = axis_text_size))

# Save plots to PDF
pdf("3_03_PCA.pdf", height = 16, width = 24) # Adjusted size to accommodate plots
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p8, p9, ncol = 3) # Changed layout to 3 columns
dev.off()

topPC = PCs[,1:nPC]
colnames(topPC) <- paste("PC",c(1:nPC),"_",(signif(varexp[c(1:nPC)],2)*100),"%",sep="")
save(topPC, file = "3_03_topPCs.rda")

#################################
#################################
#################################
#################################
#################################

#PART 2:
# 4_1_PCA_build_linear_model_APOE33AsBaselineE2andE4asDosage_includeSex_VIF2p5.R

### Get top PCs that explain total 80% variance in the CQN normalized peak count (logRPKM), identify biological and technical covariates that significantly correlate with these top PCs.

# We include read depth into covariates.
# We use apoe33 as baseline, and use E2 and E4 dosage as covariates.
# Use VIF cutoff of 2.5

##### Pipeline #####

### 1) Load data and format
### 2) Load top PCs of datExpr
### 3) Calculate correlation between covariates and topPCs of datExpr
### 4) Select covariates for linear model

rm(list = ls())
load("3_02_NormalizedGeneFilteredRNAseq.rda")

datMeta$haplotype = factor(datMeta$haplotype, levels = c("H1H1", "H1H2"))
#datMeta$sex = factor(datMeta$sex, levels = c("M", "F"))
#colnames(datMeta) = gsub("differentiation", "", colnames(datMeta))
#datMeta$batch = factor(datMeta$batch) 
#datMeta$batch = factor(datMeta$batch) 

datMeta = as.data.frame(datMeta)
rownames(datMeta) = datMeta$File.Name_1
datMeta = datMeta[,-c(1)] 
datMeta = datMeta[,-c(1)] 


rownames(datSeq_numeric) = datSeq_numeric$cell_line
datSeq_numeric = datSeq_numeric[,-1]
colnames(datSeq_numeric) = gsub("_", ".", colnames(datSeq_numeric))
# Note that I'm using the raw picard metrics values rather than scaled ones.
# no need to scale the datSeq_numeric

# Reorder datMeta based on the order of datExpr column names
# Reorder the columns of datExpr to match the rownames of datMeta, preserving data integrity
datExpr <- datExpr[, rownames(datMeta)]
# Reorder rows of datSeq_numeric to match the order of colnames(datExpr)
datSeq_numeric <- datSeq_numeric[match(colnames(datExpr), rownames(datSeq_numeric)), ]

stopifnot(rownames(datMeta) == colnames(datExpr))
stopifnot(rownames(datSeq_numeric) == colnames(datExpr))

#Covariates = cbind(datMeta[,-c(3)], datSeq_numeric) # 19 samples with 49 variables,
Covariates = cbind(datMeta, datSeq_numeric) # 19 samples with 49 variables,

save(datMeta, datSeq_numeric, Covariates, datExpr, rsem_gene_effLen, file = "4_01_Data_with_Scaled_datMeta_datSeq.rda") # no need of rsem_gene_effLen later, but save it anyway

##### 2) Load top PCs of datExpr #####
lnames = load("3_03_topPCs.rda") # topPC

##### 3) Calculate correlation between covariates and topPCs of datExpr #####
# Generate the model matrix
mod_mat_expr = paste(c(colnames(Covariates)), collapse = " + ")
mod_mat_expr = paste0("~ ", mod_mat_expr)
mod_mat = model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates)[,-1] # Remove intercept
save(mod_mat, file = "4_03_mod_mat.rda")

# Assuming topPC is already defined and contains principal components

# Combine top principal components with covariates
mod_mat_withPC = cbind(topPC, mod_mat)

# Identify and remove constant columns (columns with zero variance)
constant_cols <- sapply(mod_mat_withPC, function(x) var(x, na.rm = TRUE) == 0)
constant_cols_names <- names(constant_cols[constant_cols])
print("Constant columns identified and removed:")
print(constant_cols_names)

mod_mat_withPC_non_constant <- mod_mat_withPC[, !constant_cols]

# Calculate Pearson and Spearman correlations
Cor_pearson = cor(mod_mat_withPC_non_constant)
Cor_spearman = cor(mod_mat_withPC_non_constant, method = "spearman")

# Print the comparison of correlations
compare_cor = data.frame(
  variable = colnames(mod_mat_withPC_non_constant),
  pearson = Cor_pearson[1,],  # Example for PC1
  spearman = Cor_spearman[1,]
)
print(compare_cor)

# Correcting the covariate names
correct_covariates_to_update <- c("rin", "batchb2", "batchb3", "genotypeWT", "donord2", "treatmentDMSO", "genotype_treatmentMT_DMSO", "genotype_treatmentWT_Anast", "genotype_treatmentWT_DMSO", "donor_genotype_treatmentd1_MT_DMSO", "donor_genotype_treatmentd1_WT_Anast", "donor_genotype_treatmentd1_WT_DMSO", "donor_genotype_treatmentd2_MT_Anast", "donor_genotype_treatmentd2_MT_DMSO", "donor_genotype_treatmentd2_WT_Anast", "donor_genotype_treatmentd2_WT_DMSO")

# Match corrected covariates
covariate_indices_corrected = match(correct_covariates_to_update, colnames(mod_mat_withPC_non_constant))

# Print corrected covariate indices to check for NA values
print(covariate_indices_corrected)

# Remove missing covariates from the list
correct_covariates_to_update <- correct_covariates_to_update[!is.na(covariate_indices_corrected)]
covariate_indices_corrected <- covariate_indices_corrected[!is.na(covariate_indices_corrected)]

# Adjust for principal components
idx_spearman = ncol(topPC) + covariate_indices_corrected

# Update the correlation matrix with Spearman correlations for the chosen covariates
Cor = Cor_pearson # start with Pearson correlation matrix
Cor[idx_spearman,] = Cor_spearman[idx_spearman,]
Cor[,idx_spearman] = Cor_spearman[,idx_spearman]

# Save the updated correlation matrix
save(Cor, file = "updated_correlation_matrix.rda")

#Cor = cor(mod_mat_withPC)
#Cor_spearman = cor(mod_mat_withPC, method = "spearman")
#colnames(Cor)
#ncol(topPC) # 8, correct.
#idx_spearman = ncol(topPC) + c(1:2,4:8)
#Cor[idx_spearman,] = Cor_spearman[idx_spearman,]; Cor[,idx_spearman] = Cor_spearman[,idx_spearman]

## Find out which Covariates significantly correlate with the topPCs (like what Yuyan described in Feng 2022 BioRxiv)

# Calculate correlation significance p-values
Cor_sig = matrix(nrow = nrow(Cor), ncol = ncol(Cor))
rownames(Cor_sig) = colnames(Cor_sig) = colnames(Cor)

# Adjust the loop to use mod_mat_withPC_non_constant for correct dimensions
for (i in 1:ncol(mod_mat_withPC_non_constant)) {
  for (j in 1:ncol(mod_mat_withPC_non_constant)) {
    tmp = cor.test(mod_mat_withPC_non_constant[,i], mod_mat_withPC_non_constant[,j])
    Cor_sig[i,j] = tmp$p.value
    # if (tmp$p.value < 0.05) {print(paste(i, colnames(Cor)[i], j, colnames(Cor)[j]))}
  }
}

# Save the significance matrix
save(Cor_sig, file = "correlation_significance_matrix.rda")

# Correlation significance p-values
#Cor_sig = matrix(nrow = nrow(Cor), ncol = ncol(Cor))
#rownames(Cor_sig) = colnames(Cor_sig) = colnames(Cor)
#for (i in 1:ncol(mod_mat_withPC)) {
#  for (j in 1:ncol(mod_mat_withPC)) {
#    tmp = cor.test(mod_mat_withPC[,i], mod_mat_withPC[,j])
#    Cor_sig[i,j] = tmp$p.value
#    # if (tmp$p.value < 0.05) {print(paste(i, colnames(Cor)[i], j, colnames(Cor)[j]))}
#  }
#}
# View(Cor_sig)
# Observation: Haplotype highly correlates with PC7 and PC3 (cor_pval = 0.0025 and 0.14 respectively)


# Correlation significance FDR
Cor_fdr = matrix(p.adjust(Cor_sig, method = "fdr"), nrow = nrow(Cor_sig))
rownames(Cor_fdr) = colnames(Cor_fdr) = colnames(Cor_sig)
for (i in 1:ncol(mod_mat_withPC)) {
  Cor_fdr[i,i] = 1
}
# View(Cor_fdr)
# Observation: Haplotype significantly correlates with PC7 (cor_fdr = 0.02)

### Get candidate covariates by Cor_FDR or Cor(_coef)
# Option a) Criteria: FDR < 0.05
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.05)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx])
#"haplotypeH2H2"    "PCT.USABLE.BASES" "GC.NC.0.19" 
#YL:  #[1] "haplotypeH2H2"        "APOE4"                "PCT.INTERGENIC.BASES"  "PCT.RIBOSOMAL.BASES" 

# Option b) Criteria: |cor| > 0.4 [choose this in the end]
Cov_pool_idx = which(apply(Cor[1:ncol(topPC), (ncol(topPC) + 1) : ncol(Cor)], 2, function(x) any(abs(x) > 0.4)))
nPC = 8
(Cov_pool = colnames(Cor)[Cov_pool_idx + nPC]) # 43 columns now, 5 were removed
# Best to include as many as possible, exclude the ones that will make VIF > 2.5. Choose b) |cor| > 0.4

### For the easiness of selecting covariates, plot corrplot by ranking covs that show max sum of weighted correlation with the topPCs.

# Get the variance explained by each PC.
colnames(Cor)[1:nPC] # see the variance explained by each PC
VarExp = sapply(colnames(Cor)[1:nPC], function(x) unlist(str_split(x, "_"))[2])
VarExp = gsub("%", "", VarExp)
VarExp = as.numeric(VarExp)/100

Cov_pool
VarExp_idx = 1:nPC
Cov_select_df = data.frame(Cov_pool = Cov_pool, sumCorTimesVarExp = sapply(Cov_pool, function(x) t(abs(Cor[VarExp_idx,which(colnames(Cor) == x)])) %*% VarExp[VarExp_idx]))
Cov_select_df = Cov_select_df %>%
  arrange(desc(sumCorTimesVarExp))
# Note that haplotypeH2H2 doesn't come the first, but that's normal and fine.

# Corrplot of ranked candiate covariates (Cov_pool) with top PCs
pdf(paste0("4_03_Corrplot_topPCofLimmaTrendCQNnormedlogRPKM_Covariates_CovPoolOnly.pdf"), height = 35, width = 45)
corrplot(Cor[c(colnames(Cor)[1:ncol(topPC)], Cov_select_df$Cov_pool), c(colnames(Cor)[1:ncol(topPC)], Cov_select_df$Cov_pool)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey", tl.cex = 1.7, cl.cex = 1.7, number.cex = 1.7)
dev.off()
# Observation: See large blocks of colinear covariates

##### 4) Select covariates for linear model #####
Cov_pool
#[1] "haplotypeH2H2"                  "sexF"                           "rin"                            "batchb2"                        "batchb3"                       
#[6] "APOE2"                          "APOE3"                          "APOE4"                          "RIBOSOMAL.BASES"                "TOTAL.READS"                   
#[11] "CORRECT.STRAND.READS"           "ESTIMATED.LIBRARY.SIZE"         "INTERGENIC.BASES"               "INTRONIC.BASES"                 "MEDIAN.3PRIME.BIAS"            
#[16] "MEDIAN.5PRIME.BIAS"             "NUM.R1.TRANSCRIPT.STRAND.READS" "NUM.R2.TRANSCRIPT.STRAND.READS" "NUM.UNEXPLAINED.READS"          "PCT.CORRECT.STRAND.READS"      
#[21] "PCT.INTERGENIC.BASES"           "PCT.INTRONIC.BASES"             "PCT.MRNA.BASES"                 "PCT.R1.TRANSCRIPT.STRAND.READS" "PCT.R2.TRANSCRIPT.STRAND.READS"
#[26] "PCT.RIBOSOMAL.BASES"            "PCT.USABLE.BASES"               "PCT.UTR.BASES"                  "PERCENT.DUPLICATION"            "PF.BASES"                      
#[31] "READ.PAIRS.EXAMINED"            "READ.PAIR.DUPLICATES"           "PF.ALIGNED.BASES"               "STRAND.BALANCE"                 "PF.READS"                      
#[36] "READ.PAIR.OPTICAL.DUPLICATES"   "UNMAPPED.READS"                 "TOTAL.CLUSTERS"                 "ALIGNED.READS"                  "AT.DROPOUT"                    
#[41] "GC.NC.0.19"                     "GC.NC.20.39"                    "GC.NC.60.79"  

#YL:
#Cov_pool
# [1] "haplotypeH2H2"                "sexF"                        
# [3] "rin"                          "batchb2"                     
# [5] "batchb3"                      "APOE2"                       
# [7] "APOE3"                        "APOE4"                       
# [9] "RIBOSOMAL.BASES"              "TOTAL.READS"                 
#[11] "ESTIMATED.LIBRARY.SIZE"       "INTERGENIC.BASES"            
#[13] "PCT.INTERGENIC.BASES"         "PCT.RIBOSOMAL.BASES"         
#[15] "PERCENT.DUPLICATION"          "PF.BASES"                    
#[17] "READ.PAIRS.EXAMINED"          "READ.PAIR.DUPLICATES"        
#[19] "PF.ALIGNED.BASES"             "STRAND.BALANCE"              
#[21] "PF.READS"                     "READ.PAIR.OPTICAL.DUPLICATES"
#[23] "UNMAPPED.READS"               "TOTAL.CLUSTERS"              
#[25] "ALIGNED.READS"                "AT.DROPOUT"                  
#[27] "GC.DROPOUT"                   "GC.NC.0.19"                  
#[29] "GC.NC.20.39"                  "GC.NC.40.59"                 
#[31] "GC.NC.60.79"                  "GC.NC.80.100"   
### Test whether the full model violates rule "all VIF < 2.5"

## Full model
colnames(mod_mat)
logRPKM = datExpr[1,]
cur_data = as.data.frame(cbind(logRPKM, mod_mat))
expression_lm = paste0("lm(logRPKM ~ ", paste(Cov_pool, collapse = " + "), ", data = cur_data)")
fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction))
# All covariates get VIF = Inf, too many covariates - overfit

## Based on the PCA plot (3_03_PCA.plot), include Haplotype, sex, APOE2 and APOE4 first, and then go down the columns to ensure all VIF < 2.5

## step1. set up base model
expression_lm = "lm(logRPKM ~ haplotypeH2H2 + APOE2 + APOE4, data = cur_data)"
#haplotypeH2H2 + sexF 
#haplotypeH2H2 + sexF + APOE2 + APOE4
#haplotypeH2H2 + sexF  + BatchB

fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction)) # Google search "vif cutoff": The default VIF cutoff value is 5
#Variables Tolerance      VIF
#1 haplotypeH2H2 0.7439965 1.344092
#2          sexF 0.7079753 1.412479
#3         APOE2 0.6631577 1.507937
#4         APOE4 0.6886505 1.452115
#5 PCT.UTR.BASES 0.6839447 1.462106

#YL:
#      Variables Tolerance      VIF
#1 haplotypeH2H2 0.7440849 1.343933
#2          sexF 0.7936906 1.259937
#3         APOE2 0.8757166 1.141922
#4         APOE4 0.7837999 1.275836

## step2. iterate to add more covariates
tmp = unlist(strsplit(expression_lm, split = ",|~"))[2]
tmp2 = gsub(" ", "", tmp)
cur_covs = unlist(strsplit(tmp2, "\\+"))
all(cur_covs %in% Cov_select_df$Cov_pool) # T

for (i in 1:nrow(Cov_select_df)) {
  if (! Cov_select_df$Cov_pool[i] %in% cur_covs) {
    print(paste0("i=", i, ", Try add covariate ", Cov_select_df$Cov_pool[i]))
    exp1 = unlist(strsplit(expression_lm, ","))[1]
    exp2 = unlist(strsplit(expression_lm, ","))[2]
    expA = paste0(" + ", Cov_select_df$Cov_pool[i], ",")
    test_lm = paste0(exp1, expA, exp2)
    
    fit_infunction <- eval(parse(text = test_lm))
    vif_df_infunction = ols_vif_tol(fit_infunction)
    
    if (all(vif_df_infunction$VIF < 2.5)) {
      expression_lm = test_lm
      print("covariate added")
      print(expression_lm)
      i = i+1
    } else {
      print("covariate denied")
      i = i+1
    }
  }
}

## Final model:
expression_lm
#[1] "lm(logRPKM ~ haplotypeH2H2 + sexF + APOE2 + APOE4 + PCT.UTR.BASES + PCT.MRNA.BASES + ESTIMATED.LIBRARY.SIZE + MEDIAN.3PRIME.BIAS, data = cur_data)"
#YL: [1] "lm(logRPKM ~ haplotypeH2H2 + sexF + APOE2 + APOE4 + PERCENT.DUPLICATION + GC.NC.20.39 + GC.NC.80.100 + batchb2, data = cur_data)"
save(expression_lm, file = "4_04_final_lm_expr.rda")

##############
##############
##############

#PART 3

# 4_2_DEseq2_using_established_lm.R

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)

##### Pipeline #####

### (1) Load linear model and data
### (2) Run DESeq2
### (3) Histogram, MA plot, and Dispersion plot

##### (1) Load linear model and data #####
## a. linear model
#lnames = load("4_04_final_lm_expr.rda") # expression_lm

## b. count data
lnames = load("3_01_GeneFiltered_Data.rda") # datMeta, datSeq_numeric, rsem_gene, rsem_gene_effLen, geneAnno, geneAnno1, cpm_th

## c. meta data
lnames = load("4_01_Data_with_Scaled_datMeta_datSeq.rda") # datMeta, datSeq_numeric, Covariates, datExpr, rsem_gene_effLen

# Reorder datMeta based on the order of datExpr column names
# Reorder the columns of datExpr to match the rownames of datMeta, preserving data integrity
rsem_gene <- rsem_gene[, rownames(datMeta)]

# Reorder rows of datSeq_numeric to match the order of colnames(datExpr)
datSeq_numeric <- datSeq_numeric[match(colnames(rsem_gene), rownames(datSeq_numeric)), ]

stopifnot(rownames(datExpr) == rownames(rsem_gene))
#rsem_gene = rsem_gene[match(rownames(datExpr), rownames(rsem_gene)),]
stopifnot(rownames(datMeta) == colnames(rsem_gene))
stopifnot(rownames(datSeq_numeric) == colnames(rsem_gene))

## c. load CQN results
lnames = load("3_02_CQN_results.rda") # cqn.dat
str(cqn.dat$glm.offset)
stopifnot(rownames(cqn.dat$glm.offset) == rownames(rsem_gene))

# Check the Covariates data frame
head(Covariates)

Covariates$haplotype <- as.factor(Covariates$haplotype)
#Covariates$sex <- as.factor(Covariates$sex)
#Covariates$apoe <- as.factor(Covariates$apoe)
Covariates$Ancestry <- as.factor(Covariates$Ancestry)
Covariates$sex <- as.factor(Covariates$sex)
Covariates$apoe <- as.factor(Covariates$apoe)

# Load DESeq2 library
library(DESeq2)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = round(rsem_gene, 0),
  colData = Covariates,
  design = ~ haplotype + Ancestry + sex + apoe 
  #design = ~ haplotype + Ancestry + sex + apoe + haplotype:Ancestry
  
)



# Perform differential expression analysis
dds <- DESeq(dds)

## Integrate CQN offset
cqnNormFactors <- exp(cqn.dat$glm.offset)
normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))
normalizationFactors(dds) <- normFactors

## Run DESeq2
dds <- DESeq(dds)
#using pre-existing normalization factors
#estimating dispersions
#found already estimated dispersions, replacing these
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#1244 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest


resultsNames(dds) # lists the coefficients: [1] "Intercept"              "haplotype_H2H2_vs_H1H1" "sex_F_vs_M"             "APOE2"                 [5] "APOE4"                  "PCT_UTR_BASES"          "STRAND_BALANCE"         "GC_NC_20_39"

#now for H1H2 vs. H1H1 

res <- results(dds, name="haplotype_H1H2_vs_H1H1")
save(dds, res, file = "5_02_DESeq2_result_haplotype_ancestry_sex_apoe_res_is_haplotype_H1H2_vs_H1H1.rda") # 4_2_02_DESeq2_resultHaplotype_lm_HaplotypeSexAPOE2APOE4PCT.UTR.BASESGC.NC.60.79MEDIAN.5PRIME.TO.3PRIME.BIAS.rda


## Any significant genes?
formatC(range(res$padj, na.rm = T), digits = 2) # new 4.1e-23 to 1; old 2.8e-16 to 1
#sarahs: [1] "2.3e-36" "  1"    

summary(res)
# out of 18666 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 183, 0.98%
# LFC < 0 (down)     : 159, 0.85%
# outliers [1]       : 0, 0%
# low counts [2]     : 2895, 16%
# (mean count < 20)

# Order by pval
resOrdered <- res[order(res$pvalue),]

##### 3) Histogram, MA plot, and Dispersion plot #####
pdf("5_03_Histogram_wCQNglmoffset_haplotype_ancestry_sex_apoe_res_is_haplotype_H1H2_vs_H1H1.pdf") # 4_2_03_Histogram_wCQNglmoffset.pdf
hist(res$pvalue) # More significant than lm()
dev.off()
# very significant at p-val < 0.05

pdf("5_03_MA_and_Dispersion_plots_wCQNglmoffset_haplotype_ancestry_sex_apoe_res_is_haplotype_H1H2_vs_H1H1.pdf")  # 4_2_03_MA_and_Dispersion_plots_wCQNglmoffset.pdf
plotMA(res, ylim=c(-2,2)) # looks good
plotDispEsts(dds)
dev.off()
# looks good


# Filter for down-regulated genes in H2 vs H1
sig_down_in_H2 <- subset(res, log2FoldChange < 1 & pvalue < 0.05)
sig_gene_ids_down_in_H2 <- rownames(sig_down_in_H2)
write.csv(sig_gene_ids_down_in_H2, "sig_gene_ids_down_in_H1H2_haplotype_ancestry_sex_apoe_res_is_haplotype_H1H2_vs_H1H1_pval0.05_logfc1.csv", row.names = FALSE)

# Filter for up-regulated genes in H2 vs H1
sig_up_in_H2 <- subset(res, log2FoldChange > 1 & pvalue < 0.05)
sig_gene_ids_up_in_H2 <- rownames(sig_up_in_H2)
write.csv(sig_gene_ids_up_in_H2, "sig_gene_ids_up_in_H1H2_haplotype_ancestry_sex_apoe_res_is_haplotype_H1H2_vs_H1H1_pval0.05_logfc1.csv", row.names = FALSE)

# Get all info including p-value, baseMean, log2FC, etc. for down-regulated genes
sig_down_in_H2$gene_id <- rownames(sig_down_in_H2)
write.csv(sig_down_in_H2, "sig_down_in_H1H2_haplotype_ancestry_sex_apoe_res_is_haplotype_H1H2_vs_H1H1_pval0.05_logfc1.csv", row.names = FALSE)

# Get all info including p-value, baseMean, log2FC, etc. for up-regulated genes
sig_up_in_H2$gene_id <- rownames(sig_up_in_H2)
write.csv(sig_up_in_H2, "sig_up_in_H1H2_haplotype_ancestry_sex_apoe_res_is_haplotype_H1H2_vs_H1H1_pval0.05_logfc1.csv", row.names = FALSE)

#to print out all genes: 
# Add gene IDs as a new column
res$gene_id <- rownames(res)
# Output all data to CSV
write.csv(res, "all_genes_results_haplotype_ancestry_sex_apoe_res_is_haplotype_H1H2_vs_H1H1.csv", row.names = FALSE)


#now for EA vs AA

res <- results(dds, name="Ancestry_EA_vs_AA")
save(dds, res, file = "5_02_DESeq2_result_haplotype_ancestry_sex_apoe_res_is_Ancestry_EA_vs_AA.rda") # 4_2_02_DESeq2_resultHaplotype_lm_HaplotypeSexAPOE2APOE4PCT.UTR.BASESGC.NC.60.79MEDIAN.5PRIME.TO.3PRIME.BIAS.rda

## Any significant genes?
formatC(range(res$padj, na.rm = T), digits = 2) # new 4.1e-23 to 1; old 2.8e-16 to 1
#sarahs: [1] "2.3e-36" "  1"    

summary(res)
# out of 18666 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 183, 0.98%
# LFC < 0 (down)     : 159, 0.85%
# outliers [1]       : 0, 0%
# low counts [2]     : 2895, 16%
# (mean count < 20)

# Order by pval
resOrdered <- res[order(res$pvalue),]

##### 3) Histogram, MA plot, and Dispersion plot #####
pdf("5_03_Histogram_wCQNglmoffset_haplotype_ancestry_sex_apoe_res_is_Ancestry_EA_vs_AA.pdf") # 4_2_03_Histogram_wCQNglmoffset.pdf
hist(res$pvalue) # More significant than lm()
dev.off()
# very significant at p-val < 0.05

pdf("5_03_MA_and_Dispersion_plots_wCQNglmoffset_haplotype_ancestry_sex_apoe_res_is_Ancestry_EA_vs_AA.pdf")  # 4_2_03_MA_and_Dispersion_plots_wCQNglmoffset.pdf
plotMA(res, ylim=c(-2,2)) # looks good
plotDispEsts(dds)
dev.off()
# looks good


# Filter for down-regulated genes in H2 vs H1
sig_down_in_H2 <- subset(res, log2FoldChange < 1 & pvalue < 0.05)
sig_gene_ids_down_in_H2 <- rownames(sig_down_in_H2)
write.csv(sig_gene_ids_down_in_H2, "sig_gene_ids_down_in_EA_haplotype_ancestry_sex_apoe_res_is_Ancestry_EA_vs_AA_pval0.05_logfc1.csv", row.names = FALSE)

# Filter for up-regulated genes in H2 vs H1
sig_up_in_H2 <- subset(res, log2FoldChange > 1 & pvalue < 0.05)
sig_gene_ids_up_in_H2 <- rownames(sig_up_in_H2)
write.csv(sig_gene_ids_up_in_H2, "sig_gene_ids_up_in_EA_haplotype_ancestry_sex_apoe_res_is_Ancestry_EA_vs_AA_pval0.05_logfc1.csv", row.names = FALSE)

# Get all info including p-value, baseMean, log2FC, etc. for down-regulated genes
sig_down_in_H2$gene_id <- rownames(sig_down_in_H2)
write.csv(sig_down_in_H2, "sig_down_in_EA_haplotype_ancestry_sex_apoe_res_is_Ancestry_EA_vs_AA_pval0.05_logfc1.csv", row.names = FALSE)

# Get all info including p-value, baseMean, log2FC, etc. for up-regulated genes
sig_up_in_H2$gene_id <- rownames(sig_up_in_H2)
write.csv(sig_up_in_H2, "sig_up_in_EA_haplotype_ancestry_sex_apoe_res_is_Ancestry_EA_vs_AA_pval0.05_logfc1.csv", row.names = FALSE)

#to print out all genes: 
# Add gene IDs as a new column
res$gene_id <- rownames(res)
# Output all data to CSV
write.csv(res, "all_genes_results_haplotype_ancestry_sex_apoe_res_is_Ancestry_EA_vs_AA.csv", row.names = FALSE)


#NOW FOR EUR VS. AA

res <- results(dds, name="Ancestry_EUR_vs_AA")
save(dds, res, file = "5_02_DESeq2_result_haplotype_ancestry_sex_apoe_res_is_Ancestry_EUR_vs_AA.rda") # 4_2_02_DESeq2_resultHaplotype_lm_HaplotypeSexAPOE2APOE4PCT.UTR.BASESGC.NC.60.79MEDIAN.5PRIME.TO.3PRIME.BIAS.rda

## Any significant genes?
formatC(range(res$padj, na.rm = T), digits = 2) # new 4.1e-23 to 1; old 2.8e-16 to 1
#sarahs: [1] "2.3e-36" "  1"    

summary(res)
# out of 18666 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 183, 0.98%
# LFC < 0 (down)     : 159, 0.85%
# outliers [1]       : 0, 0%
# low counts [2]     : 2895, 16%
# (mean count < 20)

# Order by pval
resOrdered <- res[order(res$pvalue),]

##### 3) Histogram, MA plot, and Dispersion plot #####
pdf("5_03_Histogram_wCQNglmoffset_haplotype_ancestry_sex_apoe_res_is_Ancestry_EUR_vs_AA.pdf") # 4_2_03_Histogram_wCQNglmoffset.pdf
hist(res$pvalue) # More significant than lm()
dev.off()
# very significant at p-val < 0.05

pdf("5_03_MA_and_Dispersion_plots_wCQNglmoffset_haplotype_ancestry_sex_apoe_res_is_Ancestry_EUR_vs_AA.pdf")  # 4_2_03_MA_and_Dispersion_plots_wCQNglmoffset.pdf
plotMA(res, ylim=c(-2,2)) # looks good
plotDispEsts(dds)
dev.off()
# looks good


# Filter for down-regulated genes in H2 vs H1
sig_down_in_H2 <- subset(res, log2FoldChange < 1 & pvalue < 0.05)
sig_gene_ids_down_in_H2 <- rownames(sig_down_in_H2)
write.csv(sig_gene_ids_down_in_H2, "sig_gene_ids_down_in_EUR_haplotype_ancestry_sex_apoe_res_is_Ancestry_EUR_vs_AA_pval0.05_logfc1.csv", row.names = FALSE)

# Filter for up-regulated genes in H2 vs H1
sig_up_in_H2 <- subset(res, log2FoldChange > 1 & pvalue < 0.05)
sig_gene_ids_up_in_H2 <- rownames(sig_up_in_H2)
write.csv(sig_gene_ids_up_in_H2, "sig_gene_ids_up_in_EUR_haplotype_ancestry_sex_apoe_res_is_Ancestry_EUR_vs_AA_pval0.05_logfc1.csv", row.names = FALSE)

# Get all info including p-value, baseMean, log2FC, etc. for down-regulated genes
sig_down_in_H2$gene_id <- rownames(sig_down_in_H2)
write.csv(sig_down_in_H2, "sig_down_in_EUR_haplotype_ancestry_sex_apoe_res_is_Ancestry_EUR_vs_AA_pval0.05_logfc1.csv", row.names = FALSE)

# Get all info including p-value, baseMean, log2FC, etc. for up-regulated genes
sig_up_in_H2$gene_id <- rownames(sig_up_in_H2)
write.csv(sig_up_in_H2, "sig_up_in_EUR_haplotype_ancestry_sex_apoe_res_is_Ancestry_EUR_vs_AA_pval0.05_logfc1.csv", row.names = FALSE)

#to print out all genes: 
# Add gene IDs as a new column
res$gene_id <- rownames(res)
# Output all data to CSV
write.csv(res, "all_genes_results_haplotype_ancestry_sex_apoe_res_is_Ancestry_EUR_vs_AA.csv", row.names = FALSE)




