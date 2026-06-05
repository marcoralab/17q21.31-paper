# Load necessary libraries
library(tidyverse)
library(limma)
library(edgeR)
library(ggplot2)
library(biomaRt)
library(sva)
library(RColorBrewer)
library(dplyr)

rm(list = ls())

# Remove version numbers from Ensembl IDs
remove_version <- function(ensembl_ids) {
  sapply(ensembl_ids, function(x) strsplit(x, "\\.")[[1]][1])
}

# ===============================================================
# CONFIG — edit these paths for your environment.
# ===============================================================
# This study's CQN-normalized expression objects, one .rda per cell type
# (produced by the gene-expression filtering / CQN normalization step).
# Each file restores: datExpr_* (genes x samples, normalized) and
# datMeta_* (sample metadata with `File.Name_1`, `haplotype`, `Ancestry`).
ASTRO_RDA     <- "data/astro_NormalizedGeneFilteredRNAseq.rda"      # datExpr_astro,  datMeta_astro
NEURON_RDA    <- "data/neuron_NormalizedGeneFilteredRNAseq.rda"     # datExpr_neuron, datMeta_neuron
MICROGLIA_RDA <- "data/microglia_NormalizedGeneFilteredRNAseq.rda"  # datExpr,        datMeta

# Public comparison datasets (download from GEO). Counts/expression matrices:
GSE260517_FILE <- "data/GSE260517_star_gene.txt"                   # Bowles 2024
GSE190185_FILE <- "data/GSE190185_isogenicAPOE_featureCounts.txt"  # TCW 2022
GSE73721_FILE  <- "data/GSE73721_Human_and_mouse_table_with_ensembl_id.csv"  # Zhang 2016
GSE182307_FILE <- "data/GSE182307_geneLevelCounts_allSamples.csv"  # Leng 2022
GSE97904_FILE  <- "data/GSE97904_expression.tsv"

# Where plots are written.
OUTPUT_DIR <- "results/plots"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
# ===============================================================

# ASTROCYTE DATA:
# Restores datExpr_astro and datMeta_astro
load(ASTRO_RDA)
h1h1_afr_samples_astro <- datMeta_astro$File.Name_1[datMeta_astro$haplotype == "H1H1" & datMeta_astro$Ancestry == "AA"]
h1h2_afr_samples_astro <- datMeta_astro$File.Name_1[datMeta_astro$haplotype == "H1H2" & datMeta_astro$Ancestry == "AA"]
h1h1a_afr_samples_astro <- datMeta_astro$File.Name_1[datMeta_astro$haplotype == "H1H1_A" & datMeta_astro$Ancestry == "AA"]
h1h1_eur_samples_astro <- datMeta_astro$File.Name_1[datMeta_astro$haplotype == "H1H1" & datMeta_astro$Ancestry == "EUR"]
h1h2_eur_samples_astro <- datMeta_astro$File.Name_1[datMeta_astro$haplotype == "H1H2" & datMeta_astro$Ancestry == "EUR"]
rownames(datExpr_astro) <- remove_version(rownames(datExpr_astro))

# NEURON DATA:
# Restores datExpr_neuron and datMeta_neuron
load(NEURON_RDA)
h1h1_afr_samples_neuron <- datMeta_neuron$File.Name_1[datMeta_neuron$haplotype == "H1H1" & datMeta_neuron$Ancestry == "AA"]
h1h2_afr_samples_neuron <- datMeta_neuron$File.Name_1[datMeta_neuron$haplotype == "H1H2" & datMeta_neuron$Ancestry == "AA"]
h1h1a_afr_samples_neuron <- datMeta_neuron$File.Name_1[datMeta_neuron$haplotype == "H1H1_A" & datMeta_neuron$Ancestry == "AA"]
h1h1_eur_samples_neuron <- datMeta_neuron$File.Name_1[datMeta_neuron$haplotype == "H1H1" & datMeta_neuron$Ancestry == "EUR"]
h1h2_eur_samples_neuron <- datMeta_neuron$File.Name_1[datMeta_neuron$haplotype == "H1H2" & datMeta_neuron$Ancestry == "EUR"]
rownames(datExpr_neuron) <- remove_version(rownames(datExpr_neuron))

# MICROGLIA DATA:
# Restores datExpr and datMeta. NOTE: microglia samples are not plotted, but
# this dataset is intentionally kept in the common-gene intersection and the
# combined matrix below so the gene set (and therefore the PCA) matches the
# published figure.
load(MICROGLIA_RDA)
h1h1_afr_samples_microglia <- datMeta$File.Name_1[datMeta$haplotype == "H1H1" & datMeta$Ancestry == "AA"]
h1h2_afr_samples_microglia <- datMeta$File.Name_1[datMeta$haplotype == "H1H2" & datMeta$Ancestry == "AA"]
#h1h1a_afr_samples_microglia <- datMeta$File.Name_1[datMeta$haplotype == "H1H1_A" & datMeta$Ancestry == "AA"]
h1h1_eur_samples_microglia <- datMeta$File.Name_1[datMeta$haplotype == "H1H1" & datMeta$Ancestry == "EUR"]
h1h2_eur_samples_microglia <- datMeta$File.Name_1[datMeta$haplotype == "H1H2" & datMeta$Ancestry == "EUR"]
rownames(datExpr) <- remove_version(rownames(datExpr))


# Load GSE260517 data
# NEW CODE FOR GSE260517
# Load GSE260517 data
gse260517_data <- read.delim(GSE260517_FILE, header = TRUE, stringsAsFactors = FALSE)
colnames(gse260517_data)[1] <- "gene_id"
selected_samples_gse260517 <- c("gene_id", "S75_11_IW1A12_AST", "S75_11_IW1F3_AST", "S75_11_IW1G8_AST", 
                                "S300_12_NWBB7_AST", "GP1_1_WAC11_AST", "F13505_S2A3_AST", "F13505_NEA11_AST",
                                "S75_11_IW1A12_6W", "S75_11_IW1F3_6W", "S75_11_IW1G8_6W", 
                                "S300_12_NWBB7_6W", "GP1_1_WAC11_6W", "F13505_S2A3_6W", "F13505_NEA11_6W")
gse260517_data <- gse260517_data[, selected_samples_gse260517]
gse260517_data <- as_tibble(gse260517_data)

# Process gene IDs - extract Ensembl ID before the pipe symbol
gse260517_data$gene_id <- sapply(gse260517_data$gene_id, function(x) {
  parts <- strsplit(x, "\\|")[[1]]
  return(parts[1])
})

# Convert to data frame and set row names
gse260517_data <- gse260517_data %>% column_to_rownames(var = "gene_id")

# Convert to numeric
gse260517_data[] <- lapply(gse260517_data, function(x) as.numeric(as.character(x)))

# Create DGEList object and normalize
gse260517_dge <- DGEList(counts = as.matrix(gse260517_data))
gse260517_dge <- calcNormFactors(gse260517_dge)
logCPM_gse260517 <- cpm(gse260517_dge, log = TRUE, prior.count = 0.5)

# Print a check of gene IDs
print("First few gene IDs in GSE260517 dataset:")
print(head(rownames(logCPM_gse260517)))

# Define groups for GSE260517 data
Bowles_2024_astrocytes <- c("S75_11_IW1A12_AST", "S75_11_IW1F3_AST", "S75_11_IW1G8_AST", 
                           "S300_12_NWBB7_AST", "GP1_1_WAC11_AST", "F13505_S2A3_AST", "F13505_NEA11_AST")
Bowles_2024_6w_neurons <- c("S75_11_IW1A12_6W", "S75_11_IW1F3_6W", "S75_11_IW1G8_6W", 
                          "S300_12_NWBB7_6W", "GP1_1_WAC11_6W", "F13505_S2A3_6W", "F13505_NEA11_6W")



# Load GSE190185 data
gse190185_data <- read.delim(GSE190185_FILE, header = TRUE, stringsAsFactors = FALSE)
colnames(gse190185_data)[1] <- "gene_id"
selected_samples_gse190185 <- c("gene_id", "LA_TCW1_33A", "LA_TCW2_33A", "LA_MC2E_33", "LA_RC1H_33")
gse190185_data <- gse190185_data[, selected_samples_gse190185]
gse190185_data <- as_tibble(gse190185_data)
gse190185_data <- gse190185_data %>% column_to_rownames(var = "gene_id")
gse190185_dge <- DGEList(counts = as.matrix(gse190185_data))
gse190185_dge <- calcNormFactors(gse190185_dge)
logCPM_gse190185 <- cpm(gse190185_dge, log = TRUE, prior.count = 0.5)

#Load zhang data
zhang_data <- read.csv(GSE73721_FILE, header = TRUE, stringsAsFactors = FALSE)
colnames(zhang_data)[1] <- "ENSEMBL_ID"
selected_samples_zhang_data <- c("ENSEMBL_ID", "X21yo.Hippocampus.astro", "X22yo.Hippocampus.astro", "X53yo.A.Hippocampus.astro", 
                                 "X53yo.B.Hippocampus.astro", "Fetal.ctx.1.astro",           "Fetal.ctx.2.astro" ,   "Fetal.ctx.3.astro",       
                                 "Fetal.ctx.4.astro"    ,       "Fetal.ctx.5.astro"   ,       "Fetal.ctx.6.astro"      ,     
                                 "X8yo.ctx.astro"    ,          "X13yo.ctx.astro"  ,         "X16yo.ctx.astro"  ,  "X21yo.ctx.astro",
                                 "X22yo.ctx.astro"      ,       "X35yo.ctx.astro"     ,        "X47yo.ctx.astro"         ,    "X51yo.ctx.astro" ,           
                                 "X53yo.ctx.astro"       ,      "X60yo.ctx.astro"    ,     "X25yo.ctx.neuron",    "X63yo.ctx.1.astro", "X45yo.ctx.myeloid", "X51yo.ctx.myeloid", "X63.yo.ctx.myeloid")
zhang_data <- zhang_data[, selected_samples_zhang_data]
zhang_data <- as_tibble(zhang_data)
# Remove rows where ENSEMBL_ID is NA or empty
cleaned_zhang_data <- zhang_data %>%
  filter(!is.na(ENSEMBL_ID) & ENSEMBL_ID != "")
cleaned_zhang_data <- cleaned_zhang_data %>% column_to_rownames(var = "ENSEMBL_ID")
zhang_data_dge <- DGEList(counts = as.matrix(cleaned_zhang_data))
zhang_data_dge <- calcNormFactors(zhang_data_dge)
logCPM_zhang_data <- cpm(zhang_data_dge, log = TRUE, prior.count = 0.5)

# Load GSE182307 data
gse182307_data <- read.csv(GSE182307_FILE, header = TRUE, stringsAsFactors = FALSE)
colnames(gse182307_data)[1] <- "gene_id"
colnames(gse182307_data) <- make.names(colnames(gse182307_data))
selected_samples_gse182307 <- c("gene_id", "iAstrocyte_Vehicle_1", "Krencik.astrocyte_Vehicle_1", "Li.astrocyte_Vehicle_1", "TCW.astrocyte_Vehicle_1", "iNeuron_Vehicle_1")
gse182307_data <- gse182307_data[, selected_samples_gse182307]
gse182307_data <- as_tibble(gse182307_data)
gse182307_data <- gse182307_data %>% column_to_rownames(var = "gene_id")
gse182307_dge <- DGEList(counts = as.matrix(gse182307_data))
gse182307_dge <- calcNormFactors(gse182307_dge)
logCPM_gse182307 <- cpm(gse182307_dge, log = TRUE, prior.count = 0.5)

# Load GSE97904 data
gse97904_data <- read.delim(GSE97904_FILE, header = TRUE, stringsAsFactors = FALSE)
colnames(gse97904_data)[1] <- "gene_id"
gse97904_data[gse97904_data < 0] <- 0  # Set negative values to zero
gse97904_data <- as_tibble(gse97904_data)
gse97904_data <- gse97904_data %>% column_to_rownames(var = "gene_id")
gse97904_dge <- DGEList(counts = as.matrix(gse97904_data))
gse97904_dge <- calcNormFactors(gse97904_dge)
logCPM_gse97904 <- cpm(gse97904_dge, log = TRUE, prior.count = 0.5)

# other things that need to be labeled:
neuron_data_kamp <- c("iNeuron_Vehicle_1")
astro_data_kamp <- c("iAstrocyte_Vehicle_1", "Krencik.astrocyte_Vehicle_1", "Li.astrocyte_Vehicle_1", "TCW.astrocyte_Vehicle_1")
#tcw_astro <- c("TCW.astrocyte_Vehicle_1")

neuron_data_zhang <- c("X25yo.ctx.neuron")
myeloid_data_zhang <- c("X63yo.ctx.1.astro", "X45yo.ctx.myeloid", "X51yo.ctx.myeloid", "X63.yo.ctx.myeloid")

adult_hipp_ast <- c("X21yo.Hippocampus.astro", "X22yo.Hippocampus.astro", "X53yo.A.Hippocampus.astro", "X53yo.B.Hippocampus.astro")   
fetal_ctx_ast <- c("Fetal.ctx.1.astro",           "Fetal.ctx.2.astro" ,   "Fetal.ctx.3.astro", "Fetal.ctx.4.astro"    ,       "Fetal.ctx.5.astro"   ,       "Fetal.ctx.6.astro")
adult_ctx_ast <- c("X8yo.ctx.astro"    ,          "X13yo.ctx.astro"  ,         "X16yo.ctx.astro"  ,  "X21yo.ctx.astro", "X22yo.ctx.astro"      ,       "X35yo.ctx.astro"     ,        "X47yo.ctx.astro"         ,    "X51yo.ctx.astro" ,           
"X53yo.ctx.astro"       ,      "X60yo.ctx.astro"    ,         "X63yo.ctx.1.astro")


rownames(logCPM_gse182307) <- remove_version(rownames(logCPM_gse182307))
rownames(logCPM_gse190185) <- remove_version(rownames(logCPM_gse190185))
rownames(logCPM_gse97904) <- remove_version(rownames(logCPM_gse97904))
rownames(logCPM_zhang_data) <- remove_version(rownames(logCPM_zhang_data))
rownames(logCPM_gse260517) <- remove_version(rownames(logCPM_gse260517))

# Find common genes across all datasets
#common_genes <- Reduce(intersect, list(rownames(logCPM_gse190185), rownames(logCPM_gse182307), rownames(logCPM_gse97904), rownames(datExpr_astro), rownames(datExpr_microglia), rownames(datExpr_neurons), rownames(logCPM_zhang_data)))
common_genes <- Reduce(intersect, list(rownames(logCPM_gse190185), rownames(logCPM_gse182307), rownames(datExpr_astro), rownames(datExpr_neuron), rownames(logCPM_zhang_data)))

# Find common genes across all datasets including GSE260517
common_genes <- Reduce(intersect, list(
  rownames(logCPM_gse190185), 
  rownames(logCPM_gse182307), 
  rownames(datExpr_astro), 
  rownames(datExpr_neuron), 
  rownames(logCPM_zhang_data),
  rownames(logCPM_gse260517),
  rownames(datExpr)
))
# Check the number of common genes
length(common_genes)
print(head(common_genes))

# Combine datasets by matching common genes
# Combine datasets by matching common genes
combined_data <- cbind(
  logCPM_gse190185[common_genes, ],
  logCPM_gse182307[common_genes, ],
  datExpr_astro[common_genes, ],
  datExpr_neuron[common_genes, ],
  logCPM_zhang_data[common_genes, ],
  logCPM_gse260517[common_genes, ],
  datExpr[common_genes, ]
  # Add the new dataset
)

# Perform PCA on the combined data
combined_data <- t(scale(t(combined_data), scale = FALSE))
PCA_res_combined <- prcomp(combined_data, center = FALSE)

# Create a metadata frame for combined data
# Create a vector of all your specified samples
all_specified_samples <- c(
  colnames(logCPM_gse190185),
  #adult_ctx_ast,
   astro_data_kamp,
  neuron_data_kamp,
  #neuron_data_zhang,
 # myeloid_data_zhang,
  h1h1_afr_samples_astro,
  h1h2_afr_samples_astro,
  h1h1a_afr_samples_astro,
  h1h1_eur_samples_astro,
  h1h2_eur_samples_astro,
  h1h1_afr_samples_neuron,
  h1h1a_afr_samples_neuron,
  h1h1_eur_samples_neuron,
  h1h2_eur_samples_neuron,
  h1h2_afr_samples_neuron,
 #h1h1_afr_samples_microglia,
 #h1h1_eur_samples_microglia,
 #h1h2_eur_samples_microglia,
 #h1h2_afr_samples_microglia,
  Bowles_2024_astrocytes, 
  Bowles_2024_6w_neurons
 
)

# Filter your combined_data to include only specified samples
combined_data <- combined_data[, colnames(combined_data) %in% all_specified_samples]

# Now run your PCA on the filtered data
combined_data <- t(scale(t(combined_data), scale = FALSE))
PCA_res_combined <- prcomp(combined_data, center = FALSE)

# Create metadata for the filtered samples
combined_meta <- data.frame(
  sample = colnames(combined_data),
  source = NA  # Initialize with NA
)

# Assign categories in a way that avoids duplicates
combined_meta$source[combined_meta$sample %in% colnames(logCPM_gse190185)] <- "TCW 2022 iAstrocytes"
#combined_meta$source[combined_meta$sample %in% adult_ctx_ast] <- "Adult Cortex Astrocytes Zhang 2016"
combined_meta$source[combined_meta$sample %in% astro_data_kamp] <- "Leng 2022 iAstrocyte"
combined_meta$source[combined_meta$sample %in% neuron_data_kamp] <- "Leng 2022 iNeuron"
#combined_meta$source[combined_meta$sample %in% neuron_data_zhang] <- "Adult Hippocampal Neurons Zhang 2016"
#combined_meta$source[combined_meta$sample %in% myeloid_data_zhang] <- "Adult Hippocampal Myeloid Zhang 2016"
combined_meta$source[combined_meta$sample %in% h1h1_afr_samples_astro] <- "H1H1 AFR iPSC-derived Astrocytes"
combined_meta$source[combined_meta$sample %in% h1h2_afr_samples_astro] <- "H1H2 AFR iPSC-derived Astrocytes"
combined_meta$source[combined_meta$sample %in% h1h1a_afr_samples_astro] <- "H1H1_A AFR iPSC-derived Astrocytes"
combined_meta$source[combined_meta$sample %in% h1h1_eur_samples_astro] <- "H1H1 EUR iPSC-derived Astrocytes"
combined_meta$source[combined_meta$sample %in% h1h2_eur_samples_astro] <- "H1H2 EUR iPSC-derived Astrocytes"
combined_meta$source[combined_meta$sample %in% h1h1_afr_samples_neuron] <- "H1H1 AFR iPSC-derived Neurons"
combined_meta$source[combined_meta$sample %in% h1h1a_afr_samples_neuron] <- "H1H1_A AFR iPSC-derived Neurons"
combined_meta$source[combined_meta$sample %in% h1h1_eur_samples_neuron] <- "H1H1 EUR iPSC-derived Neurons"
combined_meta$source[combined_meta$sample %in% h1h2_eur_samples_neuron] <- "H1H2 EUR iPSC-derived Neurons"
combined_meta$source[combined_meta$sample %in% h1h2_afr_samples_neuron] <- "H1H2 AFR iPSC-derived Neurons"
combined_meta$source[combined_meta$sample %in% Bowles_2024_astrocytes] <- "Bowles 2024 iPSC-derived Astrocytes"
combined_meta$source[combined_meta$sample %in% Bowles_2024_6w_neurons] <- "Bowles 2024 6W Neurons"

#combined_meta$source[combined_meta$sample %in% h1h1_afr_samples_microglia] <- "H1H1 AFR iPSC-derived Microglia"
#combined_meta$source[combined_meta$sample %in% h1h1a_afr_samples_microglia] <- "H1H1_A AFR iPSC-derived Microglia"
#combined_meta$source[combined_meta$sample %in% h1h1_eur_samples_microglia] <- "H1H1 EUR iPSC-derived Microglia"
#combined_meta$source[combined_meta$sample %in% h1h2_eur_samples_microglia] <- "H1H2 EUR iPSC-derived Microglia"
#combined_meta$source[combined_meta$sample %in% h1h2_afr_samples_microglia] <- "H1H2 AFR iPSC-derived Microglia"

# Check if any samples still have NA as source (shouldn't happen after filtering)
na_samples <- combined_meta$sample[is.na(combined_meta$source)]
if (length(na_samples) > 0) {
  warning("Some samples still have NA as source. Check your group definitions.")
  print(na_samples)
}

# Extract PCA results and merge with metadata
PCs_combined <- as.data.frame(PCA_res_combined$rotation[, 1:2])
PCs_combined$sample <- rownames(PCs_combined)
PCs_combined <- merge(PCs_combined, combined_meta, by = "sample")

# Load the required library
library(RColorBrewer)

# Create a vector of colors based on the source
source_names <- unique(combined_meta$source)
num_sources <- length(source_names)

# If the number of sources exceeds the maximum available colors in Set3, extend the palette
if (num_sources > 12) {
  colors <- c(brewer.pal(12, "Set3"), brewer.pal(num_sources - 12, "Paired"))
} else {
  colors <- brewer.pal(num_sources, "Set3")
}

# Adjust length of colors to match number of sources
colors <- colors[1:num_sources]  
source_colors <- setNames(colors, source_names)

# Verify lengths
length(colors)
length(source_names)
    
# Save the plot to a file
png(filename = file.path(OUTPUT_DIR, "PCs_plot_with_ancestry_lines_with_leng.png"), width = 1600, height = 1200, res = 150)

# Plot PCA with consistent colors for datasets
ggplot(PCs_combined, aes(x = PC1, y = PC2, color = source)) +
  geom_point(size = 3) +
  scale_color_manual(values = source_colors) +  # Use the same colors for each dataset
  labs(title = "PCA of Combined Data", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

dev.off()






# Create a more detailed dataset group column with proper order of pattern matching
PCs_combined$dataset_group <- NA

# The order matters! Put the Bowles 2024 pattern FIRST so it gets priority
PCs_combined$dataset_group[grepl("Bowles 2024", PCs_combined$source)] <- "Bowles 2024"

# Then match the remaining patterns
PCs_combined$dataset_group[is.na(PCs_combined$dataset_group) & grepl("iPSC-derived Astrocytes", PCs_combined$source)] <- "iPSC-derived Astrocytes"
PCs_combined$dataset_group[is.na(PCs_combined$dataset_group) & grepl("iPSC-derived Neurons", PCs_combined$source)] <- "iPSC-derived Neurons"
PCs_combined$dataset_group[is.na(PCs_combined$dataset_group) & grepl("Leng 2022", PCs_combined$source)] <- "Leng 2022"
PCs_combined$dataset_group[is.na(PCs_combined$dataset_group) & grepl("TCW 2022", PCs_combined$source)] <- "TCW 2022"

# Create a mapping of shapes to dataset groups (with 5 NON-CIRCLE shapes)
dataset_groups <- unique(PCs_combined$dataset_group)
dataset_shapes <- setNames(
  c(17, 15, 18, 8, 5),  # Triangle, square, diamond, asterisk, diamond plus
  dataset_groups
)


# Manually specify colors for each source (excluding yellow and grey tones)
source_colors <- c(
  "H1H1 AFR iPSC-derived Astrocytes" = "#E41A1C",    # Red
  "H1H2 AFR iPSC-derived Astrocytes" = "#377EB8",    # Blue
  "H1H1_A AFR iPSC-derived Astrocytes" = "#4DAF4A",  # Green
  "H1H1 EUR iPSC-derived Astrocytes" = "#984EA3",    # Purple
  "H1H2 EUR iPSC-derived Astrocytes" = "#FF7F00",    # Orange
  "H1H1 AFR iPSC-derived Neurons" = "#A65628",       # Brown
  "H1H1_A AFR iPSC-derived Neurons" = "#F781BF",     # Pink
  "H1H1 EUR iPSC-derived Neurons" = "#1B9E77",       # Teal
  "H1H2 EUR iPSC-derived Neurons" = "#D95F02",       # Vermilion
  "H1H2 AFR iPSC-derived Neurons" = "#7570B3",       # Violet
  "Bowles 2024 iPSC-derived Astrocytes" = "#E7298A", # Magenta
  "Bowles 2024 6W Neurons" = "#66A61E",              # Olive
  "TCW 2022 iAstrocytes" = "#5E4FA2",                # Indigo
  "Leng 2022 iAstrocyte" = "#DC050C",                # Dark red
  "Leng 2022 iNeuron" = "#1965B0"                    # Navy blue
)

# Check if we've covered all source categories
missing_sources <- setdiff(unique(PCs_combined$source), names(source_colors))
if(length(missing_sources) > 0) {
  warning("Some sources are missing from the color mapping: ", 
          paste(missing_sources, collapse = ", "))
  # Add colors for any missing sources
  additional_colors <- c("#882E72", "#B17BA6", "#332288", "#88CCEE", "#44AA99")
  for(i in 1:length(missing_sources)) {
    if(i <= length(additional_colors)) {
      source_colors[missing_sources[i]] <- additional_colors[i]
    }
  }
}

# Save the plot to a file
png(filename = file.path(OUTPUT_DIR, "PCs_plot_five_shapes_no_circle.png"),
    width = 1600, height = 1200, res = 150)

# Plot PCA with five non-circle shapes
ggplot(PCs_combined, aes(x = PC1, y = PC2, color = source)) +
  geom_point(aes(shape = dataset_group), size = 4) +
  scale_color_manual(values = source_colors) +
  scale_shape_manual(values = dataset_shapes) +
  labs(
    title = "PCA of Combined RNA-seq Data", 
    x = "PC1", 
    y = "PC2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

dev.off()