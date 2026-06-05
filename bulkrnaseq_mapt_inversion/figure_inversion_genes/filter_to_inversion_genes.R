library(rtracklayer)

# ===============================================================
# CONFIG — edit these paths for your environment (see also `de_dir` and the
# `comparisons` list further down).
# ===============================================================
# GENCODE annotation GTF (used to map gene_id -> gene_name).
gtf_file <- "path/to/gencode.annotation.gtf.gz"
# NCBI gene list for the 17q21 inversion locus (+/- 1MB); column: Gene_symbol.
gene_list_file <- "path/to/NCBI_Genes_in_Inversion_locus_and1MB.csv"
# Where the inversion-filtered CSVs are written.
output_dir <- "inversion_filtered_csvs"
# ===============================================================

# Read the GTF file once
gtf_data <- import(gtf_file, format = "gtf")
gene_mapping <- gtf_data[gtf_data$type == "gene", c("gene_id", "gene_name")]
gene_mapping$gene_id_clean <- gsub("\\..*", "", gene_mapping$gene_id)

# Read the inversion gene list once
gene_list <- read.csv(gene_list_file, stringsAsFactors = FALSE)
gene_symbols <- gene_list$Gene_symbol

# Function to process a single file
process_file <- function(input_path, name) {
    cat("Processing file:", name, "\n")
    
    # Output directory (from CONFIG above)
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }
    
    output_name <- file.path(output_dir, paste0("inversion_", name, ".csv"))
    
    # Read the results file
    results <- read.csv(input_path)
    
    # Clean gene IDs
    results$gene_id_clean <- gsub("\\..*", "", results$gene_id)
    
    # Merge with gene mapping
    results_with_symbols <- merge(results, 
                                as.data.frame(gene_mapping), 
                                by.x = "gene_id_clean", 
                                by.y = "gene_id_clean", 
                                all.x = TRUE)
    
    # Subset to include only inversion genes
    subset_results <- results_with_symbols[results_with_symbols$gene_name %in% gene_symbols, ]
    
    # Recalculate adjusted p-values
    subset_results$padj <- p.adjust(subset_results$pvalue, method = "BH")
    
    # Save results
    write.csv(subset_results, output_name, row.names = FALSE)
    
    # Return the number of genes found
    return(list(
        file = name,
        genes_found = nrow(subset_results)
    ))
}

# Folder of annotated DE result CSVs (output of annotate_de_results.py).
# Name each file <comparison>.csv using the comparison names listed below;
# each is filtered to inversion genes and written as inversion_<comparison>.csv.
de_dir <- "path/to/de_results"

comparisons <- c(
    # flipped EUR-vs-AA ancestry contrasts
    "h1h1_aa_vs_h1h1_eur_astro_flipped_eur_vs_aa",
    "h1h2_aa_vs_h1h2_eur_astro_flipped_eur_vs_aa",
    "h1h1_aa_vs_h1h1_eur_brain_caudate_flipped_eur_vs_aa",
    "h1h2_aa_vs_h1h2_eur_brain_caudate_flipped_eur_vs_aa",
    "h1h1_aa_vs_h1h1_eur_brain_mfg_flipped_eur_vs_aa",
    "h1h2_aa_vs_h1h2_eur_brain_mfg_flipped_eur_vs_aa",
    "h1h1_aa_vs_h1h1_eur_neurons_flipped_eur_vs_aa",
    "h1h2_aa_vs_h1h2_eur_neurons_flipped_eur_vs_aa",
    # caudate (haplotype + ancestry)
    "aa_h1h1a_vs_aa_h1h1_brain_caudate",
    "aa_h1h2_vs_aa_h1h1_brain_caudate",
    "eur_h1h2_vs_eur_h1h1_brain_caudate",
    "h1h1_aa_vs_h1h1_eur_brain_caudate",
    "h1h2_aa_vs_h1h2_eur_brain_caudate",
    # MFG
    "aa_h1h1a_vs_aa_h1h1_brain_mfg",
    "aa_h1h2_vs_aa_h1h1_brain_mfg",
    "eur_h1h2_vs_eur_h1h1_brain_mfg",
    "h1h1_aa_vs_h1h1_eur_brain_mfg",
    "h1h2_aa_vs_h1h2_eur_brain_mfg",
    # astrocytes
    "h1h1_aa_vs_h1h1_eur_astro",
    "h1h2_aa_vs_h1h2_eur_astro",
    "aa_h1h1a_vs_aa_h1h1_astro",
    "aa_h1h2_vs_aa_h1h1_astro",
    "eur_h1h2_vs_eur_h1h1_astro",
    # neurons
    "h1h1_aa_vs_h1h1_eur_neurons",
    "h1h2_aa_vs_h1h2_eur_neurons",
    "aa_h1h2_vs_aa_h1h1_neuron",
    "eur_h1h2_vs_eur_h1h1_neuron",
    "aa_h1h1a_vs_aa_h1h1_neuron"
)

# Process all files
results <- list()
for (name in comparisons) {
    results[[name]] <- process_file(file.path(de_dir, paste0(name, ".csv")), name)
}

# Print summary
cat("\nProcessing Summary:\n")
for (name in names(results)) {
    result <- results[[name]]
    cat(sprintf("File: %s\nGenes found: %d\n\n", result$file, result$genes_found))
}