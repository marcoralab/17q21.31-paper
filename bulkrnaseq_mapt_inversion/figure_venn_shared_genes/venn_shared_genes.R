# Venn diagrams of differentially expressed genes shared across five
# haplotype/ancestry comparisons. Works for any cell type or brain region:
# point the five inputs at that group's DE tables and set `label` accordingly
# (e.g. label = "neurons" reproduces `Downregulated_Genes_Venn_neurons.png`).
# Also writes the up-regulated, combined, and intersection-table outputs.

library(VennDiagram)
library(dplyr)
library(grid)

# ===============================================================
# CONFIG
# Label for this run (cell type / region) — used in the output filenames.
label <- "neurons"

# Folder of annotated DE result CSVs for this dataset (from annotate_de_results.py).
# Each file must contain: Gene_Symbol, log2FoldChange, pvalue, padj.
de_dir <- "path/to/de_results"
H1H1A_AA_vs_H1H1_AA  <- read.csv(file.path(de_dir, "h1h1A_aa_vs_h1h1_aa.csv"), header = TRUE)
H1H2_AA_vs_H1H2_EUR  <- read.csv(file.path(de_dir, "h1h2_aa_vs_h1h2_eur_flipped.csv"), header = TRUE)
H1H2_EUR_vs_H1H1_EUR <- read.csv(file.path(de_dir, "h1h2_eur_vs_h1h1_eur.csv"), header = TRUE)
H1H2_AA_vs_H1H1_AA   <- read.csv(file.path(de_dir, "h1h2_aa_vs_h1h1_aa.csv"), header = TRUE)
H1H1_AA_vs_H1H1_EUR  <- read.csv(file.path(de_dir, "h1h1_aa_vs_h1h1_eur_flipped.csv"), header = TRUE)

# Where outputs are written.
output_dir <- "results/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# ===============================================================

# Keep the strongest effect per gene symbol (collapse duplicates)
handle_duplicates <- function(df) {
  df_cleaned <- df %>%
    group_by(Gene_Symbol) %>%
    slice(which.max(abs(log2FoldChange))) %>%
    ungroup()
  return(as.data.frame(df_cleaned))
}

H1H1A_AA_vs_H1H1_AA  <- handle_duplicates(H1H1A_AA_vs_H1H1_AA)
H1H2_AA_vs_H1H2_EUR  <- handle_duplicates(H1H2_AA_vs_H1H2_EUR)
H1H2_EUR_vs_H1H1_EUR <- handle_duplicates(H1H2_EUR_vs_H1H1_EUR)
H1H2_AA_vs_H1H1_AA   <- handle_duplicates(H1H2_AA_vs_H1H1_AA)
H1H1_AA_vs_H1H1_EUR  <- handle_duplicates(H1H1_AA_vs_H1H1_EUR)

# Filter (drop NA / 0 / 1 padj) and split into up/down/combined significant sets
split_sig <- function(df) {
  f <- df %>% dplyr::filter(!is.na(padj) & padj != 1 & padj != 0)
  list(
    DW   = subset(f, pvalue < 0.05 & log2FoldChange < -1),
    UP   = subset(f, pvalue < 0.05 & log2FoldChange > 1),
    COMB = subset(f, pvalue < 0.05 & abs(log2FoldChange) > 1)
  )
}

s1 <- split_sig(H1H1A_AA_vs_H1H1_AA)
s2 <- split_sig(H1H2_AA_vs_H1H2_EUR)
s3 <- split_sig(H1H2_EUR_vs_H1H1_EUR)
s4 <- split_sig(H1H2_AA_vs_H1H1_AA)
s5 <- split_sig(H1H1_AA_vs_H1H1_EUR)

DW <- list(
  "H1H1A_AA vs H1H1_AA"  = s1$DW$Gene_Symbol,
  "H1H2_AA vs H1H2_EUR"  = s2$DW$Gene_Symbol,
  "H1H2_EUR vs H1H1_EUR" = s3$DW$Gene_Symbol,
  "H1H2_AA vs H1H1_AA"   = s4$DW$Gene_Symbol,
  "H1H1_AA vs H1H1_EUR"  = s5$DW$Gene_Symbol
)
UP <- list(
  "H1H1A_AA vs H1H1_AA"  = s1$UP$Gene_Symbol,
  "H1H2_AA vs H1H2_EUR"  = s2$UP$Gene_Symbol,
  "H1H2_EUR vs H1H1_EUR" = s3$UP$Gene_Symbol,
  "H1H2_AA vs H1H1_AA"   = s4$UP$Gene_Symbol,
  "H1H1_AA vs H1H1_EUR"  = s5$UP$Gene_Symbol
)
COMB <- list(
  "H1H1A_AA vs H1H1_AA"  = s1$COMB$Gene_Symbol,
  "H1H2_AA vs H1H2_EUR"  = s2$COMB$Gene_Symbol,
  "H1H2_EUR vs H1H1_EUR" = s3$COMB$Gene_Symbol,
  "H1H2_AA vs H1H1_AA"   = s4$COMB$Gene_Symbol,
  "H1H1_AA vs H1H1_EUR"  = s5$COMB$Gene_Symbol
)

# Report the genes unique to each comparison (no overlap)
get_unique_counts <- function(gene_list, list_names) {
  unique_genes <- list()
  for (name in list_names) {
    unique_genes[[name]] <- setdiff(
      gene_list[[name]],
      unlist(gene_list[names(gene_list) != name])
    )
  }
  cat("\nNumber of unique genes (no overlap) in each group:\n")
  for (name in list_names) {
    cat(sprintf("%s: %d genes\n", name, length(unique_genes[[name]])))
  }
  return(unique_genes)
}

cat("\nDOWNREGULATED GENES:")
invisible(get_unique_counts(DW, names(DW)))
cat("\nUPREGULATED GENES:")
invisible(get_unique_counts(UP, names(UP)))
cat("\nCOMBINED GENES (UP AND DOWN):")
invisible(get_unique_counts(COMB, names(COMB)))

venn_fill <- c("#EFC000FF", "orange", "salmon", "red", "purple")

draw_venn <- function(gene_lists, title, out_png) {
  p <- venn.diagram(
    x = gene_lists,
    filename = NULL,
    fill = venn_fill,
    main = title,
    main.cex = 2,
    cex = 1.5,                                  # size of the count numbers
    euler.d = TRUE,
    scaled = TRUE,
    category.names = c("", "", "", "", ""),     # no category labels
    margin = 0.2
  )
  png(file.path(output_dir, out_png), width = 3000, height = 3000, res = 300)
  grid.draw(p)
  dev.off()
}

# The paper figure is the downregulated Venn:
draw_venn(DW,   "Downregulated Genes",             paste0("Downregulated_Genes_Venn_", label, ".png"))
draw_venn(UP,   "Upregulated Genes",               paste0("Upregulated_Genes_Venn_", label, ".png"))
draw_venn(COMB, "Combined DE Genes (Up and Down)", paste0("Combined_DE_Genes_Venn_", label, ".png"))

# Write every set-intersection (genes specific to each combination) to CSV
get_all_intersections <- function(gene_list) {
  list_names <- names(gene_list)
  n <- length(list_names)
  all_intersections <- list()
  for (i in 1:n) {
    combs <- combn(list_names, i, simplify = FALSE)
    for (comb in combs) {
      intersection_genes <- Reduce(intersect, gene_list[comb])
      other_lists <- list_names[!list_names %in% comb]
      if (length(other_lists) > 0) {
        other_genes <- unique(unlist(gene_list[other_lists]))
        intersection_genes <- setdiff(intersection_genes, other_genes)
      }
      if (length(intersection_genes) > 0) {
        all_intersections[[paste(comb, collapse = " & ")]] <- intersection_genes
      }
    }
  }
  return(all_intersections)
}

write_intersections_to_csv <- function(intersections, filename) {
  max_length <- max(sapply(intersections, length))
  df <- data.frame(matrix(NA, nrow = max_length, ncol = length(intersections)))
  colnames(df) <- names(intersections)
  for (name in names(intersections)) {
    df[1:length(intersections[[name]]), name] <- intersections[[name]]
  }
  write.csv(df, file.path(output_dir, filename), row.names = FALSE, na = "")
}

write_intersections_to_csv(get_all_intersections(DW),   paste0("downregulated_genes_intersections_", label, ".csv"))
write_intersections_to_csv(get_all_intersections(UP),   paste0("upregulated_genes_intersections_", label, ".csv"))
write_intersections_to_csv(get_all_intersections(COMB), paste0("combined_genes_intersections_", label, ".csv"))

cat("\nDone. Venn diagrams and intersection tables written to", output_dir, "\n")
