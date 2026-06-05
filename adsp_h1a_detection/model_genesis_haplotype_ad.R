# Load required packages
library(lme4)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Biobase)
library(Matrix)
library(tidyverse)
library(magrittr)
library(GENESIS)

rm(list = ls())

# -----------------------------
# 1. Load data 
# -----------------------------

# Modeling input produced by the Python merge step
# (merge_with_phenotype_file.py).
input_file <- "results/mapt_genotypes_with_covar_merged.csv"
output_dir <- "results"
model_dir <- file.path(output_dir, "models")
plot_dir <- file.path(output_dir, "plots")
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Configure analysis
# Modify these variables to match dataset structure
id_column <- "sample.id"           # ID column in data
outcome <- "diagnosis_code"        # Name of outcome variable
glmm_family <- "binomial"          # GLM family (usually binomial for case/control)
pc_cols <- paste0("V", 1:5)        # Principal component columns to include
random_covars <- "sequencing_center_platform"  # Variables to be included as random effects

# Path to the genetic relatedness matrix (GRM): a sample x sample .rds
# matrix (e.g. a GCTA/PC-Relate kinship matrix) with sample IDs as
# row/column names. Used as the random effect in the GENESIS null model.
grm_path <- "path/to/your_grm.rds"


#LOAD DATA
df <- read.csv(input_file, stringsAsFactors = FALSE)
cat("Original dataset size:", nrow(df), "samples\n")

# Filter to only include cases and controls 
df_filtered <- df %>%
  filter(diagnosis %in% c("AD", "Control"))  # Replace with case/control values

# Remove rows with missing data in key columns
# Modify this list to match the variables needed for analysis
df_complete <- df_filtered %>%
  filter(
    !is.na(Haplotype) & 
    !is.na(sex) & 
    !is.na(diagnosis) & 
    !is.na(apoe) & 
    !is.na(admixture_super_pop_max)
  )


# Generate group_info if it doesn't exist (modify according to the dataset's structure)
if (!"group_info" %in% names(df_complete)) {
  df_complete <- df_complete %>% 
    rowwise() %>% 
    mutate(group_info = {
      eth <- admixture_super_pop_max
      if (Agreement == "Agree") {
        geno <- trimws(gsub("\\|", "/", rs8070723))
        recode <- ifelse(geno == "0/0", "H1H1",
                         ifelse(geno %in% c("0/1", "1/0"), "H1H2",
                                ifelse(geno == "1/1", "H2H2", "Unknown")))
        paste0(recode, "_", eth)
      } else {
        paste0("Disagree_", eth)
      }
    }) %>% 
    ungroup()
}

# Create binary outcome variable (modify as needed for outcome)
df_complete <- df_complete %>% 
  mutate(diagnosis_code = ifelse(diagnosis == "AD", 1, 0))

# Create haplotype2 variable (modify as needed for haplotype coding)
df_complete <- df_complete %>%
  mutate(haplotype2 = case_when(
    startsWith(group_info, "Disagree") ~ "H1_A Carrier",
    startsWith(group_info, "H1H1") ~ "H1H1",
    startsWith(group_info, "H1H2") ~ "H1H2",
    startsWith(group_info, "H2H2") ~ "H2H2",
    TRUE ~ NA_character_
  ))

# Check distributions
cat("\nDiagnosis distribution in final dataset:\n")
print(table(df_complete$diagnosis))

cat("\nHaplotype2 distribution in final dataset:\n")
print(table(df_complete$haplotype2, useNA = "ifany"))

cat("\nAncestry distribution in final dataset:\n")
print(table(df_complete$admixture_super_pop_max))

cat("\nAPOE genotype distribution in final dataset:\n")
print(table(df_complete$apoe, useNA = "ifany"))

# -----------------------------
# 3. Prepare data for modeling
# -----------------------------

# Define base covariates for all models
base_covars <- c("sex", "apoe")  # Modify based on covariates

# Transform the dataframe
df_transformed <- df_complete %>%
  # Create diagnosis_code from diagnosis
  mutate(diagnosis_code = case_when(
    diagnosis == "AD" ~ 1,
    diagnosis == "Control" ~ 0,
    TRUE ~ NA_real_
  ),
  # Create combined sequencing_center_platform column
  sequencing_center_platform = paste(sequencing_center, sequencing_platform, sep = "_")) %>%
  # Convert base covariates to factors
  mutate_at(base_covars, as.factor) %>%
  # Also convert admixture_super_pop_max to factor for overall model
  mutate(admixture_super_pop_max = as.factor(admixture_super_pop_max))

# Filter for European ancestry samples (modify for ancestry designations)
df_eur <- df_transformed %>% 
  filter(admixture_super_pop_max == "EUR")

# Filter for African ancestry samples
df_afr <- df_transformed %>% 
  filter(admixture_super_pop_max == "AFR")

# Print sample counts for each ancestry group
cat("Overall samples:", nrow(df_transformed), "\n")
cat("European ancestry samples:", nrow(df_eur), "\n")
cat("African ancestry samples:", nrow(df_afr), "\n")

overall_covars <- base_covars
eur_afr_covars <- base_covars

# -----------------------------
# 4. Load Genetic Relatedness Matrix (GRM) if available
# -----------------------------

### GENESIS — the GRM configured above is used as the kinship/covariance matrix
cov_mat <- readRDS(grm_path)

cov_mat_samples <- rownames(cov_mat)

# Check how many samples are in the covariance matrix
cat("Number of samples in covariance matrix:", length(cov_mat_samples), "\n")
# Print sample counts before filtering
cat("Samples in df_transformed before filtering:", nrow(df_transformed), "\n")
cat("Samples in df_eur before filtering:", nrow(df_eur), "\n")
cat("Samples in df_afr before filtering:", nrow(df_afr), "\n")

# Filter each dataframe to only keep samples that are in the covariance matrix
df_transformed_filtered <- df_transformed %>%
  filter(!!sym(id_column) %in% cov_mat_samples)

df_eur_filtered <- df_eur %>%
  filter(!!sym(id_column) %in% cov_mat_samples)

df_afr_filtered <- df_afr %>%
  filter(!!sym(id_column) %in% cov_mat_samples)

## Filter GRM for samples in analysis
cov_mat_filtered_overall <- cov_mat[df_transformed_filtered$sample.id, df_transformed_filtered$sample.id]
cov_mat_filtered_eur <- cov_mat[df_eur_filtered$sample.id, df_eur_filtered$sample.id]
cov_mat_filtered_afr <- cov_mat[df_afr_filtered$sample.id, df_afr_filtered$sample.id]

## Add random effects to model
createCategoryMatrix <- function(indices, values) {
  num_indices <- length(indices)
  stopifnot(num_indices == length(values))
  category_matrix_i <- rep(0, num_indices^2)
  category_matrix_j <- rep(0, num_indices^2)
  ii = 0
  for (i in 1:num_indices) {
    for (j in 1:num_indices) {
      if (values[i] == values[j]) {
        ii <- ii + 1
        category_matrix_i[ii] <- i
        category_matrix_j[ii] <- j
      }
    }
  }
  category_matrix_i = category_matrix_i[category_matrix_i != 0]
  category_matrix_j = category_matrix_j[category_matrix_j != 0]
  category_matrix <- Matrix::sparseMatrix(i = category_matrix_i,
                                          j = category_matrix_j,
                                          x = 1,
                                          dims = c(num_indices, num_indices),
                                          dimnames = list(indices, indices))
  
  return(category_matrix)
}

random_cov_mat <- list()

## Generate covariance matrix for other covariates to be coded as random effects
for (i in random_covars){
  print(i)
  # Pull out the column as a vector
  values <- df_transformed_filtered[[i]]  # Use double brackets to extract vector
  random_cov_mat[[i]] <- createCategoryMatrix(df_transformed_filtered$sample.id, values)
}

# -----------------------------
# 5. Set up parameters for analysis
# -----------------------------

# Define output directories
dir.create("results/models", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# Modified two-step analysis function with corrected residual extraction
run_two_step_analysis <- function(data, population_name, cov_matrix, random_cov_matrix = NULL, covariates) {
  # First, make H1H1 explicitly the reference level for the second step
  data$haplotype2 <- as.character(data$haplotype2)  # Convert to character first
  data$haplotype2_ordered <- factor(
    ifelse(data$haplotype2 == "H1_A Carrier", "zH1_A Carrier", data$haplotype2),
    levels = c("H1H1", "H1H2", "H2H2", "zH1_A Carrier")
  )
  
  # Print haplotype distribution for this population
  cat(paste("\nHaplotype distribution for", population_name, "population:\n"))
  print(table(data$haplotype2_ordered, data$diagnosis))
  
  # Create annotated dataframe with our modified data
  sampd <- Biobase::AnnotatedDataFrame(data)
  
  # Set up covariance matrices for random effects
  if (is.null(random_cov_matrix)) {
    cov_mats <- list(grm = cov_matrix)
  } else {
    cov_mats <- c(list(grm = cov_matrix), random_cov_matrix)
  }
  
  # STEP 1: Fit the null model WITHOUT haplotype (just covariates)
  cat(paste("Step 1: Fitting GENESIS model WITHOUT haplotype for", population_name, "...\n"))
  
  # Add try-catch to handle potential errors
  nullmod <- tryCatch({
    model <- fitNullModel(sampd,
                          outcome = outcome,
                          family = glmm_family,
                          covars = covariates,  # No haplotype here!
                          cov.mat = cov_mats,
                          verbose = TRUE,
                          group.var = NULL)
    model
  }, error = function(e) {
    cat(paste("Error in fitting null model:", e$message, "\n"))
    return(NULL)
  })
  
  # Save the model if it was successfully created
  if (!is.null(nullmod)) {
    model_file <- paste0("results/models/nullmod_nohap_", tolower(population_name), ".rds")
    saveRDS(nullmod, model_file)
    
    # STEP 2: Analysis using Cholesky residuals with haplotype as sole predictor
    cat(paste("\nStep 2: Performing Cholesky residual analysis with haplotype for", population_name, "...\n"))
    
    # Extract Cholesky residuals - CORRECTED: they are in nullmod$fit$resid.cholesky
    if (is.null(nullmod$fit$resid.cholesky) || length(nullmod$fit$resid.cholesky) == 0) {
      cat("\nWarning: No Cholesky residuals found in the model. Skipping second step.\n")
      return(list(model = nullmod, lm_model = NULL, lm_summary = NULL))
    }
    
    chol_residuals <- nullmod$fit$resid.cholesky
    
    # Print diagnostic information
    cat(paste("Number of Cholesky residuals:", length(chol_residuals), "\n"))
    cat(paste("Number of samples in data:", nrow(data), "\n"))
    
    # Create a data frame with residuals and haplotype information
    resid_data <- data.frame(
      residual = chol_residuals,
      haplotype = data$haplotype2_ordered
    )
    
    # Run linear model with haplotype as sole predictor
    lm_model <- lm(residual ~ haplotype, data = resid_data)
    
    # Print summary of linear model
    cat(paste("\nLinear model of Cholesky residuals for", population_name, ":\n"))
    print(summary(lm_model))
    
    # Save the linear model
    lm_file <- paste0("results/models/lm_cholesky_", tolower(population_name), "_h1h1ref.rds")
    saveRDS(lm_model, lm_file)
    
    # Extract linear model results
    lm_coef <- summary(lm_model)$coefficients
    
    # Create results dataframe for linear model (exclude intercept)
    if (nrow(lm_coef) > 1) {
      lm_results <- data.frame(
        Population = population_name,
        Term = gsub("^haplotype", "", rownames(lm_coef)[-1]),  # Exclude intercept
        Estimate = lm_coef[-1, "Estimate"],
        StdError = lm_coef[-1, "Std. Error"],
        PValue = lm_coef[-1, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
      
      # Clean up Term names
      lm_results$Term <- gsub("^z", "", lm_results$Term)
      
      # Calculate effect size statistics
      lm_results <- lm_results %>%
        mutate(
          Label = paste0(Population, ": ", Term),
          CI_Lower = Estimate - 1.96 * StdError,
          CI_Upper = Estimate + 1.96 * StdError
        )
    } else {
      lm_results <- data.frame()
      cat("No coefficients found for linear model (possibly singular fit).\n")
    }
    
    return(list(
      model = nullmod, 
      lm_model = lm_model,
      lm_summary = lm_results
    ))
  } else {
    cat("Null model fitting failed. Cannot proceed to second step.\n")
    return(list(model = NULL, lm_model = NULL, lm_summary = NULL))
  }
}

# Modified approach: create random_cov_mat for each population separately
# For European ancestry
random_cov_mat_eur <- list()
for (i in random_covars) {
  values <- df_eur_filtered[[i]]
  random_cov_mat_eur[[i]] <- createCategoryMatrix(df_eur_filtered$sample.id, values)
}

# For African ancestry
random_cov_mat_afr <- list()
for (i in random_covars) {
  values <- df_afr_filtered[[i]]
  random_cov_mat_afr[[i]] <- createCategoryMatrix(df_afr_filtered$sample.id, values)
}

# Run the two-step analysis for each population
# 1. Overall population
overall_results <- run_two_step_analysis(
  data = df_transformed_filtered,
  population_name = "Overall",
  cov_matrix = cov_mat_filtered_overall,
  random_cov_matrix = random_cov_mat,
  covariates = c(overall_covars, pc_cols)
)
overall_lm_summary <- overall_results$lm_summary

# European ancestry
eur_results <- run_two_step_analysis(
  data = df_eur_filtered,
  population_name = "European",
  cov_matrix = cov_mat_filtered_eur,
  random_cov_matrix = random_cov_mat_eur,
  covariates = c(eur_afr_covars, pc_cols)
)
eur_lm_summary <- eur_results$lm_summary

# African ancestry
afr_results <- run_two_step_analysis(
  data = df_afr_filtered,
  population_name = "African",
  cov_matrix = cov_mat_filtered_afr,
  random_cov_matrix = random_cov_mat_afr,
  covariates = c(eur_afr_covars, pc_cols)
)
afr_lm_summary <- afr_results$lm_summary

# Combine the Cholesky residual linear model results
combined_lm_summary <- bind_rows(afr_lm_summary, eur_lm_summary, overall_lm_summary) %>%
  mutate(
    # Create ordering for the plot
    model_order = case_when(
      Population == "African" ~ 1,
      Population == "European" ~ 2,
      Population == "Overall" ~ 3,
      TRUE ~ 4
    ),
    # Define haplotype order
    hap_order = case_when(
      Term == "H1H2" ~ 1,
      Term == "H2H2" ~ 2,
      Term == "H1_A Carrier" ~ 3,
      TRUE ~ 4
    ),
    # Create a combined order variable
    plot_order = model_order * 10 + hap_order,
    # Format p-values for display
    PValue_Display = ifelse(PValue < 0.001, "p < 0.001", 
                            ifelse(PValue < 0.01, paste0("p = ", sprintf("%.3f", PValue)),
                                   paste0("p = ", sprintf("%.2f", PValue)))),
    # Format the label with p-value
    Label_With_P = paste0(Population, ": ", Term, " (", PValue_Display, ")")
  ) %>%
  arrange(plot_order)

# Save the combined summary
write.csv(combined_lm_summary, "results/haplotype_effects_cholesky_summary.csv", row.names = FALSE)

# Create a forest plot for Cholesky residual analysis
cholesky_plot <- ggplot(combined_lm_summary, aes(x = Estimate, y = Label_With_P, color = Population)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  scale_color_manual(values = c("African" = "blue", "European" = "darkgreen", "Overall" = "purple")) +
  labs(
    title = "Haplotype Effects on AD Risk",
    subtitle = "Coefficient Estimates with 95% CI (Reference = H1H1)",
    x = "Coefficient Estimate (95% CI)",
    y = "",
    color = "Population"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Save plot
ggsave("results/plots/haplotype_forest_plot_cholesky_only.png", cholesky_plot, width = 10, height = 6, dpi = 300)

# Filter out African H2H2 before creating the plot
filtered_summary <- combined_lm_summary %>%
  filter(!(Population == "African" & Term == "H2H2"))

# Create a forest plot with the filtered data
cholesky_plot <- ggplot(filtered_summary, aes(x = Estimate, y = Label_With_P, color = Population)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  scale_color_manual(values = c("African" = "blue", "European" = "darkgreen", "Overall" = "purple")) +
  labs(
    title = "Haplotype Effects on AD Risk",
    subtitle = "Coefficient Estimates with 95% CI (Reference = H1H1)",
    x = "Coefficient Estimate (95% CI)",
    y = "",
    color = "Population"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Save plot
ggsave("results/plots/haplotype_forest_plot_cholesky_filtered.png", cholesky_plot, width = 10, height = 6, dpi = 300)

# Extract residual statistics
# Function to extract residual statistics
extract_residual_stats <- function(model_path, population) {
  # Load the model
  nullmod <- readRDS(model_path)
  
  # Extract Cholesky residuals
  if (is.null(nullmod$fit$resid.cholesky) || length(nullmod$fit$resid.cholesky) == 0) {
    cat(paste("\nWarning: No Cholesky residuals found for", population, "model.\n"))
    return(NULL)
  }
  
  chol_residuals <- nullmod$fit$resid.cholesky
  
  # Calculate statistics
  stats <- data.frame(
    Population = population,
    Mean = mean(chol_residuals, na.rm = TRUE),
    SD = sd(chol_residuals, na.rm = TRUE),
    Median = median(chol_residuals, na.rm = TRUE),
    Min = min(chol_residuals, na.rm = TRUE),
    Max = max(chol_residuals, na.rm = TRUE),
    N = length(chol_residuals)
  )
  
  return(stats)
}

# Define paths to saved models
model_files <- c(
  overall = "results/models/nullmod_nohap_overall.rds",
  european = "results/models/nullmod_nohap_european.rds",
  african = "results/models/nullmod_nohap_african.rds"
)

# Extract statistics for each population
residual_stats <- lapply(names(model_files), function(pop) {
  extract_residual_stats(model_files[pop], pop)
}) %>% bind_rows()

# Print the results
print(residual_stats)