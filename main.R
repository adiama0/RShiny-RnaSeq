# Import required packages
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(tidyverse)
library(SummarizedExperiment)
library(DelayedArray)
library(dplyr)
library(tibble)
library(ggridges)
library(edgeR)
library(DESeq2)
library(DT)
library(pheatmap)
library(viridis)

# Summarized Table
# Create summary table from counts dataframe 
summarized_table <- function(dataf) {
  columns <- colnames(dataf)
  condition <- sapply(strsplit(columns, "\\."), `[`, 1)
  replicates <- sapply(strsplit(columns, "\\."), `[`, 2)
  coldata <- data.frame(condition = condition, replicate = replicates, stringsAsFactors = FALSE)
  coldata <- coldata[complete.cases(coldata), , drop = FALSE]
  
  summary_tibble <- tibble(
    `Column name` = colnames(coldata),
    Type = sapply(coldata, class),
    `Mean(sd) or Distinct values` = sapply(coldata, function(col) paste(unique(col), collapse = ", "))
  )
  
  number_of_genes_row <- tibble(
    `Column name` = "Number of genes",
    Type = class(nrow(dataf)),
    `Mean(sd) or Distinct values` = as.character(nrow(dataf))
  )
  
  bind_rows(summary_tibble, number_of_genes_row)
}

# Filter Rows with Zero Values
# Function to filter out rows with zero counts
filter_zeros <- function(dataf) {
  numeric_data <- dataf[, -ncol(dataf)]  # Exclude the last column (Genes)
  filtered_data <- dataf[apply(numeric_data, 1, function(row) all(row != 0)), ]
  return(filtered_data)
}

# Filter and count data of RAW COUNTs based on min and max filters for PLOTS
filter_pivot_data <- function(dataf, min, max) {
  
  # Filter rows where the sum of the counts across samples 
  filtered_data <- dataf %>% 
    filter(rowSums(.[, -ncol(dataf)]) >= min & rowSums(.[, -ncol(dataf)]) <= max)
  # Pivot the filtered data to a longer format
  long_data <- pivot_longer(filtered_data, 
                            cols = -Genes, 
                            names_to = "replicates", 
                            values_to = "counts")
  
  return(long_data)
}

# Get Library Size
# Create table of total counts for NORMALIZATION by CPM
get_library_size <- function(dataf) {
  # Exclude the gene column and calculate total counts for each sample
  total_counts <- colSums(dataf[, -ncol(dataf)])
  
  # Create a tibble with sample names and their total counts
  plot_data <- tibble(
    sample = names(total_counts),
    value = total_counts
  )
  
  return(plot_data)
}

# Normalize by CPM
# function for NORMALIZATION by CPM
normalize_by_cpm <- function(dataf) {
  
  # Get the library sizes using the helper function
  library_sizes <- get_library_size(dataf)$value
  # Exclude the gene column for normalization
  gene_column <- dataf[[ncol(dataf)]]
  # Normalize counts to CPM
  normalized_data <- dataf[, -ncol(dataf)] %>%
    mutate(across(everything(), ~ .x / (library_sizes[cur_column()] / 1e6)))
  # Add back the gene column
  normalized_data <- cbind(normalized_data, Genes = gene_column)
  
  return(normalized_data)
}

# Normalize by DESeq
# Function for NORMALIZATION by deseq
normalize_by_deseq <- function(dataf) {
  
  # Extract count data (excluding the "Genes" column)
  countData <- dataf[, -ncol(dataf)]
  
  # Create colData with the same number of rows as countData columns
  colData <- data.frame(
    condition = rep(c("N", "P"), each=3)
  )
  row.names(colData) <- colnames(countData)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(countData),
                                colData = colData,
                                design = ~ condition)
  
  # Normalize counts using DESeq2
  dds <- DESeq(dds)
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Convert to a tibble and add "Genes" column
  normalized_counts_df <- as_tibble(normalized_counts)
  
  # Add the Genes column as the first column
  normalized_counts_df$Genes <- dataf$Genes
  
  # Reorder columns so that "Genes" is the first column
  normalized_counts_df <- normalized_counts_df %>%
    select(Genes, everything())
  
  return(normalized_counts_df)
}



# Normalize by Limma
# function for NORMALIZATION by Limma 
normalize_by_limma <- function(dataf) {
  
  # Create colData (condition for each sample)
  colData <- data.frame(
    condition = factor(rep(c("N", "P"), length.out = ncol(dataf) - 1))  # Dynamically adjusts conditions
  )
  row.names(colData) <- colnames(dataf)[-ncol(dataf)]  # Exclude 'Genes' column
  
  # Prepare the count matrix by removing the 'Genes' column
  count_data <- as.matrix(dataf[, -ncol(dataf)])  # Convert to matrix for voom compatibility
  # Use the voom function from limma for normalization
  v <- voom(count_data, design = model.matrix(~condition, data = colData), plot = TRUE)  # Plot diagnostics
  # Extract the normalized log2 counts per million (log-CPM) from the voom object
  normalized_counts <- v$E
  # Convert to a data frame and add gene names
  normalized_counts_df <- as.data.frame(normalized_counts)
  normalized_counts_df$Genes <- dataf$Genes
  # Reorder columns to put 'Genes' first
  normalized_counts_df <- normalized_counts_df[, c("Genes", setdiff(colnames(normalized_counts_df), "Genes"))]
  # Convert back to tibble if necessary
  normalized_counts_df <- as_tibble(normalized_counts_df)
  
  return(normalized_counts_df)
}

# Normalization Summary

# Calculate Statistics
calculate_stats <- function(name, data, filtered_data) {
  total_genes <- nrow(data)
  total_samples <- ncol(data) - 1  # Exclude the "Genes" column
  passing_genes <- nrow(filtered_data)
  failing_genes <- total_genes - passing_genes
  
  passing_percentage <- round((passing_genes / total_genes) * 100, 2)
  failing_percentage <- round((failing_genes / total_genes) * 100, 2)
  
  data.frame(
    Technique = name,
    `Number of Samples` = total_samples,
    `Total Genes` = total_genes,
    `Passing Genes` = passing_genes,
    `Percent passing` = passing_percentage
  )
}

# Create summarization table for NORMALIZATION counts 
normalization_summary <- function(original_data, filtered_data_cpm, filtered_data_deseq, filtered_data_limma) {
  # Combine the summary for all techniques
  summary_cpm <- calculate_stats("CPM", original_data, filtered_data_cpm)
  summary_deseq <- calculate_stats("DESeq2", original_data, filtered_data_deseq)
  summary_limma <- calculate_stats("Limma", original_data, filtered_data_limma)
  
  # Combine results
  rbind(summary_cpm, summary_deseq, summary_limma)
}

# Plot Boxplot
plot_boxplot <- function(dataf, normalization_method) {
  data_long <- pivot_longer(
    dataf, 
    cols = -Genes, 
    names_to = "Samples", 
    values_to = "Counts"
  )
  ggplot(data_long, aes(x = Samples, y = log10(Counts + 1), fill = Samples)) +
    geom_boxplot() +
    labs(
      title = paste("Boxplot of", normalization_method, "Normalized Data"),
      x = "Samples",
      y = "Log10(Counts)"
    )
}

# Plot Variance vs Mean
plot_variance_vs_mean <- function(data, normalization_method) {
  # Calculate the mean and the variance
  variance <- data %>% select(-Genes) %>% apply(1, var)
  mean <- data %>% select(-Genes) %>% apply(1, mean)
  
  # create a tibble for the plot 
  variancevsmean <- tibble(
    variance = variance,
    mean = mean
  )
  
  # add rownames in case of labeling
  rownames(variancevsmean) <- data$Genes
  
  # plot and add smooth for average line
  plot <- variancevsmean %>% ggplot() +
    geom_point(aes(x=rank(mean), y=log10(variance))) +
    geom_smooth(aes(x = rank(mean), y = log10(variance))) +
    ggtitle(paste("Boxplot of", normalization_method, "Normalized Data"))
  
  return(plot)
}

# Generate Heatmap
generate_heatmap <- function(data) {
  # Prepare the normalized data matrix
  normalized_matrix <- data %>%
    as.data.frame() %>%  # Ensure it's a data frame
    select(-Genes) %>%   # Remove the Genes column
    `rownames<-`(data$Genes) %>%  # Set row names
    as.matrix()  # Convert to matrix
  
  # Log-transform the data
  log_matrix <- log2(normalized_matrix + 1)
  
  # Generate the heatmap
  heatmap(
    log_matrix,
    Rowv = TRUE,  # Hierarchical clustering for rows
    Colv = TRUE,  # Hierarchical clustering for columns
    col = colorRampPalette(c("blue", "white", "red"))(100),  # Color scale
    scale = "none",  # Data is already normalized/log-transformed
    labRow = FALSE,  # Remove gene names from the y-axis
    margins = c(5, 5),  # Adjust margins
    main = "Clustered Heatmap of Normalized Counts",
    xlab = "Samples",
    ylab = "Genes"
  )
  
  # Add a legend manually
  par(xpd = TRUE)  # Allow plotting outside plot area
  legend(
    "top",
    legend = c("Low", "Median", "High"),
    fill = c("blue", "white", "red"),
    title = "Expression Level",
    border = "black",
    box.col = "black"
  )
}

# Plot PCA Variance
plot_pca_variance <- function(data) {
  # Transpose the data to have samples as rows and genes as columns
  data_matrix <- as.matrix(data %>% select(-Genes))  # Remove Genes column
  
  # Perform PCA
  pca_results <- prcomp(scale(t(data_matrix), center = TRUE, scale = TRUE))
  
  # Calculate variance explained
  variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2)
  
  # Create a tibble for plotting
  variance_df <- tibble(
    PC = factor(paste0("PC", seq_along(variance_explained))),
    Variance = variance_explained
  )
  
  # Plot the variance explained by each PC
  ggplot(variance_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(
      title = "Variance Explained by Principal Components",
      x = "Principal Component",
      y = "Variance Explained"
    ) +
    theme_minimal()
}

# Plot PCA Scatter
plot_pca_scatter <- function(data) {
  # Create matrix for data
  data_matrix <- as.matrix(data %>% select(-Genes))
  
  # Create metadata table
  sample_names <- colnames(data_matrix)
  timepoint <- substr(sample_names, 1, 1)
  replicate <- substr(sample_names, 3, 4)
  meta <- tibble(sample = sample_names,
                 timepoint = timepoint,
                 replicate = replicate)
  
  # Perform PCA on the data
  pca_results <- prcomp(scale(t(data_matrix), center = TRUE, scale = TRUE))
  pc_scores <- as_tibble(pca_results$x[, 1:2])
  pc_scores <- pc_scores %>%
    mutate(sample = meta$timepoint[match(row.names(pca_results$x), meta$sample)])
  
  # Calculate variance explained
  variance_explained <- summary(pca_results)$importance[2, ] * 100
  
  # Create the PCA plot
  plot <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = sample)) +
    geom_point(size = 4, alpha = 0.7) +
    labs(
      title = "PCA Plot",
      x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
      y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
      color = "Sample"
    ) +
    theme_minimal()
  
  return(plot)
}


# Differential Expression Functions

# Helper function to run differential expression based on the selected method
run_difex <- function(method, dataf) {
  if (method == "DESeq") {
    return(run_deseq(dataf))
  } else if (method == "EdgeR") {
    return(run_edger(dataf))
  } else if (method == "Limma") {
    return(run_limma(dataf))
  } else {
    stop("Invalid differential expression method selected.")
  }
}

# DESeq2 differential expression
run_deseq <- function(dataf) {
  colData <- data.frame(condition = rep(c("N", "P"), each = 3))
  rownames(colData) <- colnames(dataf)[-ncol(dataf)]
  countData <- dataf %>% select(-Genes)
  dds <- DESeqDataSetFromMatrix(countData = countData, 
                                colData = colData, 
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, name = "condition_P_vs_N")
  # Assign Genes as row names
  res_df <- as.data.frame(res)
  rownames(res_df) <- dataf$Genes
  return(res_df)
}

# EdgeR differential expression
run_edger <- function(dataf) {
  colData <- data.frame(group = rep(c("N", "P"), each = 3))
  countData <- dataf %>% select(-Genes)
  y <- DGEList(counts = countData, group = colData$group)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y)
  et <- exactTest(y)
  res <- topTags(et, n = nrow(y))$table
  return(res)
}

# Limma differential expression
run_limma <- function(dataf) {
  group <- rep(c("N", "P"), each = 3)
  design <- model.matrix(~ group)
  countData <- dataf %>% select(-Genes)
  dge <- DGEList(counts = countData, group = group)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE)
  res <- topTable(fit, coef = ncol(design), n = dim(fit)[1], sort.by = "P")
  return(res)
}

# Create Facets
# Function to combine differential expression results into a long format
create_facets <- function(deseq, edger, limma) {
  # Extract p-values and log fold changes
  deseq_df <- data.frame(
    logFC = deseq$log2FoldChange,
    pval = deseq$pvalue,
    package = "DESeq2"
  )
  edger_df <- data.frame(
    logFC = edger$logFC,
    pval = edger$PValue,
    package = "EdgeR"
  )
  limma_df <- data.frame(
    logFC = limma$logFC,
    pval = limma$P.Value,
    package = "Limma"
  )
  
  # Combine all three data frames
  combined_df <- dplyr::bind_rows(deseq_df, edger_df, limma_df)
  return(combined_df)
}

# Custom Theme for Plot
# Function to create a styled volcano plot
theme_plot <- function(volcano_data, threshold) {
  volcano_data <- volcano_data %>% 
    mutate(status = ifelse(-log10(pval) > threshold, "Significant", "Not Significant")) %>%
    na.omit()
  
  ggplot(volcano_data, aes(x = logFC, y = -log10(pval), color = status)) +
    geom_point(alpha = 0.7) +
    facet_wrap(~ package, scales = "free") +
    labs(
      x = "Log2 Fold Change",
      y = "-Log10(P-value)",
      title = "Volcano Plots for Differential Expression Methods"
    ) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
    theme_minimal() +
    scale_color_manual(values = c("Significant" = "#DB6389", "Not Significant" = "#85E7F2"))
}
