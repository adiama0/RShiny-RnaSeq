# test_app.R
source("main.R") # Ensure the main.R file is sourced
source("app.R")
library(testthat)
library(shinytest)

# Test summarized_table function
describe("summarized_table()", {
  summary <- summarized_table(mock_data)
  
  it("should return a tibble", {
    expect_true(is_tibble(summary))
  })
  
  it("should include 'Number of genes'", {
    expect_true(any(summary$`Column name` == "Number of genes"))
  })
  
  it("should correctly calculate the number of genes", {
    number_of_genes_row <- summary %>% filter(`Column name` == "Number of genes")
    expect_equal(as.numeric(number_of_genes_row$`Mean(sd) or Distinct values`), nrow(mock_data))
  })
})

# Test filter_zeros function
describe("filter_zeros()", {
  filtered_data <- filter_zeros(mock_data)
  
  it("should remove rows with zero counts in all numeric columns", {
    expect_equal(nrow(filtered_data), 5)  # Gene A/C/E/G/H remain
  })
  
  it("should retain only rows with non-zero counts in all samples", {
    expect_true(all(rowSums(filtered_data[, -ncol(filtered_data)]) > 0))
  })
})

# Test filter_pivot_data function
describe("filter_pivot_data()", {
  pivoted_data <- filter_pivot_data(filter_zeros(mock_data), min = 10, max = 200)
  
  it("should return a tibble in long format", {
    expect_true(is_tibble(pivoted_data))
    expect_true(all(c("replicates", "counts") %in% colnames(pivoted_data)))
  })
  
  it("should filter rows based on min and max thresholds", {
    expect_equal(nrow(pivoted_data), 30)  # 5 genes x 6 replicates remain
  })
})

# Test get_library_size function
describe("get_library_size()", {
  library_sizes <- get_library_size(mock_data)
  
  it("should return a tibble with sample names and total counts", {
    expect_true(is_tibble(library_sizes))
    expect_equal(colnames(library_sizes), c("sample", "value"))
  })
  
  it("should correctly calculate the total counts for each sample", {
    expected_counts <- colSums(mock_data[, -ncol(mock_data)])
    expect_equal(library_sizes$value, expected_counts)
  })
})

# Filter mock data for normalization functions!
mock_filtered_data <- function(dataf, count_filter) {
  # Separate the gene column
  gene_column <- dataf[[ncol(dataf)]]  # Extract the gene column
  numeric_data <- dataf[, -ncol(dataf)]  # Exclude the gene column for filtering
  
  # Filter rows where the sum of numeric data is greater than or equal to the count filter
  keep_rows <- rowSums(numeric_data) >= count_filter # automatically selected filter
  filtered_numeric_data <- numeric_data[keep_rows, , drop = FALSE]
  
  # Add back the gene column for the filtered rows
  filtered_genes <- gene_column[keep_rows]
  filtered_data <- cbind(filtered_numeric_data, Genes = filtered_genes)
  
  return(filtered_data)
}

# Test normalize_by_cpm function
describe("normalize_by_cpm()", {
  # Filter data using a mock count_filter (e.g., 20)
  filtered_data <- mock_filtered_data(mock_data, count_filter = 10)
  normalized_data <- normalize_by_cpm(filtered_data)
  
  it("should return a data frame with normalized counts", {
    expect_true(is.data.frame(normalized_data)) 
  })
  
  it("should include the original 'Genes' column", {
    expect_true("Genes" %in% colnames(normalized_data))
  })
  
  it("should correctly normalize counts to CPM", {
    library_sizes <- colSums(filtered_data[, -ncol(filtered_data)])  # Total counts per sample
    expected_cpm <- filtered_data$Sample1 / (library_sizes[1] / 1e6)
    expect_equal(normalized_data$Sample1, expected_cpm)
  })
})

# Test normalize_by_deseq function
describe("normalize_by_deseq()", {
  # Filter data using a mock count_filter (e.g., 10)
  filtered_data <- mock_filtered_data(mock_data, count_filter = 10)
  normalized_deseq <- normalize_by_deseq(filtered_data)
  
  it("should return a data frame with normalized counts", {
    expect_true(is_tibble(normalized_deseq))
  })
})

# Test normalize_by_limma function
describe("normalize_by_limma()", {
  # Filter data using a mock count_filter (e.g., 20)
  filtered_data <- mock_filtered_data(mock_data, count_filter = 10)
  normalized_limma <- normalize_by_limma(filtered_data)
  
  it("should return a tibble with normalized counts", {
    expect_true(is_tibble(normalized_limma))
  })
  
  it("should include the original 'Genes' column", {
    expect_true("Genes" %in% colnames(normalized_limma))
  })
})

# Test calculate_stats function
describe("calculate_stats()", {
  # Use mock data and filtered data with count filter = 10
  filtered_data <- mock_filtered_data(mock_data, count_filter = 10)
  stats <- calculate_stats("Test Technique", mock_data, filtered_data)
  
  it("should return a data frame", {
    expect_true(is.data.frame(stats))
  })
  
  it("should correctly calculate the total number of genes", {
    expect_equal(stats$Total.Genes, nrow(mock_data))
  })
  
  it("should correctly calculate the number of passing genes", {
    passing_genes <- nrow(filtered_data)
    passing_percentage <- round((passing_genes / nrow(mock_data)) * 100, 2)
    expected_passing <- passing_genes
    expect_equal(stats$Passing.Genes, expected_passing)
  })
})

# Test normalization_summary function
describe("normalization_summary()", {
  # Step 1: Use mock_filtered_data to filter the mock data
  filtered <- mock_filtered_data(mock_data, count_filter = 10)  # Adjust `count_filter` as needed
  
  # Step 2: Apply normalization functions on the filtered data
  filtered_cpm <- normalize_by_cpm(filtered)
  filtered_deseq <- normalize_by_deseq(filtered)
  filtered_limma <- normalize_by_limma(filtered)
  
  # Step 3: Generate the summary table
  summary <- normalization_summary(
    original_data = mock_data,
    filtered_data_cpm = filtered_cpm,
    filtered_data_deseq = filtered_deseq,
    filtered_data_limma = filtered_limma
  )
  
  # Test cases
  it("should return a data frame summarizing normalization results", {
    expect_true(is.data.frame(summary))
  })
})

describe("run_deseq()", {
  
  # Run the function
  res_df <- run_deseq(mock_data)
  
  # Test if the result is a data frame
  it("should return a data frame", {
    expect_true(is.data.frame(res_df))
  })
  
  # Test if the number of rows matches the number of genes
  it("should have the same number of rows as genes in the input data", {
    expect_equal(nrow(res_df), nrow(mock_data))
  })
  
  # Test if the key columns are present in the result
  it("should contain key columns like log2FoldChange, pvalue, and padj", {
    expect_true(all(c("log2FoldChange", "pvalue", "padj") %in% colnames(res_df)))
  })
  
  # Test if the gene names are assigned as row names
  it("should have gene names as row names", {
    expect_equal(rownames(res_df), mock_data$Genes)
  })
})

describe("run_edger()", {
  # Run the function
  res <- run_edger(mock_data)
  
  # Test if the result is a data frame
  it("should return a data frame", {
    expect_true(is.data.frame(res))
  })
  
  # Test if the number of rows is less than or equal to the number of rows in the input data
  it("should have rows <= the number of rows in the input data due to filtering", {
    expect_lte(nrow(res), nrow(mock_data))
  })
  
  # Test if the result contains key columns
  it("should contain key columns like logFC, PValue, and FDR", {
    expect_true(all(c("logFC", "PValue", "FDR") %in% colnames(res)))
  })
  
  # Test if the filtering step removes some genes
  it("should filter out genes with low expression", {
    countData <- mock_data %>% select(-Genes)
    colData <- data.frame(group = rep(c("N", "P"), each = 3))
    y <- DGEList(counts = countData, group = colData$group)
    keep <- filterByExpr(y)
    filtered_genes <- sum(keep)
    expect_equal(nrow(res), filtered_genes)
  })
})

describe("run_limma()", {
  # Run the function
  res <- run_limma(mock_data)
  
  # Test if the result is a data frame
  it("should return a data frame", {
    expect_true(is.data.frame(res))
  })
  
  # Test if the number of rows is less than or equal to the number of rows in the input data
  it("should have rows <= the number of rows in the input data due to filtering", {
    expect_lte(nrow(res), nrow(mock_data))
  })
  
  # Test if the result contains key columns
  it("should contain key columns like logFC, P.Value, and adj.P.Val", {
    expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(res)))
  })
  
  # Test if the filtering step removes some genes
  it("should filter out genes with low expression", {
    group <- rep(c("N", "P"), each = 3)
    design <- model.matrix(~ group)
    countData <- mock_data %>% select(-Genes)
    dge <- DGEList(counts = countData, group = group)
    keep <- filterByExpr(dge, design)
    filtered_genes <- sum(keep)
    expect_equal(nrow(res), filtered_genes)
  })
})

describe("create_facets()", {
  filtered <- mock_filtered_data(mock_data, count_filter = 10)
  
  deseq_mock <- run_deseq(filtered)
  edger_mock <- run_edger(filtered)
  limma_mock <- run_limma(filtered)
  
  # Run the create_facets function
  result <- create_facets(deseq_mock, edger_mock, limma_mock)
  
  it("should return a data frame", {
    expect_true(is.data.frame(result))
  })
  
  it("Check that the number of rows is correct", {
    expect_equal(nrow(result), nrow(deseq_mock) + nrow(edger_mock) + nrow(limma_mock)) 
  })
})
