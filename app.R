## Author: Aron Diamond
## aronad@bu.edu
## BU BF591
## Final Project

# Required packages
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(tidyverse)
library(bigPint) # for dataset download
library(SummarizedExperiment)
library(DelayedArray)
library(dplyr)
library(tibble)
library(ggridges)
library(edgeR)
library(DESeq2)
library(ggbeeswarm)
library(DT)


ui <- fluidPage(
  titlePanel("Soybean gene expression after iron exposure"),
  p("This experiment assayed gene expression from soybean leaves that were exposed to iron-rich (iron -postive) soil conditions versus leaves that were exposed to iron-poor (iron-negative) soil conditions. The data was collected 120 minutes after iron conditions were initiated"),
  # create 4 to 5 tab panels
  tabsetPanel(
    tabPanel("Sample Information",
      # create sidebar layout
      sidebarLayout(
        sidebarPanel(
          fileInput("upload", "Load counts data"),
          sliderInput("filter_raw", "Raw count filter", min = 10, max = 16000, value = c(10,16000), step = 1)
        ),
        mainPanel(
          # create 3 tab panels
          tabsetPanel(
            tabPanel("Sample Summary",
              tableOutput("table_summary")
            ),
            tabPanel("Data Table",
              tableOutput("table_unsorted"),
              div(style = "text-align: center;", 
                  actionButton("sort_button", "Sort Table"),
                  actionButton("unsort_button", "Unsort Table"),
                  actionButton("remove_button", "Remove 0 counts")),
              div(style = "text-align: center;", 
                  textOutput("sort_message"),
                  textOutput("unsort_message"),
                  textOutput("remove_message"))
            ),
            tabPanel("Plots",
              plotOutput("histogram_raw"),
              plotOutput("violin_raw"),
              plotOutput("ridged")
            )
          )
        )
      )
    ),
    tabPanel("Count Matrix",
      # create sidebar for input filters and output display
      sidebarLayout(
        sidebarPanel(
          h3("Normalization"),
          numericInput("count_filter", "Choose minimum count filter", 10),
          radioButtons("normalize_method", 
            "Choose your normalization method", 
            c("CPM", "DESeq", "Limma")
            ),
          numericInput("num_pcs", 
            "Select the number of principal components (PCs) to plot:", 
            3,
            min = 1,
            max = 6)
        ),
        mainPanel(
          # create tabset for 3 outputs
          tabsetPanel(
            tabPanel("Summary",
              tableOutput("normalization_summary"),
              tableOutput("normalized_table")
            ),
            tabPanel("Diagnostic Plots",
              plotOutput("boxplot_distribution"),
              plotOutput("variance_vs_mean")
            ),
            tabPanel("Heatmap",
              plotOutput("heatmap")
            ),
            tabPanel("PCA",
              plotOutput("variance_explained"),
              plotOutput("PCA_plot"),
              plotOutput("beeswarm")
            )
          )
        )
      )
    ),
    tabPanel("Differential Expression Analysis",
      # create sidebar layout, See from assignment 8
      sidebarLayout(
        sidebarPanel(
          radioButtons("difex_method",
            "Differential Expression Method",
            c("DESeq", "Limma", "EdgeR")
          ),
          actionButton("run_difex", 
                      "Run Differential Expression"
          ),
          sliderInput("pval", 
                      "Set P-value Threshold (-log10 scale):", 
                      min = 1, max = 100, value = 25, step = 1
          )
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Table",
              p('Differentially expressed genes'),
              DT::DTOutput("difex_genes")
            ),
            tabPanel("Plot",
              plotOutput("volcanoes")
            )
          )
        )
      )
    )
  )
)



server <- function(input, output, session) {
  
  ################### Sample Information ###################
  
  # Load data from uploaded file, including gene names as a column
  load_data <- reactive({
    req(input$upload)
    data <- read.csv(input$upload$datapath, row.names = 1)
    # Convert rownames to a column called "Genes"
    data$Genes <- rownames(data)
    # Rename any column that might be incorrectly labeled as 'X' to 'Genes'
    colnames(data)[colnames(data) == 'X'] <- 'Genes'
    return(data)
  })
  
  #' Create summary table from counts dataframe 
  #' @param dataf The loaded data frame
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
    
    summary_tibble <- bind_rows(summary_tibble, number_of_genes_row)
    return(summary_tibble)
  }
  
  #' Sort the loaded dataframe of RAW COUNTS
  #' @param dataf The loaded data frame
  # Extract gene names from rownames and add them as a column
  sort_data <- function(dataf) {
    
    # Extract gene names from rownames and add them as a column
    gene_names <- rownames(dataf)
    dataf$Genes <- gene_names
    # Exclude non-numeric columns when calculating row sums
    numeric_data <- dataf[sapply(dataf, is.numeric)]
    total_counts <- rowSums(numeric_data)
    # Sort the data frame by the total counts in descending order
    sorted_dataf <- dataf[order(total_counts, decreasing = TRUE), ]
    
    # Return the sorted data frame
    return(sorted_dataf)
  }
  
  # Filter loaded dataframe of RAW COUNTS 
  filter_zeros <- function(dataf) {
    
    numeric_data <- dataf[, -ncol(dataf)]  # Exclude the last column (Genes)
    filtered_data <- dataf[apply(numeric_data, 1, function(row) all(row != 0)), ]
    
    return(filtered_data)
  }
  
  # Filter and picot data of RAW COUNTs based on min and max filters for plots
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
  
  # filter RAW COUNTs based on interactive filters 
  filtered_long_data <- reactive({
    req(input$upload)
    data <- load_data()
    
    # Get slider values
    min_count <- input$filter_raw[1]
    max_count <- input$filter_raw[2]
    
    # Apply filtering and pivoting
    filter_pivot_data(data, min_count, max_count)
  })
  
  ################### Count Matrix ###################
  
  # Filter data based on count_filter input for NORMALIZATION 
  filtered_data <- reactive({
    req(input$count_filter)
    data <- load_data()
    
    # Separate the gene column
    gene_column <- data[[ncol(data)]]  # Extract the gene column
    numeric_data <- data[, -ncol(data)]  # Exclude the gene column for filtering
    
    # Filter rows where the sum of numeric data is greater than or equal to the count filter
    keep_rows <- rowSums(numeric_data) >= input$count_filter
    filtered_numeric_data <- numeric_data[keep_rows, , drop = FALSE]
    
    # Add back the gene column for the filtered rows
    filtered_genes <- gene_column[keep_rows]
    filtered_data <- cbind(filtered_numeric_data, Genes = filtered_genes)
    
    return(filtered_data)
  })
  
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
    
  # Function for NORMALIZATION by deseq
  normalize_by_deseq <- function(dataf) {
    # create coldata 
    colData <- data.frame( 
      condition = rep(c("N", "P"), each=3)
    )
    row.names(colData) <- colnames(dataf)[-ncol(dataf)]
    
    # create dds object (dataf already filtered by input)
    dds <- DESeqDataSetFromMatrix(countData = dataf[,-ncol(dataf)],
                                  colData = colData,
                                  design= ~ condition)
    # run DESeq on dds object
    dds <- DESeq(dds)
    
    # retreive results
    res <- results(dds, name="condition_P_vs_N")
    # Extract normalized count data
    normalized_counts <- counts(dds, normalized = TRUE)
    # Convert to a data frame and add the gene names
    normalized_counts_df <- as_tibble(normalized_counts)
    # Extract the gene names from the 'Genes' column
    normalized_counts_df$Genes <- dataf$Genes
    
    return(normalized_counts_df)
  }
  
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
  
  # Reactive function for NORMALIZATION 
  normalized_data <- reactive({
    req(input$normalize_method)
    data <- filtered_data()
    
    if (input$normalize_method == "CPM") {
      return(normalize_by_cpm(data))
    } else if (input$normalize_method == "DESeq") {
      return(normalize_by_deseq(data))
    } else if (input$normalize_method == "Limma") {
      return(normalize_by_limma(data))
    } else {
      return(data)
    }
  })
  
  # Create summarization table for NORMALIZATION counts 
  normalization_summary <- function(original_data, filtered_data_cpm, filtered_data_deseq, filtered_data_limma) {
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
        `Genes Passing Filter` = paste(passing_genes, "(", passing_percentage, "%)"),
        `Genes Not Passing Filter` = paste(failing_genes, "(", failing_percentage, "%)")
      )
    }
    
    # Combine the summary for all techniques
    summary_cpm <- calculate_stats("CPM", original_data, filtered_data_cpm)
    summary_deseq <- calculate_stats("DESeq2", original_data, filtered_data_deseq)
    summary_limma <- calculate_stats("Limma", original_data, filtered_data_limma)
    
    rbind(summary_cpm, summary_deseq, summary_limma)
  }
  
  # create reactive event for NORMALIZATION summary table 
  normalized_summaries <- reactive({
    req(input$upload)
    data <- load_data()
    filtered <- filtered_data()
    
    summary_table <- normalization_summary(
      original_data = data,
      filtered_data_cpm = normalize_by_cpm(filtered),
      filtered_data_deseq = normalize_by_deseq(filtered),
      filtered_data_limma = normalize_by_limma(filtered)
    )
    
    summary_table
  })
  
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
  
  plot_pca_variance <- function(data) {
    # Transpose the data to have samples as rows and genes as columns
    data_matrix <- as.matrix(data %>% select(-Genes))  # Remove Genes column
    
    # Perform PCA
    pca_results <- prcomp(scale(t(data_matrix), center=TRUE, scale=TRUE))
    
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
  
  plot_pca_scatter <- function(data) {
    
    # create matrix for data
    data_matrix <- as.matrix(data %>% select(-Genes))
    
    # create metadata table
    sample_names <- colnames(data_matrix)
    timepoint <- substr(sample_names, 1, 1)
    replicate <- substr(sample_names, 3, 4)
    meta <- tibble(sample = sample_names,
                   timepoint = timepoint,
                   replicate = replicate)
    
    # perform PCA on the data
    pca_results <- prcomp(scale(t(data_matrix), center=TRUE, scale=TRUE))
    pc_scores <- as_tibble(pca_results$x[, 1:2])
    pc_scores <- pc_scores %>%
      mutate(sample = meta$timepoint[match(row.names(pca_results$x), meta$sample)])
    
    variance_explained <- summary(pca_results)$importance[2, ] * 100
    
    # Create the PCA plot
    plot <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = sample)) +
      geom_point(size = 4, alpha = 0.7)+
      labs(
        title = "PCA Plot",
        x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
        y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
        color = "Sample"
      ) +
      theme_minimal()
   
   return(plot)
  }
  
  pca_scores <- reactive({
    req(normalized_data())  # Ensure normalized data is available
    data <- as.matrix(normalized_data() %>% select(-Genes))
    pca_result <- prcomp(t(data), center = TRUE, scale. = TRUE)
    return(pca_result)  # Return PCA scores as a data frame
  })
  
  plot_beeswarm <- function(pca_data, num_pcs) {
    
    # Extract the first `num_pcs` principal components as a tibble
    pc_scores <- as_tibble(pca_data$x[, 1:num_pcs])
    
    # Pivot longer to reshape the data
    pcs_longer <- pivot_longer(
      pc_scores, 
      names_to = "PC", 
      values_to = "values", 
      cols = everything() # Selects all columns in pc_scores
    )
    
    # Generate beeswarm plot
    ggplot(pcs_longer, aes(x = PC, y = values, color = PC)) +
      geom_beeswarm(cex = 1.5) +
      theme_minimal() +
      labs(
        title = paste("Beeswarm Plot for Top", num_pcs, "Principal Components"),
        x = "Principal Components",
        y = "Value"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  ################### DIFFERENTIAL EXPRESSION ############
  
  # Differential Expression Functions
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
  
  # Reactive for differential expression
  difex_results <- eventReactive(input$run_difex, {
    req(normalized_data())
    data <- filtered_data()
    if (input$difex_method == "DESeq") {
      return(run_deseq(data))
    } else if (input$difex_method == "EdgeR") {
      return(run_edger(data))
    } else if (input$difex_method == "Limma") {
      return(run_limma(data))
    }
  })
  
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
  
  # Function to create a styled volcano plot
  theme_plot <- function(volcano_data, threshold) {
    volcano_data <- volcano_data %>% 
      mutate(status = ifelse(-log10(pval) > threshold, "Significant", "Not Significant"))
    
    ggplot(volcano_data, aes(x = logFC, y = -log10(pval), color = status)) +
      geom_point(alpha = 0.7) +
      facet_wrap(~ package, scales = "free") +
      labs(
        x = "Log2 Fold Change",
        y = "-Log10(P-value)",
        title = "Volcano Plots for Differential Expression Methods"
      ) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      theme_minimal() +
      scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "gray"))
  }
  
  # Reactive to create combined differential expression results
  combined_volcano_data <- reactive({
    req(difex_results())
    
    # Assuming the three methods' results are available
    deseq_res <- run_deseq(filtered_data())
    edger_res <- run_edger(filtered_data())
    limma_res <- run_limma(filtered_data())
    
    # Combine the results
    create_facets(deseq_res, edger_res, limma_res)
  })
  
  
  ################### RENDER FUNCTIONS ################### 
  # Render summary table
  output$table_summary <- renderTable({
    req(input$upload)
    summarized_table(load_data())
  })
  
  # Render unsorted table initially
  output$table_unsorted <- renderTable({
    req(input$upload)
    head(load_data(), 10)
  })
  
  # Display sorted table after button is pressed
  observeEvent(input$sort_button, {
    output$table_unsorted <- renderTable({
      req(input$upload)
      head(sort_data(load_data()), 10)
    })
    
    # Display sorting message after table is sorted
    output$sort_message <- renderText({
      "The table is now sorted."
    })
  })
  
  # Display original unsorted table 
  observeEvent(input$unsort_button, {
    output$table_unsorted <- renderTable({
      req(input$upload)
      head(load_data(), 10)
    })
    
    # Display sorting message after table is sorted
    output$unsort_message <- renderText({
      "You are viewing the original unsorted table"
    })
  })
  
  # Display non-zero table
  observeEvent(input$remove_button, {
    output$table_unsorted <- renderTable({
      req(input$upload)
      head(filter_zeros(load_data()), 10)
    })
    
    # Display filtering message after table is updated
    output$remove_message <- renderText({
      "The table is now filtered to remove rows containing zeros."
    })
  })
  
  output$histogram_raw <- renderPlot({
    req(filtered_long_data())
    ggplot(filtered_long_data(), aes(x = counts)) +
      geom_histogram() +
      labs(title = "Filtered Data Distribution", x = "Counts", y = "Frequency")
  })
  
  output$violin_raw <- renderPlot({
    req(filtered_long_data())
    ggplot(filtered_long_data(), aes(x = replicates, y = counts)) +
      geom_violin() +
      labs(title = "Violin Plot of Counts Across Replicates", x = "Replicates", y = "Counts") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  })
  
  output$ridged <- renderPlot({
    req(filtered_long_data())
    ggplot(filtered_long_data(), aes(x = counts, y = replicates, fill = replicates)) +
      geom_density_ridges(alpha = 0.7, scale = 1.5) +  # Create ridgeline plot
      scale_fill_viridis_d() +  # Optional: Add a color scale
      labs(title = "Ridgeline Density Plot of Counts Across Replicates", 
           x = "Counts", y = "Replicates") +
      theme_minimal() +
      theme(legend.position = "none") 
  })
  
  output$normalized_table <- renderTable({
    req(normalized_data())
    head(normalized_data(), 10)  
  })
  
  output$normalization_summary <- renderTable({
    req(normalized_summaries())
    normalized_summaries()
  })
  
  output$boxplot_distribution <- renderPlot({
    req(input$normalize_method)
    plot_boxplot(normalized_data(), input$normalize_method)
  })
  
  output$variance_vs_mean <- renderPlot({
    req(input$normalize_method)
    plot_variance_vs_mean(normalized_data(), input$normalize_method)
  })
  
  output$heatmap <- renderPlot({
    req(normalized_data())
    # Pass the normalized data to the heatmap function
    generate_heatmap(normalized_data())
  }, height = 900, width = 600)
  
  output$variance_explained <- renderPlot({
    req(normalized_data())
    plot_pca_variance(normalized_data())
  })
  
  output$PCA_plot <- renderPlot({
    req(normalized_data())   
    plot_pca_scatter(normalized_data())
  })
  
  output$beeswarm <- renderPlot({
    req(pca_scores())  # Ensure PCA scores are available
    num_pcs <- input$num_pcs  # Get the selected number of PCs
    plot_beeswarm(pca_scores(), num_pcs)  # Call the updated beeswarm function
  })
  
  # Output: Differential Expression Table
  output$difex_genes <- DT::renderDT({difex_results()})
  
  # Output: Volcano Plot
  output$volcanoes <- renderPlot({
    req(combined_volcano_data(), input$pval)  # Ensure data and slider input are available
    theme_plot(combined_volcano_data(), input$pval)  # Pass the slider value
  })
}

# Run the application
shinyApp(ui = ui, server = server)
