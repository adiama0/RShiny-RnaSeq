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

ui <- fluidPage(
  titlePanel("Soybean cotyledon gene expression across development"),
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
              div(style = "text-align: center;", actionButton("sort_button", "Sort Table"),
                  actionButton("unsort_button", "Unsort Table"),
                  actionButton("remove_button", "Remove 0 counts")),
              div(style = "text-align: center;", textOutput("sort_message"),
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
          radioButtons("normalize_method", "Choose your normalization method", c("CPM", "DESeq", "Limma")),
          p("Slider to include genes with at least X percentile of variance"),
          p("Slider to include genes with at least X samples that are non-zero")
        ),
        mainPanel(
          # create tabset for 3 outputs
          tabsetPanel(
            tabPanel("Summary",
              tableOutput("normalization_summary"),
              tableOutput("CPM")
            ),
            tabPanel("Diagnostic Plots"),
            tabPanel("Heatmap"),
            tabPanel("PCA")
          )
        )
      )
    ),
    tabPanel("Differential Expression Analysis",
      # create sidebar layout, See from assignment 8
      sidebarLayout(
        sidebarPanel(
          p("Sortable filters displaying differential expression results")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Table",
              p("sorted table")
            ),
            tabPanel("Plot",
              p("Volcano Plot")
            )
          )
        )
      )
    )
  )
)


server <- function(input, output, session) {
  
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
  
  normalize_by_limma <- function(dataf) {
    # Create colData (condition for each sample)
    colData <- data.frame(
      condition = rep(c("N", "P"), each = 3)  # Adjust based on your experimental design
    )
    row.names(colData) <- colnames(dataf)[-ncol(dataf)]  # Exclude 'Genes' column
    
    # Prepare the count matrix by removing the 'Genes' column
    count_data <- dataf[, -ncol(dataf)]
    
    # Use the voom function from limma for normalization
    # voom() normalizes the counts to log2 counts per million (log-CPM)
    v <- voom(count_data, design = model.matrix(~condition, data = colData), plot = FALSE)
    # Extract the normalized log2 counts per million (log-CPM) from the voom object
    normalized_counts <- (v$E)*10 # multiply by 10 to maintain same x axis
    # Convert to a data frame
    normalized_counts_df <- as_tibble(normalized_counts)
    # Add the gene names from the 'Genes' column
    normalized_counts_df$Genes <- dataf$Genes
    
    return(normalized_counts_df)
  }
  
  # Reactive function for NORMALIZATION 
  normalized_data <- reactive({
    req(input$normalize_method)
    data <- filtered_data()
    
    if (input$normalize_method == "CPM") {
      return(head(normalize_by_cpm(data), 10))
    } else if (input$normalize_method == "DESeq") {
      return(head(normalize_by_deseq(data), 10))
    } else if (input$normalize_method == "Limma") {
      return(head(normalize_by_limma(data), 10))
    } else {
      return(head(data), 10)
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
  
  output$CPM <- renderTable({
    req(normalized_data())
    head(normalized_data(), 10)  
  })
  
  output$normalization_summary <- renderTable({
    req(normalized_summaries())
    normalized_summaries()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
