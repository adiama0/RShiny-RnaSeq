## Author: Aron Diamond
## aronad@bu.edu
## BU BF591
## Final Project

source("main.R")

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
                 sliderInput("filter_raw", "Raw count filter", min = 10, max = 16000, value = c(10,16000), step = 10)
               ),
               mainPanel(
                 # create 3 tab panels
                 tabsetPanel(
                   tabPanel("Sample Summary",
                            DT::DTOutput("table_summary") 
                   ),
                   tabPanel("Data Table",
                            DT::DTOutput("table_unsorted"),  # Use DTOutput
                            div(style = "text-align: center;", 
                                actionButton("remove_button", "Remove 0 counts")),
                            div(style = "text-align: center;", 
                                textOutput("remove_message"))
                   ),
                   tabPanel("Plots",
                            fluidRow(
                              column(6, plotOutput("histogram_raw")),  # Left plot
                              column(6, plotOutput("violin_raw"))     # Right plot
                            ),
                            fluidRow(
                              column(12, plotOutput("ridged"))        # Bottom plot spanning full width
                            )
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
                            DT::DTOutput("normalization_summary"),
                            DT::DTOutput("normalized_table")
                   ),
                   tabPanel("Diagnostic Plots",
                            fluidRow(
                              column(6, plotOutput("boxplot_distribution")),  # Left plot
                              column(6, plotOutput("variance_vs_mean"))     # Right plot
                            )    
                   ),
                   tabPanel("Heatmap",
                            plotOutput("heatmap")
                   ),
                   tabPanel("PCA",
                            fluidRow(
                              column(6,  plotOutput("beeswarm")),  # Left plot
                              column(6, plotOutput("PCA_plot"))     # Right plot
                            ),
                            fluidRow(
                              column(12, plotOutput("variance_explained"))        # Bottom plot spanning full width
                            )
                   )
                 )
               )
             )
    ),
    tabPanel("Differential Expression Analysis",
             # create sidebar layout, See from assignment 8
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 radioButtons("difex_method",
                              "Differential Expression Method",
                              c("DESeq", "Limma", "EdgeR")
                 ),
                 actionButton("run_difex", 
                              "Run Differential Expression"
                 ),
                 sliderInput("pval", 
                             "Set P-value Threshold (-log10 scale):", 
                             min = 1, max = 50, value = 25, step = 1
                 )
               ),
               mainPanel(
                 width = 9,
                 tabsetPanel(
                   tabPanel("Table",
                            DT::DTOutput("difex_genes")
                   ),
                   tabPanel("Plot",
                            fluidRow(
                              column(12,  plotOutput("volcanoes"))
                            ),
                            fluidRow(
                              column(6,
                                     p("DESeq Results"),
                                     DT::DTOutput("difex_genes_deseq")),
                              column(6, 
                                     p("EdgeR Results"),
                                     DT::DTOutput("difex_genes_edger"))
                            )
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
  
  pca_scores <- reactive({
    req(normalized_data())  # Ensure normalized data is available
    data <- as.matrix(normalized_data() %>% select(-Genes))
    pca_result <- prcomp(t(data), center = TRUE, scale. = TRUE)
    return(pca_result)  # Return PCA scores as a data frame
  })
  
  ################### DIFFERENTIAL EXPRESSION ############
  
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
  
  difex_genes_deseq <- reactive({
    req(normalized_data())  # Ensure normalized data is available
    data <- filtered_data()  # Get the filtered data
    data <- run_deseq(data)  # Run DESeq to get results
    return(data)
  })
  
  difex_genes_edger <- reactive({
    req(normalized_data())  # Ensure normalized data is available
    data <- filtered_data()  # Get the filtered data
    data <- run_edger(data)  # Run DESeq to get results
    return(data)
  })
  
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
  output$table_summary <- DT::renderDT({
    req(input$upload)
    data <- summarized_table(load_data())
    
    # Format the table using DT for better aesthetics
    datatable(
      data,
      options = list(
        pageLength = 5,  # Show 5 rows per page
        autoWidth = TRUE,  # Automatically adjust column widths
        dom = "t",  # Disable search, pagination, and info elements
        columnDefs = list(list(className = 'dt-center', targets = "_all"))  # Center-align text
      ),
      rownames = FALSE,  # Hide row names
      caption = "Summary of Sample Data"  # Add a caption to the table
    )
  })
  
  # Render unsorted table initially
  output$table_unsorted <- DT::renderDT({
    req(input$upload)
    data <- load_data()
    data <- data[, !colnames(data) %in% "Genes"]  # Exclude the Genes column
    data  # Display the loaded data as an interactive datatable
  }, options = list(pageLength = 10))  # Show 10 rows per page
  
  # Update the table to display non-zero rows when "Remove 0 Counts" button is clicked
  observeEvent(input$remove_button, {
    output$table_unsorted <- DT::renderDT({
      req(input$upload)
      data <- filter_zeros(load_data())# Filter the table to remove rows with zeros
      data <- data[, !colnames(data) %in% "Genes"]  # Exclude the Genes column
      data
    }, options = list(pageLength = 10))  # Show 10 rows per page
    
    # Display a message indicating the table has been filtered
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
  
  output$normalization_summary <- DT::renderDT({
    req(normalized_summaries())
    datatable(normalized_summaries(), options = list(pageLength = 10))
  })
  
  output$normalized_table <- DT::renderDT({
    req(normalized_data())
    data <- normalized_data() 
    data <- data %>% select(-Genes)
    datatable(
      data,
      options = list(pageLength = 10, scrollX = TRUE)  # Add horizontal scrolling for small width
    )
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
  
  output$difex_genes_deseq <- DT::renderDT({
    data <- difex_genes_deseq()  # Get the DESeq results
    
    # Convert the slider input to a p-value threshold
    pval_threshold <- 10^(-input$pval)  # Convert slider value to p-value
    
    # Filter the DESeq results based on the p-value threshold
    filtered_data <- data[data$padj < pval_threshold, ]  # Filter for p-value > threshold
    filtered_data <- na.omit(filtered_data)
    
    # Display the filtered data in a datatable
    datatable(filtered_data %>%
                select(padj, log2FoldChange) 
    )
  })
  
  output$difex_genes_edger <- DT::renderDT({
    data <- difex_genes_edger()  # Get the DESeq results
    
    # Convert the slider input to a p-value threshold
    pval_threshold <- 10^(-input$pval)  # Convert slider value to p-value
    
    # Filter the DESeq results based on the p-value threshold
    filtered_data <- data[data$PValue < pval_threshold, ]  # Filter for p-value > threshold
    filtered_data <- na.omit(filtered_data)
    
    # Display the filtered data in a datatable
    datatable(filtered_data %>%
                select(PValue, logFC)
    )
  })
  
  # Output: Volcano Plot
  output$volcanoes <- renderPlot({
    req(combined_volcano_data(), input$pval)  # Ensure data and slider input are available
    theme_plot(combined_volcano_data(), input$pval)  # Pass the slider value
  })
}

# Run the application
shinyApp(ui = ui, server = server)
