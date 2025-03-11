## Author: Aron Diamond
## aronad@bu.edu
## Developed at BU - BF591

source("main.R")

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = 'flatly'), 
  titlePanel("Soybean Gene Expression after Iron Exposure"),
  p("This experiment assayed gene expression from soybean leaves that were exposed to iron-rich (iron-positive) soil conditions versus leaves that were exposed to iron-poor (iron-negative) soil conditions. The data was collected 120 minutes after iron conditions were initiated."),
  tabsetPanel(
    tabPanel("Sample Information",
             sidebarLayout(
               sidebarPanel(
                 fileInput("upload", "Load counts data"),
                 actionButton("remove_button", "Remove 0 counts"),
                 textOutput("remove_message")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Sample Summary",
                            DTOutput("table_summary"),
                            uiOutput("download_summary_ui")
                   ),
                   tabPanel("Data Table",
                            DTOutput("table_unsorted")
                   ),
                   tabPanel("Plots",
                            sliderInput("filter_raw", "Raw count filter", min = 10, max = 16000, value = c(10, 16000), step = 10),
                            helpText("Adjust the slider to filter raw counts."),
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
             sidebarLayout(
               sidebarPanel(
                 h3("Normalization"),
                 radioButtons("normalize_method", 
                              "Choose your normalization method", 
                              c("CPM", "DESeq", "Limma")
                 )
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            fluidRow(
                              column(6, numericInput("count_filter", "Choose minimum count filter", 10)),
                              column(6, textOutput("normalization_summary_text"))  # Text output
                            ),
                            DTOutput("normalized_table")  # Table output for normalized data
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
                              column(6, plotOutput("PCA_plot")),     # Left plot
                              column(6, plotOutput("variance_explained"))  # Right plot
                            )
                   )
                 )
               )
             )
    ),
    tabPanel("Differential Expression Analysis",
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 radioButtons("difex_method",
                              "Differential Expression Method",
                              c("DESeq", "Limma", "EdgeR")
                 ),
                 actionButton("run_difex", 
                              "Run Differential Expression"
                 )
               ),
               mainPanel(
                 width = 9,
                 tabsetPanel(
                   tabPanel("Table",
                            DT::DTOutput("difex_genes") %>% shinycssloaders::withSpinner()
                   ),
                   tabPanel("Plot",
                            fluidRow(
                              column(12, sliderInput("pval", 
                                                     "Set P-value Threshold (-log10 scale):", 
                                                     min = 1, max = 50, value = 25, step = 1
                              ))
                            ),
                            fluidRow(
                              column(12, plotOutput("volcanoes") %>% shinycssloaders::withSpinner())
                            ),
                            fluidRow(
                              column(6,
                                     p("DESeq Results"),
                                     DT::DTOutput("difex_genes_deseq") %>% shinycssloaders::withSpinner()),
                              column(6, 
                                     p("EdgeR Results"),
                                     DT::DTOutput("difex_genes_edger") %>% shinycssloaders::withSpinner())
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
  
  # Load data from uploaded file
  load_data <- reactive({
    req(input$upload)
    data <- read.csv(input$upload$datapath, row.names = 1)
    data$Genes <- rownames(data)
    colnames(data)[colnames(data) == 'X'] <- 'Genes'
    return(data)
  })
  
  # Dynamic download button
  output$download_summary_ui <- renderUI({
    req(input$upload)  # Ensure data is uploaded
    downloadButton("download_summary", "Download Summary Table")
  })
  
  # Download handler for summary table
  output$download_summary <- downloadHandler(
    filename = function() {
      paste("summary_table_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(summarized_table(load_data()), file, row.names = FALSE)
    }
  )
  
  # Filter data based on raw count filter
  filtered_long_data <- reactive({
    req(input$upload)
    data <- load_data()
    min_count <- input$filter_raw[1]
    max_count <- input$filter_raw[2]
    filter_pivot_data(data, min_count, max_count)
  })
  
  # Remove rows with zero counts
  observeEvent(input$remove_button, {
    data <- load_data()
    filtered_data <- filter_zeros(data)
    removed_rows <- nrow(data) - nrow(filtered_data)
    output$remove_message <- renderText({
      paste("Removed", removed_rows, "rows with zero counts.")
    })
    output$table_unsorted <- renderDT({
      datatable(filtered_data[, -ncol(filtered_data)],  # Exclude the Genes column
                rownames = TRUE,  # Show gene names as row names
                options = list(pageLength = 10))
    })
  })
  
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
  
  # Render raw count plots
  output$histogram_raw <- renderPlot({
    req(filtered_long_data())
    ggplot(filtered_long_data(), aes(x = counts)) +
      geom_histogram(binwidth = 50) +
      labs(title = "Histogram of Raw Counts", x = "Counts", y = "Frequency")
  })
  
  output$violin_raw <- renderPlot({
    req(filtered_long_data())
    ggplot(filtered_long_data(), aes(x = replicates, y = counts)) +
      geom_violin() +
      labs(title = "Violin Plot of Raw Counts", x = "Replicates", y = "Counts")
  })
  
  output$ridged <- renderPlot({
    req(filtered_long_data())
    ggplot(filtered_long_data(), aes(x = counts, y = replicates, fill = replicates)) +
      geom_density_ridges(alpha = 0.7) +
      labs(title = "Ridgeline Plot of Raw Counts", x = "Counts", y = "Replicates")
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
  # Example of normalized_data() structure
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
  
  # Create normalization summary text
  output$normalization_summary_text <- renderText({
    req(input$count_filter)
    data <- load_data()
    filtered <- filtered_data()
    
    total_genes <- nrow(data)
    passing_genes <- nrow(filtered)
    passing_percentage <- round((passing_genes / total_genes) * 100, 2)
    
    paste("Based on the minimum count filter of", input$count_filter, 
          ",", passing_percentage, "% of genes passed (", passing_genes, 
          "out of", total_genes, "genes).")
  })
  
  # Render normalization table
  output$normalized_table <- DT::renderDT({
    req(normalized_data())
    data <- normalized_data()
    
    # Display the normalized data as an interactive datatable
    datatable(
      data,
      options = list(pageLength = 10, scrollX = TRUE)  # Add horizontal scrolling for small width
    )
  })
  
  # plot boxplot
  output$boxplot_distribution <- renderPlot({
    req(normalized_data(), input$normalize_method)
    plot_boxplot(normalized_data(), input$normalize_method)
  })
  
  # plot scatterplot 
  output$variance_vs_mean <- renderPlot({
    req(normalized_data(), input$normalize_method)
    plot_variance_vs_mean(normalized_data(), input$normalize_method)
  })
  
  # plot heatmap
  output$heatmap <- renderPlot({
    req(normalized_data())
    generate_heatmap(normalized_data())
  }, , height = 900, width = 600)
  
  # Render PCA plots
  output$PCA_plot <- renderPlot({
    req(normalized_data())   
    plot_pca_scatter(normalized_data())
  })
  
  output$variance_explained <- renderPlot({
    req(normalized_data())
    plot_pca_variance(normalized_data())
  })
  
  ################### DIFFERENTIAL EXPRESSION ############
  
  # Reactive values for caching differential expression results
  cache <- reactiveValues(
    deseq_results = NULL,
    edger_results = NULL,
    limma_results = NULL
  )
  
  # Reactive for differential expression
  difex_results <- eventReactive(input$run_difex, {
    req(normalized_data())
    data <- filtered_data()
    
    tryCatch({
      method <- input$difex_method
      if (is.null(cache[[paste0(method, "_results")]])) {
        cache[[paste0(method, "_results")]] <- run_difex(method, data)
      }
      return(cache[[paste0(method, "_results")]])
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Reactive for DESeq results
  difex_genes_deseq <- reactive({
    req(normalized_data())  # Ensure normalized data is available
    data <- filtered_data()  # Get the filtered data
    data <- run_deseq(data)  # Run DESeq to get results
    return(data)
  })
  
  # Reactive for EdgeR results
  difex_genes_edger <- reactive({
    req(normalized_data())  # Ensure normalized data is available
    data <- filtered_data()  # Get the filtered data
    data <- run_edger(data)  # Run EdgeR to get results
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
  
  # Provide user feedback
  observeEvent(input$run_difex, {
    showNotification("Differential expression analysis complete!", type = "message")
  })
  
  # Output: Differential Expression Table
  output$difex_genes <- DT::renderDT({
    difex_results()
  }) 
  
  # Output: DESeq Results Table
  output$difex_genes_deseq <- DT::renderDT({
    data <- difex_genes_deseq()
    pval_threshold <- 10^(-input$pval)
    filtered_data <- data[data$padj < pval_threshold, ]
    filtered_data <- na.omit(filtered_data)
    datatable(filtered_data %>% select(padj, log2FoldChange))
  })
  
  # Output: EdgeR Results Table
  output$difex_genes_edger <- DT::renderDT({
    data <- difex_genes_edger()
    pval_threshold <- 10^(-input$pval)
    filtered_data <- data[data$PValue < pval_threshold, ]
    filtered_data <- na.omit(filtered_data)
    datatable(filtered_data %>% select(PValue, logFC))
  })
  
  # Output: Volcano Plot
  output$volcanoes <- renderPlot({
    req(combined_volcano_data(), input$pval)
    theme_plot(combined_volcano_data(), input$pval)
  })

}

# Run the application
shinyApp(ui = ui, server = server)
