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


ui <- fluidPage(
  titlePanel("Soybean cotyledon gene expression across development"),
  p("This experiment assayed gene expression in the cotyledons, or embryonic leaves, for soybean plants at different ages after planting."),
  # create 4 to 5 tab panels
  tabsetPanel(
    tabPanel("Sample Information",
      # create sidebar layout
      sidebarLayout(
        sidebarPanel(
          fileInput("upload", "Load counts data")
        ),
        mainPanel(
          # create 3 tab panels
          tabsetPanel(
            tabPanel("Sample Summary",
              tableOutput("table_summary")
            ),
            tabPanel("Data Table"),
            tabPanel("Plots")
          )
        )
      )
    ),
    tabPanel("Count Matrix",
      # create sidebar for input filters and output display
      sidebarLayout(
        sidebarPanel(
          p("Filters"),
          p("Slider to include genes with at least X percentile of variance"),
          p("Slider to include genes with at least X samples that are non-zero")
        ),
        mainPanel(
          # create tabset for 3 outputs
          tabsetPanel(
            tabPanel("Summary"),
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
  
  load_data <- reactive({
    req(input$upload)  # Require file input
    data <- read.csv(input$upload$datapath)
    return(data)
  })
  
  
  #' create summary table from counts dataframe 
  #' @param dataf The loaded data frame
  summarized_table <- function(dataf) {
    
    # Extract the column names (e.g., S1.1, S1.2, ...)
    columns <- colnames(dataf)
    
    # Separate the stages and replicates
    stages <- sapply(strsplit(columns, "\\."), `[`, 1)      # Extracts S1, S2, S3
    replicates <- sapply(strsplit(columns, "\\."), `[`, 2)  # Extracts 1, 2, 3
    
    # Create the coldata for the summary tibble
    coldata <- data.frame(
      stage = stages,
      replicate = replicates,
      stringsAsFactors = FALSE
      
    )
    
    # Ensure coldata doesn't have any unwanted rows (e.g., X, NA)
    coldata <- coldata[complete.cases(coldata), , drop = FALSE]
    
    # Create the summary tibble
    summary_tibble <- tibble(
      `Column name` = colnames(coldata),
      Type = sapply(coldata, class),
      `Mean(sd) or Distinct values` = sapply(coldata, function(col) paste(unique(col), collapse = ", "))
    )
    
    # Add a row for the total number of genes
    number_of_genes_row <- tibble(
      `Column name` = "Number of genes",
      Type = "numeric",
      `Mean(sd) or Distinct values` = as.character(nrow(dataf))
    )
    
    # Append the new row to the summary_tibble
    summary_tibble <- bind_rows(summary_tibble, number_of_genes_row)
    
    # Return the summary tibble
    return(summary_tibble)
  }
  
  output$table_summary <- renderTable({
    req(input$upload)
    summarized_table(load_data())
  })
}

# Run the application
shinyApp(ui = ui, server = server)
