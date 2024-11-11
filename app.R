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


ui <- fluidPage(
  titlePanel("Soybean cotyledon gene expression across development"),
  p("This experiment assayed gene expression in the cotyledons, or embryonic leaves, for soybean plants at different ages after planting."),
  # create 4 to 5 tab panels
  tabsetPanel(
    tabPanel("Sample Information",
      # create sidebar layout
      sidebarLayout(
        sidebarPanel(
          p("input data set")
        ),
        mainPanel(
          # create 3 tab panels
          tabsetPanel(
            tabPanel("Sample Summary"),
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
  
}

# Run the application
shinyApp(ui = ui, server = server)
