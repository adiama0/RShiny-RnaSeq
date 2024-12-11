## BF591 Rshiny Final Project

### Description
The application combines Shiny for interactive data visualization with advanced statistical analysis tools for RNA-seq data. It allows users to upload gene expression data, perform PCA (Principal Component Analysis), and differential expression analysis using DESeq2, EdgeR, and Limma. The app includes features for generating various plots such as PCA scatter plots, beeswarm plots for principal components, and volcano plots for differential expression results. The app also enables users to interactively select which differential expression method to apply and visualize the results of multiple methods in a single interface.

The app has three main components:
1. **Data Upload and Preprocessing**: Users upload a gene expression dataset in CSV format, which is processed to extract relevant data and metadata.
2. **PCA Plotting**: The app performs PCA on the uploaded data and visualizes the first two principal components. The user can choose to generate a beeswarm plot to view the distribution of the top principal components.
3. **Differential Expression Analysis**: Users can select from three different differential expression methods (DESeq2, EdgeR, Limma). The app displays the results in the form of volcano plots, showing log fold changes and p-values, with significance highlighted based on a user-defined threshold.

### Features
- **PCA Visualization**: Displays a scatter plot of the first two principal components and a beeswarm plot for the selected number of principal components.
- **Differential Expression Analysis**: Implements DESeq2, EdgeR, and Limma for RNA-seq data analysis. Users can compare the results of different methods side by side.
- **Volcano Plot**: Visualizes the results of differential expression analyses with customizable significance thresholds.
- **Interactivity**: The user can interact with the app by uploading files, selecting analysis methods, and adjusting visualization parameters (such as the p-value threshold for volcano plots).

### Components of the Application

1. **`app.R` (Shiny Application UI and Server)**: 
   - This file handles the user interface (UI) and the server logic of the Shiny application. The UI provides components for file upload, data display, plot rendering, and selection of analysis methods. The server-side code manages data processing, computation, and plot generation based on user inputs.
   - Main UI components include:
     - File upload for RNA-seq data (CSV format).
     - Input controls for selecting the number of principal components to visualize.
     - Dropdowns to select differential expression methods.
     - Interactive plots for PCA, beeswarm, and volcano plots.

2. **`main.R` (Differential Expression Functions)**:
   - This file contains functions for running differential expression analysis using DESeq2, EdgeR, and Limma. It also includes functions for performing PCA and generating plots.
   - Main functions in `main.R`:
     - `plot_pca_scatter()`: Generates a PCA scatter plot.
     - `plot_beeswarm()`: Creates a beeswarm plot for principal components.
     - `run_deseq()`, `run_edger()`, `run_limma()`: Perform differential expression analysis using DESeq2, EdgeR, and Limma, respectively.
     - `create_facets()`: Combines the results of differential expression methods into a long format for easier comparison.
     - `theme_plot()`: Customizes volcano plots with coloring based on significance.

3. **`main2.R` (Helper Functions and Plotting)**:
   - This file contains additional helper functions for formatting and visualizing the data, particularly for plotting volcano plots and customizing the appearance of plots.

### Usage Example

1. **Start the Shiny App**: 
   - Load the `app.R` file and run the Shiny application. This opens the interactive app in your browser.
   
   ```r
   shiny::runApp('path_to_your_app_folder')
   ```

2. **Upload Data**: 
   - Once the app is open, upload your RNA-seq data in CSV format. The file should have genes as rows and samples as columns.

3. **Select PCA Parameters**:
   - Choose the number of principal components (PCs) to visualize. The app will generate the corresponding PCA scatter plot.
   
4. **Run Differential Expression Analysis**:
   - Choose one of the three differential expression methods (DESeq2, EdgeR, Limma) and click to run the analysis.
   - The app will generate the results, including log fold changes, p-values, and volcano plots.
   
5. **View Results**:
   - The app will display the volcano plot and allow you to adjust the threshold for significance. You can also compare the results from different differential expression methods side by side.

### Example Workflow in the App

1. **Upload Gene Expression Data**: Upload a CSV file containing gene expression counts. The app will process this data.
2. **Visualize PCA**: Select the number of principal components to visualize and view the PCA scatter plot and beeswarm plot.
3. **Perform Differential Expression Analysis**: Select a differential expression method (DESeq2, EdgeR, or Limma) and run the analysis.
4. **View Volcano Plot**: The app will generate a volcano plot showing the results of differential expression analysis, allowing you to customize the significance threshold and explore the results.

### Requirements
- R (>= 4.0.0)
- Shiny
- DESeq2
- EdgeR
- Limma
- ggplot2
- dplyr
- tidyr
- tibble

### Example Code to Run the App

```r
library(shiny)
runApp('path_to_your_app_folder')
```

This setup allows for an interactive, user-friendly experience in RNA-seq data analysis, enabling users to perform PCA, differential expression analysis, and visualize the results directly within a web interface.
