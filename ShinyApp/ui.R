library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)

# Define UI for application
shinyUI(navbarPage("Single-Cell Analysis Dashboard",

    # --- QC Tab ---
    tabPanel("Quality Control",
        sidebarLayout(
            sidebarPanel(
                h3("QC Parameters"),
                # Input: Specify data path --- Consider making this dynamic later
                textInput("data_path", "Data Path:", value = "Data"), # Relative to ShinyApp dir (assuming Data is moved inside)

                # Input: QC Thresholds based on README
                numericInput("max_rna_count", "Max nCount_RNA:", value = 75000, min = 1),
                numericInput("max_features", "Max nFeature_RNA:", value = 5000, min = 1),
                numericInput("max_mt_percent", "Max Mitochondrial %:", value = 15, min = 0, max = 100),
                numericInput("max_rb_percent", "Max Ribosomal %:", value = 50, min = 0, max = 100),

                actionButton("run_qc", "Load Data & Run QC")
            ),
            mainPanel(
                h3("QC Metrics Visualization"),
                h4("Before Filtering"),
                plotOutput("qcVlnPlotPre"),
                plotOutput("qcScatterPlotsPre"),
                hr(),
                h4("After Filtering"),
                plotOutput("qcVlnPlotPost"),
                verbatimTextOutput("qcSummary") # To show cell counts before/after
            )
        )
    ),

    # --- Normalization & Dim Reduc Tab ---
    tabPanel("Normalization & Dimensionality Reduction",
        sidebarLayout(
            sidebarPanel(
                h3("Parameters"),
                # Normalization is standard LogNormalize
                numericInput("n_variable_features", "Number of Variable Features:", value = 3000, min = 100),
                # PCA/Clustering Params
                numericInput("pca_dims", "Number of PCA Dimensions for Clustering:", value = 10, min = 1),
                numericInput("cluster_resolution", "Clustering Resolution:", value = 0.5, min = 0.1, step = 0.1),

                actionButton("run_norm_dimreduc", "Run Normalization, PCA & Clustering")
            ),
            mainPanel(
                h3("Outputs"),
                plotOutput("variableFeaturePlot"),
                plotOutput("pcaElbowPlot"),
                plotOutput("umapPlot"),
                plotOutput("tsnePlot"),
                verbatimTextOutput("clusterInfo") # Show cluster sizes
            )
        )
    ),

    # --- Marker Genes Tab ---
    tabPanel("Marker Genes",
        sidebarLayout(
            sidebarPanel(
                h3("Parameters"),
                numericInput("marker_cluster_id", "Find Markers for Cluster:", value = 1, min = 0), # Assuming clusters are 0-indexed or need adjustment
                numericInput("marker_logfc", "LogFC Threshold:", value = 0.5, min = 0, step = 0.1),
                numericInput("marker_min_pct", "Min Pct Expression:", value = 0.25, min = 0, max = 1, step = 0.05),
                actionButton("find_markers", "Find Cluster Markers"),
                hr(),
                textInput("feature_plot_gene", "Gene(s) for FeaturePlot (comma-separated):", value = "MS4A1,CD79A,CD8A"),
                actionButton("run_featureplot", "Generate FeaturePlot")

            ),
            mainPanel(
                h3("Marker Gene Outputs"),
                h4("Top Markers for Selected Cluster"),
                tableOutput("markerTable"),
                plotOutput("markerVlnPlot"),
                hr(),
                h4("Feature Plots"),
                plotOutput("featurePlot")
            )
        )
    )

    ), # Close Marker Genes tabPanel

    # --- Cell Type Annotation Tab ---
     tabPanel("Cell Type Annotation",
         sidebarLayout(
             sidebarPanel(
                 h3("Annotation Parameters"),
                 # For simplicity, start with pre-defined references from README
                 checkboxGroupInput("references_to_use", "Select Reference Datasets:",
                                    choices = list("Human Primary Cell Atlas" = "primary",
                                                   "Blueprint/ENCODE" = "blueprint"),
                                    selected = c("primary", "blueprint")), # Default to both for consensus
                 actionButton("run_annotation", "Run Cell Type Annotation (SingleR)")
             ),
             mainPanel(
                 h3("Annotation Results"),
                 p("Note: Annotation can take several minutes, especially the first time references are downloaded."),
                 h4("Consensus Annotations (UMAP)"),
                 plotOutput("umapAnnotationConsensus"),
                 h4("Human Primary Cell Atlas Annotations (UMAP)"),
                 plotOutput("umapAnnotationPrimary"),
                 h4("Blueprint/ENCODE Annotations (UMAP)"),
                 plotOutput("umapAnnotationBlueprint"),
                 h4("Annotation Scores (Primary Atlas)"),
                 plotOutput("scoreHeatmapPrimary"),
                 h4("Annotation Scores (Blueprint/ENCODE)"),
                 plotOutput("scoreHeatmapBlueprint"),
                 h4("Cell Type Counts (Consensus)"),
                 verbatimTextOutput("annotationSummary")
             )
         )
     ), # Close Subcluster Analysis tabPanel

    # --- Differential Expression Tab ---
     tabPanel("Differential Expression",
         sidebarLayout(
             sidebarPanel(
                 h3("Differential Expression Parameters"),
                 # Allow user to select cell types for comparison
                 selectInput("diffex_cell_types", "Select Cell Types to Compare:",
                             choices = NULL, # Will be updated dynamically
                             multiple = TRUE),
                 # Add options for FindAllMarkers parameters if needed
                 numericInput("diffex_logfc", "LogFC Threshold:", value = 0.5, min = 0, step = 0.1),
                 numericInput("diffex_min_pct", "Min Pct Expression:", value = 0.5, min = 0, max = 1, step = 0.05), # README used 0.5 here
                 # Option to use SCTransform (as per README step 9) - might be slow
                 # checkboxInput("use_sctransform_diffex", "Use SCTransform for Normalization (Slower)", value = FALSE),
                 actionButton("run_diffex", "Run Differential Expression")
             ),
             mainPanel(
                 h3("Differential Expression Results"),
                 p("Compares selected cell types to find differentially expressed genes between them."),
                 h4("Top Differentially Expressed Genes (All Groups)"),
                 tableOutput("diffexMarkerTable"),
                 h4("Top Marker Violin Plots"),
                 plotOutput("diffexVlnPlot")
             )
         )
     )
    ), # Close Annotation tabPanel

    # --- Subcluster Analysis Tab ---
     tabPanel("Subcluster Analysis",
         sidebarLayout(
             sidebarPanel(
                 h3("Subsetting Parameters"),
                 # Allow user to select cell types based on consensus annotation
                 selectInput("subset_cell_types", "Select Cell Types to Subset:",
                             choices = NULL, # Will be updated dynamically
                             multiple = TRUE),
                 # Add options for re-analysis parameters if needed
                 # numericInput("subset_pca_dims", "PCA Dims for Subset:", value = 10, min = 1),
                 # numericInput("subset_resolution", "Resolution for Subset:", value = 0.5, min = 0.1, step = 0.1),
                 actionButton("run_subset_analysis", "Run Subcluster Analysis")
             ),
             mainPanel(
                 h3("Subcluster Results"),
                 p("This performs re-clustering and marker finding on the selected subset of cells."),
                 h4("Subset UMAP"),
                 plotOutput("umapSubset"),
                 h4("Top Markers for Subset Clusters"),
                 tableOutput("subsetMarkerTable"),
                 h4("Top Marker Feature Plots (Subset)"),
                 plotOutput("subsetFeaturePlot")
                 # Potentially add fine-grained annotation outputs here later
             )
         )
     )
    # Add more tabs for further analysis steps if needed...
))