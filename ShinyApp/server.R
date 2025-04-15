# Load necessary libraries
library(shiny)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(magrittr) # For pipe operator if used beyond Seurat/dplyr
library(RColorBrewer) # Potentially for plots later
library(ggrepel) # For non-overlapping labels
library(SingleR) # For cell type annotation
library(celldex) # For reference datasets
library(SingleCellExperiment) # Required for SingleR input

# Define server logic
shinyServer(function(input, output, session) {

    # Reactive value to store the Seurat object and intermediate results
    vals <- reactiveValues(
        seurat_obj = NULL,
        seurat_obj_qc = NULL, # Store object after QC
        seurat_obj_processed = NULL, # Store object after clustering/dim reduc
        markers = NULL,
        qc_summary_text = "Data not loaded yet.",
        # Annotation results
        singler_primary_results = NULL,
        singler_blueprint_results = NULL,
        annotation_summary_text = "Annotation not run yet.",
        # Subset analysis results
        seurat_obj_subset = NULL,
        subset_markers = NULL,
        # Differential Expression results
        diffex_markers = NULL
    )

    # --- 1. Quality Control ---
    observeEvent(input$run_qc, {
        req(input$data_path) # Ensure data path is provided

        # Show progress to the user
        withProgress(message = 'Loading Data & Running QC...', value = 0, {

            incProgress(0.1, detail = "Reading 10X data...")
            # Step 1: Load the Dataset
            data <- tryCatch({
                Read10X(data.dir = input$data_path)
            }, error = function(e) {
                showNotification(paste("Error reading 10X data:", e$message), type = "error")
                return(NULL)
            })

            if (is.null(data)) return()

            # Replace underscores with dashes if necessary (as per README)
            rownames(data) <- gsub("_", "-", rownames(data))

            incProgress(0.2, detail = "Creating Seurat object...")
            # Step 2: Initialize Seurat Object
            df <- tryCatch({
                 CreateSeuratObject(counts = data, project = "sc_dashboard", min.cells = 3, min.features = 200)
            }, error = function(e) {
                showNotification(paste("Error creating Seurat object:", e$message), type = "error")
                return(NULL)
            })

            if (is.null(df)) return()
            rm(data) # Free memory

            n_cells_pre_qc <- ncol(df) # Store initial cell count

            incProgress(0.3, detail = "Calculating QC metrics...")
            # Step 3: Calculate QC Metrics
            # Use tryCatch for PercentageFeatureSet as patterns might not match all datasets
             tryCatch({
                 df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
             }, error = function(e) {
                 showNotification("Could not calculate MT%, setting to 0. Check pattern '^MT-' matches gene names.", type = "warning")
                 df[["percent.mt"]] <- 0
             })
             tryCatch({
                 df[["percent.rb"]] <- PercentageFeatureSet(df, pattern = "^RP[SL]")
             }, error = function(e) {
                 showNotification("Could not calculate RB%, setting to 0. Check pattern '^RP[SL]' matches gene names.", type = "warning")
                 df[["percent.rb"]] <- 0
             })


            vals$seurat_obj <- df # Store pre-filtered object temporarily for plotting

            incProgress(0.5, detail = "Generating pre-filter plots...")
            # Step 4: Visualize Initial QC Metrics
            output$qcVlnPlotPre <- renderPlot({
                req(vals$seurat_obj)
                VlnPlot(vals$seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
                        ncol = 4, pt.size = 0.1) + ggtitle("QC Metrics Before Filtering")
            })

            output$qcScatterPlotsPre <- renderPlot({
                 req(vals$seurat_obj)
                 plot1 <- FeatureScatter(vals$seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("RNA Count vs Mito %")
                 plot2 <- FeatureScatter(vals$seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("RNA Count vs Feature Count")
                 plot3 <- FeatureScatter(vals$seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.rb") + ggtitle("RNA Count vs Ribo %")
                 plot1 + plot2 + plot3
            })

            incProgress(0.7, detail = "Filtering cells...")
            # Step 5: Filter Based on QC Metrics
            df_filtered <- tryCatch({
                subset(df, subset = nCount_RNA < input$max_rna_count &
                                  nFeature_RNA < input$max_features &
                                  percent.mt < input$max_mt_percent &
                                  percent.rb < input$max_rb_percent)
            }, error = function(e) {
                showNotification(paste("Error during subsetting:", e$message), type = "error")
                return(NULL)
            })

            if (is.null(df_filtered)) return()

            n_cells_post_qc <- ncol(df_filtered)
            vals$seurat_obj_qc <- df_filtered # Store filtered object

            incProgress(0.9, detail = "Generating post-filter plots...")
            # Step 6: Visualize QC Metrics After Filtering
            output$qcVlnPlotPost <- renderPlot({
                req(vals$seurat_obj_qc)
                VlnPlot(vals$seurat_obj_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
                        ncol = 4, pt.size = 0.1) + ggtitle("QC Metrics After Filtering")
            })

            # Update summary text
            vals$qc_summary_text <- paste(
                "Cells before QC:", n_cells_pre_qc, "\n",
                "Cells after QC:", n_cells_post_qc
            )
            showNotification("QC step completed.", type = "message")
        }) # End withProgress
    })

    # Output for QC Summary
    output$qcSummary <- renderText({
        vals$qc_summary_text
    })


    # --- 2. Normalization & Dimensionality Reduction ---
    observeEvent(input$run_norm_dimreduc, {
        req(vals$seurat_obj_qc) # Need the QC'd object

        df_processed <- vals$seurat_obj_qc # Start with the QC'd object

        withProgress(message = 'Running Analysis...', value = 0, {

            incProgress(0.1, detail = "Normalizing data...")
            # Step 1: Normalize Data
            df_processed <- NormalizeData(df_processed, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

            incProgress(0.2, detail = "Finding variable features...")
            # Step 2: Find Variable Features
            df_processed <- FindVariableFeatures(df_processed, selection.method = "vst", nfeatures = input$n_variable_features, verbose = FALSE)
            top10 <- head(VariableFeatures(df_processed), 10)

            incProgress(0.3, detail = "Generating variable feature plot...")
            # Step 3: Visualize Highly Variable Features
            output$variableFeaturePlot <- renderPlot({
                plot1 <- VariableFeaturePlot(df_processed) + scale_x_continuous(trans = "log1p") + theme_minimal()
                LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0) +
                    ggtitle("Top Variable Features")
            })

            incProgress(0.4, detail = "Scaling data...")
            # Step 4: Scale Data (using variable features)
            all.genes <- rownames(df_processed) # Scale all genes for PCA visualization later if needed, but PCA uses HVGs
            df_processed <- ScaleData(df_processed, features = all.genes, verbose = FALSE) # Scaling all genes can be slow

            incProgress(0.5, detail = "Running PCA...")
            # Step 5: Perform PCA (on variable features)
            df_processed <- RunPCA(df_processed, features = VariableFeatures(object = df_processed), verbose = FALSE)

            incProgress(0.6, detail = "Generating Elbow plot...")
            # Step 6: Determine Optimal Number of PCs (Elbow Plot)
            output$pcaElbowPlot <- renderPlot({
                ElbowPlot(df_processed) + ggtitle("PCA Elbow Plot")
            })

            incProgress(0.7, detail = "Finding neighbors and clusters...")
            # Step 7: Find Neighbors and Identify Clusters
            dims_to_use <- 1:input$pca_dims
            df_processed <- FindNeighbors(df_processed, dims = dims_to_use, verbose = FALSE)
            df_processed <- FindClusters(df_processed, resolution = input$cluster_resolution, verbose = FALSE)

            # Display cluster info
            output$clusterInfo <- renderPrint({
                req(df_processed)
                print("Cells per cluster:")
                table(Idents(df_processed))
            })

            incProgress(0.8, detail = "Running UMAP and t-SNE...")
            # Step 8: Run Dimensional Reduction (UMAP and t-SNE)
            # Wrap in tryCatch as they can sometimes fail with certain parameters/data
             tryCatch({
                 df_processed <- RunUMAP(df_processed, dims = dims_to_use, verbose = FALSE)
             }, error = function(e) {
                 showNotification(paste("Error running UMAP:", e$message), type = "error")
             })
             tryCatch({
                 df_processed <- RunTSNE(df_processed, dims = dims_to_use, verbose = FALSE) # Note: t-SNE can be slow
             }, error = function(e) {
                 showNotification(paste("Error running t-SNE:", e$message), type = "error")
             })


            incProgress(0.9, detail = "Generating DimPlots...")
            # Step 9: Visualize Clusters
            output$umapPlot <- renderPlot({
                req(df_processed[["umap"]]) # Check if UMAP ran successfully
                DimPlot(df_processed, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("UMAP of Clusters")
            })
            output$tsnePlot <- renderPlot({
                 req(df_processed[["tsne"]]) # Check if tSNE ran successfully
                DimPlot(df_processed, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("t-SNE of Clusters")
            })

            vals$seurat_obj_processed <- df_processed # Store the fully processed object
            showNotification("Normalization, Dimensionality Reduction, and Clustering complete.", type = "message")

        }) # End withProgress
    })


    # --- 3. Marker Genes ---
    observeEvent(input$find_markers, {
        req(vals$seurat_obj_processed) # Need processed object
        req(input$marker_cluster_id)

        df_markers_obj <- vals$seurat_obj_processed

        # Ensure cluster ID is valid
        cluster_id_input <- input$marker_cluster_id
        available_clusters <- levels(Idents(df_markers_obj))

        if (!as.character(cluster_id_input) %in% available_clusters) {
             showNotification(paste("Cluster ID", cluster_id_input, "not found. Available clusters:", paste(available_clusters, collapse=", ")), type = "warning")
             return()
        }

        withProgress(message = paste('Finding markers for cluster', cluster_id_input), value = 0, {

            incProgress(0.1, detail = "Running FindMarkers...")
            # Step 1: Find Markers for Selected Cluster
            markers <- tryCatch({
                FindMarkers(df_markers_obj, ident.1 = cluster_id_input,
                            min.pct = input$marker_min_pct,
                            logfc.threshold = input$marker_logfc,
                            verbose = FALSE)
            }, error = function(e) {
                showNotification(paste("Error finding markers:", e$message), type = "error")
                return(NULL)
            })

            if (is.null(markers) || nrow(markers) == 0) {
                 showNotification(paste("No markers found for cluster", cluster_id_input, "with current settings."), type = "warning")
                 vals$markers <- NULL # Clear previous markers
                 # Clear table and plot outputs
                 output$markerTable <- renderTable(NULL)
                 output$markerVlnPlot <- renderPlot(NULL)
                 return()
            }

            # Add gene symbol column for easier reading
            markers$gene <- rownames(markers)
            markers <- markers %>% arrange(desc(avg_log2FC)) # Sort by Log2FC

            vals$markers <- markers # Store markers

            incProgress(0.5, detail = "Generating marker table...")
            # Step 2: Display Top Markers Table
            output$markerTable <- renderTable({
                head(vals$markers, n = 20) # Show top 20 markers
            }, rownames = FALSE)

            incProgress(0.8, detail = "Generating Violin plot...")
            # Step 3: Visualize Top Markers (Violin Plot)
            output$markerVlnPlot <- renderPlot({
                req(vals$markers)
                top_markers_for_plot <- head(rownames(vals$markers), 2) # Get top 2 genes
                if (length(top_markers_for_plot) > 0) {
                    VlnPlot(df_markers_obj, features = top_markers_for_plot, idents = cluster_id_input, pt.size = 0.1) + # Show only selected cluster for clarity
                        ggtitle(paste("Top 2 Markers for Cluster", cluster_id_input))
                } else {
                    # Should not happen due to check above, but as fallback:
                    plot.new()
                    title("No significant markers found.")
                }
            })
            showNotification(paste("Marker finding for cluster", cluster_id_input, "complete."), type = "message")
        }) # End withProgress
    })

    # Feature Plot Generation
    observeEvent(input$run_featureplot, {
        req(vals$seurat_obj_processed)
        req(input$feature_plot_gene)

        genes_to_plot <- trimws(unlist(strsplit(input$feature_plot_gene, ",")))
        genes_to_plot <- genes_to_plot[genes_to_plot != ""] # Remove empty strings

        if (length(genes_to_plot) == 0) {
            showNotification("Please enter at least one gene.", type = "warning")
            return()
        }

        # Check if genes exist in the object
        genes_exist <- genes_to_plot %in% rownames(vals$seurat_obj_processed)
        if (!all(genes_exist)) {
             showNotification(paste("Gene(s) not found:", paste(genes_to_plot[!genes_exist], collapse=", ")), type = "warning")
             genes_to_plot <- genes_to_plot[genes_exist] # Plot only existing genes
             if (length(genes_to_plot) == 0) {
                 output$featurePlot <- renderPlot(NULL) # Clear plot if no valid genes
                 return()
             }
        }


        output$featurePlot <- renderPlot({
             req(vals$seurat_obj_processed[["umap"]]) # Ensure UMAP exists
             FeaturePlot(vals$seurat_obj_processed, features = genes_to_plot, reduction = "umap") &
                theme(plot.title = element_text(size=10)) # Adjust title size if needed
                # plot_layout(ncol = ceiling(sqrt(length(genes_to_plot)))) # Arrange plots nicely - FeaturePlot handles this
        })
         showNotification("FeaturePlot generated.", type = "message")
    })

    # --- 4. Cell Type Annotation ---
    observeEvent(input$run_annotation, {
        req(vals$seurat_obj_processed) # Need processed object
        req(input$references_to_use)   # Need at least one reference selected

        seurat_obj <- vals$seurat_obj_processed
        selected_refs <- input$references_to_use

        withProgress(message = 'Running Cell Type Annotation...', value = 0, {

            incProgress(0.1, detail = "Converting to SingleCellExperiment...")
            # Convert Seurat object to SingleCellExperiment (use DietSeurat for efficiency)
            sce <- tryCatch({
                as.SingleCellExperiment(DietSeurat(seurat_obj))
            }, error = function(e) {
                showNotification(paste("Error converting to SCE:", e$message), type = "error")
                return(NULL)
            })
            if (is.null(sce)) return()

            # Initialize results storage
            results_primary <- NULL
            results_blueprint <- NULL
            ref_primary <- NULL
            ref_blueprint <- NULL

            # Run SingleR for selected references
            if ("primary" %in% selected_refs) {
                incProgress(0.2, detail = "Loading Human Primary Cell Atlas...")
                ref_primary <- tryCatch({
                    celldex::HumanPrimaryCellAtlasData()
                }, error = function(e) {
                    showNotification(paste("Error loading Primary Atlas:", e$message), type = "error")
                    return(NULL)
                })
                if (!is.null(ref_primary)) {
                    incProgress(0.3, detail = "Running SingleR (Primary Atlas)...")
                    results_primary <- tryCatch({
                        SingleR(test = sce, assay.type.test = 1, ref = ref_primary, labels = ref_primary$label.main)
                    }, error = function(e) {
                        showNotification(paste("Error running SingleR (Primary):", e$message), type = "error")
                        return(NULL)
                    })
                    vals$singler_primary_results <- results_primary # Store full results for heatmap
                }
            }

            if ("blueprint" %in% selected_refs) {
                incProgress(0.5, detail = "Loading Blueprint/ENCODE reference...")
                 ref_blueprint <- tryCatch({
                    celldex::BlueprintEncodeData()
                }, error = function(e) {
                    showNotification(paste("Error loading Blueprint/ENCODE:", e$message), type = "error")
                    return(NULL)
                })
                 if (!is.null(ref_blueprint)) {
                    incProgress(0.6, detail = "Running SingleR (Blueprint/ENCODE)...")
                    results_blueprint <- tryCatch({
                        SingleR(test = sce, assay.type.test = 1, ref = ref_blueprint, labels = ref_blueprint$label.main)
                    }, error = function(e) {
                        showNotification(paste("Error running SingleR (Blueprint):", e$message), type = "error")
                        return(NULL)
                    })
                    vals$singler_blueprint_results <- results_blueprint # Store full results
                 }
            }

            incProgress(0.8, detail = "Adding annotations to Seurat object...")
            # Add annotations to metadata
            seurat_obj$celltype_primary <- NA # Initialize columns
            seurat_obj$celltype_blueprint <- NA
            seurat_obj$celltype_consensus <- NA

            if (!is.null(results_primary)) {
                seurat_obj$celltype_primary <- results_primary$pruned.labels
            }
            if (!is.null(results_blueprint)) {
                seurat_obj$celltype_blueprint <- results_blueprint$pruned.labels
            }

            # Create consensus annotation (prioritize primary if available)
            if (!is.null(results_primary) && !is.null(results_blueprint)) {
                 seurat_obj$celltype_consensus <- ifelse(!is.na(seurat_obj$celltype_primary),
                                                      seurat_obj$celltype_primary,
                                                      seurat_obj$celltype_blueprint)
            } else if (!is.null(results_primary)) {
                 seurat_obj$celltype_consensus <- seurat_obj$celltype_primary
            } else if (!is.null(results_blueprint)) {
                 seurat_obj$celltype_consensus <- seurat_obj$celltype_blueprint
            }

            vals$seurat_obj_processed <- seurat_obj # Update the stored object

            incProgress(0.9, detail = "Generating annotation plots...")
            # Generate plots
            output$umapAnnotationConsensus <- renderPlot({
                req(vals$seurat_obj_processed$celltype_consensus)
                DimPlot(vals$seurat_obj_processed, reduction = "umap", group.by = "celltype_consensus", label = TRUE, repel = TRUE) +
                    NoLegend() + ggtitle("Consensus Cell Type Annotations")
            })
            output$umapAnnotationPrimary <- renderPlot({
                req(vals$seurat_obj_processed$celltype_primary)
                DimPlot(vals$seurat_obj_processed, reduction = "umap", group.by = "celltype_primary", label = TRUE, repel = TRUE) +
                    NoLegend() + ggtitle("Human Primary Cell Atlas Annotations")
            })
            output$umapAnnotationBlueprint <- renderPlot({
                 req(vals$seurat_obj_processed$celltype_blueprint)
                 DimPlot(vals$seurat_obj_processed, reduction = "umap", group.by = "celltype_blueprint", label = TRUE, repel = TRUE) +
                    NoLegend() + ggtitle("Blueprint/ENCODE Annotations")
            })
            output$scoreHeatmapPrimary <- renderPlot({
                req(vals$singler_primary_results)
                plotScoreHeatmap(vals$singler_primary_results) + ggtitle("Annotation Scores - Primary Atlas")
            })
             output$scoreHeatmapBlueprint <- renderPlot({
                req(vals$singler_blueprint_results)
                plotScoreHeatmap(vals$singler_blueprint_results) + ggtitle("Annotation Scores - Blueprint/ENCODE")
            })

            # Generate summary text
            vals$annotation_summary_text <- capture.output(
                 print("Consensus Cell Type Counts:"),
                 print(table(vals$seurat_obj_processed$celltype_consensus))
            )

            showNotification("Cell type annotation complete.", type = "message")
        }) # End withProgress
    })

    # Output for Annotation Summary
    output$annotationSummary <- renderText({
        vals$annotation_summary_text
    })

    })

    # --- Dynamically update subset choices based on annotation ---
    observe({
        req(vals$seurat_obj_processed)
        # Check if the consensus annotation column exists
        if ("celltype_consensus" %in% colnames(vals$seurat_obj_processed@meta.data)) {
            available_cell_types <- na.omit(unique(vals$seurat_obj_processed$celltype_consensus))
            # Sort alphabetically for better UI
            available_cell_types <- sort(available_cell_types)
            updateSelectInput(session, "subset_cell_types",
                              choices = available_cell_types,
                              selected = NULL) # Start with none selected
        } else {
             # If no annotation, keep choices empty
             updateSelectInput(session, "subset_cell_types", choices = list(), selected = NULL)
        }
    })

    # --- 5. Subcluster Analysis ---
    observeEvent(input$run_subset_analysis, {
        req(vals$seurat_obj_processed)
        req(input$subset_cell_types)
        req("celltype_consensus" %in% colnames(vals$seurat_obj_processed@meta.data)) # Ensure annotation ran

        selected_types <- input$subset_cell_types
        seurat_full <- vals$seurat_obj_processed

        # Clear previous subset results before starting
        vals$seurat_obj_subset <- NULL
        vals$subset_markers <- NULL
        output$umapSubset <- renderPlot(NULL)
        output$subsetMarkerTable <- renderTable(NULL)
        output$subsetFeaturePlot <- renderPlot(NULL)


        withProgress(message = 'Running Subcluster Analysis...', value = 0, {

            incProgress(0.1, detail = "Subsetting data...")
            # Step 1: Subset based on selected cell types
            # Ensure the identity is set correctly before subsetting
            Idents(seurat_full) <- "celltype_consensus"
            df_subset <- tryCatch({
                subset(seurat_full, idents = selected_types)
            }, error = function(e) {
                showNotification(paste("Error subsetting:", e$message), type = "error")
                return(NULL)
            })

            if (is.null(df_subset) || ncol(df_subset) == 0) {
                showNotification("Subsetting resulted in zero cells.", type = "warning")
                return()
            }

            # Re-run analysis on the subset for finer resolution
            # Resetting the default assay might be necessary if SCTransform was used previously
            # DefaultAssay(df_subset) <- "RNA" # Assuming RNA assay for consistency with earlier steps

            incProgress(0.2, detail = "Re-normalizing subset...")
            # Re-normalize (optional but often good practice for subsets)
            df_subset <- NormalizeData(df_subset, verbose = FALSE)

            incProgress(0.3, detail = "Finding variable features for subset...")
            # Re-find variable features
            df_subset <- FindVariableFeatures(df_subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE) # Use a reasonable default or add input

            incProgress(0.4, detail = "Scaling subset data...")
            # Re-scale
            all.genes.subset <- rownames(df_subset)
            df_subset <- ScaleData(df_subset, features = all.genes.subset, verbose = FALSE)

            incProgress(0.5, detail = "Running PCA on subset...")
            # Re-run PCA
            df_subset <- RunPCA(df_subset, features = VariableFeatures(object = df_subset), verbose = FALSE)

            incProgress(0.6, detail = "Clustering subset...")
            # Re-cluster (use default dims or add input)
            dims_subset <- 1:10 # Or use input$subset_pca_dims if added
            res_subset <- 0.5 # Or use input$subset_resolution if added
            df_subset <- FindNeighbors(df_subset, dims = dims_subset, verbose = FALSE)
            df_subset <- FindClusters(df_subset, resolution = res_subset, verbose = FALSE)

            incProgress(0.7, detail = "Running UMAP on subset...")
            # Re-run UMAP
             tryCatch({
                 df_subset <- RunUMAP(df_subset, dims = dims_subset, verbose = FALSE)
             }, error = function(e) {
                 showNotification(paste("Error running UMAP on subset:", e$message), type = "error")
                 # Still store the object even if UMAP fails, clustering might be useful
                 vals$seurat_obj_subset <- df_subset
                 return() # Stop further processing in this block if UMAP fails
                 })
             
                 # --- Dynamically update DiffEx choices based on annotation ---
                 observe({
                     req(vals$seurat_obj_processed)
                     # Check if the consensus annotation column exists
                     if ("celltype_consensus" %in% colnames(vals$seurat_obj_processed@meta.data)) {
                         available_cell_types <- na.omit(unique(vals$seurat_obj_processed$celltype_consensus))
                         # Sort alphabetically for better UI
                         available_cell_types <- sort(available_cell_types)
                         updateSelectInput(session, "diffex_cell_types",
                                           choices = available_cell_types,
                                           selected = NULL) # Start with none selected
                     } else {
                          # If no annotation, keep choices empty
                          updateSelectInput(session, "diffex_cell_types", choices = list(), selected = NULL)
                     }
                 })
             
             
                 # --- 6. Differential Expression ---
                 observeEvent(input$run_diffex, {
                     req(vals$seurat_obj_processed)
                     req(input$diffex_cell_types)
                     req("celltype_consensus" %in% colnames(vals$seurat_obj_processed@meta.data)) # Ensure annotation ran
             
                     selected_types_diffex <- input$diffex_cell_types
                     if (length(selected_types_diffex) < 2) {
                         showNotification("Please select at least two cell types for comparison.", type = "warning")
                         return()
                     }
             
                     seurat_full <- vals$seurat_obj_processed
             
                     # Clear previous results
                     vals$diffex_markers <- NULL
                     output$diffexMarkerTable <- renderTable(NULL)
                     output$diffexVlnPlot <- renderPlot(NULL)
             
                     withProgress(message = 'Running Differential Expression...', value = 0, {
             
                         incProgress(0.1, detail = "Subsetting data...")
                         # Step 1: Subset based on selected cell types
                         Idents(seurat_full) <- "celltype_consensus"
                         df_diffex_subset <- tryCatch({
                             subset(seurat_full, idents = selected_types_diffex)
                         }, error = function(e) {
                             showNotification(paste("Error subsetting for DiffEx:", e$message), type = "error")
                             return(NULL)
                         })
             
                         if (is.null(df_diffex_subset) || ncol(df_diffex_subset) == 0) {
                             showNotification("Subsetting resulted in zero cells.", type = "warning")
                             return()
                         }
             
                         # Step 2: Normalize and Prep (Using SCTransform as per README Step 9 example)
                         # Note: This re-runs SCTransform on the subset, which can be slow.
                         # An alternative would be to use the existing RNA assay data if SCT wasn't desired.
                         assay_to_use <- "RNA" # Default to RNA assay
                         tryCatch({
                             incProgress(0.3, detail = "Running SCTransform on subset (can be slow)...")
                             # Ensure variable features are found *before* SCTransform if using it
                             # DefaultAssay(df_diffex_subset) <- "RNA" # Make sure RNA is default before finding features
                             # df_diffex_subset <- FindVariableFeatures(df_diffex_subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE) # Find features on RNA assay
                             df_diffex_subset <- SCTransform(df_diffex_subset, verbose = FALSE)
                             assay_to_use <- "SCT" # Use SCT assay for marker finding
                             incProgress(0.6, detail = "SCTransform complete.")
                         }, error = function(e) {
                              showNotification(paste("Error during SCTransform:", e$message, ". Proceeding with RNA assay."), type = "warning")
                              # If SCT fails, ensure RNA assay is default and scaled (might need re-scaling if not done before)
                              DefaultAssay(df_diffex_subset) <- "RNA"
                              if (!"scale.data" %in% slotNames(df_diffex_subset@assays$RNA)) {
                                 incProgress(0.4, detail = "Scaling RNA data...")
                                 all.genes.subset <- rownames(df_diffex_subset)
                                 df_diffex_subset <- ScaleData(df_diffex_subset, features = all.genes.subset, verbose = FALSE)
                              }
                              assay_to_use <- "RNA"
                         })
             
             
                         incProgress(0.7, detail = "Finding differential markers...")
                         # Step 3: Differential Expression Analysis using FindAllMarkers
                         # Ensure Idents are set to the cell types being compared
                         Idents(df_diffex_subset) <- "celltype_consensus"
                         diffex_markers_all <- tryCatch({
                             FindAllMarkers(
                               df_diffex_subset,
                               only.pos = TRUE, # Find positive markers for each group vs all others
                               min.pct = input$diffex_min_pct,
                               logfc.threshold = input$diffex_logfc,
                               assay = assay_to_use, # Use SCT or RNA assay
                               verbose = FALSE
                             )
                         }, error = function(e) {
                              showNotification(paste("Error finding differential markers:", e$message), type = "error")
                              return(NULL)
                         })
             
                          if (!is.null(diffex_markers_all) && nrow(diffex_markers_all) > 0) {
                              diffex_markers_all <- diffex_markers_all %>% arrange(cluster, desc(avg_log2FC))
                              vals$diffex_markers <- diffex_markers_all
                         } else {
                              vals$diffex_markers <- NULL
                              showNotification("No differential markers found with current settings.", type = "warning")
                         }
             
                         incProgress(0.9, detail = "Generating plots...")
                         # Generate Outputs
                         output$diffexMarkerTable <- renderTable({
                              req(vals$diffex_markers)
                              # Show top N markers per cluster (comparison group)
                              vals$diffex_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
                         }, rownames = FALSE)
             
                         output$diffexVlnPlot <- renderPlot({
                              req(vals$diffex_markers)
                              req(df_diffex_subset) # Need the subsetted object
                              # Get top 1 marker per cluster for plotting
                              top_diffex_markers <- vals$diffex_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
                              if (nrow(top_diffex_markers) > 0) {
                                  # Use the assay determined earlier (SCT or RNA)
                                  DefaultAssay(df_diffex_subset) <- assay_to_use
                                  VlnPlot(df_diffex_subset, features = unique(top_diffex_markers$gene), pt.size = 0.1) +
                                     ggtitle("Top Differentially Expressed Gene per Group (Violin)") +
                                     theme(plot.title = element_text(size=10))
                              } else {
                                  plot.new()
                                  title("No differential markers to plot.")
                              }
                         })
             
                         showNotification("Differential expression analysis complete.", type = "message")
             
                     }) # End withProgress
                 })
             
             })


            vals$seurat_obj_subset <- df_subset # Store the subsetted and re-analyzed object

            incProgress(0.8, detail = "Finding markers for subset clusters...")
            # Find markers for the new subset clusters
            subset_markers_all <- tryCatch({
                # Ensure Idents are set to the new clusters found in the subset
                Idents(df_subset) <- df_subset@meta.data[[paste0("RNA_snn_res.", res_subset)]] # Adjust if using different assay/resolution name
                FindAllMarkers(df_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE) # Use default thresholds or add inputs
            }, error = function(e) {
                 showNotification(paste("Error finding subset markers:", e$message), type = "error")
                 return(NULL)
            })

            if (!is.null(subset_markers_all) && nrow(subset_markers_all) > 0) {
                 subset_markers_all <- subset_markers_all %>% arrange(cluster, desc(avg_log2FC))
                 vals$subset_markers <- subset_markers_all
            } else {
                 vals$subset_markers <- NULL
                 showNotification("No markers found for subset clusters.", type = "warning")
            }


            incProgress(0.9, detail = "Generating subset plots...")
            # Generate subset plots
            output$umapSubset <- renderPlot({
                req(vals$seurat_obj_subset[["umap"]])
                DimPlot(vals$seurat_obj_subset, reduction = "umap", label = TRUE, repel = TRUE) +
                    ggtitle(paste("UMAP of Subset:", paste(selected_types, collapse=", ")))
            })

            output$subsetMarkerTable <- renderTable({
                 req(vals$subset_markers)
                 # Show top N markers per cluster
                 vals$subset_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
            }, rownames = FALSE)

            output$subsetFeaturePlot <- renderPlot({
                 req(vals$seurat_obj_subset[["umap"]])
                 req(vals$subset_markers)
                 # Get top 1 marker per cluster for plotting
                 top_subset_markers <- vals$subset_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
                 if (nrow(top_subset_markers) > 0) {
                     FeaturePlot(vals$seurat_obj_subset, features = unique(top_subset_markers$gene), reduction = "umap") &
                        theme(plot.title = element_text(size=10))
                 } else {
                     plot.new()
                     title("No subset markers to plot.")
                 }
            })

            showNotification("Subcluster analysis complete.", type = "message")
        }) # End withProgress
    })

})