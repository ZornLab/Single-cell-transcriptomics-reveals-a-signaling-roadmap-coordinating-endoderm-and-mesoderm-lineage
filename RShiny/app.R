if (!require("shiny")) install.packages("shiny",repos="http://cran.us.r-project.org")
if (!require("servr")) install.packages("servr",repos="http://cran.us.r-project.org")
if (!require("DT")) install.packages("DT",repos="http://cran.us.r-project.org")
if (!require("Seurat")) install.packages("Seurat",repos="http://cran.us.r-project.org")
if (!require("ggplot2")) install.packages("ggplot2",repos="http://cran.us.r-project.org")
if (!require("dplyr")) install.packages("dplyr",repos="http://cran.us.r-project.org")
if (!require("ggcorrplot")) install.packages("ggcorrplot",repos="http://cran.us.r-project.org")
if (!require("reticulate")) install.packages("reticulate",repos="http://cran.us.r-project.org")

library(shiny)
library(reticulate)
library(servr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(DT)
library(ggcorrplot)


source("helper.r")


ui <- fluidPage(
  
  # App title ----
  titlePanel("Single cell transcriptomics reveals a signaling roadmap coordinating endoderm and mesoderm lineage diversification during foregut organogenesis"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of Dataset ----
      selectInput("Dataset",
                  label = "Choose a Lineage",
                  choices = c("Foregut Endoderm Cells", "Foregut Mesoderm Cells", "Foregut All Cells"),
                  selected = "Foregut Endoderm Cells"),
                  
      selectInput("DR",
                  label = "Choose a Dimension Reduction",
                  choices = c("tsne", "umap"),
                  selected = "tsne"),
                  
      selectInput("Colorby",
                  label = "Group Cells By",
                  choices = c("Stages", "Stage Specific Lineage Annotations", "Integration based seurat clusters"),
                  selected = "Stage Specific Lineage Annotations"),
                  
                  
      selectizeInput("Gene",
                  label = "Choose a Gene of Interest",
                  choices = genelist,
                  selected = "0610007P14Rik", options = list(placeholder = "Type a gene", maxOptions = 19111)),
                  
      fileInput("File",
                  label = "Optional: Upload Gene List for Dot Plot or Correlation Plot (Max Limit: DotPlot (40 Genes) and Correlation Plot (100 Genes) (Accepted Format: .txt)",
                  accept = c("text", "text", ".txt")),
                  
      selectInput("Correlation",
                  label = "Correlation Method",
                  choices = c("pearson", "kendall", "spearman"),
                  selected = "pearson"),
      
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("About", fluidRow(
                    p(strong("Authors:"), "Lu Han, Hiroyuki Koike, Praneet Chaturvedi, Keishi Kishimoto, Kentaro Iwasawa, Kirsten Giesbrecht, Phillip C Witcher, Alexandra Eicher, Talia Nasr, Lauren Haines, John M Shannon, Mitsuru Morimoto, James M Wells, Takanori Takebe, Aaron M Zorn",style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Abstract:"), "Visceral organs, such as the lungs, stomach, liver and pancreas, derive from the fetal foregut through a series of inductive interactions between the definitive endoderm (DE) epithelium and the surrounding splanchnic mesoderm (SM). This foregut patterning, which occurs between embryonic day (E) 8.5 and E9.5 in the mouse embryo, equivalent to 17-23 days of human gestation, defines the landscape of the thoracic cavity and disruptions in this process can lead to severe congenital defects. While patterning of the endoderm lineages has been fairly well studies, the SM which is known to provide many paracrine factors required for organogeneis is virtually unstudied. In particular we lack a comprehensive understanding of the molecular nature of SM regional identity, the mechanisms by which SM signaling boundaries are established, the role of the epithelium in SM patterning and how SM and DE lineages are dynamically coordinated during organogenesis. Here we used single cell transcriptomics to generate a high-resolution expression map of the embryonic mouse foregut. This uncovered an unexpected diversity in SM progenitors that developed in close register with the organ-specific epithelium. These data allowed us to infer a spatial and temporal signaling roadmap of the combinatorial endoderm-mesoderm interactions that orchestrate foregut organogenesis. We validated key predictions with mouse genetics, showing importance of epithelial signaling in mesoderm patterning. Finally, we leveraged the signaling road map to generate different SM subtypes from human pluripotent stem cells (hPSCs), which previously have been elusive.",style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Affiliation:"), "Center for Stem Cell & Organoid Medicine (CuSTOM), Perinatal Institute, Divisions of Developmental Biology, Gastroenterology and Pulmonary Biology, Cincinnati Children’s Hospital, Cincinnati OH, USA 45229. Department of Pediatrics University of Cincinnati College of Medicine. Laboratory for Lung Development, RIKEN Center for Biosystems Dynamics Research (BDR), Kobe, 650-0047, Japan, CuSTOM-RIKEN BDR collaborative laboratory, Cincinnati OH, USA 45229.",style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Explore Manuscript at bioRxiv:"), a('Manuscript', href="https://www.biorxiv.org/content/10.1101/756825v1.full", target="_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Data and Code Availability:"), "The scRNA-seq and bulk RNA-seq data are available at Gene Expression Omnibus (GEO): GSE136689 and GSE136687. All the code (scripts, R-packages and software) and their documentation has been uploaded to GitHub. All the deposited code is available to use with GPLv3.0. Access", a("GitHub", href="https://github.com/ZornLab/Single-cell-transcriptomics-reveals-a-signaling-roadmap-coordinating-endoderm-and-mesoderm-lineage", target="_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("References:"), a("Seurat", href="https://www.nature.com/articles/nbt.4096", target="_blank"), "&", a("SPRING", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6030950/", target="_blank"),style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Developed By:"), "Praneet Chaturvedi. To view tool contributions visit", a("GitHub", href="https://github.com/praneet1988", target="_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"))),
                  tabPanel("SPRING KNN GRAPH", uiOutput("Link"), textOutput("Display")),
                  tabPanel("Dimension Reduction Plot", downloadButton('downloaddimpplot', 'SavePlot'), plotOutput("DimPlot")),
                  tabPanel("Gene Expression", downloadButton('downloadexpplot', 'SavePlot'), plotOutput("Exp")),
                  tabPanel("Markers", downloadButton('downloadData', 'Markers'), DT::dataTableOutput("Markers")),
                  tabPanel("DotPlot", downloadButton('downloaddotplot', 'SavePlot'), plotOutput("DotPlot")),
                  tabPanel("Correlation Plot", downloadButton('downloadcorrplot', 'SavePlot'), plotOutput("CorrPlot")),
                  tabPanel("Tutorial", fluidRow(
                    p("This RShiny app is developed to host single cell RNA-Seq data reported in the paper Single cell transcriptomics reveals a signaling roadmap coordinating endoderm and mesoderm lineage diversification during foregut organogenesis. Mouse embryos were sequenced at three timepoints E8.5, E9.0 and E9.5 to study the foregut lineage diversification. Please see tab (About)",style="text-align:left;color:black;background-color:white;padding:20px;border-radius:15px"),
                    p(strong("Datasets:"), "Currently there are three datasets which users can explore from Choose a Lineage section in the app. 1) Foregut All Cells - all cells sequenced at each stage were integrated using Seurat v3.0 integration algorithm by finding integration anchors (cell pairs) and then integrating the data. Please read more on", a("Integration using Seurat", href="https://www.nature.com/articles/nbt.4096", target="_blank"), "2) Foregut Mesoderm Cells - all splanchnic mesoderm cells were sub-setted from all cells at each stage using expression of following markers Foxf1/Vim/Pdgfr. All mesoderm cells at each stage were then integrated using Seurat. 3) Foregut Endoderm Cells - all endoderm cells were sub-setted from all cells at each stage using expression of following markers Foxa1/Cdh1/Epcam. All endoderm cells at each stage were then integrated using Seurat.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Visualize Cells:"), "Users can visualize cells using two non-linear dimension reduction approaches 1) tSNE (t-stochastic neighbor embedding approach) 2) UMAP (Uniform manifold approximation projection) in Seurat and also in SPRING which uses a knn force directed layout. Please see Tabs (SPRING KNN GRAPH and Dimension Reduction Plot)", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Group Cells By:"), "Cells can be grouped and colored by three different categories 1) Stages (E8.5, E9.0 and E9.5) 2) Stage Specific Lineage Annotations (at each stage cells were annotated based on apriori knowledge of cell-type markers and using in-situ images in Mouse Genome informatics", a("MGI", href="http://www.informatics.jax.org/", target="_blank"), "3) Integration based Seurat clusters (clusters which were generated using graph based KNN algorithm after integration using Seurat). For simplifying annotations E8.5 is set to a, E9.0 is set to b and E9.5 is set to c in Stage Specific Lineage Annotations", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Visualize Gene Expression:"), "Expression of a Gene can be visualized by selecting a gene of interest from gene list menu. FeaturePlot function in Seurat is used to generation gene expression plot. To visualize the expression please see tab (Gene Expression)", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("DotPlot and Correlation Plot:"), "User defined list of genes can be uploaded to generate 1) DotPlot (shows expression and % of cells expressed for genes across Category selected in Group cells by) 2) Correlation Plot (shows correlation among genes based on average expression of categories selected in Group Cells By). For correlation users can select different correlation method 1) Pearson 2) Spearman 3) Kendall", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Download Plots:"), "Each plot generated can be saved to your computer using SavePlot button. User can provide a name for a plot (by default: app defines one)", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Download Markers:"), "Markers tab shows markers for each category in Group cells by and Dataset selected. All marker files can be downloaded using Markers download button on tab (Markers)", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                    p(strong("Improvements/Bugs/Suggestions:"), "Please use", a("GitHub", href="https://github.com/ZornLab/Single-cell-transcriptomics-reveals-a-signaling-roadmap-coordinating-endoderm-and-mesoderm-lineage", target="_blank"), "for reporting issues, bugs or improvements", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px")))
      )
    )
  )
)
server <- function(input, output) {

  Group_cells <- reactive({
    if(input$Colorby == "Stages")
      {
        Group_cells = "Stages"
      }
    else if(input$Colorby == "Stage Specific Lineage Annotations")
      {
        Group_cells = "LineageAnnotations"
      }
    else if(input$Colorby == "Integration based seurat clusters")
      {
        Group_cells = "seurat_clusters"
      }
    })
  
  output$Link <- renderUI(
    if(input$Dataset == "Foregut Endoderm Cells") {
      source_python("python.py")
      url <- a(href=paste0('http://localhost:8000/springViewer.html?datasets/Foregut_Endoderm'), "ForegutEndoderm", target="_blank")
    }
    
    else if(input$Dataset == "Foregut Mesoderm Cells") {
      source_python("python.py")
      url <- a(href=paste0('http://localhost:8000/springViewer.html?datasets/Foregut_Mesoderm'), "ForegutMesoderm", target="_blank")
    }
    
    else if(input$Dataset == "Foregut All Cells") {
      url <- a(href=paste0('https://kleintools.hms.harvard.edu/tools/spring.html'), "NOT AVAILABLE", target="_blank")
    }
  )
  
  DimPlotinput <- reactive({
     if(input$Dataset == "Foregut Endoderm Cells") {
       if(Group_cells() == "Stages") {
          my_cols_stages <- c("#ff9400", "#7cff79", "#0080ff", "#00007f")
          DimPlot(Integration.Endoderm_Reduced, reduction = input$DR, pt.size=1, group.by = Group_cells(), cols = my_cols_stages) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5))) 
          }
        else {
          DimPlot(Integration.Endoderm_Reduced, reduction = input$DR, pt.size=1, group.by = Group_cells()) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
         }
     }
     else if(input$Dataset == "Foregut Mesoderm Cells") {
       if(Group_cells() == "Stages") {
         my_cols_stages <- c("#ff9400", "#7cff79", "#0080ff", "#00007f")
         DimPlot(Integration.Mesoderm_Reduced, reduction = input$DR, pt.size=1, group.by = Group_cells(), cols = my_cols_stages) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
         }
        else {
          DimPlot(Integration.Mesoderm_Reduced, reduction = input$DR, pt.size=1, group.by = Group_cells()) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
          }
     }
     else if(input$Dataset == "Foregut All Cells") {
       if(Group_cells() == "Stages") {
         my_cols_stages <- c("#ff9400", "#7cff79", "#0080ff", "#00007f")
         DimPlot(Integration.Total_Reduced, reduction = input$DR, pt.size=1, cols = my_cols_stages, group.by = Group_cells()) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
         }
       else if(Group_cells() == "LineageAnnotations") {
         my_cols_lineages <- c("#8B4513", "#FF4500", "#87CEEB", "#FFD700", "#6A5ACD", "#A0A0A0", "#D8BFD8", "#6B8E23", "#FFA500")
         DimPlot(Integration.Total_Reduced, reduction = input$DR, pt.size=1, cols = my_cols_lineages, group.by = Group_cells()) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
         }
       else if(Group_cells() == "seurat_clusters") {
         DimPlot(Integration.Total_Reduced, reduction = input$DR, pt.size=1, group.by = Group_cells()) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
         }
     }
  })
  
  Exp <- reactive({ 
     if(input$Dataset == "Foregut Endoderm Cells") {
       FeaturePlot(Integration.Endoderm_Reduced, features = input$Gene, min.cutoff = "q9", reduction = input$DR, pt.size=1)
     }
     else if(input$Dataset == "Foregut Mesoderm Cells") {
       FeaturePlot(Integration.Mesoderm_Reduced, features = input$Gene, min.cutoff = "q9", reduction = input$DR, pt.size=1)
     }
     else if(input$Dataset == "Foregut All Cells") {
       FeaturePlot(Integration.Total_Reduced, features = input$Gene, min.cutoff = "q9", reduction = input$DR, pt.size=1)
     }
 })
  
 output$Markers <- DT::renderDataTable({
     if(input$Dataset == "Foregut Endoderm Cells") {
       if(Group_cells() == "Stages") {
         file_markers <- 'data/Markers_Endoderm_Integration_Stages.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
       else if(Group_cells() == "LineageAnnotations") {
         file_markers <- 'data/Markers_Endoderm_Integration_LineageAnnotations.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
       else if(Group_cells() == "seurat_clusters") {
         file_markers <- 'data/Markers_Endoderm_Integration_seurat_clusters.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
     }
     else if(input$Dataset == "Foregut Mesoderm Cells") {
       if(Group_cells() == "Stages") {
         file_markers <- 'data/Markers_Mesoderm_Integration_Stages.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
       else if(Group_cells() == "LineageAnnotations") {
         file_markers <- 'data/Markers_Mesoderm_Integration_LineageAnnotations.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
       else if(Group_cells() == "seurat_clusters") {
         file_markers <- 'data/Markers_Mesoderm_Integration_seurat_clusters.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
     }
     else if(input$Dataset == "Foregut All Cells") {
       if(Group_cells() == "Stages") {
         file_markers <- 'data/Markers_Total_Integration_Stages.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
       else if(Group_cells() == "LineageAnnotations") {
         file_markers <- 'data/Markers_Total_Integration_LineageAnnotations.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
       else if(Group_cells() == "seurat_clusters") {
         file_markers <- 'data/Markers_Total_Integration_seurat_clusters.txt'
         markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
         markers <- data.frame(markers_temp)
         DT::datatable(markers)
       }
     }
  })
  
  DotPlotinput <- reactive({ 
      inFile <- input$File
     if (is.null(inFile))
       return(NULL)
      if(input$Dataset == "Foregut Endoderm Cells") {
       if(Group_cells() == "Stages") {
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Endoderm_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
       else if(Group_cells() == "LineageAnnotations") {
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Endoderm_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
       else if(Group_cells() == "seurat_clusters") {
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Endoderm_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
     }
     else if(input$Dataset == "Foregut Mesoderm Cells") {
       if(Group_cells() == "Stages") {
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Mesoderm_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
       else if(Group_cells() == "LineageAnnotations") {
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Mesoderm_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
       else if(Group_cells() == "seurat_clusters") {
        GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Mesoderm_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
     }
     else if(input$Dataset == "Foregut All Cells") {
       if(Group_cells() == "Stages") {
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Total_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
       else if(Group_cells() == "LineageAnnotations") {
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Total_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
       else if(Group_cells() == "seurat_clusters") {
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         DotPlot(Integration.Total_Reduced, features = unique(GenesUpload), group.by = Group_cells()) + RotatedAxis()
       }
     }
  })
  
  CorrPlot <- reactive({
     inFile <- input$File
     if (is.null(inFile))
       return(NULL)
     if(input$Dataset == "Foregut Endoderm Cells") {
       if(Group_cells() == "Stages") {
         file_expression <- 'data/AverageExpression_IntegratedEndoderm_Stages.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
       else if(Group_cells() == "LineageAnnotations") {
         file_expression <- 'data/AverageExpression_IntegratedEndoderm_LineageAnnotations.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
       else if(Group_cells() == "seurat_clusters") {
         file_expression <- 'data/AverageExpression_IntegratedEndoderm_SeuratClusters.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
     }
     else if(input$Dataset == "Foregut Mesoderm Cells") {
       if(Group_cells() == "Stages") {
         file_expression <- 'data/AverageExpression_IntegratedMesoderm_Stages.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
       else if(Group_cells() == "LineageAnnotations") {
         file_expression <- 'data/AverageExpression_IntegratedMesoderm_LineageAnnotations.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
       else if(Group_cells() == "seurat_clusters") {
         file_expression <- 'data/AverageExpression_IntegratedMesoderm_SeuratClusters.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
     }
     else if(input$Dataset == "Foregut All Cells") {
       if(Group_cells() == "Stages") {
         file_expression <- 'data/AverageExpression_IntegratedTotal_Stages.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
       else if(Group_cells() == "LineageAnnotations") {
         file_expression <- 'data/AverageExpression_IntegratedTotal_LineageAnnotations.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
       else if(Group_cells() == "seurat_clusters") {
         file_expression <- 'data/AverageExpression_IntegratedTotal_SeuratClusters.txt'
         data <- as.matrix(read.table(file_expression, sep="\t", header=T, row.names=1, check.names=F))
         GenesUpload <- as.matrix(read.table(inFile$datapath, sep="\n", header=F))
         GenesUpload <- data.frame(GenesUpload)
         GenesUpload <- GenesUpload$V1
         datause <- subset(data, rownames(data) %in% GenesUpload)
         datause <- t(datause)
         cormat <- cor(datause, method=input$Correlation)
         titleuse <- paste0("Displaying Correlation Plot based on Average Expression of ", input$Dataset, " ",  Group_cells())
         ggcorrplot(cormat, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
       }
     }
  })
  
  output$DimPlot <- renderPlot({
     DimPlotinput()
  }, height = 850, width = 1200)
  
  output$Exp <-renderPlot({
     Exp()
  }, height = 850, width = 1000)
  
  output$DotPlot <- renderPlot({
     DotPlotinput()
     }, height = 1200, width = 1200)
     
  output$CorrPlot <- renderPlot({
     CorrPlot()
     }, height = 1000, width = 1000)
 
 output$Display <- renderText({"Please be patient, loading SPRING can take upto 3-5 mins"})
 
 output$downloaddimpplot <- downloadHandler(
      filename = function() {
        paste("DimPlot", "_", input$Dataset, "_", input$DR, "_", input$Colorby, "-", Sys.Date(), ".png")
      },
      content = function(file) {
        DimPlotinput()
        ggsave(file, width=15, height=15)
      })
      
 output$downloadexpplot <- downloadHandler(
      filename = function() {
        paste("GeneExpression", "_", input$Gene, "_", input$Dataset, "_", input$DR, "_", input$Colorby, "-", Sys.Date(), ".png")
      },
      content = function(file) {
        Exp()
        ggsave(file, width=15, height=15)
      })
      
 output$downloaddotplot <- downloadHandler(
      filename = function() {
        paste("DotPlot", "_", input$Dataset, "_", input$DR, "_", input$Colorby, "-", Sys.Date(), ".png")
      },
      content = function(file) {
        DotPlotinput()
        ggsave(file, width=15, height=15)
      })
  
  output$downloadcorrplot <- downloadHandler(
      filename = function() {
        paste("CorrelationPlot", "_", input$Dataset, "_", input$DR, "_", input$Colorby, "-", Sys.Date(), ".png")
      },
      content = function(file) {
        CorrPlot()
        ggsave(file, width=15, height=15)
      })
      
  output$downloadData <- downloadHandler(
     filename = function() {
        paste("Markers", "zip", sep=".")
      },
      content = function(file) {
        fs <- c()
        dir <- "Markers"
        path <- paste0(dir, ".zip")
        fs <- c(fs, path)
        zip(zipfile=file, files=fs)
        },
      contentType = "application/zip"
    )
}
scRNA_app <- shinyApp(ui = ui, server = server)
runApp(scRNA_app, launch.browser = TRUE, display.mode = "normal")