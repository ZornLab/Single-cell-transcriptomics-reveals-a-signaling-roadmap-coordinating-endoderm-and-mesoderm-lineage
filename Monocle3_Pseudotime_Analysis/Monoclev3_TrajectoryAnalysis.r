############ This script processes Pseudotime Analysis using Monocle3 [v3.0 alpha]###############
############ Please note only Markers obtained from Seurat were used to drive Pseudotime analysis
############ Please Refer to Monocle [v3.0 alpha] manual for details on parameters and functions

############## Load R packages ########
library(reticulate)
import("louvain")
library(rmarkdown)
library(plotly)
rm(list = ls())
options(warn=-1)
library(monocle)
######################################



########## Reading Data and Creating Monocle Object ##################
setwd(path) 						#########-----> Please set the path to the directory where files are present
file <- CountsFile                                   #### CountsMatrix
data <- as.matrix(read.table(file,sep="\t",header=T,check.names=F,row.names=1))
file1 <- MetaFile                   #########------> MetaFile
col_data <- as.matrix(read.table(file1,sep="\t",header=T,row.names=1,check.names=1))
col_d <- data.frame(col_data)
pd <- new("AnnotatedDataFrame", data = col_d)
rownames(pd) <- pd$Cells
file2 <- GeneInfo                   #########------> Gene Info File [Seurat Markers]
feature_data <- as.matrix(read.table(file2,sep="\t",header=T,row.names=1,check.names=F))
feature_data1 <- data.frame(feature_data)
fd <- new("AnnotatedDataFrame", data = feature_data1)
rownames(fd) <- fd$gene_short_name
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = VGAM::negbinomial.size())
######################################################################



############### Normalize and Preprocess the Data #######################
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- preprocessCDS(cds, num_dim = 20)
########################################################################


############### Dimensionality Reduction using tSNE ####################
cds <- reduceDimension(cds, reduction_method = 'tSNE')
########################################################################


###############	Partioning Cells ######################################
cds <- partitionCells(cds)
#######################################################################


############### Learning the Principal Graph ##########################
cds <- learnGraph(cds,  RGE_method = 'SimplePPT')
#######################################################################


############## Visualizing the learnt trajectory #############################
plot_cell_trajectory(cds, color_by = "Stages")
ggsave('Monoclev3_Trajectory_coloredByStages_withTSNE.png')

plot_cell_trajectory(cds, color_by = "Clusters")
ggsave('Monoclev3_Trajectory_coloredByClusters_withTSNE.png')
##############################################################################


############## Ordering of Cells ###########################################
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
 closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
root_IDs = get_correct_root_state(cds,
                                      cell_phenotype =
                                        'Clusters', "a0")      ###############------> setting the root as cluster a0
cds <- orderCells(cds, root_pr_nodes = root_IDs)
plot_cell_trajectory(cds)
ggsave('Monoclev3_Trajectory_coloredByPseudotime_withTSNE.png')
############################################################################


############## Learning and Visualizing trajectories in 3-D ################
cds <- reduceDimension(cds, max_components = 3,
                       reduction_method = 'tSNE',
                       metric="cosine",
                       verbose = F)

cds <- partitionCells(cds)

cds <- learnGraph(cds,
                  max_components = 3,
                  RGE_method = 'SimplePPT',
                  partition_component = T,
                  verbose = F)
                  
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
 closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

root_IDs = get_correct_root_state(cds, cell_phenotype = 'Clusters', "a0")
cds <- orderCells(cds, root_pr_nodes = root_IDs)

plot_3d_cell_trajectory(cds,
                        color_by="Pseudotime",
                        webGL_filename=
                          paste(getwd(), "/Monoclev3_3D_trajectory_coloredbyPseudotime_withTSNE.html", sep=""),
                        show_backbone=TRUE,
                        useNULL_GLdev=TRUE) + scale_color_gradient2(midpoint=mid_Cited1, low="blue", mid="white", high="red", space ="Lab" )
                        
plot_3d_cell_trajectory(cds,
                        color_by="Clusters",
                        webGL_filename=
                          paste(getwd(), "/Monoclev3_3D_trajectory_coloredbyClusters_withTSNE.html", sep=""),
                        show_backbone=TRUE,
                        useNULL_GLdev=TRUE)
                        
plot_3d_cell_trajectory(cds,
                        color_by="Stages",
                        webGL_filename=
                          paste(getwd(), "/Monoclev3_3D_trajectory_coloredbyStages_withTSNE.html", sep=""),
                        show_backbone=TRUE,
                        useNULL_GLdev=TRUE)
##############################################################################



################ Identifying Variable Genes across Trajectory #################
pr_graph_test <- principalGraphTest(cds, k=3, cores=1)
dplyr::add_rownames(pr_graph_test) %>% dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-qval)) %>% head(3)
pr_graph_test[fData(cds)$gene_short_name %in% c("Cdx2"),]         ##############------> Check for favourite Marker
pr_graph_test[fData(cds)$gene_short_name %in% c("Foxf1"),]        ##############------> Check for favourite Marker
nrow(subset(pr_graph_test, qval < 0.01))
DiffGenesAcrossTraj <- subset(pr_graph_test, qval < 0.01)
write.table(DiffGenesAcrossTraj, "DifferentialGenes_AcrossTrajectory.txt", sep="\t", quote=F)    
###############################################################################                   

################# Saving Robject ######################################
save(cds, file="Monoclev3_Pseudotime-based_TrajectoryAnalysis.Robj")
#######################################################################