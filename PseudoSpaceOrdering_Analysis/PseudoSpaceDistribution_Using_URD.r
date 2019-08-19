############ This script processes Pseudospace Ordering of cells using URD [v1.0]###############
############ Please Refer to URD [v1.0] manual for details on parameters and functions

#########Load R packages #########
library(URD)
library(Seurat)
library(plotly)
library(ggplot2)
library(destiny)
##################################


############### Read data and create URD Object ##########
setwd('/Volumes/Praneet_HD/Hiro_Takebe_Data/New_SetOfAnalysis_Lu_Praneet/Pseudotime_Density_Distribution_Endoderm_Mesoderm/Endoderm/E8.5')
file <- 'Expression_Endoderm_E8.5_subset.txt'
count.total <- as.matrix(read.table("Expression_Endoderm_E8.5_subset.txt", header=T, sep="\t", row.names=1, check.names=F)) 
file <- strsplit(file,".txt")
meta.total <- read.table("MetaInformation_URD_EndodermAaron_E8.5.txt", header=T, sep="\t", row.names=1, check.names=F)
total <- createURD(count.data = count.total, meta = meta.total, min.cells=3, min.counts=3)
total@group.ids$stage <- as.character(total@meta[rownames(total@group.ids),"Stages"])
total@group.ids$scluster <- as.character(total@meta[rownames(total@group.ids),"Clusters"])
############################################################

#######################Reodering Cells according to Clusters############
my_cluster_order <- c('a5', 'a1', 'a2', 'a3', 'a4', 'a0')       ############----> OPTIONAL
total@group.ids$scluster <- factor(x = total@group.ids$scluster, levels = my_cluster_order)     ##########----> OPTIONAL
########################################################################

################ Get Seurat Markers ##########################
stages <- sort(unique(total@group.ids$stage))
seuratclusters <- sort(unique(total@group.ids$scluster))

file_markers <- 'Seurat_Markers.txt'          ##############------> Using Seurat Markers
markers_temp <- as.matrix(read.table(file_markers, sep="\t", header=T))
markers <- data.frame(markers_temp)
var.genes <- sort(unique(markers$gene))
total@var.genes <- var.genes
##############################################################


############### Run PCA, DiffusionMap & calculate Pseudotime ######################
total <- calcPCA(total, mp.factor = 2)
set.seed(19)
total <- calcTsne(object = total)
total <- calcDM(total, knn = 50, sigma='local')
root.cells <- cellsInCluster(total, "scluster", "a5")           #####################------> setting root cells from cluster a5 [One can set accordingly]
total.floods <- floodPseudotime(total, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)
total <- floodPseudotimeProcess(total, total.floods, floods.name="pseudotime")
########################################################################################


################## Pseudospace ordering ####################
plotDists(total, "pseudotime", "scluster", plot.title="Pseudotime by Clusters")
ggsave('Pseudospace_Ordering_ofClusters.png')
htmlwidgets::saveWidget(as.widget(ggplotly()), 'Pseudospace_Ordering_ofClusters.html')
###########################################################
