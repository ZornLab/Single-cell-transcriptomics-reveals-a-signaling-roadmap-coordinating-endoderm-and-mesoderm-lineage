############ This script carries out Transcription Factor Enrichment Analysis using Seurat [v3.0]###############
############ Please Refer to Seurat [v3.0] manual for details on parameters and functions


#########Load R packages #########
library(Biobase)
library(Seurat)
library(plotly)
library(dplyr)
require(scales)
##################################


########## Reading Data and Creating Seurat Object ##################
setwd(path)														###############-----> Set Working Directory
file <- 'TFCounts.txt'							                ###############-----> TF counts across cells of Interest
data <- as.matrix(read.table(file,sep="\t",header=T,check.names=F,row.names=1))
file1 <- 'Metafile.txt'										    ###############------> Metafile [Please see github for example files]
col_data <- as.matrix(read.table(file1,sep="\t",header=T,row.names=1,check.names=F))
col_d <- data.frame(col_data)
pd <- new("AnnotatedDataFrame", data = col_d)
rownames(pd) <- pd$Cells
file2 <- 'GeneInfo.txt'										    ###############------> Gene Info
feature_data <- as.matrix(read.table(file2,sep="\t",header=T,row.names=1,check.names=F))
feature_data1 <- data.frame(feature_data)
fd <- new("AnnotatedDataFrame", data = feature_data1)
rownames(fd) <- fd$gene_short_name
M <- CreateSeuratObject(data, min.cells =0, min.features = 0, project= 'TFEnrichment')
######################################################################


############### Load Cell Annotations ###################
CellsMeta = M@meta.data
head(CellsMeta)
CellsMeta["Clusters"] <- pd$Clusters
head(CellsMeta)
CellsMetaTrim <- subset(CellsMeta, select = c("Clusters"))
head(CellsMetaTrim)
M <- AddMetaData(M, CellsMetaTrim)
head(M@meta.data)
##########################################################


############### Normalize and Scale Seurat Object ############
M <- NormalizeData(object = M, normalization.method = "LogNormalize", scale.factor = 10000)
M <- ScaleData(object = M)
##############################################################


################## Set Identity #############################
Idents(object = M) <- M@meta.data$Clusters   							############----> One can set this identity to any class of interest from meta data loaded
############################################################


################## TF Enrichment Analysis Using FindAllMarkers ###########
M.markers <- FindAllMarkers(object = M, only.pos = TRUE)
write.table(M.markers,"EnrichedTFs_acrossClusters.txt",quote=F,sep="\t")
##########################################################################


################### Heatmap Visualization ###########################
top10 <- M.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = M, features = top10$gene, size = 3, group.by="Clusters") + NoLegend()
ggsave('Heatmap_Top10_EnrichedTFs_acrossClutsers.png')
####################################################################

################### Save Robject ######################
save(M,file='Seuratv3.0_TFEnrichmentAnalysis.Robj')
#######################################################

