############ This script processes MetaProfiles of GeneSets to create DotPlots using Seurat [v3.0]###############
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
file <- 'MetaProfile_GeneSets.txt'							    ###############-----> Output from GenerateMetaprofile_ForGeneSets.pl
data <- as.matrix(read.table(file,sep="\t",header=T,check.names=F,row.names=1))
file1 <- 'Metafile.txt'										    ###############------> Metafile [Please github for example files]
col_data <- as.matrix(read.table(file1,sep="\t",header=T,row.names=1,check.names=1))
col_d <- data.frame(col_data)
pd <- new("AnnotatedDataFrame", data = col_d)
rownames(pd) <- pd$Cells
file2 <- 'GeneInfo.txt'										    ###############------> Gene Info
feature_data <- as.matrix(read.table(file2,sep="\t",header=T,row.names=1,check.names=F))
feature_data1 <- data.frame(feature_data)
fd <- new("AnnotatedDataFrame", data = feature_data1)
rownames(fd) <- fd$gene_short_name
M <- CreateSeuratObject(data, min.cells =3, min.features = 1, project= 'MetaPathwayProfile')
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


################# Set the Cluster Order [Optional]###################
my_cluster_order <- c("M_c14", "M_c13", "M_c9", "M_c2", "M_c10", "M_c4", "M_c12", "M_c16", "M_c0", "M_c8", "M_c1", "M_c6", "M_c15", "M_c5", "M_c7", "M_c11", "M_c3", "M_b11", "M_b8", "M_b10", "M_b2", "M_b7", "M_b5", "M_b1", "M_b6", "M_b4", "M_b3", "M_b0", "M_b9", "M_a2", "M_a1", "M_a6", "M_a0", "M_a3", "M_a5", "M_a4", "E_c10", "E_c1", "E_c11", "E_c0", "E_c9", "E_c6", "E_c2", "E_c7", "E_c5", "E_c4", "E_c8", "E_c3", "E_b5", "E_b2", "E_b7", "E_b1", "E_b6", "E_b0", "E_b4", "E_b3", "E_a0", "E_a4", "E_a3", "E_a2", "E_a1", "E_a5")

M@meta.data$Clusters <- factor(x = M@meta.data$Clusters, levels = my_cluster_order)
##########################################################


################ Generate DotPlot #########################
DotPlot(M,rownames(M), group.by="Clusters", cols=c("green", "red"), dot.scale=5, assay="RNA") + RotatedAxis()
ggsave('MetaProfile_GeneSets_DotPlot.png')
###########################################################

################ Save Robject #############################
save(M, file='Seuratv3.0_MetaProfile_DotPlot.Robj')
##########################################################