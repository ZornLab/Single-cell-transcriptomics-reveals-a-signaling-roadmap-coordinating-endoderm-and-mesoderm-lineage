############ This script processes mm10 folder outputted from CellRanger using Seurat [v2.3.4]###############
############ Script writes out counts matrix
############ Please Refer to Seurat [v2.3.4] manual for details on parameters and functions


#########Load R packages #########
library(Seurat)
library(plotly)
library(dplyr)
##################################

########## Reading Data, Creating Seurat Object and Writing Counts Matrix ##################
setwd(path) 				######### Please set the path to the directory where files are present
M.data <- Read10X(data.dir = 'pathto_mm10')
M <- CreateSeuratObject(raw.data = M.data, min.cells =0, min.genes = 0, project= 'ProjectName')
counts <- as.matrix(M@data)
write.table(counts, "Counts_FileName.txt", sep="\t", quote=F)
############################################################################################