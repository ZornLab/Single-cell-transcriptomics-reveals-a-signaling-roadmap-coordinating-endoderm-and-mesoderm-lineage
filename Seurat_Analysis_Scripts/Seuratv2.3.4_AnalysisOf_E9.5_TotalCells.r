############ This script processes E9.5 Total Cells using Seurat [v2.3.4]###############
############ Script carries out basic filtering, global scaling based normalization and scaling using Seurat Functions
############ Script regresses out Cellcycle difference between G2m and S phase using ScaleData Function
############ Please Refer to Seurat [v2.3.4] manual for details on parameters and functions

#########Load R packages #########
library(Biobase)
library(Seurat)
library(plotly)
library(dplyr)
require(scales)
##################################

########## Reading Data and Creating Seurat Object ##################
setwd(path) ######### Please set the path to the directory where files are present
file <- CountsFile                                   #### CountsMatrix
data <- as.matrix(read.table(file,sep="\t",header=T,check.names=F,row.names=1))
file1 <- MetaFile                            #### Sample Info File
col_data <- as.matrix(read.table(file1,sep="\t",header=T,row.names=1,check.names=1))
col_d <- data.frame(col_data)
pd <- new("AnnotatedDataFrame", data = col_d)
rownames(pd) <- pd$Cells
file2 <- GeneInfo                            ###### Gene Info File
feature_data <- as.matrix(read.table(file2,sep="\t",header=T,row.names=1,check.names=F))
feature_data1 <- data.frame(feature_data)
fd <- new("AnnotatedDataFrame", data = feature_data1)
rownames(fd) <- fd$gene_short_name
M <- CreateSeuratObject(raw.data = data, min.cells =3, min.genes = 100, project= 'Total_E9.5')
######################################################################

############### Load Cell Annotations ###################
CellsMeta = M@meta.data
head(CellsMeta)
CellsMeta["Stages"] <- pd$Stages
head(CellsMeta)
CellsMetaTrim <- subset(CellsMeta, select = c("Stages"))
M <- AddMetaData(M, CellsMetaTrim)

CellsMeta = M@meta.data
CellsMeta["LineageAnnotations"] <- pd$LineageAnnotations
head(CellsMeta)
CellsMetaTrim <- subset(CellsMeta, select = c("LineageAnnotations"))
M <- AddMetaData(M, CellsMetaTrim)
head(M@meta.data)
##########################################################

############## Filtering, Normalizing and Scaling #########
file <- strsplit(file,".txt")                              #### Gets the prefix from Counts Matrix filename
mito.genes <- grep(pattern = "^mt-", x = rownames(x = M@data), value = TRUE)
percent.mito <- colSums(M@raw.data[mito.genes, ])/colSums(M@raw.data)
png("QCstats.png", width = 8, height = 15, units = 'in', res = 600)
M <- AddMetaData(M, percent.mito, "percent.mito")
VlnPlot(M, c("nGene", "nUMI", "percent.mito"), nCol=3)
dev.off()
png("Gene-plot-1.png", width = 8, height = 15, units = 'in', res = 600)
par(mfrow = c(1,2))
GenePlot(M, "nUMI", "percent.mito")
dev.off()
png("Gene-plot-2.png", width = 8, height = 15, units = 'in', res = 600)
GenePlot(M, "nUMI", "nGene")
dev.off()
M <- FilterCells(object = M, subset.names = c("nGene", "percent.mito"), low.thresholds = c(-Inf, -Inf), high.thresholds = c(7500, 0.15))
M <- NormalizeData(object = M, normalization.method = "LogNormalize", scale.factor = 10000)
outfile8 <- paste0(file, "_VariableGenePlot", ".png") ############ Please change the names if rerun with new genes
png(outfile8)
M <- FindVariableGenes(object = M, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3.5, y.cutoff = 0.5)
dev.off()
length(M@var.genes)
M <- ScaleData(object = M, vars.to.regress = c("nUMI", "percent.mito"))
#############################################################

##################### Cell Cycle Estimation #################
cc.genes <- readLines(con = pathToCellCycleGene) ####Please provide celll cycle genes [See files provided on Github for reference]
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
M <- CellCycleScoring(object = M, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
head(x = M@meta.data)
M@meta.data$CC.Difference <- M@meta.data$S.Score - M@meta.data$G2M.Score
M <- RunPCA(object = M, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
png("CellCycleEffect_Before_Removing.png")
PCAPlot(object = M)
dev.off()
png("DistributionCellCycle_Markers.png")
RidgePlot(object = M, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67"),  nCol = 2)
dev.off()
##############################################################


############################### Define Colors ###############
my_cols_lineages <- c("#8B4513", "#FF4500", "#87CEEB", "#FFD700", "#D8BFD8", "#6A5ACD", "#000000", "#FFA500", "#6B8E23")
my_cols_stages <- c("#00007f", "#0080ff")
#############################################################



################ PCA, Clustering, tSNE reduction and Visualization Before CC removal ############
M <- RunPCA(object = M, pc.genes = M@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
M <- ProjectPCA(object = M, do.print = FALSE)
M <- JackStraw(object = M, num.replicate = 100)
M <- FindClusters(object = M, reduction.type = "pca", dims.use = 1:20, resolution = 0.4, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = M)
outfile1 <- paste0(file, "_static", "_TSNE_plot_BeforeCellCycle_Removal", ".png")
M <- RunTSNE(object = M, dims.use = 1:20 , do.fast = TRUE)
TSNEPlot(object = M)
ggsave(outfile1)
outfile <- paste0(file, "_interactive", "_TSNE_plot_BeforeCellCycle_Removal", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile)
dev.off()

M <- StashIdent(object = M, save.name = "CellType")
outfile99 <- paste0(file, "_static", "_TSNE_plot_with_CellCyclePhase_BeforeCellCycle_Removal", ".png")
M <- SetAllIdent(object = M, id = "Phase")
TSNEPlot(object = M, do.identify=TRUE, do.hover=TRUE)
ggsave(outfile99)
outfile98 <- paste0(file, "_interactive", "_TSNE_plot_with_CellCyclePhase_BeforeCellCycle_Removal", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile98)

M <- StashIdent(object = M, save.name = "CellType")
outfile99 <- paste0(file, "_static", "_TSNE_plot_with_Stages_BeforeCellCycle_Removal", ".png")
M <- SetAllIdent(object = M, id = "Stages")
TSNEPlot(object = M, colors.use=my_cols_lineages)
ggsave(outfile99)
outfile98 <- paste0(file, "_interactive", "_TSNE_plot_with_Stages_BeforeCellCycle_Removal", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile98)

M <- StashIdent(object = M, save.name = "CellType")
outfile99 <- paste0(file, "_static", "_TSNE_plot_with_LineageAnnotations_BeforeCellCycle_Removal", ".png")
M <- SetAllIdent(object = M, id = "LineageAnnotations")
TSNEPlot(object = M, colors.use=my_cols_lineages)
ggsave(outfile99)
outfile98 <- paste0(file, "_interactive", "_TSNE_plot_with_LineageAnnotations_BeforeCellCycle_Removal", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile98)
######################################################################################################


################ PCA, Clustering, tSNE reduction and Visualization After CC removal ############
M <- ScaleData(object = M, vars.to.regress = "CC.Difference", display.progress = FALSE)
M <- RunPCA(object = M, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
png("CellCycleEffect_After_Removing.png")
PCAPlot(object = M)
dev.off()
M <- RunPCA(object = M, pc.genes = M@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
M <- ProjectPCA(object = M, do.print = FALSE)M <- JackStraw(object = M, num.replicate = 100)
M <- FindClusters(object = M, reduction.type = "pca", dims.use = 1:20, resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
PrintFindClustersParams(object = M)
outfile1 <- paste0(file, "_static", "_TSNE_plot_AfterCellCycle_Removal", ".png")
M <- RunTSNE(object = M, dims.use = 1:20 , do.fast = TRUE)
TSNEPlot(object = M)
ggsave(outfile1)
outfile <- paste0(file, "_interactive", "_TSNE_plot_AfterCellCycle_Removal", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile)
dev.off()
M.markers <- FindAllMarkers(object = M, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(M.markers,file="Markers_All_Clusters_AfterCellCycle_Removal.txt",sep="\t",col.names= NA,quote=F)
write.table(M@ident,"Samples_Clusters_Association_AfterCellCycle_Removal.txt",sep="\t",col.names= NA,quote=F)

M <- StashIdent(object = M, save.name = "CellType")
outfile99 <- paste0(file, "_static", "_TSNE_plot_with_CellCyclePhase_AfterCellCycle_Removal", ".png")
M <- SetAllIdent(object = M, id = "Phase")
TSNEPlot(object = M, do.identify=TRUE, do.hover=TRUE)
ggsave(outfile99)
outfile98 <- paste0(file, "_interactive", "_TSNE_plot_with_CellCyclePhase_AfterCellCycle_Removal", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile98)

M <- StashIdent(object = M, save.name = "CellType")
outfile99 <- paste0(file, "_static", "_TSNE_plot_with_Stages_AfterCellCycle_Removal", ".png")
M <- SetAllIdent(object = M, id = "Stages")
TSNEPlot(object = M, colors.use=my_cols_lineages)
ggsave(outfile99)
outfile98 <- paste0(file, "_interactive", "_TSNE_plot_with_Stages_AfterCellCycle_Removal", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile98)

M <- StashIdent(object = M, save.name = "CellType")
outfile99 <- paste0(file, "_static", "_TSNE_plot_with_LineageAnnotations_AfterCellCycle_Removal", ".png")
M <- SetAllIdent(object = M, id = "LineageAnnotations")
TSNEPlot(object = M, colors.use=my_cols_lineages)
ggsave(outfile99)
outfile98 <- paste0(file, "_interactive", "_TSNE_plot_with_LineageAnnotations_AfterCellCycle_Removal", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile98)

M <- SetAllIdent(object = M, id = "res.1")
######################################################################################################

######################### Save Robject ###############################################################
save(M, file = 'Seurat_v2.3.4_AnalysisOf_E9.5_TotalCells.Robj')
######################################################################################################