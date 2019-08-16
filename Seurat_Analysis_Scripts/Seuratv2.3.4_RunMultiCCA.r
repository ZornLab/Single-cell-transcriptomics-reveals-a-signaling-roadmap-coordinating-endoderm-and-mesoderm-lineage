############ This script runs Canonical Correlation Analysis[CCA] using Seurat [v2.3.4]###############
############ Script carries out basic filtering, global scaling based normalization and scaling using Seurat Functions
############ Script regresses out Cellcycle difference between G2m and S phase using ScaleData Function
############ Please Refer to Seurat [v2.3.4] manual for details on parameters and functions
############ Please remove Ribosomal genes from counts matrices before running this Script

#########Load R packages #########
library(Biobase)
library(Seurat)
library(plotly)
library(dplyr)
require(scales)
##################################


########## Reading Data and Creating Seurat Objects for each dataset ##################
setwd(path)          #############-------> Set Working Directory
file <- 'CCA_Ananlysis.txt'               #############-----> File Name to use in naming Output files
file33 <- 'CountsMatrix1.txt' 			  #############-----> Counts Matrix 1
data <- as.matrix(read.table(file33,sep="\t",header=T,check.names=F,row.names=1))
file1 <- 'MetaFile1.txt'                  #########----> MetaFile [please see github for exmaple files]
col_data <- as.matrix(read.table(file1,sep="\t",header=T,row.names=1,check.names=F))
col_d <- data.frame(col_data)
pd <- new("AnnotatedDataFrame", data = col_d)
rownames(pd) <- pd$Cells
file2 <- 'GeneInfo.txt'                   ##########----> Gene Info [please see github for exmaple files]  
feature_data <- as.matrix(read.table(file2,sep="\t",header=T,row.names=1,check.names=F))
feature_data1 <- data.frame(feature_data)
fd <- new("AnnotatedDataFrame", data = feature_data1)
rownames(fd) <- fd$gene_short_name
HSMM <- ExpressionSet(data, phenoData = pd, featureData = fd)
ctrl <- CreateSeuratObject(raw.data = data, min.cells =3, min.genes = 100, project= 'Sample1')

file3 <- 'CountsMatrix1.txt' 			 #############-----> Counts Matrix 2
data1 <- as.matrix(read.table(file3,sep="\t",header=T,check.names=F,row.names=1))
file4 <- 'MetaFile2.txt'                 #########----> MetaFile [please see github for exmaple files]
col_data1 <- as.matrix(read.table(file4,sep="\t",header=T,row.names=1,check.names=1))
col_d1 <- data.frame(col_data1)
pd2 <- new("AnnotatedDataFrame", data = col_d1)
rownames(pd2) <- pd2$Cells
file5 <- 'GeneInfo.txt'                  ##########----> Gene Info [please see github for exmaple files]
feature_data2 <- as.matrix(read.table(file5,sep="\t",header=T,row.names=1,check.names=F))
feature_data3 <- data.frame(feature_data2)
fd2 <- new("AnnotatedDataFrame", data = feature_data3)
rownames(fd2) <- fd2$gene_short_name
HSMM2 <- ExpressionSet(data1, phenoData = pd2, featureData = fd2)
stim1 <- CreateSeuratObject(raw.data = data1, min.cells =3, min.genes = 100, project= 'Sample2')

file3 <- 'CountsMatrix3.txt' 			 #############-----> Counts Matrix 3
data1 <- as.matrix(read.table(file3,sep="\t",header=T,check.names=F,row.names=1))
file4 <- 'MetaFile3.txt'                 #########----> MetaFile [please see github for exmaple files]
col_data1 <- as.matrix(read.table(file4,sep="\t",header=T,row.names=1,check.names=1))
col_d1 <- data.frame(col_data1)
pd2 <- new("AnnotatedDataFrame", data = col_d1)
rownames(pd2) <- pd2$Cells
file5 <- 'GeneInfo.txt'                  ##########----> Gene Info [please see github for exmaple files]
feature_data2 <- as.matrix(read.table(file5,sep="\t",header=T,row.names=1,check.names=F))
feature_data3 <- data.frame(feature_data2)
fd2 <- new("AnnotatedDataFrame", data = feature_data3)
rownames(fd2) <- fd2$gene_short_name
HSMM2 <- ExpressionSet(data1, phenoData = pd2, featureData = fd2)
stim1 <- CreateSeuratObject(raw.data = data1, min.cells =3, min.genes = 100, project= 'Sample3')
########################################################################################


############### Load Cell Annotations ###################
CellsMeta = ctrl@meta.data
head(CellsMeta)
CellsMeta["Stages"] <- pd$Stages
head(CellsMeta)
CellsMetaTrim <- subset(CellsMeta, select = c("Stages"))
head(CellsMetaTrim)
ctrl <- AddMetaData(ctrl, CellsMetaTrim)
head(ctrl@meta.data)

CellsMeta = ctrl@meta.data
head(CellsMeta)
CellsMeta["LineageAnnotations"] <- pd$LineageAnnotations
head(CellsMeta)
CellsMetaTrim <- subset(CellsMeta, select = c("LineageAnnotations"))
head(CellsMetaTrim)
ctrl <- AddMetaData(ctrl, CellsMetaTrim)
head(ctrl@meta.data)
ctrl@meta.data$stim <- "Sample1"

CellsMeta1 = stim1@meta.data
head(CellsMeta1)
CellsMeta1["Stages"] <- pd2$Stages
head(CellsMeta1)
CellsMetaTrim1 <- subset(CellsMeta1, select = c("Stages"))
head(CellsMetaTrim1)
stim1 <- AddMetaData(stim1, CellsMetaTrim1)
head(stim1@meta.data)

CellsMeta1 = stim1@meta.data
head(CellsMeta1)
CellsMeta1["LineageAnnotations"] <- pd2$LineageAnnotations
head(CellsMeta1)
CellsMetaTrim1 <- subset(CellsMeta1, select = c("LineageAnnotations"))
head(CellsMetaTrim1)
stim1 <- AddMetaData(stim1, CellsMetaTrim1)
head(stim1@meta.data)
stim1@meta.data$stim <- "Sample2"

CellsMeta1 = stim@meta.data
head(CellsMeta1)
CellsMeta1["Stages"] <- pd1$Stages
head(CellsMeta1)
CellsMetaTrim1 <- subset(CellsMeta1, select = c("Stages"))
head(CellsMetaTrim1)
stim <- AddMetaData(stim, CellsMetaTrim1)
head(stim@meta.data)

CellsMeta1 = stim@meta.data
head(CellsMeta1)
CellsMeta1["LineageAnnotations"] <- pd1$LineageAnnotations
head(CellsMeta1)
CellsMetaTrim1 <- subset(CellsMeta1, select = c("LineageAnnotations"))
head(CellsMetaTrim1)
stim <- AddMetaData(stim, CellsMetaTrim1)
head(stim@meta.data)
stim@meta.data$stim <- "Sample4"
###################################################################


############## Filtering, Normalizing and Scaling #########
mito.genes <- grep(pattern = "^mt-", x = rownames(x = ctrl@data), value = TRUE)
percent.mito <- colSums(ctrl@raw.data[mito.genes, ])/colSums(ctrl@raw.data)
file <- strsplit(file,".txt")
png("QCstats.png", width = 8, height = 15, units = 'in', res = 600)
ctrl <- AddMetaData(ctrl, percent.mito, "percent.mito")
VlnPlot(ctrl, c("nGene", "nUMI", "percent.mito"), nCol=3)
dev.off()
png("Gene-plot-1.png", width = 8, height = 15, units = 'in', res = 600)
par(mfrow = c(1,2))
GenePlot(ctrl, "nUMI", "percent.mito")
dev.off()
png("Gene-plot-2.png", width = 8, height = 15, units = 'in', res = 600)
GenePlot(ctrl, "nUMI", "nGene")
dev.off()
ctrl <- FilterCells(object = ctrl, subset.names = c("nGene", "percent.mito"), low.thresholds = c(-Inf, -Inf), high.thresholds = c(7500, 0.15))
ctrl <- NormalizeData(object = ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
ctrl <- ScaleData(object = ctrl, vars.to.regress = c("nUMI", "percent.mito"))

mito.genes <- grep(pattern = "^mt-", x = rownames(x = stim@data), value = TRUE)
percent.mito <- colSums(stim@raw.data[mito.genes, ])/colSums(stim@raw.data)
png("QCstats.png", width = 8, height = 15, units = 'in', res = 600)
stim <- AddMetaData(stim, percent.mito, "percent.mito")
VlnPlot(stim, c("nGene", "nUMI", "percent.mito"), nCol=3)
dev.off()
png("Gene-plot-1.png", width = 8, height = 15, units = 'in', res = 600)
par(mfrow = c(1,2))
GenePlot(stim, "nUMI", "percent.mito")
dev.off()
png("Gene-plot-2.png", width = 8, height = 15, units = 'in', res = 600)
GenePlot(stim, "nUMI", "nGene")
dev.off()
stim <- FilterCells(object = stim, subset.names = c("nGene", "percent.mito"), low.thresholds = c(-Inf, -Inf), high.thresholds = c(7500, 0.15))
stim <- NormalizeData(object = stim, normalization.method = "LogNormalize", scale.factor = 10000)
stim <- ScaleData(object = stim, vars.to.regress = c("nUMI", "percent.mito"))

mito.genes <- grep(pattern = "^mt-", x = rownames(x = stim1@data), value = TRUE)
percent.mito <- colSums(stim1@raw.data[mito.genes, ])/colSums(stim1@raw.data)
png("QCstats.png", width = 8, height = 15, units = 'in', res = 600)
stim1 <- AddMetaData(stim1, percent.mito, "percent.mito")
VlnPlot(stim1, c("nGene", "nUMI", "percent.mito"), nCol=3)
dev.off()
png("Gene-plot-1.png", width = 8, height = 15, units = 'in', res = 600)
par(mfrow = c(1,2))
GenePlot(stim1, "nUMI", "percent.mito")
dev.off()
png("Gene-plot-2.png", width = 8, height = 15, units = 'in', res = 600)
GenePlot(stim1, "nUMI", "nGene")
dev.off()
stim1 <- FilterCells(object = stim1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(-Inf, -Inf), high.thresholds = c(7000, 0.15))
stim1 <- NormalizeData(object = stim1, normalization.method = "LogNormalize", scale.factor = 10000)
stim1 <- ScaleData(object = stim1, vars.to.regress = c("nUMI", "percent.mito"))
##########################################################################


################ Cell Cycle Estimation and Removal #######################
cc.genes <- readLines(con = pathToCellCycleGene) ####Please provide celll cycle genes [See files provided on Github for reference]
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]


ctrl <- CellCycleScoring(object = ctrl, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
head(x = ctrl@meta.data)
ctrl@meta.data$CC.Difference <- ctrl@meta.data$S.Score - ctrl@meta.data$G2M.Score
ctrl <- RunPCA(object = ctrl, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)

stim <- CellCycleScoring(object = stim, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
head(x = stim@meta.data)
stim@meta.data$CC.Difference1 <- stim@meta.data$S.Score - stim@meta.data$G2M.Score
stim <- RunPCA(object = stim, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)

stim1 <- CellCycleScoring(object = stim1, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
head(x = stim1@meta.data)
stim1@meta.data$CC.Difference2 <- stim1@meta.data$S.Score - stim1@meta.data$G2M.Score
stim1 <- RunPCA(object = stim1, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)

ctrl <- ScaleData(object = ctrl, vars.to.regress = "CC.Difference", display.progress = FALSE)
stim <- ScaleData(object = stim, vars.to.regress = "CC.Difference1", display.progress = FALSE)
stim1 <- ScaleData(object = stim1, vars.to.regress = "CC.Difference2", display.progress = FALSE)
#########################################################################

###################### Estimating Variable Genes and Selecting Variable Geneset #################
ctrl <- FindVariableGenes(object = ctrl, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3.5, y.cutoff = 0.5, do.plot=F)
length(ctrl@var.genes)
stim1 <- FindVariableGenes(object = stim1, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3.5, y.cutoff = 0.5, do.plot=F)
length(stim1@var.genes)
stim <- FindVariableGenes(object = stim, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3.5, y.cutoff = 0.5, do.plot=F)
length(stim@var.genes)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
g.3 <- head(rownames(stim1@hvg.info), 1000)
genes.use <- g.2
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))
genes.use <- intersect(genes.use, rownames(stim1@scale.data))
################################################################################################



#################### Define Colors ###############
my_cols_lineages <- c("#8B4513", "#FF4500", "#87CEEB", "#FFD700", "#D8BFD8", "#A9A9A9", "#6A5ACD", "#FFA500", "#6B8E23")
my_cols_stages <- c("#ff9400", "#7cff79", "#00007f", "#0080ff")
####################################################


################ CCA, Clustering, tSNE reduction and Visualization ############
Total.list <- list(ctrl, stim1, stim)
M_Main <- RunMultiCCA(object.list = Total.list, genes.use = genes.use, num.ccs = 30)  #### Change Object Name
PrintDim(M_Main,reduction.type = 'cca')
p1 <- DimPlot(object = M_Main, reduction.use = "cca", group.by = "stim", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = M_Main, features.plot = "CC1", group.by = "stim", do.return = TRUE)
png("CCA_DimPlot_and_VlnPlot.png", width = 20, height = 20, units = 'in', res = 600)
plot_grid(p1, p2)
dev.off()
PrintDim(object = M_Main, reduction.type = "cca", dims.print = 1:2, genes.print = 10)
png("CorrelationStrentgh_CCA_analysis.png", width = 20, height = 20, units = 'in', res = 600)
p3 <- MetageneBicorPlot(M_Main, grouping.var = "stim", dims.eval = 1:30, display.progress = FALSE)
dev.off()
png("Heatmaps_Using_10CC.png", width = 20, height = 20, units = 'in', res = 600)
DimHeatmap(object = M_Main, reduction.type = "cca", cells.use = 500, dim.use = 1:10, do.balanced = TRUE)
dev.off()
M_Main <- AlignSubspace(M_Main, reduction.type = "cca", grouping.var = "stim", dims.align = 1:20)
M_Main <- RunTSNE(M_Main, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
M_Main <- FindClusters(M_Main, reduction.type = "cca.aligned", resolution = 1, dims.use = 1:20)
p1 <- TSNEPlot(M_Main, do.return = T, pt.size = 0.8, group.by = "stim")
p2 <- TSNEPlot(M_Main, do.label = T, do.return = T, pt.size = 0.8)
outfile2 <- paste0("CCA_TotalCells", ".png")
plot_grid(p1, p2)
ggsave(outfile2)
outfile2 <- paste0("CCA_ANALYSIS_TSNE_PLOT_static", ".png")
TSNEPlot(M_Main, do.label = T, pt.size = 0.5)
ggsave(outfile2)
outfile3 <- paste0("CCA_ANALYSIS_TSNE_PLOT_interactive", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile3)

write.table(M_Main@ident,"Samples_Clusters_Association_CCA.txt",sep="\t",col.names= NA,quote=F)

M_Main.markers <- FindAllMarkers(object = M_Main, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(M_Main.markers,file="Markers_All_CCA_Clusters.txt",sep="\t",col.names= NA,quote=F)

outfile2 <- paste0('CCA_Analysis', "_static", "_TSNE_plot_with_stages", ".png")
M_Main <- SetAllIdent(object = M_Main, id = "Stages")
TSNEPlot(object = M_Main, colors.use=my_cols_stages)
ggsave(outfile2)
outfile3 <- paste0('CCA_Analysis', "_static", "_TSNE_plot_with_stages", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile3)

M_Main <- StashIdent(M_Main, save.name = "celltype")
outfile2 <- paste0('CCA_Analysis', "_static", "_TSNE_plot_with_LineageAnnotations", ".png")
M_Main <- SetAllIdent(object = M_Main, id = "LineageAnnotations")
TSNEPlot(object = M_Main, colors.use=my_cols_lineages)
ggsave(outfile2)
outfile3 <- paste0('CCA_Analysis', "_static", "_TSNE_plot_with_LineageAnnotations", ".html")
htmlwidgets::saveWidget(as.widget(ggplotly()), outfile3)

M_Main <- SetAllIdent(object = M_Main, id = "res.1")
#############################################################################################


######################### Save Robject ###############################################################
save(M_Main, file = 'Seurat_v2.3.4_MultiCCA_Analysis_Sample1_Sample2_Smaple3.Robj')
######################################################################################################