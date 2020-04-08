load('data/Reduced_IntegrationEndoderm.Robj')
load('data/Reduced_IntegrationMesoderm.Robj')
load('data/Integration_Total_Reduced_03-30-2020.Robj')
Genes <- as.matrix(read.table('data/GenesSingleCell.txt', sep="\n", header=F))
Genes <- data.frame(Genes)
genelist <- Genes$V1