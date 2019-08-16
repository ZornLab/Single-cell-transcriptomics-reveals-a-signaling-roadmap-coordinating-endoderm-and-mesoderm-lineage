#################### Single Cell Voting using KNN classification algorithm
#################### Please see Github for sample files
#################### Files required : TrainSet.txt and TestSet.txt
#################### Please use >= R/3.4.4
#################### Some Issues that can occur : problems with installation of KODAMA and knnflex packages and if using > R/3.5.0 then please follow Part2 of the script



############## Part1 : use when using R/3.4.4 and packages [KODAMA, knnflex] are installed

############## Load R Packages ###########
library(class)
library(KODAMA)
library(knnflex)
##########################################


############ Read Train/ Test datasets and Run KNN ##################
setwd(path)

file1 <- TrainSet.txt
data1 <- read.table(file1,sep="\t",header=T)
file2 <- TestSet.txt
data2 <- read.table(file2,sep="\t",header=T)
train_data <- data.frame(data1)
test_data <- data.frame(data2)
train_data_num <- train_data[,2:1199]     #########---->1199 here is number of Markers/Features
test_data_num <- test_data[,2:1199]		  #########---->1199 here is number of Markers/Features

cl <- factor(train_data$Cluster)
x <- rbind(train_data_num,test_data_num)
kdist <- knn.dist(x)
nn3 <- knn(train_data_num,test_data_num,cl,prob=TRUE,k=3)        ###########---> k=3 means 3-nearest neighbors [Please adjust accordingly]

## The resulting confusion matrix
conf <- as.matrix(table(test_data[,'Cluster'],nn3))
write.table(conf,"ConfusionMatrix.txt",sep="\t",quote=F)

probs <- knn.probability(1:803, 804:1852, cl, kdist, k=3, ties.meth="min")  ##########-----> 803 is number of cells in Train set and 1852 is total number of cells including Train and Test
Probablity <- head(t(probs),1049)`                                          ##########-----> 1049 is number of cells in Test set
Probablity <- cbind(data2$Cells, data2$Cluster, Probablity)
write.table(Probablity,"ProbabilityMatrix.txt",sep="\t",quote=F)

#######################################################################

#################### End Part1 #################################################################


############## Part2 : use when using > R/3.5.0 and have problems installing KODAMA and knnflex

############## Load R Packages ###########
library(class)
##########################################


############ Read Train/ Test datasets and Run KNN ##################
setwd(path)

file1 <- TrainSet.txt
data1 <- read.table(file1,sep="\t",header=T)
file2 <- TestSet.txt
data2 <- read.table(file2,sep="\t",header=T)
train_data <- data.frame(data1)
test_data <- data.frame(data2)
train_data_num <- train_data[,2:1199]     #########---->1199 here is number of Markers/Features
test_data_num <- test_data[,2:1199]		  #########---->1199 here is number of Markers/Features

cl <- factor(train_data$Cluster)
x <- rbind(train_data_num,test_data_num)
kdist <- knn.dist(x)
nn3 <- knn(train_data_num,test_data_num,cl,prob=TRUE,k=3)      ###########---> k=3 means 3-nearest neighbors [Please adjust accordingly]

## The resulting confusion matrix
conf <- as.matrix(table(test_data[,'Cluster'],nn3))
write.table(conf,"ConfusionMatrix.txt",sep="\t",quote=F)
#####################################################################
######################## End Part2 #############################################################