################################# Please Download SPRING [v1.0] from (https://github.com/AllonKleinLab/SPRING)
################################# This script is a modified version of the regular SPRING processing script
################################# Script is designed to learn Principle Component Space from most complex dataset [in terms of lineage diversification] and then transform the whole dataset using learnt PC space to enhance lineage diversification through time
################################# Please see SPRING github page on information on file formats and functions
################################# Please use python v2.7


print 'Loading Libraries'
import pickle, numpy as np
from preprocessing_python import *
from numpy import genfromtxt
import collections
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA

print 'Loading expression matrix'
E = genfromtxt('CompeleteDataset.txt', delimiter='\t')  ###### Data Matrix with cells as rows and genes as columns [Firstly we have 2994 E9.5 cells then 1550 E8.5 cells and lastly 941 E8.5 cells] --> Complete data
E1 = genfromtxt('ComplexDataset.txt', delimiter='\t')            ###### Data Matrix with cells as rows and genes as columns [we have 2994 E9.5 cells] ---> part of the data

print 'Filtering cells'
E,cell_filter = filter_cells(E,1000)								 ###### Filtering cells on the complete data


print 'Row-normalizing'
E = row_normalize(E)												 ###### Normalizing cells of the complete data


print 'Filtering genes'
_,gene_filter = filter_genes(E,0.1,3)								 ###### Filtering Genes from the complete data

print 'Filtering Data'
Filtered_data=Zscore(E[:,gene_filter])

print 'Creating data to learn PCs'
data_to_learnPC=Filtered_data[:E1.shape[0],:]  						###### Creating a learning set data [E9.5 2994 cells]

print 'Learning PC space from complex dataset and using PC space to transform whole dataset'
pca = PCA(n_components=40)
pca.fit(data_to_learnPC) 											###### Learnt PC space from E9.5 [2994 cells]
Etrans = pca.transform(Filtered_data)                               ###### Transforming Complete data using dimensions learnt from E9.5 [2994 cells]


print 'Getting distance matrix'
D = get_distance_matrix(Etrans)										###### Calculating Distance Matrix


print 'Reading Cell Groupings'
cell_groupings = defaultdict(list)
with open("CellGroupings_CompleteDataset.txt") as f:
    for line in f:
       (key, val) = line.split()
       cell_groupings[key].append(val)

print 'Reading Gene List'
with open('Geneslist.txt') as f:
    gene_list = f.read().splitlines()
    
print 'Reading Custom Colors'
custom_colors = defaultdict(list)
with open("GeneExpressionAcrossCells_ForCustomColors.txt") as f:
    for line in f:
       (key, val) = line.split()
       custom_colors[key].append(float(val))
    
print 'Saving SPRING plot'
save_spring_dir(E,D,5,gene_list,'datasets/SPRINGDIR_DATASET', cell_groupings=cell_groupings, custom_colors=custom_colors, coarse_grain_X=1)  ######### SAVING SPRING PLOT
