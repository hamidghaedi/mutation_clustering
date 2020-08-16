# Somatic Mutation Clustering
Somatic Mutation data Clustering
I am dealing with similar task. Here are my findings that may be of help for somebody else.

Somatic mutations are said to be spare and heterogeneous. So using them for clustering is not going to be straightforward task. Before jumping to clustering methods there are suggestions on how you may go to de-sparsify your data. For instance, knowledge on gene-gene network are usually considered for data de-sparsification . A detailed discussion could be find here. I am not-covering desparsification methods in this answer.
It is tricky to cluster categorical data, because it could lead to non-sense and wrong conclusion!
In contrast to classic clustering, your matrix here is not numerical. So you DONT allow to use common algorithms like k-mean clustering. If you apply , it wont complain about your data type and provide you result!
    Mutational matrix is binary or categorical. These are steps needed for clustering: 
Dissimilarity matrix calculation

In the first step you should calculate a dissimilarity matrix for clustering . Again there is difficulty regarding to math calculation on categorical/binary data. To do this, you would go for something called Gower distance. this method is available in cluster R base package. Also there are methods available in vegan R package appropriate to be applied on binary data: binomial, raup and jaccard . It depends on your data and your decision to chose what method. 

Choosing Clustering algorithms

Choosing the clustering algorithm is the next step. For categorical data you would go for hierarchical clustering (either agglomerative or divisive approach). The final steps would be assessing the clustering result. Below I am providing what I used for my case in short.

You did not provide details on your input data, so exact code is not possible to post here. But the following are the general steps you can follow to cluster your samples.

1- Making mutation count/binary matrix: In my case, I am dealing with TCGA data, and so there are maf files and could be converted to the matrix by maftools package by mutCountMatrix function.This will provide a count matrix. you may need to convert it to binary (0,1) code. 

<library(maftools)
mtx <- mutCountMatrix(maf, includeSyn = FALSE, countOnly = NULL, removeNonMutated = FALSE) #maf file contains mutation infor
#transpose mtx to have genes in columns and samples in row
mtx <- t(mtx)
#Convert counts to binary
mtx.b <- apply(mtx, 2, function(x) ifelse(x > 0, 1, x)) # So 0 = no, 1 =yes>
