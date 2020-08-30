# Somatic Mutation Clustering

*This repository is a reposting  of my answer on* [*Biostar*](https://www.biostars.org/p/288370/#444130).

The research question was: How would TCGA bladder cancer samples cluster together regarding to have somatic mutations?

Somatic mutations are said to be spare and heterogeneous. So using them for clustering is not going to be straightforward task. Before diving to clustering methods there are suggestions on how you may go to de-sparsify your data. For instance, knowledge on gene-gene network are usually considered for data de-sparsification . A detailed discussion could be find [here](https://www.nature.com/articles/s41416-018-0109-7). I am not-covering desparsification methods in this thread. A would strongly recommend if you are going to cluster samples based on mutation in a large list of genes do not use my approach , instead consult with packages like [SAMBAR](https://github.com/mararie/SAMBAR) or any other package that deals with data sparsity.

It is tricky to cluster categorical data, because it could lead to non-sense and wrong conclusion! In the case of somatc mutation in contrast to rotine clustering, your matrix here is not numerical. So you DONT allow to use common algorithms like k-mean clustering. If you apply , it wont complain about your data type and provide you result! Mutational matrix is binary or frequency count. These are steps needed for clustering:

## 1.Dissimilarity matrix calculation

In the first step you should calculate a dissimilarity matrix for clustering . Again there is difficulty regarding to math calculation on categorical/binary data. 
In a dissimilarity matrix distances between all possible pair would be calculated and then subsequent algorithm will use these stats to make clusters.
For doing this there are diffrent methods: Gower distance from ```cluster``` R base package, methods available in ```vegan``` R package a: ```binomial```, ```chi-square```,  ```raup``` and ```jaccard``` . It depends on your data and your decision to chose which method. 

## 2.Choosing Clustering algorithms

Choosing the clustering algorithm is the next step. For categorical data you would go for hierarchical clustering (either agglomerative or divisive approach). 

## 3.Assessing the clustering result
The final steps would be assessing the clustering result. We make clusters to see how each data point come together in one cluster and how clusters are distinct from each other. Personally, I preferer to inspect clusters by looking at the cluster dendrograms . However there are some methods helps to inspect clustering experiment. Two common approach are *Elbow method* and *Silhouette method*. The Elbow method provides details on the clustering compactness and the Silhouette mainly concerns about cluster separation.  

Below is my script I am using to make clusters out of somatic mutation in 10 epi-genetic related genes believed to be driver for bladder cancer (10.1038/s41588-019-0572-y). 

```R
# loading packages
require(maftools)
require(vegan)
require(cluster)
require(fpc)
require("pheatmap")

# preparing input matrix for clustering
## downloading and merging data from TCGA into one file
library(TCGAbiolinks)
muse <- GDCquery_Maf("BLCA", pipelines = "muse")
mutect  <- GDCquery_Maf("BLCA", pipelines = "mutect")
somaticsniper  <- GDCquery_Maf("BLCA", pipelines = "somaticsniper")
varscan2 <- GDCquery_Maf("BLCA", pipelines = "varscan2")
all_maf <- rbind(muse,mutect,somaticsniper,varscan2) # this file contains all known variant for TCGA-BLCA

# Filteration
cols <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
all_maf$unique <-apply( all_maf[ , cols ] , 1 , paste , collapse = "-" )
#Retain variant with unique record in unique
all_maf <- all_maf[!duplicated(all_maf$unique), ]
# retaing epigenetic related genes 
epigen <- data.frame(gene = c("KMT2D", "KDM6A", "ARID1A", "EP300", "CREBBP", "KMT2A", "NCOR1", "ARID2", "BAP1", "DHX30"))
all_maf <- all_maf[which (all_maf$Hugo_Symbol %in% epigen$gene), ]
# Filtering out variants
all_maf <- all_maf[all_maf$MC3_Overlap == "TRUE", ] ## select variants overlapped with pan_cancer (MC3) vaiants
##Remove silent and IMPACT$LOW and IMPACT$MODIFIER
all_maf <- all_maf[-which (all_maf$IMPACT == "LOW"), ]
all_maf <- all_maf[-which (all_maf$IMPACT == "MODIFIER"), ]

# convert plain maf file to count matrix
maf = read.maf(maf = all_maf)
count_mtx <- mutCountMatrix(maf,includeSyn = FALSE,countOnly = NULL,removeNonMutated = FALSE)
count_mtx <- t(count_mtx)
#make a binary matrix with 0 and 1. One hot coding is also possibe here, howver final result would not change dramatically
cout_mtx.b <- as.data.frame(apply(count_mtx, 2, function(x) ifelse(x > 0, 1, x)))
# Clustering steps:
#----- Dissimilarity Matrix -----#
dis.gower <- daisy(cout_mtx.b, metric = c("gower")) # this shows you a warning message about your data type that is not factors. Its OK.
#---AGGLOMERATIVE CLUSTERING ----#
ag.clust.com <- hclust(dis.gower, method = "complete") # there are different method for HC here we use coplete. 
plot(ag.clust.com, cex = 0.7, main = "complete linkages agglomerative HC") # looking at dendrogram 4 clusters is imaginable
```
![alt text]https://github.com/hamid-gen/mutation_clustering/blob/master/dendrogram.PNG
To show clusters on denderogram.
```R
rect.hclust(ag.clust.com, k = 4, border = 2:5)
```
![alt text]https://github.com/hamid-gen/mutation_clustering/blob/master/dendrogram_2.PNG

Cluster stats are comming out as list which is not convenient  to read. I borrowed the *cstats.table* function from @Anastasia Reusova from [Towards Data Science](https://towardsdatascience.com/hierarchical-clustering-on-categorical-data-in-r-a27e578f2995?source=user_profile---------2-----------------------) to convet stats into a data frame.
```R
#----ASSESSING THE CLUSTERS-----#
cstats.table <- function(dist, tree, k) {
  clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between",
                    "wb.ratio","dunn2","avg.silwidth")
  clust.size <- c("cluster.size")
  stats.names <- c()
  row.clust <- c()
  output.stats <- matrix(ncol = k, nrow = length(clust.assess))
  cluster.sizes <- matrix(ncol = k, nrow = k)
  for(i in c(1:k)){
    row.clust[i] <- paste("Cluster-", i, " size")
  }
  for(i in c(2:k)){
    stats.names[i] <- paste("Test", i-1)
    
    for(j in seq_along(clust.assess)){
      output.stats[j, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.assess])[j]
      
    }
    for(d in 1:k) {
      cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
      dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
      cluster.sizes[d, i]
      
    }
  }
  output.stats.df <- data.frame(output.stats)
  cluster.sizes <- data.frame(cluster.sizes)
  cluster.sizes[is.na(cluster.sizes)] <- 0
  rows.all <- c(clust.assess, row.clust)
  # rownames(output.stats.df) <- clust.assess
  output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
  colnames(output) <- stats.names[2:k]
  rownames(output) <- rows.all
  is.num <- sapply(output, is.numeric)
  output[is.num] <- lapply(output[is.num], round, 2)
  output}

#
stats.aggl.clustering <-cstats.table(dis.gower, ag.clust.com, 8) 
stats.aggl.clustering # in our case we perefred to select 4 clusters with Silhouette coefficient as 0.22. However there still ks with higher Silhouette coefficient. 

## To find what samples belong to which cluster
sub_gp <- data.frame(clust = cutree(ag.clust.com, k = 4)) #this return a table with sample name in row name and a column with clust.
table(sub_gp$clust) # how many sample in different clusters

## heatmap Visualization

df <- cout_mtx.b
row.names(df) <- NULL
pheatmap(df, cutree_rows = 4, border_color = NA, cellwidth = 10)
```
![alt text]https://github.com/hamid-gen/mutation_clustering/blob/master/heatmap.PNG
