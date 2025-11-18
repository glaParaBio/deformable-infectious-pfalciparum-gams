library(dendextend)
library(DescTools)

#Read all PCL files and create lists for time course datasets

young <- read.table("ref_data/microarrays/filteredInput/young.pcl", header=T, sep="\t", quote='')
hu <- read.table("ref_data/microarrays/filteredInput/hu.pcl", header=T, sep="\t", quote='')
leRoch <- read.table("ref_data/microarrays/filteredInput/leRoch.pcl", header=T, sep="\t", quote='')
llinas3D7 <- read.table("ref_data/microarrays/filteredInput/llinas3D7.pcl", header=T, sep="\t", quote='')
llinasDD2 <- read.table("ref_data/microarrays/filteredInput/llinasDD2.pcl", header=T, sep="\t", quote='')
llinasHB3 <- read.table("ref_data/microarrays/filteredInput/llinasHB3.pcl", header=T, sep="\t", quote='')
voss <- read.table("ref_data/microarrays/filteredInput/voss.pcl", header=T, sep="\t", quote="")

timecourse <- list(young, hu, leRoch, llinas3D7, llinasDD2, llinasHB3, voss)

#Read all PCL files and create lists for in vivo datasets

daily <- read.table("ref_data/microarrays/filteredInput/daily.pcl", header=T, sep="\t", quote='')
leRoux <- read.table("ref_data/microarrays/filteredInput/leRoux.pcl", header=T, sep="\t", quote='')
milner <- read.table("ref_data/microarrays/filteredInput/milner.pcl", header=T, sep="\t", quote='')

inVivo <- list(daily, leRoux, milner)

#Preprocessing steps per dataset:
#  * Average rows for duplicated genes (arising as a result of merging of gene models between the annotation originally used for processing and the current annotation)
#  * Remove unnecessary columns
#  * Scale and centre
#  * Impute missing values (KNN, k=10)
#  * Refactor for simple conversion to matrix

#Define function
preProcess <- function (df) {
  colnames(df)[1] <- "X"
  df <- aggregate.data.frame(df, by=list(df$X), FUN=mean)
  row.names(df) <- df$Group.1
  df <- df[, 5:ncol(df)]
  mat <- scale(df)
  mat <- ImputeKnn(mat, k=10, scale=F)
  return (mat)
}

#Process timecourse datasets
timecourse <- lapply(timecourse, preProcess)
names(timecourse) <- c("young", "hu", "leRoch", "llinas3D7", "llinasDD2", "llinasHB3", "voss")

#Process in vivo datasets
inVivo <- lapply(inVivo, preProcess)
names(inVivo) <- c("daily", "leRoux", "milner")

#Calculate similarity
#Per dataset:
#  * Calculate a matrix of pairwise Pearson correlation R values between genes
#  * Carry out Fisher's Z transform on all R values
#  * Make a list of matrices

#Define functions
makeDistMatrix <- function(df) {
  pearCor <- cor(t(df), method="pearson")
  pearCor <-FisherZ(pearCor)
  return(pearCor)
}

makeDistMatrixList <- function(datasetList) {
  distList <- list()
  for (i in 1: length(datasetList)) {
    df <- datasetList[[i]]
    pearCor <- makeDistMatrix(df)
    distList[[i]] <- pearCor
  }
  return (distList)
}

#Process timecourse datasets
timecourseDist <- makeDistMatrixList(timecourse)
names(timecourseDist) <- c("young", "hu", "leRoch", "llinas3D7", "llinasDD2", "llinasHB3", "voss")

#Process in vivo datasets
inVivoDist <- makeDistMatrixList(inVivo)
names(inVivoDist) <- c("daily", "leRoux", "milner")


#Average datasets

#Average all timecourse datasets
all <- rbind(as.data.frame.table(timecourseDist[[1]]), as.data.frame.table(timecourseDist[[2]]), as.data.frame.table(timecourseDist[[3]]), as.data.frame.table(timecourseDist[[4]]), as.data.frame.table(timecourseDist[[5]]), as.data.frame.table(timecourseDist[[6]]), as.data.frame.table(timecourseDist[[7]]))
timeCourseMergedWithVoss <- as.matrix(tapply(all[[3]], all[-3], mean, default=NA))

#Average all in vivo datasets
all <- rbind(as.data.frame.table(inVivoDist[[1]]), as.data.frame.table(inVivoDist[[2]]), as.data.frame.table(inVivoDist[[3]]))
inVivoMerged <- as.matrix(tapply(all[[3]], all[-3], mean, default=NA))

#Average merged sets together
all <- rbind(as.data.frame.table(timeCourseMergedWithVoss), as.data.frame.table(inVivoMerged))
finalMergedWithVoss <- as.matrix(tapply(all[[3]], all[-3], mean, default=NA))


#Impute missing values in distance matrices
#Missing values arise when gene pairs don't occur in the same dataset.  E.g., if gene 1 only occurs in dataset 1 1 but not 2, and gene 2 occurs in dataset 2 but not 1, there can be no pairwise comparison.

#These values are imputed by:
#  * Finding the 10 genes with the strongest connections to gene 1 and extracting their connections to gene 2
#  * Finding the 10 genes with the strongest connections to gene 2 and extracting their connections to gene 1
#  * Averaging these extracted weights


#Define function
matrixImpute <- function(matrix) {
  naVals <- which (is.na(matrix), arr.ind=T)
  
  imputed <- matrix
  
  for (i in 1:nrow(naVals)) {
    naVal <- naVals[i,]
    
    top10G1 <- names(sort(matrix[naVal["Var1"],], decreasing=T))[2:11] #1st item is always same gene so exclude this
    top10G2 <- names(sort(matrix[naVal["Var2"],], decreasing=T))[2:11]
    
    weights <- c(matrix[top10G1, naVal["Var2"]], matrix[top10G2, naVal["Var1"]])
    
    imputed[naVal["Var1"], naVal["Var2"]] <- mean(weights, na.rm=T)
  }
  return(imputed)
}

#Carry out imputation
finalMergedWithVoss.imputed <- matrixImpute(finalMergedWithVoss)


#Hierarchical clustering
hrWithVoss <- hclust(as.dist(1 - finalMergedWithVoss.imputed))
hrWithVoss$height <- hrWithVoss$height+4

hrdWithVoss <- as.dendrogram(hrWithVoss, hang=0.01) %>%  set("labels_to_character") %>% color_branches(h=quantile(hrWithVoss$height, .95))

#Cut tree to make clusters
clusterWithVoss <- cutree(hrdWithVoss, h=quantile(hrWithVoss$height, 0.95))

genes <- as.data.frame(names(clusterWithVoss))
colnames(genes) <- "GeneIds"

hrdWithVoss <- set_labels(hrdWithVoss, rep(NA, nleaves(hrdWithVoss)))
pdf('output/geneDendrogram.pdf', width=18/2.54, height=12/2.54)
plot(hrdWithVoss)
dev.off()

getOldIds <- function(x) {
  oldIdList <- tolower(strsplit(strsplit(x, ': ')[[1]][2], ';')[[1]])
  return (oldIdList)
}

invertList <- function (i, x, n) {
  lapply(x[[i]], function(y, n){newList <- n; names(newList) <- y; return(newList)}, n=n[[i]])
}

idMappings <- read.table("ref_data/microarrays/auxiliary_files/allGenesOldIds.txt", header=T, sep="\t", quote="")
oldIdList <- lapply(as.character(idMappings$Previous.ID.s.), getOldIds)
names(oldIdList) <- idMappings$Gene.ID

oldIdList <- unlist(lapply(seq_along(oldIdList), invertList, x=oldIdList, n=names(oldIdList)))

#Gene Products

allGeneProducts <- read.table("ref_data/microarrays/auxiliary_files/GeneProducts.txt", sep="\t", header=T, quote="")

geneData <- merge(genes, allGeneProducts, by.x="GeneIds", by.y="X.Gene.ID.", no.dups=F)
geneData <- as.data.frame(cbind(as.character(geneData$GeneIds), as.character(geneData$X.Product.Description.)))

geneData$Cluster <- unlist(lapply(as.character(geneData$V1), function(x, clusters) {return (clusters[x])}, clusters=clusterWithVoss))
geneData <- unique(geneData)
write.table(geneData, "ref_data/geneData.tsv", sep='\t', quote=F, row.names=F)

