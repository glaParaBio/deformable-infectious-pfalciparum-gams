
# Following protocol from Pelle *et al* to filter microarray data prior to
# analysis with Sleipnir.  This includes:
# 1. Remove paralogs of *var*, *rifin* and *stevor* gene families
# 2. Remove any gene that does not have data in more than half of the *in
# vitro* timecourse studies or more than half of the *in vivo* studies
#
# Paralogs of variable gene families were downloaded from PlasmoDB.  Searches
# were carried out for genes annotated as *var*, *rifin* or *stevor* in the
# Gene Product Description. An orthology transform was used to find un- or
# mis-annotated paralogs of these genes.  The results of the transform were
# unioned with the results of the original text search.  As the microarray
# datasets all use old gene ids, previous ids were included in the downloads.
#
# In addition, the Voss dataset uses current gene ids while most of the older
# datasets use old gene ids. To build a network incorporating the Voss data,
# these must all be the same.  We will therefore need to convert all gene ids
# to current ids.

#Obtain Lists of *var*, *rifin* and *stevor* genes

getOldIds <- function(x) {
  oldIdList <- tolower(strsplit(strsplit(x, ': ')[[1]][2], ';')[[1]])
  return (oldIdList)
}

varGenes <- read.table("ref_data/microarrays/rawInput/varWithOldIds.txt", header=T, sep="\t", quote="")
oldVarIds <- unlist(lapply(as.character(varGenes$Previous.ID.s.), getOldIds))

rifinGenes <- read.table("ref_data/microarrays/rawInput/rifinWithOldIds.txt", header=T, sep="\t", quote="")
oldVarIds <- c(oldVarIds, unlist(lapply(as.character(rifinGenes$Previous.ID.s.), getOldIds)))

stevorGenes <- read.table("ref_data/microarrays/rawInput/stevorWithOldIds.txt", header=T, sep="\t", quote="")
oldVarIds <- c(oldVarIds, unlist(lapply(as.character(stevorGenes$Previous.ID.s.), getOldIds)))

oldVarIds <- unique(oldVarIds)
oldVarIds <- oldVarIds[!is.na(oldVarIds)]

"%ni%" <- Negate("%in%")

datasetFilter <- function(dataset, oldIds) {
  filtered <- dataset[which(tolower(dataset[,1]) %ni% oldIds),]
  return(filtered)
}

getGeneList <- function(dataset) {
  return(dataset[,1])
}

filterDataByCount <- function(dataset, keep) {
  filtered <- dataset[which(dataset[,1] %in% keep),]
  return(filtered)
}

getKeepList <- function(datasetList, x) {
  geneList <- unlist(lapply(datasetList, getGeneList), use.names=F)
  counts <- tabulate(match(geneList, unique(geneList)))
  names(counts) <- unique(geneList)
  keep <- names(counts)[which(counts >= x)]
  return(keep)
}

idConvert <- function(dataset, oldIdList) {
  dataset[,1] <- unlist(lapply(as.character(dataset[,1]), function(x, oldIdList) {x <- tolower(x); geneList <- unique(oldIdList[which(names(oldIdList) == x)]); geneList <- (paste(geneList,collapse="; "))}, oldIdList=oldIdList))
  dataset <- dataset[which(dataset[,1] != ""),]
  return (dataset)
}

invertList <- function (i, x, n) {
  lapply(x[[i]], function(y, n){newList <- n; names(newList) <- y; return(newList)}, n=n[[i]])
}

idMappings <- read.table("ref_data/microarrays/rawInput/allGenesOldIds.txt", header=T, sep="\t", quote="")
oldIdList <- lapply(as.character(idMappings$Previous.ID.s.), getOldIds)
names(oldIdList) <- idMappings$Gene.ID

oldIdList <- unlist(lapply(seq_along(oldIdList), invertList, x=oldIdList, n=names(oldIdList)))

inVivo <- list(daily=read.table("ref_data/microarrays/rawInput/Daily.pcl.gz", sep="\t", quote="", header=T),
               leRoux=read.table("ref_data/microarrays/rawInput/LeRoux.pcl.gz", sep="\t", quote="", header=T), 
               milner=read.table("ref_data/microarrays/rawInput/Malawi.pcl.gz", sep="\t", quote="", header=T))

inVivo <- lapply(inVivo, datasetFilter, oldIds=oldVarIds)

keep <- getKeepList(inVivo, 2)

inVivo <- lapply(inVivo, filterDataByCount, keep=keep)

inVivo <- lapply(inVivo, idConvert, oldIdList=oldIdList)

timecourse <- list(llinas3D7=read.table("ref_data/microarrays/rawInput/3D7_Llinas.pcl.gz", sep="\t", quote="", header=T),
                   llinasDD2=read.table("ref_data/microarrays/rawInput/DD2_Llinas.pcl.gz", sep="\t", quote="", header=T),
                   llinasHB3=read.table("ref_data/microarrays/rawInput/HB3_Llinas.pcl.gz", sep="\t", quote="", header=T),
                   hu=read.table("ref_data/microarrays/rawInput/Hu.pcl.gz", sep="\t", quote="", header=T),
                   leRoch=read.table("ref_data/microarrays/rawInput/LeRoch.pcl.gz", sep="\t", quote="", header=T),
                   young=read.table("ref_data/microarrays/rawInput/young2005.pcl.gz", sep="\t", quote="", header=T))

timecourse <- lapply(timecourse, datasetFilter, oldIds=oldVarIds)

keep <- getKeepList(timecourse, 3)

timecourse <- lapply(timecourse, filterDataByCount, keep=keep)

timecourse <- lapply(timecourse, idConvert, oldIdList=oldIdList)

voss <- read.table("ref_data/microarrays/rawInput/vossMicroarray.txt.gz", sep="\t", quote="", header=T)
voss <- cbind(voss[,c(1,3:13)])

voss <- voss[which(voss$geneid %ni% varGenes$Gene.ID),]
voss <- voss[which(voss$geneid %ni% rifinGenes$Gene.ID),]
voss <- voss[which(voss$geneid %ni% stevorGenes$Gene.ID),]

keep <- getKeepList(timecourse, 3)
voss <- filterDataByCount(voss, keep=keep)

dir.create('ref_data/microarrays/filteredInput', recursive=TRUE, showWarnings=FALSE)
saveTable <- function(i, datasets, names) {
  dataset <- datasets[[i]]
  name = names(datasets)[i]
  cols <- colnames(dataset)
  dataset <-cbind(dataset[,1], NA, rep(1, nrow(dataset)), dataset[,2:ncol(dataset)])
  colnames(dataset) <- c(cols[1], "name", "GWeight", cols[2:length(cols)])
  write.table(dataset, paste(c("ref_data/microarrays/filteredInput/", name, ".pcl"), collapse=""), sep="\t", quote=F, col.names=T, row.names=F, na="")
}

invisible(lapply(seq_along(inVivo), saveTable, datasets=inVivo, names=names(inVivo)))
invisible(lapply(seq_along(timecourse), saveTable, datasets=timecourse, names=names(timecourse)))

cols <- colnames(voss)
voss <- cbind(voss[,1], NA, rep(1, nrow(voss)), voss[,2:ncol(voss)])
colnames(voss) <- c(cols[1], "name", "GWeight", cols[2:length(cols)])
write.table(voss, "ref_data/microarrays/filteredInput/voss.pcl", sep="\t", quote=F, col.names=T, row.names=F, na="")
