library(data.table)
library(DescTools)
library(qgraph)

preProcess <- function (df) {
    df <- aggregate.data.frame(df, by=list(rownames(df)), FUN=mean)
    row.names(df) <- df$Group.1
    df$Group.1 <- NULL
    mat <- scale(df)
    nna <- sum(is.na(mat))
    N <- nrow(mat) * ncol(mat)
    if(any(is.na(mat)) == TRUE) {
        mat <- ImputeKnn(mat, k=10, scale=F)
    }
    return(mat)
}

makeQgraph <- function(mat, annotation, annotate='Cluster', min_cor=0.75) {
    stopifnot(c(annotate, 'Classification') %in% names(annotation))
    cormat <- cor(mat, use='pairwise.complete.obs')

    is_single <- apply(cormat, 1, function(x) sort(abs(x), decreasing=TRUE)[2] < min_cor)
    cormat <- cormat[!is_single, !is_single]

    keep_ann <- annotation[match(colnames(cormat), annotation[[annotate]])]
    stopifnot(identical(keep_ann[[annotate]], colnames(cormat)))
    groups <- list()
    for(x in unique(keep_ann$Classification)) {
        grp <- keep_ann[Classification == x][[annotate]]
        idx <- which(colnames(cormat) %in% grp)
        groups[[x]] <- idx
    }
    Q <- qgraph(cormat, groups=groups, minimum=min_cor, layout='spring', posCol='blue', negCol='red')
    return(Q)
}

pclDir <- snakemake@input[['pclDir']] 
annotation <- fread(snakemake@input[['annotation']])
annotation[, Cluster := as.character(Cluster)]
annotation <- unique(annotation[, list(Cluster, Classification)])

clstData <- fread(snakemake@input[['clstData']], )
setnames(clstData, 'V1', 'gene_id')
clstData <- clstData[, list(gene_id, Cluster=as.character(Cluster))]

annotation <- merge(annotation, clstData, by='Cluster')

pclFiles <- c('daily.pcl.gz',
    'hu.pcl.gz',
    'leRoch.pcl.gz',
    'leRoux.pcl.gz',
    'llinas3D7.pcl.gz',
    'llinasDD2.pcl.gz',
    'llinasHB3.pcl.gz',
    'milner.pcl.gz',
    'voss.pcl.gz',
    'young.pcl.gz')


pcl <- list()
for(fin in pclFiles) {
    dat <- fread(file.path(pclDir, fin))
    dat[, name := NULL]
    dat[, GWeight := NULL]
    setnames(dat, 1, 'gene_id')
    dat <- as.matrix(dat, rownames='gene_id')
    dat <- preProcess(dat)
    dat <- as.data.table(dat, keep.rownames='gene_id')
    dat <- melt(data=dat, id.vars= 'gene_id', value.name= 'gex', variable.name= 'array_id')
    dat[, dataset := sub('\\.pcl$', '', sub('\\.gz$', '', fin))]
    pcl[[length(pcl) + 1]] <- dat
}
pcl <- rbindlist(pcl)
pcl[, array_id := sprintf('%s:%s', dataset, array_id)]

clst_pcl <- merge(annotation, pcl, by='gene_id')

clst_size <- annotation[, .N, Cluster]
clst_avg <- clst_pcl[Cluster %in% clst_size[N > 5]$Cluster, list(gex=mean(gex, na.rm=TRUE), n_genes=.N), by=list(array_id, Cluster)]

mat <- as.matrix(dcast(clst_avg, array_id ~ Cluster, value.var= 'gex'), rownames='array_id')
pdf(snakemake@output[['network_clusters']], width=18/2.54, height=14/2.54)
Q <- makeQgraph(mat, annotation, annotate='Cluster')
title('Correlation network between cluster profiles')
dev.off()

keep_genes <- annotation[!is.na(Classification)]$gene_id
mat <- as.matrix(dcast(clst_pcl[keep_genes], array_id ~ gene_id, value.var= 'gex'), rownames='array_id')
pdf(snakemake@output[['network_genes']], width=18/2.54, height=14/2.54)
Q <- makeQgraph(mat, annotation, annotate='gene_id', min_cor=0.95)
title('Correlation network between gene profiles')
dev.off()
