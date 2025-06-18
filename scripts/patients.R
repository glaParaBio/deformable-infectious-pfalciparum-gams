library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(ggplot2)
library(ggbeeswarm)

oldToNewGeneIDs <- function(allGenesOldIdsFile, pcl) {
    ids <- fread(allGenesOldIdsFile, header=TRUE)
    ids[, previous_id := sub('Previous IDs: ', '', `Previous ID(s)`)]
    previous_id <- strsplit(ids$previous_id, ';')
    n_ids <- sapply(previous_id, length)
    gene_id <- rep(ids$`Gene ID`, n_ids)
    old2new <- unique(data.table(previous_id=tolower(unlist(previous_id)), gene_id))
    
    gene_col <- names(pcl)[1]
    setnames(pcl, gene_col, 'previous_id')
    pcl <- melt(data=pcl, id.vars= 'previous_id', value.name= 'gex', variable.name= 'array_id')
    pcl[, previous_id := tolower(previous_id)]

    pcl <- merge(pcl, old2new, by='previous_id', all.x=TRUE)
    # Average expression when two old IDs map to the same new ID.  There are also
    # cases where the old ID maps to more than 1 new IDs. The merge will
    # effectively duplicate these genes. This is not ideal but it's ok.
    pcl <- pcl[, list(gex=mean(gex)), list(array_id, gene_id)]
    pcl <- dcast(pcl, gene_id ~ array_id, value.var='gex')
    setnames(pcl, 'gene_id', gene_col)
    return(pcl)
}

read_pcl <- function(pcl, allGenesOldIdsFile) {
    dataset <- sub('\\.pcl$', '', sub('\\.gz$', '', basename(pcl)))
    dat <- fread(pcl)
    if(dataset == 'Joiceetaldata') {
        dat <- oldToNewGeneIDs(allGenesOldIdsFile, dat)
    }

    if('name' %in% names(dat)) {dat[, name := NULL]}
    if('GWeight' %in% names(dat)) {dat[, GWeight := NULL]}

    setnames(dat, 1, 'gene_id')
    dat <- melt(data=dat, id.vars= 'gene_id', value.name= 'gex', variable.name= 'array_id')
    return(dat)
}

options(warn=2)

peakAll <- fread(snakemake@input[['peakAll']]) # archive/allDatasets/tablesWithVoss/peaksofMeanAll.tsv
peakAll <- peakAll[MinMax == 'Max', list(Cluster, Llinas.Timepoint, Llinas.Value, Voss.Timepoint, Voss.Value)]
peakAll[, Cluster := as.character(Cluster)]
peakAll[, Voss.Timepoint := as.numeric(sub('day', '', Voss.Timepoint))]
peakAll[, Llinas.Timepoint := as.numeric(sub('Timepoint\\.', '', Llinas.Timepoint))]
# xord <- unique(peakAll[order(vtp), list(Voss.Timepoint, vtp)])
# peakAll[, Voss.Timepoint := factor(Voss.Timepoint, xord$Voss.Timepoint)]

allGenesOldIdsFile <- snakemake@input[['allGenesOldIds']] # archive/inputFiles/allGenesOldIds.txt 

# Patient dataset
pcl <- read_pcl(snakemake@input[['pcl']], allGenesOldIdsFile) # 'archive/inputFiles/filteredInput/daily.pcl.gz'

annotation <- fread(snakemake@input[['annotation']]) # output/cluster_annotation/cluster_enrichment.tsv
annotation[, Cluster := as.character(Cluster)]
annotation <- unique(annotation[, list(Cluster, Classification)])

clstData <- fread(snakemake@input[['clstData']]) # archive/allDatasets/tablesWithVoss/geneData.tsv
setnames(clstData, 'V1', 'gene_id')
clstData <- clstData[, list(gene_id, Cluster=as.character(Cluster))]

annotation <- merge(annotation, clstData, by='Cluster')

clst_pcl <- merge(annotation, pcl, by='gene_id')
clst_pcl[, zscore := scale(gex), by=array_id]
clst_avg <- clst_pcl[, list(zscore=mean(zscore)), by= list(array_id, Cluster, Classification)]

peakAll <- merge(peakAll, unique(annotation[, list(Cluster, Classification)]), by='Cluster')

orderby <- data.table(Classification=c('asexual', 'gam', 'shared', 'shared'), by=c('Llinas.Timepoint', 'Voss.Timepoint', 'Llinas.Timepoint', 'Voss.Timepoint'))
clst_avg <- merge(clst_avg, orderby, by='Classification', allow.cartesian=TRUE)
clst_avg[, ClusterClass := sprintf('%s_%s', Cluster, sub('\\.Timepoint', '', by))]

peak_order <- c()
for(x in c('asexual', 'gam', 'shared')) {
    if(x == 'asexual') {
        xord <- peakAll[Classification == x][order(Llinas.Timepoint, Llinas.Value)]
        xord[, ClusterClass := sprintf('%s_%s', Cluster, 'Llinas')]
        peak_order <- c(peak_order, xord$ClusterClass)
    } else if(x == 'gam') {
        xord <- peakAll[Classification == x][order(Voss.Timepoint, Voss.Value)]
        xord[, ClusterClass := sprintf('%s_%s', Cluster, 'Voss')]
        peak_order <- c(peak_order, xord$ClusterClass)
    } else if(x == 'shared') {
        xord <- peakAll[Classification == x][order(Llinas.Timepoint, Llinas.Value)]
        xord[, ClusterClass := sprintf('%s_%s', Cluster, 'Llinas')]
        peak_order <- c(peak_order, xord$ClusterClass)
        xord <- peakAll[Classification == x][order(Voss.Timepoint, Voss.Value)]
        xord[, ClusterClass := sprintf('%s_%s', Cluster, 'Voss')]
        peak_order <- c(peak_order, xord$ClusterClass)
    } else {
        stop()
    }
}
clst_avg[, ClusterClass := factor(ClusterClass, peak_order)]

# All clusters
# ============

mat <- dcast(data=clst_avg, array_id ~ ClusterClass, value.var='zscore')
mat <- as.matrix(mat, rownames='array_id')
mat <- mat[, match(peak_order, colnames(mat))]

clst_ann <- unique(clst_avg[, list(Classification=sprintf('%s_%s', Classification, sub('\\.Timepoint', '', by)), ClusterClass=as.character(ClusterClass))])
clst_ann <- clst_ann[match(colnames(mat), ClusterClass)]
stopifnot(identical(clst_ann$ClusterClass, colnames(mat)))
ds <- unique(clst_ann$Classification)
cols <- brewer.pal(n=length(ds), 'Dark2')
names(cols) <- ds
cls <- HeatmapAnnotation(df=clst_ann[, list(Classification=Classification)], col=list(Classification=cols), show_legend=FALSE, annotation_label='',
    cluster=anno_text(sub('_.*', '', clst_ann$ClusterClass), gp=gpar(fontsize=6)))

col_labels <- melt(data=peakAll[, list(Cluster, Classification, Llinas.Timepoint, Voss.Timepoint)], id.vars=c('Cluster', 'Classification'),
    value.name='Timepoint')
col_labels[, ClusterClass := sprintf('%s_%s', Cluster, sub('\\.Timepoint', '', variable))]
col_labels <- col_labels[match(colnames(mat), ClusterClass)]
stopifnot(identical(col_labels$ClusterClass, colnames(mat)))

tp <- rep(NA, nrow(col_labels))
tp[1] <- col_labels$Timepoint[1]
cur <- tp[1]
for(i in 2:nrow(col_labels)) {
    if(cur == col_labels$Timepoint[i]) {
        tp[i] <- ''
    } else {
        tp[i] <- col_labels$Timepoint[i]
        cur <- tp[i]
    }
}
col_labels[, tp := tp]

cluster_rows <- ladderize(as.dendrogram(hclust(dist(mat))))

dt <- snakemake@wildcards[['dataset']]

pdf(snakemake@output[['hm']], height=((0.3 * nrow(mat)) + 4)/2.54, width=60/2.54)
Heatmap(mat, name='Z-score\nof gene\nexpression', 
    column_split=clst_ann$Classification,
    column_gap=unit(4, 'mm'),
    column_title=unique(clst_ann$Classification),
    row_labels=rep('', nrow(mat)), 
    column_labels=col_labels$tp, 
    column_names_rot=90,
    # column_title='Cluster', 
    row_title=sprintf('Array from %s', dt), 
    top_annotation=cls, 
    cluster_columns=FALSE, cluster_rows=cluster_rows)
dev.off()

# Line plot
# =========

class_avg <- clst_pcl[, list(zscore=mean(zscore), sd=sd(zscore)), list(Cluster, Classification)]
class_avg <- merge(class_avg, orderby, by='Classification', allow.cartesian=TRUE)
class_avg[, ClusterClass := sprintf('%s_%s', Cluster, sub('\\.Timepoint', '', by))]
class_avg[, ClusterClass := factor(ClusterClass, peak_order)]
class_avg[, Classification := sprintf('%s_%s', Classification, sub('\\.Timepoint', '', by))]

gg <- ggplot(data=class_avg, aes(x=ClusterClass, y=zscore, group=Classification)) +
    geom_segment(aes(xend=ClusterClass, y=zscore-sd, yend=zscore+sd), col='grey80') +
    geom_line() +
    facet_wrap(~Classification, nrow=1, scales='free_x') +
    xlab('Cluster') +
    ylab('Gene expression z-score') +
    ggtitle(sprintf('Average and SD across samples from %s', dt)) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(), axis.text.x=element_blank(), strip.text=element_text(colour='black'))
ggsave(snakemake@output[['line']], width=24, height=10, units='cm')

## ## Clusters as boxplots
## # ---------------------
## 
## boxdata <- clst_pcl[Cluster %in% keep_cluster]
## boxdata <- boxdata[, list(gex=mean(gex)), by=list(Classification, Cluster, dataset)]
## 
## gg <- ggplot(data=boxdata, aes(x=dataset, y=gex)) +
##     geom_quasirandom() +
##     geom_point(data=boxdata[, list(gex=median(gex)), list(dataset, Classification)], pch='-', size=20, colour='orange') +
##     facet_wrap(~Classification) +
##     theme_light() +
##     ylab('Mean cluster expression') +
##     theme(strip.text=element_text(colour='black'))
## ggsave(snakemake@output[['clst_expr_by_classification']], width=16, height=16, units='cm')
## 
## gg <- ggplot(data=boxdata, aes(x=Classification, y=gex)) +
##     geom_quasirandom(size=0.5) +
##     geom_point(data=boxdata[, list(gex=median(gex)), list(dataset, Classification)], pch='-', size=15, colour='orange') +
##     facet_wrap(~dataset) +
##     theme_light() +
##     ylab('Mean cluster expression') +
##     theme(strip.text=element_text(colour='black'))
## ggsave(snakemake@output[['clst_expr_by_dataset']], width=18, height=8, units='cm')
## 
## 
## # Selected clusters
## # =================
## 

## gsea_clusters <- fread(snakemake@input[['gsea_clusters']]) # 'output/gsea/gsea_clusters.tsv'
## gsea_clusters[, Cluster := as.character(Cluster)]

## gsea_trait <- unique(gsea_clusters[padj < 0.01 & trait == 'Deformability', list(Cluster, NES)])
## keep_cluster <- unique(gsea_trait$Cluster)
## mat <- dcast(data=clst_avg[Cluster %in% keep_cluster], Cluster ~ array_id, value.var='zscore')
## mat <- as.matrix(mat, rownames='Cluster')
## mat <- mat[match(cluster_order, rownames(mat)),]
## mat <- mat[!is.na(rownames(mat)),]
## 
## ann_df <- unique(clst_avg[, list(array_id, dataset)])
## ann_df <- ann_df[match(colnames(mat), array_id)]
## ds <- unique(ann_df$dataset)
## cols <- brewer.pal(n=length(ds), 'Dark2')
## names(cols) <- ds
## ha <- HeatmapAnnotation(df= ann_df[, list(dataset=dataset)], col=list(dataset=cols))
## 
## clst_ann <- unique(annotation[!is.na(Classification), list(Cluster, Classification)])
## clst_ann <- merge(clst_ann, gsea_trait, by='Cluster')
## clst_ann <- clst_ann[order(NES)]
## mat <- mat[match(clst_ann$Cluster, rownames(mat)),]
## stopifnot(identical(clst_ann$Cluster, rownames(mat)))
## ds <- unique(clst_ann$Classification)
## class_cols <- brewer.pal(n=length(ds), 'Dark2')
## names(class_cols) <- ds
##  
## cls <- rowAnnotation(
##     `GSEA\nenrich.\ndeform.`=anno_barplot(clst_ann$NES, border=FALSE, width=unit(2, 'cm')),
##     df=clst_ann[, list(Classification=Classification)], col=list(Classification=class_cols)
##     )
## 
## # cluster_rows <- ladderize(as.dendrogram(hclust(dist(mat))))
## cluster_columns <- ladderize(as.dendrogram(hclust(dist(t(mat)))))
## 
## pdf(snakemake@output[['hm_gsea']], height=16/2.54, width=18/2.54)
## Heatmap(mat, name='Z-score\nof gene\nexpression', column_labels=rep('', ncol(mat)), 
##     row_title='Cluster', column_title='Array', bottom_annotation=ha, left_annotation=cls,
##     cluster_rows=FALSE, cluster_columns=cluster_columns)
## dev.off()
