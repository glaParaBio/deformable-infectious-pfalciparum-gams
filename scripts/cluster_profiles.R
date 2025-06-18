library(data.table)
library(ggplot2)
library(dendextend)

hclust_order <- function(mat) {
    dd <- dist(mat)
    hc <- as.dendrogram(hclust(dd))
    hc <- ladderize(hc)    
    xord <- labels(hc)
    return(xord)
}

plot_clusters <- function(data, cluster='Cluster', idvar='gene_id', xvar='Time', yvar='zscore', minSize=15) {
    dat <- copy(data)
    setnames(dat, c(idvar, cluster, xvar, yvar), c('idvar', 'cluster', 'xvar', 'yvar'))
    clst_avg <- dat[!is.na(yvar), list(yvar=mean(yvar), std=sd(yvar), .N), by=list(xvar, cluster)]
    clst_avg[, dist := abs(yvar) - (std * 2)]
    clst_avg[, big_dist := abs(yvar) - (std * 3)]
    keep <- clst_avg[N >= minSize, list(not_zero=sum((dist) > 0), big=sum(big_dist > 0)), cluster][not_zero >= 2 | big > 0]$cluster
    clst_avg <- clst_avg[cluster %in% keep]
    dat <- dat[cluster %in% keep]
    mat <- dcast(clst_avg, cluster ~ xvar, value.var= 'yvar')
    mat <- as.matrix(mat, rownames= 'cluster')
    mat <- t(scale(t(mat)))
    xord <- hclust_order(mat)
    dat <- dat[cluster %in% xord]
    dat[, cluster := factor(cluster, xord)]
    clst_avg[, cluster := factor(cluster, xord)]

    gg <- ggplot(data=dat, aes(x=xvar, y=yvar)) +
        geom_line(aes(group=idvar), colour='grey30', alpha=0.2) +
        geom_line(data=clst_avg, aes(group=cluster), colour='orange', size=0.5) +
        geom_segment(data=clst_avg, aes(xend=xvar, y=yvar - std, yend=yvar + std), colour='orange', size=0.25) +
        geom_hline(yintercept=0, colour='blue', linetype='dashed') +
        theme_light() +
        xlab('Day | Voss time course') +
        ylab('Z-score of gene expression') +
        facet_wrap(~cluster) +
        scale_x_discrete(breaks=sort(unique(dat$xvar)), labels=sub('^0', '', sort(unique(dat$xvar)))) +
        theme(strip.text=element_text(colour='black'))
    
    #setnames(dat, c('idvar', 'cluster', 'xvar', 'yvar'), c(idvar, cluster, xvar, yvar))
    #setnames(clst_avg, c('cluster', 'xvar', 'yvar'), c(cluster, xvar, yvar))

    return(list(data=dat, data_avg=clst_avg, gg=gg))
}

clstData <- fread(snakemake@input[['clstData']], select=c('V1', 'Cluster'))
setnames(clstData, 'V1', 'gene_id')
vossPcl <- fread(snakemake@input[['vossPcl']]) 
vossPcl[, name := NULL]
vossPcl[, GWeight := NULL]
setnames(vossPcl, 'geneid', 'gene_id')
vossPcl <- melt(data= vossPcl, id.vars= 'gene_id', value.name= 'gex', variable.name= 'Time')

vossPcl[, Time := as.character(Time)]
vossPcl[, Time := sprintf('%02d', as.integer(sub('^day', '', Time)))]
vossPcl[, zscore := scale(gex), by= gene_id]

clst_pcl <- merge(clstData, vossPcl, by='gene_id')

clst <- plot_clusters(clst_pcl)
ggsave(snakemake@output[['clusters_voss_time_course']], width=30, height=24, units='cm')

# ===============
# 
# hm <- plot_hdbscan_profiles(voss_pcl, 'gene_id', 'Time', 'zscore', minPts=5)
# hm$gg <- hm$gg + xlab('Day') + scale_x_discrete(breaks=unique(voss_pcl$Time))
# 
# cv <- function(x) {
#     x <- x - min(x)
#     return(sd(x) / mean(x))
# }
# 
# xcv <- voss_pcl[!is.na(gex), list(cv=cv(gex)), gene_id]
# hm$clusters <- merge(hm$clusters, xcv, by='gene_id')
# 
# noise <- voss_pcl[gene_id %in% hm$clusters[clst_id == 'noise']$gene_id]
# hm_noise <- plot_hdbscan_profiles(noise, 'gene_id', 'Time', 'zscore', minPts=5)
# 
# noise2 <- voss_pcl[gene_id %in% hm_noise$clusters[clst_id == 'noise']$gene_id]
# hm_noise2 <- plot_hdbscan_profiles(noise2, 'gene_id', 'Time', 'zscore', minPts=5)
# 
# noise3 <- voss_pcl[gene_id %in% hm_noise2$clusters[clst_id == 'noise']$gene_id]
# hm_noise3 <- plot_hdbscan_profiles(noise3, 'gene_id', 'Time', 'zscore', minPts=5)
# 
# noise4 <- voss_pcl[gene_id %in% hm_noise3$clusters[clst_id == 'noise']$gene_id]
# hm_noise4 <- plot_hdbscan_profiles(noise4, 'gene_id', 'Time', 'zscore', minPts=5)
# 
# hm$clusters[clst_id == 1]
# hm_noise$clusters[clst_id == 17]
# 
# 
# keep <- c(hm$clusters[clst_id != 'noise']$gene_id, 
#         hm_noise$clusters[clst_id != 'noise']$gene_id,
#         hm_noise2$clusters[clst_id != 'noise']$gene_id,
#         hm_noise3$clusters[clst_id != 'noise']$gene_id,
#         hm_noise4$clusters[clst_id != 'noise']$gene_id)
# 
# final <- voss_pcl[gene_id %in% keep]
# hm_final <- plot_hdbscan_profiles(final, 'gene_id', 'Time', 'zscore', minPts=5)
# 
# # Load *in vivo* datasets
# load('/export/projects/III-data/wcmp_bioinformatics/db291g/data/20210128_matt_infection/scripts/mattNetwork.RData')
# 
# inVivoLong <- list()
# for(dataset in names(inVivo)) {
#     ll <- as.data.table(inVivo[[dataset]], keep.rownames='gene_id')
#     ll <- melt(ll, id.vars='gene_id', variable.name='time_point', value.name='gex')
#     ll[, dataset := dataset]
#     inVivoLong[[length(inVivoLong) + 1]] <- ll
# }
# inVivoLong <- rbindlist(inVivoLong)
# 
# # Choose clusters to plot
# gsea <- fread('gsea/gsea_clusters.tsv')
# keep <- unique(gsea[padj < 0.01, list(Cluster, NES)])[, .SD[which.max(abs(NES))], by=Cluster][order(-NES)]
# keep <- gsea[Cluster %in% keep$Cluster]
# 
# mat <- inVivoLong[gene_id %in% keep[Cluster == 95]$gene_id]
# mat <- dcast(mat, gene_id ~ time_point, value.var='gex')
# mat <- as.matrix(mat, rownames='gene_id')
# 
# hm <- Heatmap(mat, name= 'Some name', cluster_columns=FALSE)
# pdf('tmp.pdf', width= 24/2.54, height= 18/2.54)
# draw(hm, heatmap_legend_side = "left")
# dev.off()
