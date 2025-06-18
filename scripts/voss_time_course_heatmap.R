library(data.table)
library(ComplexHeatmap)
library(ggplot2)

whichMax <- function(x, which=1) {{
    xord <- order(-x)[which]
    return(xord)
}}

order_clusters <- function(clst_avg) {
    dat <- clst_avg[, list(Time, Cluster, zscore)]

    # Order clusters by those that peak at T0, break ties peak time
    keycol <- c()
    for(i in 1:length(unique(dat$Time))) {{
        peaks <- dat[, .SD[whichMax(zscore, i)], Cluster]
        dat <- merge(dat, peaks[, list(Cluster, peakTime=Time)], by='Cluster')
        k <- sprintf('peakTime%s', i)
        setnames(dat, 'peakTime', k)
        keycol <- c(keycol, k)
    }}
    setorderv(dat, keycol)
    cluster_order <- unique(dat$Cluster)
    return(cluster_order)
}

plot_heatmap <- function(clst_pcl, show_clusters) {
    clst_avg <- clst_pcl[, list(zscore=mean(zscore), sd=sd(zscore), .N), by=list(Cluster, Time)][N > 5]
    cluster_order <- order_clusters(clst_avg)

    # stopifnot(show_clusters %in% cluster_order)
    row_labels <- ifelse(cluster_order %in% show_clusters, cluster_order, '')

    mat <- as.matrix(dcast(clst_avg, Cluster ~ Time, value.var='zscore'), rownames='Cluster')
    mat <- mat[match(cluster_order, rownames(mat)), ]

    #pdf(snakemake@output[['hm']], width=14/2.54, height=16/2.54)
    hm <- Heatmap(mat,
        row_labels=row_labels,
        row_names_side='left',
        row_title='Cluster',
        row_names_gp = gpar(fontsize = 10),
        column_names_rot=0,
        column_title='Time point (Voss dataset)',
        column_names_side='top',
        cluster_columns=FALSE, 
        cluster_rows=FALSE,
        name='Z-score\nof gene\nexpression')
    # dev.off()
    return(hm)
}

kc_class <- fread(snakemake@input[['kc_class']])
gsea <- fread(snakemake@input[['gsea']])
vossPcl  <- fread(snakemake@input[['vossPcl']])
clstData <- fread(snakemake@input[['clstData']])
setnames(clstData, 'V1', 'gene_id')

vossPcl[, name := NULL]
vossPcl[, GWeight := NULL]
setnames(vossPcl, 'geneid', 'gene_id')
vossPcl <- as.matrix(vossPcl, rownames='gene_id')
vossPcl <- limma::normalizeQuantiles(vossPcl)
vossPcl <- melt(data= as.data.table(vossPcl, keep.rownames='gene_id'), id.vars= 'gene_id', value.name= 'gex', variable.name= 'Time')

vossPcl[, Time := as.character(Time)]
vossPcl[, Time := as.integer(sub('^day', '', Time))]
vossPcl[, zscore := scale(gex), by= gene_id]

show_clusters <- unique(gsea[padj < 0.01 & trait == 'deformable_vs_mixed']$Cluster)

clst_pcl <- merge(clstData, vossPcl, by='gene_id')
hh <- plot_heatmap(clst_pcl, show_clusters)
pdf(snakemake@output[['hm']])
print(hh)
dev.off()

use_class <- kc_class[Consensus == 'gam']$Cluster
clst_sub <- merge(clstData, vossPcl, by='gene_id')[Cluster %in% use_class]
hh <- plot_heatmap(clst_sub, show_clusters)
pdf(snakemake@output[['hm_gam']])
print(hh)
dev.off()

use_class <- kc_class[Consensus == 'shared']$Cluster
clst_sub <- merge(clstData, vossPcl, by='gene_id')[Cluster %in% use_class]
hh <- plot_heatmap(clst_sub, show_clusters)
pdf(snakemake@output[['hm_shared']])
print(hh)
dev.off()

# Line plots
# ==========

clst_avg <- clst_pcl[, list(zscore=mean(zscore), sd=sd(zscore), .N), by=list(Cluster, Time)][N > 5]
sel_clst <- merge(clst_pcl[Cluster %in% show_clusters], clst_avg[, list(Cluster, Time, avg_zscore=zscore)], by=c('Cluster', 'Time'))

cluster_order <- order_clusters(clst_avg)
clst_pcl[, Cluster := factor(Cluster, cluster_order)]
clst_avg[, Cluster := factor(Cluster, cluster_order)]

tp <- unique(clst_pcl$Time)
tp <- tp[tp %% 2 == 0]
gg <- ggplot(data=clst_pcl[Cluster %in% show_clusters], aes(x=Time, y=zscore, group=gene_id)) +
    geom_hline(yintercept=0, linewidth=0.05) +
    geom_line(colour='grey60', linewidth=0.25, alpha=0.5) +
    geom_line(data=clst_avg[Cluster %in% show_clusters], colour='orange', aes(group=Cluster)) +
    geom_segment(data=clst_avg[Cluster %in% show_clusters], aes(xend=Time, y=zscore - sd, yend=zscore + sd, group=NULL), colour='orange') +
    scale_x_continuous(breaks=tp) +
    facet_wrap(~Cluster, nrow=2) +
    xlab('Time point (Voss dataset)') +
    ylab('Z-score of gene expression') +
    theme_light() +
    theme(strip.text=element_text(colour='black'))
ggsave(snakemake@output[['line']], width=20, height=10, units='cm')
