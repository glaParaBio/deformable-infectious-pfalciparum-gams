library(data.table)
library(DESeq2)
library(ggplot2)
library(sva)
library(ggrepel)
library(parallel)
library(dplyr)
library(ggbeeswarm)

source(file.path(snakemake@scriptdir, 'utils.R'))
ss <- fread(snakemake@input[['ss']]) 
lss <- parse_sample_sheet(ss) 
ss <- dcast(lss, sample_ID + plate ~ trait, value.var='value')

# genes <- fread('ref/Pfalciparum3D7.genes.tsv')
genes <- fread(snakemake@input[['genes']])

cntfiles <- snakemake@input[['cntfiles']] 

counts <- list()
for(x in cntfiles) {
    sample_id <- sub('\\.count\\.gz$', '', basename(x))
    if(sample_id %in% ss$sample_ID) {
        cnt <- fread(x, header= FALSE, col.names= c('gene_id', 'count'))
        smry <- cnt[grepl('^__', gene_id)]
        cnt <- cnt[!gene_id %in% smry$gene_id]
        tot <- (sum(cnt$count) + sum(smry$count))
        pct_assigned <- 100 * (sum(smry$count) / tot)
        # cat(sprintf('%s | tot = %s | %.1f%% assigned\n', sample_id, tot, pct_assigned))
        cnt[, sample_id := sample_id]
        counts[[length(counts) + 1]] <- cnt
    } 
}
counts <- rbindlist(counts)
counts[, gene_id := sub('\\..*|-.*', '', sub('rna_', '', gene_id))]
counts <- counts[, list(count= sum(count)), by= list(gene_id, sample_id)]

id_not_found <- unique(counts[!gene_id %in% genes$gene_id][count > 0]$gene_id)
stopifnot(length(id_not_found) < 20)
counts <- counts[gene_id %in% genes$gene_id]

ss <- ss[sample_ID %in% counts$sample_id]
lss <- lss[sample_ID %in% counts$sample_id]

mat <- as.matrix(dcast(data= counts, gene_id ~ sample_id, value.var= 'count'), rownames= 'gene_id')
mat <- mat[, match(ss$sample_ID, colnames(mat))]
stopifnot(identical(colnames(mat), ss$sample_ID)) 

fwrite(as.data.table(mat, keep.rownames='gene_id'), snakemake@output[['cnt']], sep='\t')

keep <- edgeR::filterByExpr(mat, min.prop= 0.1)
mat <- mat[keep, ]

vsd <- vst(mat, fitType= 'local')
select <- order(apply(vsd, 1, var), decreasing=TRUE)[1:500]
pca <- prcomp(t(vsd[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
d <- data.table(PC1= pca$x[,1], PC2= pca$x[,2], sample_ID= colnames(vsd))
pcaout <- merge(d, ss[, list(sample_ID, plate)], by= 'sample_ID', all= TRUE)

xy <- range(d[, list(PC1, PC2)])
gg1 <- ggplot(data= pcaout, aes(x= PC1, y= PC2, colour= plate, label= sample_ID)) +
    geom_point() +
    scale_color_brewer(palette= 'Dark2') +
    ylim(xy[1], xy[2]) +
    xlim(xy[1], xy[2]) +
    xlab(sprintf('PC1: %.1f%%', percentVar[1] * 100)) +
    ylab(sprintf('PC2: %.1f%%', percentVar[2] * 100)) +
    ggtitle('Raw counts') +
    theme_light() 

bat <- ComBat_seq(mat, batch= ss$plate, group= NULL)
vsd <- vst(bat, fitType= 'local')
# rm(bat) # Remove, so you don't get confused between raw and adjusted count matrices
select <- order(apply(vsd, 1, var), decreasing=TRUE)[1:500]
pca <- prcomp(t(vsd[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
d <- data.table(PC1= pca$x[,1], PC2= pca$x[,2], sample_ID= colnames(vsd))
pcaout <- merge(d, ss[, list(sample_ID, plate)], by= 'sample_ID', all= TRUE)

xy <- range(d[, list(PC1, PC2)])
gg2 <- ggplot(data= pcaout, aes(x= PC1, y= PC2, colour= plate, label= sample_ID)) +
    geom_point() +
    # geom_text_repel(show.legend = FALSE) +
    scale_color_brewer(palette= 'Dark2') +
    ylim(xy[1], xy[2]) +
    xlim(xy[1], xy[2]) +
    xlab(sprintf('PC1: %.1f%%', percentVar[1] * 100)) +
    ylab(sprintf('PC2: %.1f%%', percentVar[2] * 100)) +
    ggtitle('Counts corrected for plate') +
    theme_light() 

gg <- gridExtra::arrangeGrob(grobs = list(gg1, gg2), nrow= 1)
ggsave(snakemake@output[['gex_pca_combatseq']], gg, width= 24, height= 10, units= 'cm')

pcaout <- merge(d, lss, by= 'sample_ID', all= TRUE)

gg <- ggplot(data= pcaout, aes(x= value, y= PC1, colour=day_offset)) +
    geom_smooth(se= TRUE, method= 'lm') +
    geom_point(size= 0.5) +
    scale_color_brewer(palette='Dark2') +
    facet_wrap(~trait, scales= 'free_x', strip.position= 'bottom') +
    xlab(NULL) +
    theme_light() +
    theme(strip.background= element_blank(), strip.placement= 'outside', strip.text= element_text(colour= 'black', size= 12))
ggsave(snakemake@output[['gex_pca_traits']], width= 18, height= 16, units= 'cm')

vsd <- melt(as.data.table(vsd, keep.rownames='gene_id'), id.vars='gene_id', variable.name='sample_ID', value.name='gex')
fwrite(vsd, snakemake@output[['vsd']], sep='\t')

# DGE
# ===
lss[, dicho_value2 := paste0(trait, dicho_value)]
lss[, dicho_value2 := sub('log2(Exflagellation/ml)', 'Exfl_', dicho_value2, fixed=TRUE)]
lss[, dicho_value2 := sub('Deformability', '', dicho_value2, fixed=TRUE)]
lss[, dicho_value2 := sub('% Infected mosquitoes', 'Mosq_', dicho_value2, fixed=TRUE)]
lss[, dicho_value2 := sub('Oocysts per mosquito', 'Oocy_', dicho_value2, fixed=TRUE)]
lss[, dicho_value2 := gsub(' ', '_', dicho_value2)]

dat <- dcast(lss[trait %in% c('log2(Exflagellation/ml)', '% Infected mosquitoes', 'Oocysts per mosquito', 'Deformability')], sample_ID ~ trait, value.var='dicho_value2')
setnames(dat, c('log2(Exflagellation/ml)', '% Infected mosquitoes', 'Oocysts per mosquito', 'Deformability'), c('Exfl', 'Mosq', 'Oocy', 'Defo'))
for(x in names(dat)) {
    dat[[x]] <- ifelse(is.na(dat[[x]]), sprintf('%s_NA', x), dat[[x]])
}
dat <- data.table(dat)
dat[, group := paste(Exfl, Mosq, Oocy, Defo, sep= '.')]

design <- model.matrix(~0 + group, dat)
rownames(design) <- dat$sample_ID
colnames(design) <- sub('^group', '', colnames(design))

contrasts <- limma::makeContrasts(
    # "deformable" aka "down"; "mixed" aka "up" 
    deformable_vs_mixed=Exfl_NA.Mosq_NA.Oocy_NA.Deformable - 
                      ((Exfl_Not_detected.Mosq_NA.Oocy_NA.Mixed + 
                        Exfl_Detected.Mosq_NA.Oocy_NA.Mixed)/2),

    infected_mosq_high_vs_low=((Exfl_Detected.Mosq_High_infection.Oocy_Detected.Defo_NA + 
                                Exfl_NA.Mosq_High_infection.Oocy_Detected.Defo_NA)/2) - 
                              ((Exfl_Detected.Mosq_Low_infection.Oocy_Detected.Defo_NA +
                                Exfl_Detected.Mosq_Low_infection.Oocy_Not_detected.Defo_NA + 
                                Exfl_Not_detected.Mosq_Low_infection.Oocy_Not_detected.Defo_NA + 
                                Exfl_NA.Mosq_Low_infection.Oocy_Detected.Defo_NA)/4),

    exfl_detected_vs_not_detected=((Exfl_Detected.Mosq_Low_infection.Oocy_Detected.Defo_NA + 
                                    Exfl_Detected.Mosq_NA.Oocy_NA.Defo_NA + 
                                    Exfl_Detected.Mosq_High_infection.Oocy_Detected.Defo_NA + 
                                    Exfl_Detected.Mosq_Low_infection.Oocy_Not_detected.Defo_NA + 
                                    Exfl_Detected.Mosq_NA.Oocy_NA.Mixed)/5) - 
                                  ((Exfl_Not_detected.Mosq_Low_infection.Oocy_Not_detected.Defo_NA + 
                                      Exfl_Not_detected.Mosq_NA.Oocy_NA.Mixed)/2),

    oocy_detected_vs_not_detected=((Exfl_Detected.Mosq_Low_infection.Oocy_Detected.Defo_NA + 
                                    Exfl_Detected.Mosq_High_infection.Oocy_Detected.Defo_NA + 
                                    Exfl_NA.Mosq_Low_infection.Oocy_Detected.Defo_NA + Exfl_NA.Mosq_High_infection.Oocy_Detected.Defo_NA)/4) - 
                                  ((Exfl_Detected.Mosq_Low_infection.Oocy_Not_detected.Defo_NA + 
                                    Exfl_Not_detected.Mosq_Low_infection.Oocy_Not_detected.Defo_NA)/2),
    levels=design
)
stopifnot(identical(rownames(contrasts), colnames(design)))

deseq <- DESeqDataSetFromMatrix(bat[, match(rownames(design), colnames(bat))], colData=dat, design=design)
deseq <- DESeq(deseq)

dge <- list()
for(cntr in colnames(contrasts)) {
    print(cntr)
    res <- lfcShrink(deseq, type= 'ashr', contrast= contrasts[, cntr], quiet= TRUE)
    res <- as.data.table(as.data.frame(res), keep.rownames= 'gene_id')
    res[, trait := cntr]
    dge[[length(dge) + 1]] <- res
}
dge <- rbindlist(dge)
dge <- merge(dge, genes, by= 'gene_id', sort= FALSE)

fwrite(dge, snakemake@output[['dge']], sep='\t')

dat <- copy(dge)
logfc_cut <- 0.2
dat <- dat[abs(log2FoldChange) > logfc_cut & pvalue < 0.05]
dat[, trait := gsub('_', ' ', sub('_vs_', '\nvs\n', trait))]
xord <- dat[, mean(abs(log2FoldChange)), trait][order(-V1)]
dat[, trait := factor(trait, xord$trait)]
dat[, padj := ifelse(padj < 0.0001, 0.0001, padj)]
cnt <- dat[, list(change=sign(log2FoldChange), trait)][, .N, by=list(change, trait)]

gg <- ggplot(data=dat, aes(x=trait, y=log2FoldChange, colour=-log10(padj))) +
    geom_hline(yintercept=logfc_cut, colour='orange') +
    geom_hline(yintercept=-logfc_cut, colour='orange') +
    geom_quasirandom(size=0.5, width=0.4) +
    geom_text(data=cnt[change == 1], aes(y=Inf, colour=NULL, label=N), vjust=2, hjust=1.5) +
    geom_text(data=cnt[change == -1], aes(y=-Inf, colour=NULL, label=N), vjust=-0.5, hjust=1.5) +
    xlab('') +
    ylab('Log2(fold change)') +
    theme_light()
ggsave(snakemake@output[['gex_logfc']], width=16, height=10, units='cm')
