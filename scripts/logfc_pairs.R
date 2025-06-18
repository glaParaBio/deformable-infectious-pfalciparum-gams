library(data.table)
library(ggplot2)

source(file.path(snakemake@scriptdir, 'utils.R'))

dge <- fread(snakemake@input[['dge']])
dge[, trait := ifelse(
  trait == 'deformable_vs_mixed', 'Deformability', ifelse(
  trait == 'exfl_detected_vs_not_detected', 'Exflagellation', ifelse(
  trait == 'infected_mosq_high_vs_low', 'Infected mosquitoes', ifelse(
  trait == 'oocy_detected_vs_not_detected', 'Oocysts in mosquitoes', NA))))]
stopifnot(is.na(dge$trait) == FALSE)

cmat <- as.matrix(dcast(dge, gene_id ~ trait, value.var='log2FoldChange'), rownames='gene_id')
M <- cor(cmat, method='spearman')
M[upper.tri(M, diag=TRUE)] <- NA
M <- melt(as.data.table(M, keep.rownames='t1'), id.vars='t1', variable.name='t2', value.name='cor')
M <- M[!is.na(cor)]
M[, label := sprintf('%s vs\n%s', t1, t2)]

## Same order as trait correlations
xord <- rev(c('Oocysts in mosquitoes vs\nInfected mosquitoes',
'Infected mosquitoes vs\nExflagellation',
'Oocysts in mosquitoes vs\nExflagellation',
'Exflagellation vs\nDeformability',
'Infected mosquitoes vs\nDeformability',
'Oocysts in mosquitoes vs\nDeformability'))

#xord <- M[order(cor)]$label
M[, label := factor(label, xord)]

ci_low <- rep(NA, nrow(M))
ci_high <- rep(NA, nrow(M))
for(i in 1:nrow(M)) {
    x <- M$t1[i]
    y <- M$t2[i]
    lh <- spearman_CI(cmat[, x], cmat[, y], use='pairwise.complete.obs')
    ci_low[i] <- lh[1]
    ci_high[i] <- lh[2]
}
M[, ci_low := ci_low]
M[, ci_high := ci_high]

gg <- ggplot(data=M, aes(x=cor, y= label)) +
    geom_segment(aes(x=0, xend=cor, yend=label), colour='grey60') +
    geom_segment(aes(x=ci_low, xend=ci_high, yend=label), colour='black') +
    geom_point() +
    xlab('Correlation in differential expression\n(i.e. log2(fold change)) with 95% CI') +
    ylab('') +
    xlim(0, 1) +
    theme_light()
ggsave(snakemake@output[['corr']], width=12, height=9, units='cm')

ldge <- dcast(dge, gene_id ~ trait, value.var='log2FoldChange')
ldge[, gene_id := NULL]
tt <- combn(names(ldge), 2)
dat <- list()
for(i in seq(1, ncol(tt))) {
    t1 <- tt[1, i]
    t2 <- tt[2, i]
    dat[[length(dat) + 1]] <- data.table(t1=t1, t2=t2, lfc1=unlist(ldge[, t1, with=FALSE]), lfc2=unlist(ldge[, t2, with=FALSE]))
}
dat <- rbindlist(dat)
dat[, label := sprintf('%s vs\n%s', t1, t2)]
xord <- dat[, list(cor=cor(lfc1, lfc2)), by=label][order(-cor)]$label
dat[, label := factor(label, xord)]

dat[, coldens := densCols(lfc1, lfc2, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))), by=label]
glist <- list()
for(l in sort(unique(dat$label))) {
    gdat <- dat[label == l]
    xlab <- unique(gdat$t1)
    ylab <- unique(gdat$t2)
    gg <- ggplot(data=gdat, aes(x=lfc1, y=lfc2, colour=coldens)) +
        geom_point(pch='.') +
        scale_color_identity() +
        xlab(xlab) +
        ylab(ylab) +
        theme_light() 
    glist[[length(glist) + 1]] <- gg
}
gg <- gridExtra::arrangeGrob(grobs = glist, nrow= 2)
ggsave(snakemake@output[['logfc_pairs']], gg, width= 18, height= 12, units= 'cm')
