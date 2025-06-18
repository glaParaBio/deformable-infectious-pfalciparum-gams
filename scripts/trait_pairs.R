library(data.table)
library(ggplot2)
library(ggbeeswarm)

source(file.path(snakemake@scriptdir, 'utils.R'))

ss <- fread(snakemake@input[['ss']]) 

lss <- parse_sample_sheet(ss) 

defo_key <- unique(lss[trait == 'Deformability', list(value, dicho_value)])
stopifnot(nrow(defo_key) == 2)

defo_legend <- sprintf('Deformability: %s: %s; %s: %s', 
                       defo_key$value[1], defo_key$dicho_value[1],
                       defo_key$value[2], defo_key$dicho_value[2])

# NB: These cutoff are only for plotting. The actual cutoffs are in
# parse_sample_sheet
dicho_cutoff <- fread('trait\tcutoff
                       log2(Exflagellation/ml)\t0.75
                       Oocysts per mosquito\t0.75
                       % Infected mosquitoes\t50
                       Deformability\t0.5')

gg <- ggplot(data= lss, aes(x= value)) +
    geom_histogram(bins= 20) +
    geom_vline(data=dicho_cutoff, aes(xintercept=cutoff), colour='blue', linetype='dashed') + 
    facet_wrap(~trait, scales= 'free', strip.position = 'bottom') +
    xlab(NULL) +
    ylab('N of samples') +
    ggtitle(defo_legend) +
    theme_light() +
    theme(strip.background= element_blank(), strip.placement= 'outside', strip.text= element_text(colour= 'black', size= 12), plot.title = element_text(size=10))
ggsave(snakemake@output[['hist']], width= 16, height= 16, units= 'cm')

#dicho_cutoff <- fread('trait\tcutoff
#                       log2(Exflagellation/ml)\t0.25
#                       Oocysts per mosquito\t0.25
#                       % Infected mosquitoes\t50
#                       Deformability\t1.5')


pairs <- combn(unique(lss[trait != 'Deformability']$trait), 2)
# We manually change the order of pairs to be in order of correlation. We could
# do it programmatically, but it would be long winded. It only affects the
# order of panels so it's not crucial to have solid code for it
pairs <- pairs[, c(3, 1, 2)]

gglist <- list()
for(i in 1:ncol(pairs)) {
    x <- lss[trait == pairs[1, i]]
    # x <- merge(x, dicho_cutoff, by='trait')
    y <- lss[trait == pairs[2, i]]
    # y <- merge(y, dicho_cutoff, by='trait')
    xy <- merge(x, y, by= c('sample_ID', 'Day'))
    if(nrow(xy) > 0) {
        xcor <- cor.test(xy$value.x, xy$value.y, method='spearman', use='pairwise.complete.obs')
        xcor <- sprintf('r = %.2f; p = %.3f', xcor$estimate, xcor$p.value)
        gg <- ggplot(data= xy, aes(x= value.x, y= value.y)) +
            annotate('text', x=-Inf, y=-Inf, label=xcor, hjust=-0.1, vjust=-4, size=3, colour='grey30') +
            geom_smooth(formula='y ~ x', se=FALSE, colour='grey40', size=0.5, method='loess') +
            geom_quasirandom(width=0.15, groupOnX=TRUE) +
            xlab(sprintf('%s\nN = %s', pairs[1, i], nrow(xy))) +
            ylab(pairs[2, i]) +
            theme_light() +
            theme(axis.title= element_text(size= 10), legend.position='none')
        gglist[[length(gglist) + 1]] <- gg
    }
}
ggsave(file = snakemake@output[['pairs']], gridExtra::arrangeGrob(grobs = gglist, ncol = 3), width= 22, height= 8, units= 'cm')

# Rank sum test of deformability vs exflagellation
# x <- lss[trait == 'Deformability']
# y <- lss[trait == 'log2(Exflagellation/ml)']
# xy <- merge(x, y, by= c('sample_ID', 'Day'))
# wilcox.test(xy[dicho_value.x == 'Mixed']$value.y, xy[dicho_value.x == 'Deformable']$value.y)

# Correlation between pairs of traits
# ============================================

cmat <- as.matrix(dcast(lss[trait != 'Deformability'], sample_ID ~ trait, value.var='value'), rownames='sample_ID')
M <- cor(cmat, method='spearman', use='pairwise.complete.obs')
M[upper.tri(M, diag=TRUE)] <- NA
M <- melt(as.data.table(M, keep.rownames='t1'), id.vars='t1', variable.name='t2', value.name='cor')
M <- M[!is.na(cor)]
M[, label := sprintf('%s vs\n%s', t1, t2)]
xord <- M[order(cor)]$label
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
    geom_vline(xintercept=0, colour='grey60', linetype='dotted') +
    geom_segment(aes(x=0, xend=cor, yend=label), colour='grey60') +
    geom_segment(aes(x=ci_low, xend=ci_high, yend=label), colour='black') +
    geom_point() +
    xlab('Spearman correlation with 95% CI') +
    ylab('') +
    xlim(-1, 1) +
    theme_light()
ggsave(snakemake@output[['corr']], width=12, height=6, units='cm')

