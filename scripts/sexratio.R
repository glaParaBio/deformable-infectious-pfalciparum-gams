library(data.table)
library(ggplot2)
library(lme4)
library(emmeans)

tpm <- function(counts, gene.length, log= TRUE, prior.count= 2) {
    if(log == TRUE){
        counts <- counts + prior.count
    }
    rate <- counts / gene.length
    x <- rate / sum(rate) * 1e6
    if(log == TRUE){
        return(log2(x))
    } else {
        return(x)
    }
}

ss <- fread(snakemake@config[['rna_ss']])
tx <- fread(cmd=sprintf('grep ">" %s', snakemake@input[['fa']]), header=FALSE)
tx <- tx[, list(transcript_id=sub('>', '', V1), gene_id=sub('gene=', '', V2), len=sub('length=', '', V7))]
tx[, len := as.numeric(len)]
glen <- tx[, list(.N, maxlen=max(len)), gene_id]

counts <- fread(snakemake@input[['deseq_counts']])
counts <- as.matrix(counts, rownames='gene_id')
xtpm <- data.table(gene_id= rownames(counts), apply(counts, 2, tpm, glen[match(rownames(counts), gene_id),]$maxlen, log= TRUE, prior.count= 2))
xtpm <- melt(xtpm, id.vars='gene_id', variable.name='sample_ID', value.name='tpm')
xtpm <- merge(xtpm, ss[, list(sample_ID, Day, deformability, day_offset)], by='sample_ID')
marks <- data.table(
    gene_id = c('PF3D7_1311100', 'PF3D7_1469900', 'PF3D7_0208900', 'PF3D7_1031000', 'PF3D7_0630000', 'PF3D7_1351600', 'PF3D7_0903800'),
    marker = c('Male', 'Male', 'Male', 'Female', 'Female', 'Female', 'Female')
)
marks <- merge(marks, glen, by='gene_id')

xtpm <- merge(xtpm, marks, by='gene_id')
xtpm[, Day := as.numeric(Day)]

avg <- xtpm[, list(tpm=mean(tpm)), by=list(gene_id, marker, Day, deformability, day_offset)]
xord <- avg[, list(tpm=mean(tpm)), by=list(gene_id, marker)][order(marker, tpm)]
avg[, gene_id := factor(gene_id, levels=xord$gene_id)]
xtpm[, gene_id := factor(gene_id, levels=xord$gene_id)]

gg <- ggplot(data=avg[day_offset=='deform_exfl'], aes(x=Day, y=tpm, colour=marker, group=paste(gene_id, deformability), shape=deformability)) +
    geom_point(data=xtpm[day_offset=='deform_exfl'], colour='grey80', size=2.5) +
    geom_line() +
    facet_wrap(~gene_id) +
    ylab('Gene expression') +
    scale_x_continuous(breaks = unique(avg$Day)) +
    theme_light() +
    theme(panel.grid.minor.x = element_blank()) +
    ggtitle('Samples with deformability data')
ggsave(snakemake@output[['markers_sex_defo']], width=24, height=24, units='cm')

gg <- ggplot(data=avg[day_offset=='only_infect'], aes(x=Day, y=tpm, colour=marker, group=gene_id)) +
    geom_point(data=xtpm[day_offset=='only_infect'], colour='grey80', size=2) +
    geom_line() +
    facet_wrap(~gene_id) +
    ylab('Gene expression') +
    scale_x_continuous(breaks = unique(avg$Day)) +
    theme_light() +
    theme(panel.grid.minor.x = element_blank()) +
    ggtitle('Samples with infectivity data')
ggsave(snakemake@output[['markers_sex_infect']], width=24, height=24, units='cm')

avg <- xtpm[, list(tpm=mean(tpm)), by=list(sample_ID, Day, deformability, day_offset, marker)]
sr <- dcast(data=avg, sample_ID + Day + deformability + day_offset ~ marker, value.var='tpm')
sr[, ratio := Female - Male]

gg <- ggplot(data=sr, aes(x=Day, y=ratio, colour=deformability)) +
    geom_hline(yintercept=0, colour='orange', linetype='dotted') +
    geom_point() +
    facet_wrap(~day_offset, ncol=1) +
    ylab('Expression of female/male markers') +
    scale_x_continuous(breaks = unique(sr$Day))
ggsave(snakemake@output[['markers_sex_ratio']], width=18, height=16, units='cm')

## Difference deformability
daykeep <- xtpm[!is.na(deformability), list(N=length(unique(deformability))), list(Day)][N == 2]$Day
xdat <- xtpm[!is.na(deformability) & Day %in% daykeep]
xdat[, gene_id := sprintf('%s (%s)', marker, gene_id)]
xdat[, Day := as.factor(Day)]

fit <- lm(tpm ~ Day + gene_id * deformability, data=xdat)
emm <- emmeans(fit, ~deformability | gene_id)
contrast(emm, 'pairwise')

fit <- lm(tpm ~ Day + marker + gene_id * deformability, data=xdat)
emm <- emmeans(fit, ~deformability | marker, nesting='gene_id %in% marker')
contrast(emm, 'pairwise')

emm <- emmeans(fit, ~deformability | gene_id, nesting='gene_id %in% marker')
contrast(emm, 'pairwise')

fiter <- lmer(tpm ~ Day + gene_id*deformability + (1|sample_ID), data=xdat, REML=FALSE)
emm <- emmeans(fiter, ~deformability | gene_id)
contrast(emm, 'pairwise')

## Piecemeal
for (gid in sort(unique(xdat$gene_id))) {
    fit <- lm(tpm ~ Day + deformability, xdat[gene_id == gid])
    emm <- emmeans(fit, ~deformability)
    cat(sprintf('%s\\n', gid))
    print(contrast(emm, 'pairwise'))
}

vossPcl  <- fread(snakemake@input[['vossPcl']])
vossPcl[, name := NULL]
vossPcl[, GWeight := NULL]
setnames(vossPcl, 'geneid', 'gene_id')
vossPcl <- as.matrix(vossPcl, rownames='gene_id')
# vossPcl <- limma::normalizeQuantiles(vossPcl)
vossPcl <- melt(data= as.data.table(vossPcl, keep.rownames='gene_id'), id.vars= 'gene_id', value.name= 'gex', variable.name= 'Time')
vossPcl[, Time := as.character(Time)]
vossPcl[, Time := as.integer(sub('^day', '', Time))]
vossPcl <- merge(vossPcl, marks, by='gene_id')
xord <- vossPcl[, list(gex=mean(gex)), by=list(gene_id, marker)][order(marker, gex)]
vossPcl[, gene_id := factor(gene_id, xord$gene_id)]

gg <- ggplot(data=vossPcl, aes(x=Time, y=gex, colour=marker)) +
    geom_line() +
    facet_wrap(~gene_id, nrow=1) +
    ylab('Gene expression') +
    theme(legend.position='bottom') +
    scale_x_continuous(breaks = unique(vossPcl$Time)) +
    theme_light() +
    theme(panel.grid.minor.x = element_blank(), legend.position='bottom') +
    ggtitle('Voss time course')
ggsave(snakemake@output[['markers_sex_voss']], width=24, height=9, units='cm')

avg <- vossPcl[, list(gex=mean(gex)), by=list(gene_id, marker, Time)]
xord <- avg[, list(gex=mean(gex)), by=list(gene_id, marker)][order(marker, gex)]
avg[, gene_id := factor(gene_id, levels=xord$gene_id)]
vossPcl[, gene_id := factor(gene_id, levels=xord$gene_id)]

avg <- vossPcl[, list(gex=mean(gex)), by=list(Time, marker)]
sr <- dcast(data=avg, Time ~ marker, value.var='gex')
sr[, ratio := Female - Male]
gg <- ggplot(data=sr, aes(x=Time, y=ratio)) +
    geom_hline(yintercept=0, colour='orange', linetype='dotted') +
    geom_point() +
    geom_line() +
    ylab('Expression of female/male markers') +
    xlab('Day (Voss time course)') +
    scale_x_continuous(breaks = unique(sr$Time)) +
    theme_light() +
    theme(panel.grid.minor.x = element_blank())
ggsave(snakemake@output[['markers_sex_ratio_voss']], width=18, height=8, units='cm')


