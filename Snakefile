import glob
import re
import pandas
import tabulate

os.makedirs('slurm', exist_ok=True)

PLASMODB_RELEASE = '55'

rule all:
    input:
        ## Figure 1
        os.path.join(workflow.basedir, 'results/histogram_traits.pdf'),
        os.path.join(workflow.basedir, 'results/upset_traits.pdf'),

        ## Figure 2
        os.path.join(workflow.basedir, 'results/infection_time_course.deform.pdf'),
        os.path.join(workflow.basedir, 'results/infection_time_course.only_infect.pdf'),
        os.path.join(workflow.basedir, 'results/trait_pairs.pdf'), # Fig 2

        ## Figure 3
        os.path.join(workflow.basedir, 'results/gex_logfc.pdf'),
        os.path.join(workflow.basedir, 'results/dge_trait_corr.pdf'),
        'gsea/go/infection_deformability.pdf',

        ## Figure 4
        os.path.join(workflow.basedir, 'results/gsea_clusters.pdf'),

        ## Figure 5
        os.path.join(workflow.basedir, 'results/selected_clusters_patients.pdf'),
        os.path.join(workflow.basedir, 'results/candidate_genes_patients.pdf'),

        ## Supp Table and Figures
        'deseq/infect_dge.tsv.gz',
        'deseq/overlapping_dge.tsv.gz',
        os.path.join(workflow.basedir, 'results/deformability_infection_data.tsv'),
        'gsea/go/infection_deformability.tsv.gz',
        os.path.join(workflow.basedir, 'results/gex_pca_combatseq.pdf'),
        os.path.join(workflow.basedir, 'results/gex_pca_traits.pdf'),
        'deseq/overlapping_dge.tsv.gz',
        os.path.join(workflow.basedir, 'results/voss_time_course_heatmap_shared.pdf'),
        os.path.join(workflow.basedir, 'results/voss_time_course_heatmap_gam.pdf'),
        os.path.join(workflow.basedir, 'results/voss_time_course_selected.pdf'),


rule patients:
    input:
       os.path.join(workflow.basedir, 'results/patient_heatmap.pdf'),
       os.path.join(workflow.basedir, 'results/patient_lineplot.pdf'),


rule infection_time_course:
    input:
        ss= config['rna_ss'], 
        utils= os.path.join(workflow.basedir, 'scripts/utils.R'),
    output:
        only_infect=os.path.join(workflow.basedir, 'results/infection_time_course.only_infect.pdf'),
        deform=os.path.join(workflow.basedir, 'results/infection_time_course.deform.pdf'),
        tsv=os.path.join(workflow.basedir, 'results/deformability_infection_data.tsv'),
    script:
        os.path.join(workflow.basedir, 'scripts/infection_time_course.R')


rule trait_pairs:
    input:
        ss= config['rna_ss'],
        utils= os.path.join(workflow.basedir, 'scripts/utils.R'),
    output:
        hist= os.path.join(workflow.basedir, 'results/histogram_traits.pdf'),
        pairs= os.path.join(workflow.basedir, 'results/trait_pairs.pdf'),
        corr= os.path.join(workflow.basedir, 'results/trait_pairs_corr.pdf'),
    script:
        os.path.join(workflow.basedir, 'scripts/trait_pairs.R')


rule upset_traits:
    input:
        ss= config['rna_ss'],
        utils= os.path.join(workflow.basedir, 'scripts/utils.R'),
    output:
        upset= os.path.join(workflow.basedir, 'results/upset_traits.pdf'),
    script:
        os.path.join(workflow.basedir, 'scripts/upset_traits.R')
        

rule gff:
    output:
        gff= f'ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gff',
    shell:
        r"""
        curl -s https://plasmodb.org/common/downloads/release-{PLASMODB_RELEASE}/Pfalciparum3D7/gff/data/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gff > {output.gff}
        """

rule gaf:
    output:
        gaf= f'ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gaf',
    shell:
        r"""
        curl -s https://plasmodb.org/common/downloads/release-{PLASMODB_RELEASE}/Pfalciparum3D7/gaf/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7_GO.gaf > {output.gaf}
        """

rule gene_descriptions:
    input:
        gff= f'ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gff',
    output:
        gd= 'ref/Pfalciparum3D7.genes.tsv',
    run:
        from urllib.parse import unquote
        fout = open(output.gd, 'w')
        fout.write('\t'.join(['gene_id', 'description', 'name']) + '\n')

        with open(input.gff) as fin:
            for line in fin:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                if line[2] in ['gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene']:
                    attr = line[8].split(';')
                    gid = [re.sub('^ID=', '', x) for x in attr if x.startswith('ID=')]
                    assert len(gid) == 1
                    desc = [unquote(re.sub('^description=', '', x)) for x in attr if x.startswith('description=')]
                    assert len(desc) == 1
                    name = [unquote(re.sub('^Name=', '', x)) for x in attr if x.startswith('Name=')]
                    assert len(name) == 0 or len(name) == 1
                    if name == []:
                        name = ['N/A']
                    outline = gid + desc + name
                    fout.write('\t'.join(outline) + '\n')
        fout.close()


rule deseq:
    input:
        ss= config['rna_ss'],
        cntfiles= glob.glob(f"{config['counts']}/*.count.gz"),
        genes= 'ref/Pfalciparum3D7.genes.tsv',
        utils= os.path.join(workflow.basedir, 'scripts/utils.R'),
    output:
        cnt='deseq/counts.tsv.gz',
        vsd='deseq/vsd.tsv.gz',
        dge= 'deseq/infect_dge.tsv.gz',
        gex_pca_combatseq=os.path.join(workflow.basedir, 'results/gex_pca_combatseq.pdf'),
        gex_pca_traits=os.path.join(workflow.basedir, 'results/gex_pca_traits.pdf'),
        gex_logfc=os.path.join(workflow.basedir, 'results/gex_logfc.pdf'),
    script:
        os.path.join(workflow.basedir, 'scripts/deseq.R')
        

rule volcano:
    input:
        dge= 'deseq/infect_dge.tsv.gz',
    output:
        volcano= os.path.join(workflow.basedir, 'results/volcano.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

dge <- fread('{input.dge}')
dge[, trait := gsub('_', ' ', trait)]
xord <- dge[, mean(abs(log2FoldChange)), trait][order(-V1)]
dge[, trait := factor(trait, xord$trait)]

gg <- ggplot(data=dge, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(size=0.1) +
    facet_wrap(~trait) +
    theme_light() +
    theme(strip.text=element_text(colour='black'))
ggsave('{output.volcano}', width=16, height=16, units='cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


#Atlas uses v52 or v56 or v59 which are the same as v55
#grep -w 'protein_coding_gene' PlasmoDB-52_Pfalciparum3D7.gff | sed 's/;.*//' | cut -f 1,4,5,9 | sort | md5sum
#grep -w 'protein_coding_gene' PlasmoDB-55_Pfalciparum3D7.gff | sed 's/;.*//' | cut -f 1,4,5,9 | sort | md5sum
rule dge_vs_atlas:
    input:
        dge='deseq/infect_dge.tsv.gz',
        atlas= os.path.join(workflow.basedir, 'ref_data/A_single_cell_atlas_of_sexual_development_in_Plasmodium_falciparum.S11.tsv.gz'),
    output:
        pdf=os.path.join(workflow.basedir, 'results/dge_vs_atlas.pdf'),
        tsv='deseq/dge_vs_atlas.tsv.gz',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

dge <- fread('{input.dge}')
atlas <- fread('{input.atlas}')
atlas[, gene_id := sub('-', '_', gene_id)]

tsv <- list()
gglist <- list()
for(cntr in unique(dge$trait)) {{
    agex <- melt(id.vars=c('gene_id'), value.name='gex', variable.name='cluster', data=atlas[, list(gene_id, FH_exp, FL_exp, MH_exp, ML_exp)])
    agex[, cluster := sub('_exp', '', cluster)]
    agex[, cluster := sub('F', 'Female', cluster)]
    agex[, cluster := sub('M', 'Male', cluster)]
    agex[, cluster := sub('H', ' high', cluster)]
    agex[, cluster := sub('L', ' low', cluster)]
    agex <- merge(agex, dge[trait == cntr, list(gene_id, log2FoldChange, padj, trait, name, description)], by='gene_id')
    tsv[[length(tsv) + 1]] <- agex
    gg <- ggplot(data=agex, aes(x=gex, y=log2FoldChange)) +
        geom_hline(yintercept=0, colour='orange') +
        geom_point(pch='.') +
        geom_smooth() +
        xlab('Gene expression in cell atlas') +
        ylab('Log2 fold change') +
        ylim(quantile(agex$log2FoldChange, 0.001), quantile(agex$log2FoldChange, 0.999)) +
        ggtitle(sprintf('Expression in cell atlas vs comparison "%s"', gsub('_', ' ', cntr))) +
        facet_wrap(~cluster) +
        theme_light() +
        theme(strip.text=element_text(color='black'))
    gglist[[length(gglist) + 1]] <- gg
}}
gg <- gridExtra::arrangeGrob(grobs = gglist, nrow= 2)
ggsave('{output.pdf}', gg, width=35, height=30, units='cm')
tsv <- rbindlist(tsv)
fwrite(tsv, '{output.tsv}', sep='\t')
EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule maplot:
    input:
        dge='deseq/infect_dge.tsv.gz',
    output:
        maplot=os.path.join(workflow.basedir, 'results/maplot.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

dge <- fread('{input.dge}')

xord <- dge[, list(lfc=mean(abs(log2FoldChange))), trait][order(-lfc)]
dge[, trait := factor(trait, xord$trait)]

gg <- ggplot(data=dge, aes(x=baseMean, y=log2FoldChange)) +
    geom_point(size=0.1, colour='grey60', alpha=0.5) +
    geom_point(data=dge[padj < 0.01], colour='red', size=0.1) +
    geom_hline(yintercept=0, colour='black') +
    scale_x_log10(labels = scales::label_comma()) +
    xlab('Average expression') +
    ylab('log2(fold-change)') +
    facet_wrap(~trait) +
    theme_light() +
    theme(strip.text=element_text(color='black'))
ggsave('{output.maplot}', width=16, height=10, units='cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule logfc_pairs:
    input:
        dge= 'deseq/infect_dge.tsv.gz',
        utils= os.path.join(workflow.basedir, 'scripts/utils.R'),
    output:
        corr=os.path.join(workflow.basedir, 'results/dge_trait_corr.pdf'),
        logfc_pairs=os.path.join(workflow.basedir, 'results/dge_logfc_pairs.pdf'),
    script:
        os.path.join(workflow.basedir, 'scripts/logfc_pairs.R')
        

rule makeBioconductorAnnotationDbi:
    input:
        gff= f'ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gff',
        gaf= f'ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gaf',
    output:
        dbi= directory('ref/org.Pfalciparum3D7.eg.db'),
    params:
        ver= f'0.{PLASMODB_RELEASE}.0'
    shell:
        r"""
        {workflow.basedir}/scripts/makeBioconductorAnnotationDbi.r --gff {input.gff} \
            --gaf {input.gaf} --genus Plasmodium --species falciparum3D7 --taxid 36329 \
            --outdir `dirname {output.dbi}` -m 'Dario Beraldi <dario.beraldi@glasgow.ac.uk>' \
            -a 'Dario Beraldi' --pckg-version {params.ver} --install
        """


rule clusterprofiler:
    input:
        dbi= 'ref/org.Pfalciparum3D7.eg.db',
        dge= 'deseq/infect_dge.tsv.gz',
    output:
        go_tsv= 'gsea/go/infection_deformability.tsv.gz',
        kegg_tsv= 'gsea/kegg/infection_deformability.tsv.gz',
    shell:
        r"""
        {workflow.basedir}/scripts/clusterprofiler.r \
            --dge {input.dge} \
            --orgdb `basename {input.dbi}` \
            --gsea-output-tsv {output.go_tsv} \
            --kegg-output-tsv {output.kegg_tsv} \
            --kegg-organism pfa
        """

rule plot_clusterprofiler:
    input:
        tsv= 'gsea/{annotation}/infection_deformability.tsv.gz',
    output:
        pdf= 'gsea/{annotation}/infection_deformability.pdf',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

tsv <- fread('{input.tsv}')
keep <- tsv[pvalue < 1]
minus <- keep[order(trait, NES)][, rank := 1:nrow(.SD), trait]
top <- keep[order(trait, -NES)][, rank := 1:nrow(.SD), trait]
keep <- rbind(minus[rank <= 5], top[rank <= 5])[order(trait, NES)]
keep[, rank := 1:nrow(keep)]
keep[, rank := factor(rank)]
keep[, colour := ifelse(p.adjust < 0.05, 'black', 'grey60')]

gg <- ggplot(data=keep, aes(x=NES, y=rank, colour=colour)) +
    scale_color_identity() + 
    geom_vline(xintercept=0, colour='grey40', linetype='dotted') +
    geom_segment(aes(x=0, xend=NES, y=rank, yend=rank)) +
    geom_point() +
    scale_y_discrete(label=NULL) +
    geom_text(aes(label=Description), x=min(keep$NES) * 1.15, hjust=1, size=3) +
    coord_cartesian(clip='off') +
    xlab('Normalised enrichment score') +
    ylab('') +
    facet_wrap(~trait, scales='free_y', ncol=1) +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'), plot.margin = margin(0.1, 0.1, 0.5, 10, "cm"))
ggsave('{output.pdf}', width=18, height=18, units='cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule gsea_clusters:
    input:
        dge= 'deseq/infect_dge.tsv.gz',
        clstData= config['clstData'],
        cluster_class=os.path.join(workflow.basedir, 'ref_data/crouch_cluster_classification.tsv'),
    output:
        gsea_long= 'gsea/gsea_clusters.tsv',
        plot= os.path.join(workflow.basedir, 'results/gsea_clusters.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(fgsea)
library(ggplot2)

stars <- function(p) {{
    xstars <- ifelse(p < 0.001, '***', 
        ifelse(p < 0.01, '**', 
            ifelse(p < 0.05, '+', '')))
    return(xstars)
}}

dge <- fread('{input.dge}')

clusters <- fread('{input.clstData}')
setnames(clusters, c('V1', 'V2'), c('gene_id', 'description'))

cluster_class <- fread('{input.cluster_class}')

gseaOut <- list()
for(cf in unique(dge$trait)) {{
    full <- dge[trait == cf]
    full[, score := log2FoldChange]
    all_genes <- full$score
    names(all_genes) <- full$gene_id

    testable <- clusters[gene_id %in% names(all_genes)]
    clst <- list()
    knames <- as.character(unique(testable$Cluster))
    for(k in knames) {{
        clst[[k]] <- testable[Cluster == k]$gene_id
    }}
    set.seed(1234)
    gsea <- fgseaMultilevel(clst, stats= all_genes, nPermSimple= 100000, minSize= 6, scoreType= 'std', nproc= 24, eps= 0)
    gsea[, trait := cf]
    gseaOut[[length(gseaOut) + 1]] <- gsea
}}
gsea <- rbindlist(gseaOut)
setnames(gsea, 'pathway', 'Cluster')
gsea[, Cluster := as.numeric(Cluster)]

lEdge <- rep(NA, nrow(gsea))
for(i in 1:nrow(gsea)) {{
    lEdge[i] <- paste(gsea[i,]$leadingEdge[[1]], collapse= ',')
}}
gsea[, leadingEdge := lEdge]
gsea <- merge(gsea, cluster_class, by='Cluster')

gseaOut <- merge(clusters[, list(Cluster, gene_id, description)], gsea, by= 'Cluster', allow.cartesian=TRUE)
gseaOut <- gseaOut[order(trait, pval, Cluster, description)]
write.table(gseaOut, '{output.gsea_long}', row.names= FALSE, quote= FALSE, sep= '\t')

# traits <- c('log2(Exflagellation/ml)', '% infected mosquitoes', 'Oocysts per mosquito)', 'deformability')
# stopifnot(traits %in% gsea$trait)

keep <- gsea[padj < 0.01]$Cluster
tgsea <- gsea[Cluster %in% keep]
tgsea[, Cluster := sprintf('c%s %s', Cluster, Consensus)]
xord <- tgsea[, .SD[which.max(abs(NES))], Cluster][order(-NES)]$Cluster
tgsea[, Cluster := factor(Cluster, xord)]
tgsea[, stars := stars(padj)]
tgsea[, fill := ifelse(padj < 0.01, 'grey30', 
    ifelse(padj < 0.05, 'grey60', 'grey80'))]
tgsea[, trait := gsub('_', ' ', trait)]

nclst <- length(unique(tgsea$Cluster))
nc <- 3
gg <- ggplot(tgsea, aes(x= trait, y= NES, fill= fill, label= stars)) +
    geom_col() +
    geom_text() +
    scale_fill_identity() +
    geom_hline(yintercept= 0, colour= 'orange', linetype= 'dashed') +
    coord_flip() +
    facet_wrap(~Cluster, ncol= nc) +
    xlab(NULL) +
    ylab('Normalised enrichment score') +
    ggtitle('Clusters enriched in differentialy expressed genes') +
    theme_minimal() +
    theme(strip.text= element_text(colour= 'black', margin = margin(0.1, 0, 0, 0, "cm")), strip.background= element_rect(fill= 'white', colour= 'white')) 
ggsave('{output.plot}', width= 16, height= 1 + (nclst/nc) * 2.75, units= 'cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule logfc_dge_by_cluster:
    input:
        dge= 'deseq/infect_dge.tsv.gz',
        gsea= 'gsea/gsea_clusters.tsv',
    output:
        logfc= os.path.join(workflow.basedir, 'results/logfc_by_cluster.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(scales)

source('{workflow.basedir}/scripts/utils.R')

dge <- fread('{input.dge}')
gsea <- fread('{input.gsea}')

gsea[, Cluster := as.character(Cluster)]

dge[, logfc_zscore := scale(log2FoldChange), by= trait]

# Include clusters with smallish enrichment pvalue
pvalue <- 0.01
keep_clst <- unique(gsea[padj < pvalue, Cluster])
keep_gsea <- unique(gsea[Cluster %in% keep_clst, list(Cluster, gene_id)])

dge <- merge(dge, keep_gsea[, list(Cluster, gene_id)], by= 'gene_id')

xord <- gsea[padj < pvalue, list(padj= min(padj)), by= Cluster][order(padj)]$Cluster
dge[, Cluster := factor(Cluster, xord)]

#traits <- c('log2(Exflagellation/ml)', '% infected mosquitoes', 'log2(oocysts/mosquito)', 'deformability')
#stopifnot(traits %in% dge$trait)
#dge[, trait := factor(trait, traits)]

nclst <- length(unique(dge$Cluster))
nc <- 3
gg <- ggplot(data= dge, aes(y= logfc_zscore, x= trait, group= trait)) +
    geom_boxplot(outlier.shape= NA, width= 0.5, fill= NA, colour= 'grey30') +
    geom_quasirandom(width= 0.2, alpha= 0.8, size= 0.35, stroke= 0) +
    geom_hline(yintercept= 0, linetype= 'dashed', colour= 'blue') +
    facet_wrap(~Cluster, ncol= nc) +
    coord_trans_flip(y=symlog_trans(base=2)) +
    xlab('') +
    ylab('z-score of log2(fold-change)') +
    ggtitle(sprintf('Differential expression of genes in clusters with GSEA p-value < %s in at least one trait', pvalue)) +
    theme_minimal() +
    theme(strip.text= element_text(colour= 'black', margin = margin(0.1, 0, 0, 0, "cm")), strip.background= element_rect(fill= 'white', colour= 'white'),
        plot.title= element_text(size= 8)) 
ggsave('{output.logfc}', width= 20, height= (nclst/nc) * 3, units= 'cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule cluster_profile_gsea_voss:
    input:
        vossPcl= config['vossPcl'],
        gsea= 'gsea/gsea_clusters.tsv',
    output:
        clst= os.path.join(workflow.basedir, 'results/voss_enriched_clustered.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)
library(ggplot2)
library(dendextend)

voss_pcl_f <- '{input.vossPcl}'
gsea_f <- '{input.gsea}'

voss_pcl <- fread(voss_pcl_f)
voss_pcl[, name := NULL]
voss_pcl[, GWeight := NULL]
setnames(voss_pcl, 'geneid', 'gene_id')
voss_pcl <- melt(data= voss_pcl, id.vars= 'gene_id', value.name= 'gex', variable.name= 'Time')
voss_pcl[, zscore := scale(gex), by= Time]

gsea <- fread(gsea_f)

gsea <- merge(gsea[, list(gene_id, Cluster, fdr_gsea= padj, trait)], voss_pcl, by= 'gene_id')

gsea <- gsea[fdr_gsea < 0.01, list(trait= paste(trait, collapse= '\n')), by= list(gene_id, Cluster, Time, zscore)]
gsea[, trait := paste(Cluster, trait, sep= ' | ')]
gsea[, Time := as.numeric(sub('day', '', Time))]
clst_avg <- gsea[, list(zscore= mean(zscore)), by= list(Cluster, trait, Time)]

hh <- hclust(dist(as.matrix(dcast(clst_avg, trait ~ Time, value.var= 'zscore'), rownames= 'trait')))
hh <- ladderize(as.dendrogram(hh))
xord <- labels(hh)
gsea[, trait := factor(trait, levels= xord)]
clst_avg[, trait := factor(trait, levels= xord)]

gg <- ggplot(data= gsea, aes(x= Time, y= zscore, group= gene_id)) +
    geom_line(color= 'grey80') +
    geom_line(data= clst_avg, colour= 'red', aes(group= NULL)) +
    facet_wrap(~trait) +
    scale_x_continuous(breaks= unique(gsea$Time)) +
    xlab('Time (day) from Voss dataset') +
    ggtitle('Clusters enriched in genes affected by the given traits') +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'))
ggsave('{output.clst}', width= 24, height= 18, units= 'cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


# rule cluster_pct_dge_vs_pvalue:
#     input:
#         clst= os.path.join(workflow.basedir, 'results/cluster_dge_enrichment.tsv'),
#     output:
#         plot= os.path.join(workflow.basedir, 'results/cluster_dge_enrichment.pdf'),
#     shell:
#         r"""
# cat <<'EOF' > {rule}.$$.tmp.R
# 
# library(ggplot2)
# library(data.table)
# library(ggrepel)
# 
# clst <- fread('{input.clst}')
# clst <- clst[direction_change != 'either']
# fdr <- 0.01
# clst[, colour := ifelse(fisher.fdr < fdr, 'firebrick4', 'grey60')]
# 
# gg <- ggplot(data= clst, aes(x= pct_genes_de, y= -log10(fisher.pvalue), colour= colour)) +
#     geom_point(size= 0.5) +
#     geom_text_repel(data= clst[fisher.fdr < fdr], aes(label= sprintf('%s-%s', Cluster, direction_change))) + 
#     scale_color_identity() +
#     facet_wrap(~trait) +
#     xlab('% genes in cluster differentially expressed') +
#     ylab('Enrichment\n-log10(pvalue)') +
#     theme_light() +
#     theme(legend.position="none", strip.text= element_text(colour= 'black', size= 12))
# ggsave('{output.plot}', height= 14, width= 14, units= 'cm')
# 
# EOF
# Rscript {rule}.$$.tmp.R
# rm {rule}.$$.tmp.R
#         """

rule order_of_peaks:
    input:
        peak= lambda wc: config['peakAll'] if wc.dataset == 'allDatasets' else config['peakVoss'] if wc.dataset == 'vossOnly' else None,
    output:
        peak= 'array_peak/{dataset}.peaksofMeanAll.tsv',
    run:
        peak = pandas.read_csv(input.peak, sep= '\t')
        pmax = peak[peak.MinMax == 'Max']
         
        for dataset in ['Llinas', 'Voss']:
            if dataset == 'Llinas':
                prefix = 'Timepoint.'
            elif dataset == 'Voss':
                prefix = 'day'
            else:
                raise Exception()

            pmax= pmax.assign(time = [int(x.replace(prefix, '')) for x in pmax[f'{dataset}.Timepoint']])
            pmax.sort_values(['time', f'{dataset}.Value'], ascending= [True, False], inplace= True)
            pmax = pmax.assign(ClusterOrder= range(1, len(pmax) + 1))
            pmax.rename(columns= {'ClusterOrder': dataset + '.ClusterOrder'}, inplace= True)
            pmax.drop('time', axis= 1, inplace= True)

        peakOrder = pandas.merge(peak, pmax[['Cluster', 'Llinas.ClusterOrder', 'Voss.ClusterOrder']], how= 'outer', on= ['Cluster'])
        peakOrder.to_csv(output.peak, sep= '\t', index= False)


rule mcl_clustering:
    input:
        mat= 'mcl/gene_expr.mat',
    output:
        clst= temp('mcl/clusters.txt'),
    params:
        corr= 0.56
    shell:
        r"""
mcxarray -data {input.mat} -skipr 1 -skipc 1 -o mcl/expr.mcx --write-binary --pearson -co 0.2 -tf 'abs()'

# You want the "Percentage of nodes that are singletons" (S) and "Node degree"
# median (NDmed) to be low. Since they are anti-correlated, find a good
# compromise.
mcx query -imx mcl/expr.mcx --vary-correlation

# knn(k): with small k you get many small clusters. With big k you get few big
# clusters. Bigger k converges towards not using knn at all.

mcxarray -data {input.mat} -skipr 1 -skipc 1 -o mcl/corr.mci -write-tab mcl/expr.tab --pearson -co {params.corr} -tf 'abs(),add(-{params.corr})'
mcl mcl/corr.mci -o mcl/out.expr.mci -I 2 -tf '#knn(50)'
mcxdump -icl mcl/out.expr.mci -tabr mcl/expr.tab --dump-pairs -o {output.clst}

rm mcl/corr.mci mcl/expr.tab mcl/expr.mcx mcl/out.expr.mci
        """

rule gene_clusters_vs_infection:
    input:
        rna_ss= config['rna_ss'],
        mat= 'mcl/gene_expr.mat',
        mcl= 'mcl/clusters.txt',
        gdesc= 'ref/Pfalciparum3D7.genes.tsv',
    output:
        clst_smry= os.path.join(workflow.basedir, 'results/mcl_clusters.tsv'),
        corr= os.path.join(workflow.basedir, 'results/gene_clusters_vs_infection_corr.tsv'),
        mcl= 'mcl/clusters.tsv'
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

ss <- fread('{input.rna_ss}')
ss <- ss[strain == 'NF54']
gdesc <- fread('{input.gdesc}')
mat <- fread('{input.mat}')
mcl <- fread('{input.mcl}', header= FALSE, col.names= c('cluster', 'gene_id'))
mcl <- merge(mcl, gdesc, by= 'gene_id')
write.table(mcl, '{output.mcl}', row.names= FALSE, sep= '\t', quote= FALSE)

clst <- mcl[, list(cluster_size= .N), cluster][, list(n_clusters= .N), cluster_size]
clst <- clst[order(-cluster_size)]
clst[, cumsum_genes := cumsum(cluster_size * n_clusters)]
clst[, cum_pct := 100 * cumsum_genes / sum(cluster_size * n_clusters)]
clst[, cum_pct := sprintf('%.1f', cum_pct)]
write.table(clst, '{output.clst_smry}', row.names= FALSE, sep= '\t', quote= FALSE)

matLong <- melt(data= mat, id.vars= 'gene_id', variable.name= 'sample_ID', value.name= 'gex')
matLong <- merge(matLong, mcl, by= 'gene_id')

avgExpr <- matLong[, list(gex= mean(gex), .N), by= list(sample_ID, cluster)]
avgExpr <- avgExpr[N >= 5]
avgExpr <- merge(avgExpr, ss[, list(sample_ID, log2_exfl_XA_per_ml= log2(exfl_XA_per_ml + 1), inf_mosq_percent= inf_mosq_percent, log2_oocysts_per_mosq= log2(oocysts_per_mosq + 1))], by= 'sample_ID')
avgExpr <- melt(avgExpr, id.vars= c('sample_ID', 'gex', 'cluster', 'N'), variable.name= 'trait', value.name= 'value')

cor_list <- function(x, y) {{
    m <- cbind(x, y)
    m <- m[complete.cases(m),]
    r <- cor.test(m[,1], m[,2], method= 'pearson')
    out <- list(cor= r$estimate, p.value= r$p.value, n_samples= nrow(m))
    return(out)
}}

corr <- avgExpr[, cor_list(gex, value), by= list(trait, cluster)]
corr[, fdr := p.adjust(p.value, method= 'fdr'), by= trait]
sizes <- mcl[, list(cluster_size= .N), cluster]
corr <- merge(corr, sizes[, list(cluster, cluster_size)], by= 'cluster')
corr <- corr[fdr < 0.01][order(fdr)]
corr[, cor := sprintf('%.2f', cor)]
corr[, p.value := sprintf('%.2e', p.value)]
corr[, fdr := sprintf('%.2e', fdr)]
corr <- corr[, list(trait, cluster, n_samples, cluster_size, p.value, fdr)]
write.table(corr, '{output.corr}', row.names= FALSE, sep= '\t', quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule annotate_clusters:
    input:
        kc_clst= config['clstData'],
        kc_class= os.path.join(workflow.basedir, 'ref_data/crouch_cluster_classification.tsv'),
        pelle= config['pelle'],
        brancucci= os.path.join(workflow.basedir, 'ref_data/Brancucci.sexual_commitment_markers.tsv'),
        josling= os.path.join(workflow.basedir, 'ref_data/Josling.TableS4.tsv'),
        mk= os.path.join(workflow.basedir, 'ref_data/Meerstein_Kessel.TableS4.tsv'),
        dantzler= os.path.join(workflow.basedir, 'ref_data/Dantzler.tsv'),
        miao= os.path.join(workflow.basedir, 'ref_data/miao.male_female_genes.tableS5.tsv'),
    output:
        geneData= temp('cluster_annotation/geneData.tsv'),
        enrich= 'cluster_annotation/cluster_enrichment.tsv',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

kc_clst <- fread('{input.kc_clst}')
kc_class <- fread('{input.kc_class}')
pelle <- fread('{input.pelle}')
brancucci <- fread('{input.brancucci}')
josling <- fread('{input.josling}')
mk <- fread('{input.mk}')
dantzler <- fread('{input.dantzler}')
miao <- fread('{input.miao}')

setnames(kc_clst, c('V1', 'V2'), c('gene_id', 'description'))

pelle <- pelle[!is.na(stage) & !is.na(plasmodb_id), list(
    stage= paste(sort(unique(stage)), collapse= ','), 
    circulating_pvalue= paste(sort(unique(circulating_pvalue)), collapse= ','), 
    sequestered_pvalue= paste(sort(unique(sequestered_pvalue)), collapse= ',')), 
    by= plasmodb_id]

ng <- nrow(kc_clst)
clstData <- merge(kc_clst, kc_class[, list(Cluster, Classification= Consensus)], by= 'Cluster', all.x= TRUE, sort= FALSE)
stopifnot(nrow(clstData) == ng)

clstData <- merge(clstData, pelle[, list(plasmodb_id, pelle= stage, circulating_pvalue, sequestered_pvalue)], by.x= 'gene_id', by.y= 'plasmodb_id', all.x= TRUE, sort= FALSE)
stopifnot(nrow(clstData) == ng)

brancucci[, brancucci := 'commit_marker']
clstData <- merge(clstData, brancucci[, list(gene_id, brancucci)], by= 'gene_id', all.x= TRUE)
stopifnot(nrow(clstData) == ng)

josling[, josling := 'gam_ap2g_target']
clstData <- merge(clstData, unique(josling[, list(gene_id, josling)]), by= 'gene_id', all.x= TRUE)
stopifnot(nrow(clstData) == ng)

clstData <- merge(clstData, mk[, list(gene_id, meerstein_kessel= set_label)], by= 'gene_id', all.x= TRUE)
stopifnot(nrow(clstData) == ng)

clstData <- merge(clstData, dantzler[, list(gene_id, dantzler= stage)], by= 'gene_id', all.x= TRUE)
stopifnot(nrow(clstData) == ng)

clstData <- merge(clstData, miao[, list(gene_id, miao_sex= sex)], by= 'gene_id', all.x= TRUE)
stopifnot(nrow(clstData) == ng)

setcolorder(clstData, c('gene_id', 'Cluster', 'Classification'))
clstData <- clstData[order(Cluster, Classification, gene_id)]

write.table(clstData, '{output.geneData}', row.names= FALSE, sep= '\t', quote= FALSE)

ksize <- clstData[, list(k_size= .N), by= Cluster]
ann <- melt(data= clstData[, list(gene_id, Cluster, Classification, pelle, brancucci, josling, meerstein_kessel, dantzler, miao_sex)], id.vars= c('gene_id', 'Cluster', 'Classification'), value.name= 'annotation', variable.name= 'source')
ann_size <- ann[!is.na(annotation), list(annotation_size= .N), by= list(source, annotation)]
ann <- ann[, list(n_hits= .N), by= list(Cluster, Classification, source, annotation)][!is.na(annotation)]
ann <- merge(ann, ksize, by= c('Cluster'), all.y= TRUE)
ann <- merge(ann, ann_size, by= c('source', 'annotation'), all= TRUE)[order(Cluster, source, annotation)]
ann[, pct_hits := 100 * n_hits/k_size]

tot <- sum(unique(ann[, list(Cluster, k_size)])$k_size)
f.test <- rep(NA, nrow(ann))
for(i in 1:nrow(ann)) {{
    d <- ann[i,]
    m <- matrix(c(d$n_hits, d$annotation_size - d$n_hits,
                  d$k_size - d$n_hits, tot - (d$k_size - d$n_hits) - (d$annotation_size - d$n_hits) - d$n_hits), nrow= 2, byrow= TRUE)
    stopifnot(sum(m) == tot)
    ft <- fisher.test(m, alternative= 'greater')
    f.test[i] <- ft$p.value
}}
ann[, p.value := f.test]
ann[, fdr := p.adjust(p.value, method= 'fdr')]
ann <- ann[order(Cluster, fdr), list(Cluster, Classification, source, annotation, k_size, n_hits, annotation_size, pct_hits, p.value, fdr)]

write.table(ann, '{output.enrich}', row.names= FALSE, sep= '\t', quote= FALSE)

# Crouch vs Pelle supergroups
cmp <- clstData[, list(n_common= .N), by= list(Classification, pelle)]
kc <- clstData[, list(n_crouch= .N), Classification]
pl <- clstData[, list(n_pelle= .N), pelle]
cmp <- merge(cmp, kc, by= 'Classification')
cmp <- merge(cmp, pl, by= 'pelle')
setnames(cmp, c('Classification', 'pelle'), c('crouch_class', 'pelle_class'))
setcolorder(cmp, c('crouch_class', 'pelle_class', 'n_crouch', 'n_pelle', 'n_common'))
cmp[order(crouch_class, -n_common)][!grepl(',', pelle_class)]

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
    """

rule cluster_table_to_wide_format:
    input:
        kc_clst= config['clstData'],
        tlong= 'cluster_annotation/cluster_enrichment.tsv',
        peakAll= config['peakAll'],
    output:
        wide= 'cluster_annotation/cluster_enrichment_wide.tsv',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

kc_clst <- fread('{input.kc_clst}', select= c('V1', 'Cluster'), col.names= c('gene_id', 'Cluster'))
peakAll <- fread('{input.peakAll}', select= c('MinMax', 'Llinas.Timepoint', 'Voss.Timepoint', 'Cluster'))
tlong <- fread('{input.tlong}')

genes <- kc_clst[, list(gene_id= paste(sort(gene_id), collapse= ',')), by= Cluster]

keep <- c(
    'brancucci.commit_marker',
    'dantzler.Asexual',
    'dantzler.Gam',
    'dantzler.Gam_Mo',
    'dantzler.Mo',
    'dantzler.Shared',
    'josling.gam_ap2g_target',
    'meerstein_kessel.gametocyte',
    'meerstein_kessel.rest',
    'pelle.circulating',
    'pelle.gam_ring',
    'pelle.imm_gam',
    'pelle.mat_gam',
    'pelle.sequestered',
    'miao_sex.Male',
    'miao_sex.Female'
)

tlong[, source_ann := paste(source, gsub(',|/', '_', annotation), sep= '.')]

peakMax <- peakAll[MinMax == 'Max']
peakMax[, MinMax := NULL]
setnames(peakMax, c('Llinas.Timepoint', 'Voss.Timepoint'), c('Llinas.Timepoint.Max', 'Voss.Timepoint.Max'))

ann <- tlong[source_ann %in% keep]
ann <- merge(ann, peakMax, by= 'Cluster', all.x= TRUE)

wide <- dcast(data= ann, Cluster + Classification + k_size + Llinas.Timepoint.Max + Voss.Timepoint.Max ~ source_ann, value.var= c('p.value', 'fdr'), sep= '.')

wide <- merge(wide, genes, by= 'Cluster')

write.table(wide, '{output.wide}', row.names= FALSE, sep= '\t', quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule add_deformability_to_genedata:
    input:
        geneData= 'cluster_annotation/geneData.tsv',
        deformability= config['deformability'],
    output:
        geneData= 'cluster_annotation/geneDataDeformability.tsv.gz',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

defo <- fread('{input.deformability}')
geneData <- fread('{input.geneData}')

setnames(defo, 'V1', 'gene_id')
defo[, gene_id := sub('^rna_', '', gene_id)]
defo[, gene_id := sub('-\\d+$', '', gene_id)]
defo[, gene_id := sub('\\.\\d+$', '', gene_id)]

keep <- defo[, list(baseMean= max(baseMean)), by= gene_id]
defo <- merge(defo, keep)
stopifnot(length(unique(defo$gene_id)) == nrow(defo))

geneDefo <- merge(geneData, defo, by= 'gene_id', sort= FALSE, all.x= TRUE)

gz <- gzfile('{output.geneData}', 'w')
write.table(geneDefo, gz, row.names= FALSE, sep= '\t', quote= FALSE)
close(gz)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule add_infectiousness_to_genedata:
    input:
        geneData= 'cluster_annotation/geneData.tsv',
        infect= 'deseq/infect_dge.tsv.gz',
    output:
        geneData= 'cluster_annotation/geneDataInfectiousness.tsv.gz',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

infect <- fread('{input.infect}')
geneData <- fread('{input.geneData}')

geneInfect <- merge(geneData, infect[, list(gene_id, log2FoldChange, baseMean, pvalue, padj, trait)], by= 'gene_id', sort= FALSE, all.x= TRUE)

gz <- gzfile('{output.geneData}', 'w')
write.table(geneInfect, gz, row.names= FALSE, sep= '\t', quote= FALSE)
close(gz)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule infection_markers:
    input:
        dge= 'deseq/infect_dge.tsv.gz',
    output:
        markers= 'deseq/infection_markers.tsv.gz',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(metap)

dge <- fread('{input.dge}')

xdge <- as.matrix(dcast(dge, gene_id ~ trait, value.var= 'pvalue'), rownames= 'gene_id')
p <- apply(xdge, 1, function(x) {{sumz(x)$p}})
max_pvalue <- apply(xdge, 1, max)
xdge <- as.data.table(xdge, keep.rownames= 'gene_id')
xdge[, combined_pvalue := p]
xdge[, adj_combined_pvalue := p.adjust(combined_pvalue)]
xdge[, max_pvalue := max_pvalue]

ydge <- copy(dge)
ydge[, trait := sprintf('%s.log2FoldChange', trait)]
ydge <- dcast(ydge, gene_id ~ trait, value.var= 'log2FoldChange')

mrg <- merge(xdge[, list(gene_id, combined_pvalue, adj_combined_pvalue, max_pvalue)], ydge, by= 'gene_id')
mrg <- merge(mrg, unique(dge[, list(gene_id, description)]), by= 'gene_id')

gz <- gzfile('{output.markers}', 'w')
write.table(mrg, gz, row.names= FALSE, sep= '\t', quote= FALSE)
close(gz)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule voss_vs_infectivity_all_vs_all_timepoints:
    input:
        rna_ss= config['rna_ss'],
        vossPcl= config['vossPcl'],
        infect= 'deseq/logcpm.tsv.gz',
    output:
        line= os.path.join(workflow.basedir, 'results/voss_vs_infectivity_all_vs_all.pdf'),
        voss_top= os.path.join(workflow.basedir, 'results/infectivity_vs_top_voss_corr.pdf'),
        infect_top= os.path.join(workflow.basedir, 'results/voss_vs_infectivity_top_corr.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)
library(ggrepel)

ss <- fread('{input.rna_ss}')
infect <- fread('{input.infect}')
voss <- fread('{input.vossPcl}')

setnames(ss, 'sample_ID', 'sample_id')
stopifnot(ss[!is.na(inf_mosq_percent)]$sample_id == ss[!is.na(oocysts_per_mosq)]$sample_id)
ss[, inf_mosq_percent := NULL]

ss <- ss[, list(sample_id, Day, 
    `Exflagellation/ml`= as.character(exfl_XA_per_ml), 
    `Mosq infection`= as.character(oocysts_per_mosq), 
    `Deformability`= as.character(deformability))]

stopifnot(infect$sample_id %in% ss$sample_id)
infect <- merge(infect, ss, by= 'sample_id')
setnames(infect, 'logcpm', 'infect_gex')
infect[, Day := suppressWarnings(as.numeric(Day))]
infect[, Day := ifelse(is.na(Day), -1, Day)]

voss[, name := NULL]
voss[, GWeight := NULL]
setnames(voss, 'geneid', 'gene_id')
voss <- melt(voss, id.vars= 'gene_id', variable.name= 'Day', value.name= 'voss_gex')
voss[, Day := sub('day', '', Day)]
voss[, Day := as.numeric(Day)]

gex <- merge(infect, voss, by= c('gene_id'), suffixes= c('_infect', '_voss'), allow.cartesian= TRUE)
gex <- melt(data= gex, id.vars= c('gene_id', 'sample_id', 'infect_gex', 'voss_gex', 'Day_infect', 'Day_voss'), variable.name= 'trait', value.name= 'value')[!is.na(value)]
gex[, value := NULL]
gex_corr <- gex[, list(cor= cor(infect_gex, voss_gex, use= 'complete.obs'), .N), by= list(sample_id, Day_infect, Day_voss, trait)]

gex_avg <- gex_corr[, list(cor= mean(cor), sd= sd(cor)), by= list(trait, Day_infect, Day_voss)]
gex_avg[, Day_infectF := factor(Day_infect, levels= sort(unique(Day_infect)))]
gex_avg[, Day_vossF := factor(Day_voss, levels= sort(unique(Day_voss)))]

gg <- ggplot(data= gex_avg, aes(x= Day_voss, y= cor, group= Day_infect, colour= as.character(Day_infect))) +
    geom_hline(yintercept= 0, colour= 'black', size= 0.2, linetype= 'dotted') +
    geom_segment(aes(y= cor + sd, yend= cor - sd, xend= Day_voss), arrow= arrow(length = unit(0.1,"cm"), angle= 90, ends= 'both'), size= 0.2) +
    geom_line() +
    geom_point(cex= 0.2) +
    geom_text_repel(data= gex_avg[Day_voss == max(Day_voss)], aes(label= Day_infect), hjust= 1, xlim= c(max(gex_avg$Day_voss) * 1.01, NA), segment.color= 'grey', segment.size= 0.1) +
    facet_wrap(~trait, ncol= 1) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.2), add = 0), breaks= unique(gex_avg$Day_voss)) +
    xlab('Day (Voss)') +
    ylab('Mean correlation +/- 1sd') +
    ggtitle('Correlation between each Voss and each infection time points') +
    theme_light() +
    theme(legend.position= 'none', strip.text= element_text(colour= 'black'))
ggsave('{output.line}', width= 18, height= 24, units= 'cm')

# gg <- ggplot(data= gex_avg, aes(x= Day_voss, y= Day_infect)) +
#    geom_tile(aes(fill= cor)) +
#    scale_x_continuous(breaks= unique(gex_avg$Day_voss)) +
#    facet_wrap(~trait, ncol= 1, scale= 'free') +
#    xlab('Day (Voss)') +
#    ylab('Day (infection)') +
#    ggtitle('Correlation between each Voss and each infection time points') +
#    theme_light() +
#    theme(strip.text= element_text(colour= 'black'))

peak_time <- gex_avg[Day_infect > -1, .SD[which.max(cor)], list(trait, Day_infect)]

gg <- ggplot(data= peak_time, aes(x= Day_infect, y= cor)) +
    geom_segment(aes(y= cor + sd, yend= cor - sd, xend= Day_infect), arrow= arrow(length = unit(0.1,"cm"), angle= 90, ends= 'both'), size= 0.2) +
    geom_line() +
    geom_point(cex= 0.2) +
    geom_label(aes(label= Day_voss, colour= Day_vossF), vjust= -0.2, hjust= 1.2) +
    facet_wrap(~trait, ncol= 1) +
    scale_x_continuous(breaks= unique(peak_time$Day_infect), minor_breaks= NULL) +
    ggtitle('Correlation between each day of infection and the most correlated day from Voss') +
    ylab('Correlation with Voss') +
    xlab('Day (infection)') +
    theme_light() +
    theme(legend.position= 'none', strip.text= element_text(colour= 'black', size= 12))
ggsave('{output.voss_top}', width= 18, height= 24, units= 'cm')

peak_time <- gex_avg[Day_infect > -1, .SD[which.max(cor)], list(trait, Day_voss)]

gg <- ggplot(data= peak_time, aes(x= Day_voss, y= cor)) +
    geom_segment(aes(y= cor + sd, yend= cor - sd, xend= Day_voss), arrow= arrow(length = unit(0.1,"cm"), angle= 90, ends= 'both'), size= 0.2) +
    geom_line() +
    geom_point(cex= 0.2) +
    geom_label(aes(label= Day_infect, colour= Day_infectF), vjust= -0.2, hjust= 1.2) +
    facet_wrap(~trait, ncol= 1) +
    scale_x_continuous(breaks= unique(peak_time$Day_voss), minor_breaks= NULL) +
    ggtitle('Correlation between each Voss day and the most correlated infection day') +
    ylab('Correlation with infection') +
    xlab('Day (Voss)') +
    theme_light() +
    theme(legend.position= 'none', strip.text= element_text(colour= 'black', size= 12))
ggsave('{output.infect_top}', width= 18, height= 24, units= 'cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule gex_voss_vs_infectivity:
    input:
        rna_ss= config['rna_ss'],
        vossPcl= config['vossPcl'],
        infect= 'deseq/logcpm.tsv.gz',
    output:
        corr= os.path.join(workflow.basedir, 'results/voss_vs_infectivity_corr.pdf'),
        tp_corr= os.path.join(workflow.basedir, 'results/voss_vs_infectivity_timepoint_corr.pdf'),
        peak= os.path.join(workflow.basedir, 'results/voss_vs_infectivity_peak.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)

ss <- fread('{input.rna_ss}')
infect <- fread('{input.infect}')
voss <- fread('{input.vossPcl}')

setnames(ss, 'sample_ID', 'sample_id')
stopifnot(ss[!is.na(inf_mosq_percent)]$sample_id == ss[!is.na(oocysts_per_mosq)]$sample_id)
ss[, inf_mosq_percent := NULL]

ss <- ss[, list(sample_id, Day, 
    `Exflagellation/ml`= as.character(exfl_XA_per_ml), 
    `Mosq infection`= as.character(oocysts_per_mosq), 
    `Deformability`= as.character(deformability))]

stopifnot(infect$sample_id %in% ss$sample_id)
infect <- merge(infect, ss, by= 'sample_id')
setnames(infect, 'logcpm', 'infect_gex')

voss[, name := NULL]
voss[, GWeight := NULL]
setnames(voss, 'geneid', 'gene_id')
voss <- melt(voss, id.vars= 'gene_id', variable.name= 'Day', value.name= 'voss_gex')
voss[, Day := sub('day', '', Day)]

gex <- merge(infect, voss, by= c('gene_id', 'Day'))
gex <- melt(gex, id.vars= c('gene_id', 'Day', 'sample_id', 'infect_gex', 'voss_gex'), variable.name= 'trait', value.name= 'value')[!is.na(value)]
# gex <- gex[, list(infect_gex= mean(infect_gex)), by= list(Day, gene_id, voss_gex, trait)]
gex[, Day := as.integer(Day)]

cor_test <- function(x, y) {{
    cc <- cor.test(x, y, use= 'complete.obs')
    return(list(cor= cc$estimate, p.value= cc$p.value))
}}
tp_gex <- gex[, list(cor= cor(voss_gex, infect_gex, use= 'complete.obs', method= 'spearman')), by= list(trait, Day, sample_id)]
tp_avg <- tp_gex[, list(cor= mean(cor)), by= list(trait, Day)]

tp_gex[, xDay := Day + as.numeric(as.factor(trait))/10 - 0.2]

gg <- ggplot(data= tp_gex, aes(x= Day, y= cor, by= trait, colour= trait, shape= trait)) +
    geom_line(data= tp_avg) +
    geom_point(size= 1.25, aes(xDay)) +
    geom_text_repel(data= tp_avg[, .SD[which.max(Day), list(Day, cor)], by= trait], aes(label= trait), hjust= 0, xlim= c(max(tp_gex$Day) * 1.01, NA)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.3), add = 0), breaks= unique(tp_gex$Day)) +
    ylab('Gene expression correlation') +
    ggtitle('Correlation between Voss and infection at each time point') +
    theme_light() +
    theme(legend.position= 'none')
ggsave('{output.tp_corr}', width= 16, height= 10, units= 'cm')

corr_gex <- gex[, cor_test(voss_gex, infect_gex), by= list(trait, gene_id)]

xord <- corr_gex[, list(cor= mean(cor)), by= trait][order(cor)]$trait
corr_gex[, trait := factor(trait, xord)]

gg <- ggplot(data= corr_gex, aes(x= trait, y= cor)) +
    geom_hline(yintercept= 0, colour= 'orange', linetype= 'dashed') +
    geom_quasirandom(width= 0.4, pch= '.', colour= 'grey70') +
    geom_violin(fill= NA, draw_quantiles= c(0.25, 0.75), linetype= 'dashed') +
    geom_violin(fill= NA, draw_quantiles= 0.5) +
    ylab('Gene expression correlation') +
    xlab('') +
    ggtitle('Correlation between Voss and infection time courses') +
    theme_light()
ggsave('{output.corr}', width= 13, height= 10, units= 'cm')

peak_gex <- gex[, list(voss_peak= .SD$Day[which.max(voss_gex)], infect_peak= .SD$Day[which.max(infect_gex)]), by= list(trait, gene_id)]
peak_gex[, peak_day_diff := voss_peak - infect_peak]
peak_gex <- peak_gex[, .N, by= list(trait, peak_day_diff)]

xord <- peak_gex[peak_day_diff == 0][order(N)]$trait
peak_gex[, trait := factor(trait, xord)]

gg <- ggplot(data= peak_gex, aes(x= peak_day_diff, y= N)) +
    geom_col() +
    geom_vline(xintercept= 0, colour= 'orange', linetype= 'dashed') +
    facet_wrap(~trait) +
    ylab('N genes') +
    xlab('Days of difference `Voss - Trait`\nNegative means Voss peaks earlier') +
    ggtitle('Difference in peak day of gene expression') +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'))
ggsave('{output.peak}', width= 18, height= 10, units= 'cm')


EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule dge_genes_to_clusters:
    input:
        dge= 'deseq/infect_dge.tsv',
        clstData= config['clstData'],
    output:
        dge_clusters= os.path.join(workflow.basedir, 'results/dge_clusters.tsv'),
        k_summary= os.path.join(workflow.basedir, 'results/cluster_dge_enrichment.tsv'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(fgsea)

dge <- fread('{input.dge}')
clusters <- fread('{input.clstData}')
setnames(clusters, c('V1', 'V2'), c('gene_id', 'description'))

test_clusters_for_dge <- function(kdge, kdge_genes, kclusters) {{
    kdge <- kdge[gene_id %in% kclusters$gene_id]
    kclusters <- kclusters[gene_id %in% kdge$gene_id]
    kdge_genes <- kdge_genes[kdge_genes %in% kdge$gene_id]

    stopifnot(length(unique(kdge_genes)) == length(kdge_genes))
    stopifnot(length(unique(kdge$gene_id)) == length(kdge$gene_id))
    stopifnot(identical(sort(kdge$gene_id), sort(kclusters$gene_id)))
    stopifnot(kdge_genes %in% kdge$gene_id)

    kclusters[, de := ifelse(gene_id %in% kdge_genes, 'in_is_de', 'in_not_de')]
    xtab <- kclusters[, .N, list(Cluster, de)]
    ksize <- kclusters[, list(k_size= .N), by= Cluster]
    xtab <- merge(xtab, ksize, by= 'Cluster')
    xtab <- dcast(data= xtab, Cluster + k_size ~ de, value.var= 'N', fill= 0)
    xtab[, out_is_de := length(kdge_genes) - in_is_de]
    xtab[, out_not_de := length(kclusters$gene_id) - k_size - out_is_de]
    xtab[, tot := in_not_de + in_is_de + out_is_de + out_not_de]
    stopifnot(xtab$tot == length(kclusters$gene_id))
    xtab[, tot := NULL]

    xtab <- xtab[k_size > 5]
    xtab[, fisher.pvalue := NA]
    for(i in 1:length(xtab$Cluster)) {{
        k <- xtab$Cluster[i]
        m <- xtab[Cluster == k, list(in_is_de, in_not_de, out_is_de, out_not_de)]
        m <- matrix(as.numeric(m), nrow= 2, byrow= TRUE)
        ft <- fisher.test(m)
        xtab$fisher.pvalue[i] <- ft$p.value
    }}
    xtab[, fisher.fdr := p.adjust(fisher.pvalue, method= 'fdr')]
    xtab[, pct_genes_de := 100 * (in_is_de / k_size)]

    clusters_dge <- merge(kclusters[, list(gene_id, Cluster, description)], xtab, by= 'Cluster')
    clusters_dge[, is_dge := gene_id %in% kdge_genes]
    clusters_dge <- merge(clusters_dge, kdge[, list(gene_id, dge_fdr= padj, log2FoldChange)], by= 'gene_id')[order(fisher.fdr, Cluster, -is_dge, dge_fdr)]
    return(clusters_dge)
}}

write_clusters <- function(clusters_dge, outfile) {{
    clusters_dge[, log2FoldChange := sprintf('%.2f', log2FoldChange)]
    clusters_dge[, k_stars := stars(fisher.fdr)]
    clusters_dge[, dge_stars := stars(dge_fdr)]
    clusters_dge[, fisher.fdr := sprintf('%.2e %s', fisher.fdr, k_stars)]
    clusters_dge[, dge_fdr := sprintf('%.2e %s', dge_fdr, dge_stars)]
    write.table(clusters_dge[, list(trait, Cluster, gene_id, description, fisher.fdr, log2FoldChange, dge_fdr)], outfile, sep= '\t', row.names= FALSE, quote= FALSE)
}}

stars <- function(p) {{
    xstars <- ifelse(p < 0.001, '***', 
        ifelse(p < 0.01, '**', 
            ifelse(p < 0.05, '+', '')))
    return(xstars)
}}

# ---------------

k_clst <- list()

for(x in c('up', 'down', 'either')) {{
    if(x == 'up') {{
        dge2 <- dge[log2FoldChange > 0]
    }} else if(x == 'down') {{
        dge2 <- dge[log2FoldChange < 0]
    }} else if(x == 'either') {{
        dge2 <- copy(dge)
    }} else {{
        stop()
    }}
    xdge <- dcast(dge2, gene_id ~ trait, value.var= c('padj', 'log2FoldChange'))
    dge_genes <- unique(c(
        xdge[FDR_is_exfl < 0.01 & abs(log2FoldChange_is_exfl) > log2(3)]$gene_id,
        xdge[FDR_is_exfl < 0.01 & FDR_log2_exfl_XA_per_ml_only_pos < 0.01]$gene_id
    ))
    clusters_dge <- test_clusters_for_dge(dge[trait == 'is_exfl'], dge_genes, clusters)
    clusters_dge[, trait := 'Exflagellation']
    clusters_dge[, direction_change := x]
    k_clst[[length(k_clst) + 1]] <- clusters_dge

    tt <- c('% infected mosquitoes', 'is_inf_mosq')
    stopifnot(tt %in% dge2$trait)
    dge_genes <- unique(dge2[trait %in% tt & pvalue < 0.01]$gene_id)
    clusters_dge <- test_clusters_for_dge(dge[trait == '% infected mosquitoes'], dge_genes, clusters)
    clusters_dge[, trait := '% infected mosquitoes']
    clusters_dge[, direction_change := x]
    k_clst[[length(k_clst) + 1]] <- clusters_dge

    tt <- c('log2(oocysts/mosquito)', 'is_oocysts_per_mosq')
    stopifnot(tt %in% dge2$trait)
    dge_genes <- unique(dge2[trait %in% tt & pvalue < 0.01]$gene_id)
    clusters_dge <- test_clusters_for_dge(dge[trait == 'log2(oocysts/mosquito)'], dge_genes, clusters)
    clusters_dge[, trait := 'log2(oocysts/mosquito)']
    clusters_dge[, direction_change := x]
    k_clst[[length(k_clst) + 1]] <- clusters_dge

    dge_genes <- unique(dge2[trait == 'deformability' & padj < 0.01 & abs(log2FoldChange) > log2(1.5)]$gene_id)
    clusters_dge <- test_clusters_for_dge(dge[trait == 'deformability'], dge_genes, clusters)
    clusters_dge[, trait := 'deformability']
    clusters_dge[, direction_change := x]
    k_clst[[length(k_clst) + 1]] <- clusters_dge
}}
k_clst <- rbindlist(k_clst)

write_clusters(k_clst[fisher.fdr < 0.1 & !grepl('unknown function', description)][order(trait, fisher.fdr, Cluster, dge_fdr)], '{output.dge_clusters}')

# NB: k_size is the number of genes that are both in the cluster AND in the dge
# table. Therefore k_size may be smaller than the size of the original cluster.
k_summary <- unique(k_clst[, list(Cluster, k_size, in_is_de, in_not_de, out_is_de, out_not_de, fisher.pvalue, fisher.fdr, pct_genes_de, direction_change, trait)])[order(Cluster, trait, direction_change)]
write.table(k_summary, '{output.k_summary}', row.names= FALSE, sep= '\t', quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule plot_clusters_on_voss:
    input:
        vossPcl=config['vossPcl'],
        clstData=config['clstData'],
    output:
        clusters_voss_time_course=os.path.join(workflow.basedir, 'results/clusters_voss_time_course.pdf')
    script:
        os.path.join(workflow.basedir, 'scripts/cluster_profiles.R')

rule plot_network:
    input:
        pclDir=config['pclDir'],
        clstData=config['clstData'],
        annotation='cluster_annotation/cluster_enrichment.tsv',
    output:
        network_clusters=os.path.join(workflow.basedir, 'results/network_clusters.pdf'),
        network_genes=os.path.join(workflow.basedir, 'results/network_genes.pdf'),
    script:
        os.path.join(workflow.basedir, 'scripts/network.R')

rule plot_heatmap_patients:
    input:
        pcl=os.path.join(config['pclDir'], '{dataset}.pcl.gz'),
        clstData=config['clstData'],
        peakAll=config['peakAll'],
        allGenesOldIds=config['allGenesOldIds'],
        annotation='cluster_annotation/cluster_enrichment.tsv',
    output:
        hm=temp('{dataset}_patient_heatmap.pdf'),
        line=temp('{dataset}_patient_lineplot.pdf'),
    script:
        os.path.join(workflow.basedir, 'scripts/patients.R')

rule concat_patient_plots:
    input:
        pdf=expand('{dataset}_patient_{{plottype}}.pdf', dataset=['daily', 'Joiceetaldata', 'milner'])
    output:
        pdf=os.path.join(workflow.basedir, 'results/patient_{plottype}.pdf'),
    shell:
        r"""
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/printer -sOutputFile={output.pdf} {input.pdf}
        """

rule selected_clusters_in_patients:
    input:
        dge='deseq/infect_dge.tsv.gz',
        pcl=expand(os.path.join(config['pclDir'], '%s.pcl.gz' % x) for x in ['Joiceetaldata', 'milner']), # 'daily', 
        gsea='gsea/gsea_clusters.tsv',
        allGenesOldIds=config['allGenesOldIds'],
    output:
        hm=os.path.join(workflow.basedir, 'results/selected_clusters_patients.pdf'),
        patient_order=temp('patient_order.tmp.R'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)

source('{workflow.basedir}/scripts/read_pcl.R')

dge <- fread('{input.dge}')

allGenesOldIds <- '{input.allGenesOldIds}'
gsea <- fread('{input.gsea}', select=c('Cluster', 'gene_id', 'padj', 'trait', 'NES'))
# For plotting, add a cluster consisting of the selected genes

pclfiles <- strsplit('{input.pcl}', ' ')[[1]]
pcl <- list()
for(x in pclfiles) {{
    dat <- read_pcl(x)
    dataset <- sub('\\.pcl$', '', sub('\\.gz$', '', basename(x)))
    if(dataset == "Joiceetaldata") {{
        dataset <- 'Joice et al.'
    }} else if(dataset == "milner") {{
        dataset <- 'Milner et al.'
    }} else if(dataset == "daily") {{
        dataset <- 'Daily et al.'
    }} else {{
        stop(sprintf("Invalid dataset: %s", dataset))
    }}
    dat[, dataset := dataset]
    pcl[[length(pcl) + 1]] <- dat
}}
pcl <- rbindlist(pcl)
# Few genes appear to be duplicates(!?)
pcl <- pcl[, list(gex=mean(gex)), list(gene_id, array_id, dataset)]
pcl[, array_id := sprintf('%s:%s', dataset, array_id)]

A <- dcast(data=pcl, gene_id ~ array_id, value.var='gex')
A <- as.matrix(A[!is.na(gene_id)], rownames='gene_id')
A <- limma::normalizeQuantiles(A)
A <- melt(data=as.data.table(A, keep.rownames='gene_id'), variable.name='array_id', value.name='qqgex', id.vars='gene_id')
pcl <- merge(pcl, A, by=c('gene_id', 'array_id'))
pcl[, zscore := scale(qqgex), by=list(dataset, gene_id)]

keep_clusters <- unique(gsea[padj < 0.01 & trait == 'deformable_vs_mixed', .SD[which.min(padj)], Cluster])[order(NES), Cluster]

keep_genes <- unique(gsea[Cluster %in% keep_clusters & gene_id %in% pcl$gene_id, list(gene_id, Cluster)])
pcl <- merge(pcl, keep_genes, by='gene_id')

keep_dge <- merge(dge[trait == 'deformable_vs_mixed', list(log2FoldChange, gene_id)], keep_genes, by='gene_id', all.y=TRUE)

pcl[, gene_id := sprintf('%s:%s', gene_id, Cluster)]

gene_order <- c()
for(clst in keep_clusters) {{
    gene_order <- c(gene_order, keep_dge[Cluster == clst][order(-log2FoldChange), list(gene_id=sprintf('%s:%s', gene_id, Cluster))]$gene_id)
}}

patient_order <- c()
for(dt in unique(pcl$dataset)) {{
    mat <- dcast(data=pcl[dataset == dt], array_id ~ gene_id, value.var='zscore')
    mat <- as.matrix(mat, rownames='array_id')
    dd <- ladderize(as.dendrogram(hclust(dist(mat))))
    patient_order <- c(patient_order, labels(dd))
}}
sink('{output.patient_order}')
dput(patient_order)
sink()

mat <- dcast(data=pcl, gene_id ~ array_id, value.var='zscore')
mat <- as.matrix(mat, rownames='gene_id')
mat <- mat[match(gene_order, rownames(mat)), match(patient_order, colnames(mat))]

Clusters <- sub('.*:', '', rownames(mat))
Clusters <- factor(Clusters, rev(unique(Clusters)))

cols <- rep('grey30', length(Clusters))
names(cols) <- Clusters
cls <- HeatmapAnnotation(which='row',
    `Deformability log2FC`=anno_barplot(keep_dge[match(sub(':.*', '', gene_order), gene_id)]$log2FoldChange, border=FALSE, width=unit(2, 'cm'), axis_param=list(labels_rot=0)))

dataset <- sub(':.*', '', colnames(mat))
dataset <- factor(dataset, unique(dataset))
cols <- rep('grey30', length(unique(dataset)))
names(cols) <- unique(dataset)
cls_col <- HeatmapAnnotation(df=data.table(dataset=dataset), show_legend=FALSE, annotation_label='', which='column', col=list(dataset=cols), simple_anno_size = unit(1, "mm"))

pdf('{output.hm}', height=(0.05 * nrow(mat))/2.54, width=45/2.54)
Heatmap(mat, name='Z-score\nof gene\nexpression',
    row_split=Clusters,
    row_labels=rep('', nrow(mat)),
    column_split=dataset,
    column_title=unique(dataset),
    # column_labels=rep('', ncol(mat)), 
    row_title_rot=0,
    row_order=rownames(mat),
    left_annotation=cls,
    top_annotation=cls_col,
    cluster_columns=FALSE, 
    cluster_rows=FALSE)
dev.off()

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule deformability_genes_in_patients:
    input:
        dge='deseq/infect_dge.tsv.gz',
        pcl=expand(os.path.join(config['pclDir'], '%s.pcl.gz' % x) for x in ['daily', 'Joiceetaldata', 'milner']),
    output:
        hm=os.path.join(workflow.basedir, 'results/deformability_genes_in_patients.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)

source('{workflow.basedir}/scripts/read_pcl.R')

dge <- fread('{input.dge}')

pclfiles <- strsplit('{input.pcl}', ' ')[[1]]
pcl <- list()
for(x in pclfiles) {{
    dat <- read_pcl(x)
    dataset <- sub('\\.pcl$', '', sub('\\.gz$', '', basename(x)))
    dat[, dataset := dataset]
    pcl[[length(pcl) + 1]] <- dat
}}
pcl <- rbindlist(pcl)
# Few genes appear to be duplicates(!?)
pcl <- pcl[, list(gex=mean(gex)), list(gene_id, array_id, dataset)]
pcl[, array_id := sprintf('%s:%s', dataset, array_id)]

A <- dcast(data=pcl, gene_id ~ array_id, value.var='gex')
A <- as.matrix(A[!is.na(gene_id)], rownames='gene_id')
A <- limma::normalizeQuantiles(A)
A <- melt(data=as.data.table(A, keep.rownames='gene_id'), variable.name='array_id', value.name='qqgex', id.vars='gene_id')
pcl <- merge(pcl, A, by=c('gene_id', 'array_id'))

pcl[, zscore := scale(qqgex), by=list(dataset, gene_id)]

keep_genes <- dge[gene_id %in% pcl$gene_id & trait == 'deformable_vs_mixed' & padj < 0.01 & abs(log2FoldChange) > 0.53, list(gene_id, log2FoldChange)][order(-log2FoldChange)]
pcl <- merge(pcl, keep_genes, by='gene_id')

mat <- dcast(data=pcl, gene_id ~ array_id, value.var='zscore')
mat <- as.matrix(mat, rownames='gene_id')
dd <- ladderize(as.dendrogram(hclust(dist(mat))))
gene_order <- labels(dd)

patient_order <- c()
for(dt in unique(pcl$dataset)) {{
    mat <- dcast(data=pcl[dataset == dt], array_id ~ gene_id, value.var='zscore')
    mat <- as.matrix(mat, rownames='array_id')
    dd <- ladderize(as.dendrogram(hclust(dist(mat))))
    patient_order <- c(patient_order, labels(dd))
}}

mat <- dcast(data=pcl, gene_id ~ array_id, value.var='zscore')
mat <- as.matrix(mat, rownames='gene_id')
mat <- mat[match(gene_order, rownames(mat)), match(patient_order, colnames(mat))]

cls <- rowAnnotation(
    `deformable_vs_mixed \nlog2FC`=anno_barplot(keep_genes[match(gene_order, gene_id)]$log2FoldChange, border=FALSE, width=unit(2, 'cm'), axis_param=list(labels_rot=0)))

dataset <- sub(':.*', '', colnames(mat))
cols <- rep('grey30', length(unique(dataset)))
names(cols) <- unique(dataset)
cls_col <- HeatmapAnnotation(df=data.table(dataset=dataset), show_legend=FALSE, annotation_label='', which='column', col=list(dataset=cols), simple_anno_size = unit(1, "mm"))

pdf('{output.hm}', height=20/2.54, width=20/2.54)
Heatmap(mat, name='Z-score\nof gene\nexpression', 
    column_split=sub(':.*', '', colnames(mat)),
    column_title=unique(sub(':.*', '', colnames(mat))),
    row_labels=rep('', nrow(mat)),
    column_labels=rep('', ncol(mat)), 
    row_title_rot=0,
    row_order=gene_order,
    left_annotation=cls,
    top_annotation=cls_col,
    cluster_columns=FALSE, 
    cluster_rows=FALSE)
dev.off()

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule voss_time_course_heatmap:
    input:
        gsea='gsea/gsea_clusters.tsv',
        vossPcl= config['vossPcl'],
        clstData= config['clstData'],
        kc_class= os.path.join(workflow.basedir, 'ref_data/crouch_cluster_classification.tsv'),
    output:
        hm=os.path.join(workflow.basedir, 'results/voss_time_course_heatmap.pdf'),
        hm_gam=os.path.join(workflow.basedir, 'results/voss_time_course_heatmap_gam.pdf'),
        hm_shared=os.path.join(workflow.basedir, 'results/voss_time_course_heatmap_shared.pdf'),
        line=os.path.join(workflow.basedir, 'results/voss_time_course_selected.pdf'),
    script:
        os.path.join(workflow.basedir, 'scripts/voss_time_course_heatmap.R')


rule overlapping_dge:
    input:
        dge='deseq/infect_dge.tsv.gz',
    output:
        ovl='deseq/overlapping_dge.tsv.gz',
        upset=os.path.join(workflow.basedir, 'results/upset_dge.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(UpSetR)

dge <- fread('{input.dge}')
dge[, is_dge := abs(log2FoldChange) > log2(1.5) & !is.na(padj) & padj < 0.05]
dge[, is_dge := ifelse(is_dge == TRUE, 1, 0) * sign(log2FoldChange)]
ovl <- dcast(dge, gene_id + description ~ trait, value.var='is_dge')

fwrite(ovl, '{output.ovl}', sep='\t')

ovl[, description := NULL]
ovl<- as.data.table(abs(as.matrix(ovl, rownames='gene_id')), keep.rownames='gene_id')
setnames(ovl, names(ovl), gsub('_', ' ', names(ovl)))

pdf('{output.upset}', width= 24/2.54, height= 12/2.54, onefile= FALSE)
upset(ovl, 
      sets.x.label= 'N. dge genes',
      mainbar.y.label= 'Genes in intersection',
      mb.ratio = c(0.6, 0.4),
      text.scale = c(numbers_above_bars= 1.5, set_names= 2), order.by= c('freq'), decreasing= TRUE)
dev.off()

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule candidate_genes:
    input:
        ref='ref/Pfalciparum3D7.genes.tsv',
        pcl=expand(os.path.join(config['pclDir'], '%s.pcl.gz' % x) for x in ['Joiceetaldata', 'milner']), # 'daily', 
        candidate_genes=config['candidate_genes'],
        patient_order='patient_order.tmp.R',
    output:
        hm=os.path.join(workflow.basedir, 'results/candidate_genes_patients.pdf'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)

source('{workflow.basedir}/scripts/read_pcl.R')
patient_order <- dget('{input.patient_order}')
candidate_genes_genes <- fread('{input.candidate_genes}', select='gene_id')
ref <- fread('{input.ref}')
pclfiles <- strsplit('{input.pcl}', ' ')[[1]]

pcl <- list()
for(x in pclfiles) {{
    dat <- read_pcl(x)
    dataset <- sub('\\.pcl$', '', sub('\\.gz$', '', basename(x)))
    if(dataset == "Joiceetaldata") {{
        dataset <- 'Joice et al.'
    }} else if(dataset == "milner") {{
        dataset <- 'Milner et al.'
    }} else if(dataset == "daily") {{
        dataset <- 'Daily et al.'
    }} else {{
        stop(sprintf("Invalid dataset: %s", dataset))
    }}
    dat[, dataset := dataset]
    pcl[[length(pcl) + 1]] <- dat
}}
pcl <- rbindlist(pcl)
# Few genes appear to be duplicates(!?)
pcl <- pcl[, list(gex=mean(gex)), list(gene_id, array_id, dataset)]
pcl[, array_id := sprintf('%s:%s', dataset, array_id)]

A <- dcast(data=pcl, gene_id ~ array_id, value.var='gex')
A <- as.matrix(A[!is.na(gene_id)], rownames='gene_id')
A <- limma::normalizeQuantiles(A)
A <- melt(data=as.data.table(A, keep.rownames='gene_id'), variable.name='array_id', value.name='qqgex', id.vars='gene_id')
pcl <- merge(pcl, A, by=c('gene_id', 'array_id'))
pcl[, zscore := scale(qqgex), by=list(dataset, gene_id)]

candidate_genes_genes <- candidate_genes_genes[gene_id %in% pcl$gene_id]

pcl <- merge(pcl, candidate_genes_genes, by='gene_id')

mat <- dcast(data=pcl, gene_id ~ array_id, value.var='zscore')
mat <- merge(mat, ref[, list(gene_id, name)], by='gene_id', all.x=TRUE)
mat[, gene_name := ifelse(name == 'N/A', gene_id, name)]
mat[, gene_id := NULL]
mat[, name := NULL]

mat <- as.matrix(mat, rownames='gene_name')
mat <- mat[, match(patient_order, colnames(mat))]

dataset <- sub(':.*', '', colnames(mat))
dataset <- factor(dataset, unique(dataset))
cols <- rep('grey30', length(unique(dataset)))
names(cols) <- unique(dataset)
cls_col <- HeatmapAnnotation(df=data.table(dataset=dataset), show_legend=FALSE, annotation_label='', which='column', col=list(dataset=cols), simple_anno_size = unit(1, "mm"))

pdf('{output.hm}', height=16/2.54, width=45/2.54)
Heatmap(mat, name='Z-score\nof gene\nexpression', 
    column_split=dataset,
    column_title=unique(dataset),
    # column_labels=rep('', ncol(mat)), 
    row_title_rot=0,
    top_annotation=cls_col,
    cluster_columns=FALSE, 
    cluster_rows=TRUE)
dev.off()
EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


## STUB ##
rule overlapping_dge_clusters:
    input:
        ovl='deseq/overlapping_dge.tsv.gz',
        clstData=config['clstData'],
    output:
    shell:
        r"""
library(data.table)

ovl <- fread('{input.ovl}')
clstData <- fread('{input.clstData}')
clstData <- clstData[, list(gene_id=V1, Cluster=as.character(Cluster))]

ovlClst <- merge(ovl, clstData)
ovlClst[, deform_and_exfl := sprintf('%s_%s', abs(deformable_vs_mixed), abs(exfl_detected_vs_not_detected))]

clusters <- ovlClst[, .N, Cluster][N >= 10]$Cluster

for(clst in clusters) {
    dcast(ovlClst[Cluster == clst, .N, deform_and_exfl], )
}
        """
