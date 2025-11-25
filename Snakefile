import glob
import re
import pandas
import tabulate

os.makedirs("slurm", exist_ok=True)

PLASMODB_RELEASE = "55"


rule all:
    input:
        ## Figure 1
        "results/histogram_traits.pdf",
        "results/upset_traits.pdf",
        ## Figure 2
        "results/infection_time_course.deform.pdf",
        "results/infection_time_course.only_infect.pdf",
        "results/trait_pairs.pdf",
        ## Figure 3
        "results/gex_logfc.pdf",
        "results/dge_trait_corr.pdf",
        "gsea/go/infection_deformability.pdf",
        "results/upset_dge.pdf",
        ## Figure 4
        "results/gsea_clusters.pdf",
        ## Figure 5
        "results/selected_clusters_patients.pdf",
        "results/candidate_genes_patients.pdf",
        ## Supp Table and Figures
        "deseq/infect_dge.tsv.gz",
        "deseq/overlapping_dge.tsv.gz",
        "gsea/go/infection_deformability.tsv.gz",
        "results/gex_pca_combatseq.pdf",
        "results/gex_pca_traits.pdf",
        "deseq/overlapping_dge.tsv.gz",
        "results/voss_time_course_heatmap_shared.pdf",
        "results/voss_time_course_heatmap_gam.pdf",
        "results/voss_time_course_selected.pdf",
        'cluster_annotation/cluster_enrichment_wide.tsv',
        "sexratio/markers_sex_defo.pdf",
        'sexratio/markers_sex_infect.pdf',
        'sexratio/markers_sex_ratio_voss.pdf',
        'sexratio/markers_sex_voss.pdf',
        'sexratio/markers_sex_ratio.pdf',


rule microarray_filtering:
    input:
        varGenes=os.path.join(workflow.basedir, "ref_data/microarrays/auxiliary_files/varWithOldIds.txt"), 
        rifinGenes=os.path.join(workflow.basedir, "ref_data/microarrays/auxiliary_files/rifinWithOldIds.txt"),
        stevorGenes=os.path.join(workflow.basedir, "ref_data/microarrays/auxiliary_files/stevorWithOldIds.txt"),
        allGenesOldIds=os.path.join(workflow.basedir, "ref_data/microarrays/auxiliary_files/allGenesOldIds.txt"),

        daily=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/Daily.pcl.gz"),
        leRoux=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/LeRoux.pcl.gz"),
        milner=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/Malawi.pcl.gz"),
        llinas3D7=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/3D7_Llinas.pcl.gz"),
        llinasDD2=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/DD2_Llinas.pcl.gz"),
        llinasHB3=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/HB3_Llinas.pcl.gz"),
        hu=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/Hu.pcl.gz"),
        leRoch=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/LeRoch.pcl.gz"),
        young=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/young2005.pcl.gz"),
        voss=os.path.join(workflow.basedir, "ref_data/microarrays/rawInput/vossMicroarray.txt.gz"),
    output:
        pcl=['microarrays/filteredInput/daily.pcl',
        'microarrays/filteredInput/hu.pcl',
        'microarrays/filteredInput/leRoch.pcl',
        'microarrays/filteredInput/leRoux.pcl',
        'microarrays/filteredInput/llinas3D7.pcl',
        'microarrays/filteredInput/llinasDD2.pcl',
        'microarrays/filteredInput/llinasHB3.pcl',
        'microarrays/filteredInput/milner.pcl',
        'microarrays/filteredInput/voss.pcl',
        'microarrays/filteredInput/young.pcl'],
    script:
        'scripts/microarray_filtering.R'

rule rename_joice:
    input:
        pcl='microarrays/filteredInput/leRoux.pcl',
    output:
        pcl='microarrays/filteredInput/Joiceetaldata.pcl',
    shell:
        r"""
        mv {input.pcl} {output.pcl}
        """

rule makeClusters:
    input:
        allGenesOldIds=os.path.join(workflow.basedir, "ref_data/microarrays/auxiliary_files/allGenesOldIds.txt"),
        geneProducts=os.path.join(workflow.basedir, "ref_data/microarrays/auxiliary_files/GeneProducts.txt"),
        
        young="microarrays/filteredInput/young.pcl",
        hu="microarrays/filteredInput/hu.pcl",
        leRoch="microarrays/filteredInput/leRoch.pcl",
        llinas3D7="microarrays/filteredInput/llinas3D7.pcl",
        llinasDD2="microarrays/filteredInput/llinasDD2.pcl",
        llinasHB3="microarrays/filteredInput/llinasHB3.pcl",
        voss="microarrays/filteredInput/voss.pcl",
        daily="microarrays/filteredInput/daily.pcl",
        leRoux="microarrays/filteredInput/Joiceetaldata.pcl",
        milner="microarrays/filteredInput/milner.pcl",
    output:
        dendrogram='makeClusters/geneDendrogram.pdf',
        geneData='makeClusters/geneData.tsv',
    script:
        'scripts/makeClusters.R'


rule infection_time_course:
    input:
        ss=config["rna_ss"],
        utils=os.path.join(workflow.basedir, "scripts/utils.R"),
    output:
        only_infect="results/infection_time_course.only_infect.pdf",
        deform="results/infection_time_course.deform.pdf",
    script:
        os.path.join(workflow.basedir, "scripts/infection_time_course.R")


rule trait_pairs:
    input:
        ss=config["rna_ss"],
        utils=os.path.join(workflow.basedir, "scripts/utils.R"),
    output:
        hist="results/histogram_traits.pdf",
        pairs="results/trait_pairs.pdf",
        corr="results/trait_pairs_corr.pdf",
    script:
        os.path.join(workflow.basedir, "scripts/trait_pairs.R")


rule upset_traits:
    input:
        ss=config["rna_ss"],
        utils=os.path.join(workflow.basedir, "scripts/utils.R"),
    output:
        upset="results/upset_traits.pdf",
    script:
        os.path.join(workflow.basedir, "scripts/upset_traits.R")


rule gff:
    output:
        gff=f"ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gff",
    shell:
        r"""
        curl -s https://plasmodb.org/common/downloads/release-{PLASMODB_RELEASE}/Pfalciparum3D7/gff/data/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gff > {output.gff}
        """


rule gaf:
    output:
        gaf=f"ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gaf",
    shell:
        r"""
        curl -s https://plasmodb.org/common/downloads/release-{PLASMODB_RELEASE}/Pfalciparum3D7/gaf/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7_GO.gaf > {output.gaf}
        """

rule tx_fasta:
    output:
        fa=f'ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7_AnnotatedTranscripts.fasta'
    shell:
        r"""
        curl -s https://plasmodb.org/common/downloads/release-{PLASMODB_RELEASE}/Pfalciparum3D7/fasta/data/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7_AnnotatedTranscripts.fasta > {output.fa}
        """

rule sexratio:
    input:
        fa=f'ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7_AnnotatedTranscripts.fasta',
        vossPcl="microarrays/filteredInput/voss.pcl",
        deseq_counts='deseq/counts.tsv.gz',
    output:
        markers_sex_defo="sexratio/markers_sex_defo.pdf",
        markers_sex_infect='sexratio/markers_sex_infect.pdf',
        markers_sex_ratio_voss='sexratio/markers_sex_ratio_voss.pdf',
        markers_sex_voss='sexratio/markers_sex_voss.pdf',
        markers_sex_ratio='sexratio/markers_sex_ratio.pdf',
    script:
        'scripts/sexratio.R'

rule gene_descriptions:
    input:
        gff=f"ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gff",
    output:
        gd="ref/Pfalciparum3D7.genes.tsv",
    run:
        from urllib.parse import unquote

        fout = open(output.gd, "w")
        fout.write("\t".join(["gene_id", "description", "name"]) + "\n")

        with open(input.gff) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                line = line.strip().split("\t")
                if line[2] in [
                    "gene",
                    "protein_coding_gene",
                    "ncRNA_gene",
                    "pseudogene",
                ]:
                    attr = line[8].split(";")
                    gid = [re.sub("^ID=", "", x) for x in attr if x.startswith("ID=")]
                    assert len(gid) == 1
                    desc = [
                        unquote(re.sub("^description=", "", x))
                        for x in attr
                        if x.startswith("description=")
                    ]
                    assert len(desc) == 1
                    name = [
                        unquote(re.sub("^Name=", "", x))
                        for x in attr
                        if x.startswith("Name=")
                    ]
                    assert len(name) == 0 or len(name) == 1
                    if name == []:
                        name = ["N/A"]
                    outline = gid + desc + name
                    fout.write("\t".join(outline) + "\n")
        fout.close()


rule deseq:
    input:
        ss=config["rna_ss"],
        cntfiles=glob.glob(f"{config['counts']}/*.count.gz"),
        genes="ref/Pfalciparum3D7.genes.tsv",
        utils=os.path.join(workflow.basedir, "scripts/utils.R"),
    output:
        cnt="deseq/counts.tsv.gz",
        vsd="deseq/vsd.tsv.gz",
        dge="deseq/infect_dge.tsv.gz",
        gex_pca_combatseq="results/gex_pca_combatseq.pdf",
        gex_pca_traits="results/gex_pca_traits.pdf",
        gex_logfc="results/gex_logfc.pdf",
    script:
        os.path.join(workflow.basedir, "scripts/deseq.R")


rule volcano:
    input:
        dge="deseq/infect_dge.tsv.gz",
    output:
        volcano="results/volcano.pdf",
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


rule logfc_pairs:
    input:
        dge="deseq/infect_dge.tsv.gz",
        utils=os.path.join(workflow.basedir, "scripts/utils.R"),
    output:
        corr="results/dge_trait_corr.pdf",
        logfc_pairs="results/dge_logfc_pairs.pdf",
    script:
        os.path.join(workflow.basedir, "scripts/logfc_pairs.R")


rule makeBioconductorAnnotationDbi:
    input:
        gff=f"ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gff",
        gaf=f"ref/PlasmoDB-{PLASMODB_RELEASE}_Pfalciparum3D7.gaf",
    output:
        dbi=directory("ref/org.Pfalciparum3D7.eg.db"),
    params:
        ver=f"0.{PLASMODB_RELEASE}.0",
    shell:
        r"""
        {workflow.basedir}/scripts/makeBioconductorAnnotationDbi.r --gff {input.gff} \
            --gaf {input.gaf} --genus Plasmodium --species falciparum3D7 --taxid 36329 \
            --outdir `dirname {output.dbi}` -m 'Dario Beraldi <dario.beraldi@glasgow.ac.uk>' \
            -a 'Dario Beraldi' --pckg-version {params.ver} --install
        """


rule clusterprofiler:
    input:
        dbi="ref/org.Pfalciparum3D7.eg.db",
        dge="deseq/infect_dge.tsv.gz",
    output:
        go_tsv="gsea/go/infection_deformability.tsv.gz",
        kegg_tsv="gsea/kegg/infection_deformability.tsv.gz",
    shell:
        r"""
        {workflow.basedir}/scripts/clusterprofiler.R \
            --dge {input.dge} \
            --orgdb `basename {input.dbi}` \
            --gsea-output-tsv {output.go_tsv} \
            --kegg-output-tsv {output.kegg_tsv} \
            --kegg-organism pfa
        """


rule plot_clusterprofiler:
    input:
        tsv="gsea/{annotation}/infection_deformability.tsv.gz",
    output:
        pdf="gsea/{annotation}/infection_deformability.pdf",
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
        dge="deseq/infect_dge.tsv.gz",
        clstData='makeClusters/geneData.tsv',
        cluster_class=os.path.join(
            workflow.basedir, "ref_data/crouch_cluster_classification.tsv"
        ),
    output:
        gsea_long="gsea/gsea_clusters.tsv",
        plot="results/gsea_clusters.pdf",
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


rule gex_voss_vs_infectivity:
    input:
        rna_ss=config["rna_ss"],
        vossPcl="microarrays/filteredInput/voss.pcl",
        infect="deseq/logcpm.tsv.gz",
    output:
        corr="results/voss_vs_infectivity_corr.pdf",
        tp_corr=os.path.join(
            workflow.basedir, "results/voss_vs_infectivity_timepoint_corr.pdf"
        ),
        peak="results/voss_vs_infectivity_peak.pdf",
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


rule selected_clusters_in_patients:
    input:
        dge="deseq/infect_dge.tsv.gz",
        pcl=['microarrays/filteredInput/Joiceetaldata.pcl', 'microarrays/filteredInput/milner.pcl'],
        gsea="gsea/gsea_clusters.tsv",
        allGenesOldIds=config["allGenesOldIds"],
    output:
        hm="results/selected_clusters_patients.pdf",
        patient_order=temp("patient_order.tmp.R"),
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

pdf('{output.hm}', height=(0.05 * nrow(mat))/2.54, width=20/2.54)
Heatmap(mat, name='Z-score\nof gene\nexpression',
    row_split=Clusters,
    row_labels=rep('', nrow(mat)),
    column_split=dataset,
    column_title=unique(dataset),
    column_labels=rep('', ncol(mat)),
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


rule voss_time_course_heatmap:
    input:
        gsea="gsea/gsea_clusters.tsv",
        vossPcl="microarrays/filteredInput/voss.pcl",
        clstData='makeClusters/geneData.tsv',
        kc_class=os.path.join(
            workflow.basedir, "ref_data/crouch_cluster_classification.tsv"
        ),
    output:
        hm="results/voss_time_course_heatmap.pdf",
        hm_gam="results/voss_time_course_heatmap_gam.pdf",
        hm_shared="results/voss_time_course_heatmap_shared.pdf",
        line="results/voss_time_course_selected.pdf",
    script:
        os.path.join(workflow.basedir, "scripts/voss_time_course_heatmap.R")


rule overlapping_dge:
    input:
        dge="deseq/infect_dge.tsv.gz",
    output:
        ovl="deseq/overlapping_dge.tsv.gz",
        upset="results/upset_dge.pdf",
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
        ref="ref/Pfalciparum3D7.genes.tsv",
        pcl=['microarrays/filteredInput/Joiceetaldata.pcl', 'microarrays/filteredInput/milner.pcl'],
        candidate_genes=config["candidate_genes"],
        patient_order="patient_order.tmp.R",
    output:
        hm="results/candidate_genes_patients.pdf",
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

pdf('{output.hm}', height=nrow(mat)*0.4/2.54, width=18/2.54)
Heatmap(mat, name='Z-score\nof gene\nexpression',
    column_split=dataset,
    column_title=unique(dataset),
    column_labels=rep('', ncol(mat)),
    row_title_rot=0,
    top_annotation=cls_col,
    cluster_columns=FALSE,
    cluster_rows=TRUE)
dev.off()
EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule annotate_clusters:
    input:
        clstData='makeClusters/geneData.tsv',
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

kc_clst <- fread('{input.clstData}')
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
        clstData='makeClusters/geneData.tsv',
        tlong= 'cluster_annotation/cluster_enrichment.tsv',
        peakAll=os.path.join(workflow.basedir, 'ref_data/peaksofMeanAll.tsv'),
    output:
        wide= 'cluster_annotation/cluster_enrichment_wide.tsv',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

kc_clst <- fread('{input.clstData}', select= c('V1', 'Cluster'), col.names= c('gene_id', 'Cluster'))
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
