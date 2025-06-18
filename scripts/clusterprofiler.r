#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))


parser <- ArgumentParser(description='Run clusterprofiler')
parser$add_argument('--dge', help='File of differential gene expression table [required]', required=TRUE)
parser$add_argument('--orgdb', help='Organism database [required]', required=TRUE)
parser$add_argument('--kegg-organism', help='Organism abbreviation for gseKEGG. See http://www.genome.jp/kegg/catalog/org_list.html [%(default)s]', default='pfa')
parser$add_argument('--gsea-output-tsv', help='File for gsea output table', required=TRUE)
parser$add_argument('--kegg-output-tsv', help='File for kegg output table', required=TRUE)

parser$add_argument('--version', '-v', action= 'version', version='0.1.0')

xargs <- parser$parse_args()

suppressPackageStartupMessages(library(xargs$orgdb, character.only=TRUE))

dge <- fread(xargs$dge)

gseaOut <- list()
keggOut <- list()
for(cntr in unique(dge$trait)) {
    full <- dge[trait == cntr]
    geneList <- full$log2FoldChange
    names(geneList) <- full$gene_id
    geneList <- na.omit(geneList)
    geneList <- sort(geneList, decreasing = TRUE)
    
    # GSEA GO
    set.seed(1234)
    gsea <- gseGO(geneList=geneList, 
                 ont="ALL", 
                 keyType="GID", 
                 minGSSize=10, 
                 maxGSSize=1000, 
                 pvalueCutoff=1, 
                 verbose=FALSE, 
                 OrgDb=get(xargs$orgdb), # use get because you pass the object orgdb not the string 'org.Pfalciparum3D7.eg.db'
                 pAdjustMethod="BH",
                 seed=TRUE)
    
    gsea <- as.data.table(gsea@result)
    gsea[, trait := cntr]
    gseaOut[[length(gseaOut) + 1]] <- gsea

    # GSEA KEGG
    set.seed(1234)
    kegg <- gseKEGG(geneList, xargs$kegg_organism, pvalueCutoff=1,
        minGSSize=10, maxGSSize=1000, seed=TRUE)

    kegg <- as.data.table(kegg@result)
    kegg[, trait := cntr]
    keggOut[[length(keggOut) + 1]] <- kegg
}
gsea <- rbindlist(gseaOut)
gsea[, qvalues := NULL]
setcolorder(gsea, 'trait')
fwrite(gsea, xargs$gsea_output_tsv, sep='\t')

kegg <- rbindlist(keggOut)
kegg[, qvalues := NULL]
setcolorder(kegg, 'trait')
fwrite(kegg, xargs$kegg_output_tsv, sep='\t')
