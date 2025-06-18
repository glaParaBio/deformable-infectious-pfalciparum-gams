<!-- vim-markdown-toc GFM -->

* [Set up environment](#set-up-environment)
* [Run](#run)

<!-- vim-markdown-toc -->

## Set up environment

Best: Use conda/bioconda to create a separate environment and install the
project dependencies listed in requirements.txt:

```
conda create --yes -n 20210128_matt_infection
conda activate 20210128_matt_infection
conda install -n 20210128_matt_infection --yes --file requirements.txt
```

## Run

```
# DATADIR=/export/projects/III-data/wcmp_bioinformatics/db291g/data/20210128_matt_infection
DATADIR=~/sharedscratch/data/20210128_matt_infection

snakemake -p -n -j 10 \
    --cluster 'sbatch --cpus-per-task=10 --mem=5G --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    --latency-wait 60 \
    -C rna_ss=$DATADIR/deformability/RNAseq_master_samples.csv \
       counts=$DATADIR/deformability/counts/ \
       candidate_genes=$PWD/candidate_genes.tsv \
       #clstData=$DATADIR/allDatasets/tablesWithVoss/geneData.tsv \
       #peakAll=$DATADIR/allDatasets/tablesWithVoss/peaksofMeanAll.tsv \
       #peakVoss=$DATADIR/vossOnly/vossTables/peaksofMeanAll.tsv \
       #vossPcl=$DATADIR/inputFiles/filteredInput/voss.pcl.gz \
       pclDir=$DATADIR/inputFiles/filteredInput \
       #deformability=$DATADIR/deformability/outputFiles/deformabilityResults.tsv \
       #pelle=~/git_repos/glaParaBio/dariober/unsorted/20210406_annotate_pelle/output/pelle_clusters_plasmodb.tsv \
       #allGenesOldIds=$DATADIR/inputFiles/allGenesOldIds.txt \
    -d output
```
