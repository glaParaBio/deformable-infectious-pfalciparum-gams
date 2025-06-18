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
    --latency-wait 60 \
    -C rna_ss=$DATADIR/deformability/RNAseq_master_samples.csv \
       counts=$DATADIR/deformability/counts \
       candidate_genes=$PWD/candidate_genes.tsv \
       pclDir=$DATADIR/inputFiles/filteredInput \
       clstData=$DATADIR/allDatasets/tablesWithVoss/geneData.tsv \
       vossPcl=$DATADIR/inputFiles/filteredInput/voss.pcl.gz \
    -d output
```
