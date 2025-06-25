<!-- vim-markdown-toc GFM -->

* [Set up environment](#set-up-environment)
* [Run](#run)

<!-- vim-markdown-toc -->

*Work in progress*

## Set up environment

Best: Use conda/bioconda to create a separate environment and install the
project dependencies listed in requirements.txt:

```
conda create --yes -n deformable-infectious-pfalciparum-gams
conda activate deformable-infectious-pfalciparum-gams
conda install -n deformable-infectious-pfalciparum-gams --yes --file requirements.txt
```

## Run

```
DATADIR=~/sharedscratch/data/deformable-infectious-pfalciparum-gams

snakemake -p -n -j 10 \
    --latency-wait 60 \
    -C rna_ss=$PWD/ref_data/deformability_infection_data.tsv \
       counts=$DATADIR/deformability/counts \
       candidate_genes=$PWD/ref_data/candidate_genes.tsv \
       pclDir=$DATADIR/inputFiles/filteredInput \
       clstData=$DATADIR/allDatasets/tablesWithVoss/geneData.tsv \
       vossPcl=$DATADIR/inputFiles/filteredInput/voss.pcl.gz \
       allGenesOldIds=$DATADIR/inputFiles/allGenesOldIds.txt \
    -d output
```
