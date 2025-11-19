<!-- vim-markdown-toc GFM -->

* [Set up environment](#set-up-environment)
* [Data preparation and clustering](#data-preparation-and-clustering)
* [Run workflow](#run-workflow)

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

## Data preparation and clustering

```
Rscript scripts/microarray_filtering.R

cd ref_data/microarrays/filteredInput
ln -s leRoux.pcl Joiceetaldata.pcl
```

```
Rscript scripts/makeClusters.R
```

## Run workflow

```
DATADIR=~/sharedscratch/data/deformable-infectious-pfalciparum-gams

snakemake -p -n -j 10 \
    --latency-wait 60 \
    -C rna_ss=$PWD/ref_data/deformability_infection_data.tsv \
       counts=$DATADIR/deformability/counts \
       candidate_genes=$PWD/ref_data/candidate_genes.tsv \
       pclDir=$PWD/ref_data/microarrays/filteredInput \
       clstData=$PWD/ref_data/geneData.tsv \
       vossPcl=$PWD/ref_data/microarrays/filteredInput/voss.pcl \
       allGenesOldIds=$PWD/ref_data/microarrays/auxiliary_files/allGenesOldIds.txt \
    -d output
```
