<!-- vim-markdown-toc GFM -->

* [Set up environment](#set-up-environment)
* [Run workflow](#run-workflow)
* [TODO:](#todo)

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

## Run workflow

```
snakemake -p -n -j 10 \
    --latency-wait 60 \
    -C rna_ss=$PWD/ref_data/deformability_infection_data.tsv \
       counts=$PWD/ref_data/counts \
       allGenesOldIds=$PWD/ref_data/microarrays/auxiliary_files/allGenesOldIds.txt \
       pelle=$PWD/ref_data/pelle_clusters_plasmodb.tsv \
    -d output
```

For `pelle_clusters_plasmodb.tsv` see https://github.com/glaParaBio/dariober/tree/master/unsorted/20210406_annotate_pelle

## TODO:

* Include code for generating `ref_data/peaksofMeanAll.tsv`

* Include code for downloading and formatting public data in `annotate_clusters`
