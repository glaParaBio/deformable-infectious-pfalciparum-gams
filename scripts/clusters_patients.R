library(ComplexHeatmap)

clstData <- fread(snakemake@input[['clstData']])
setnames(clstData, 'V1', 'gene_id')

pcl <- fread(snakemake@input[['pcl']])
pcl[, name := NULL]
pcl[, GWeight := NULL]
setnames(pcl, 1, 'gene_id')
pcl <- melt(data= pcl, id.vars= 'gene_id', value.name= 'gex', variable.name= 'Time')
pcl[, zscore := scale(gex), by= gene_id]

