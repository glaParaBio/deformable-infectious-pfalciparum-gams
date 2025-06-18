library(data.table)

oldToNewGeneIDs <- function(allGenesOldIdsFile, pcl) {
    ids <- fread(allGenesOldIdsFile, header=TRUE)
    ids[, previous_id := sub('Previous IDs: ', '', `Previous ID(s)`)]
    previous_id <- strsplit(ids$previous_id, ';')
    n_ids <- sapply(previous_id, length)
    gene_id <- rep(ids$`Gene ID`, n_ids)
    old2new <- unique(data.table(previous_id=tolower(unlist(previous_id)), gene_id))
    
    gene_col <- names(pcl)[1]
    setnames(pcl, gene_col, 'previous_id')
    pcl <- melt(data=pcl, id.vars= 'previous_id', value.name= 'gex', variable.name= 'array_id')
    pcl[, previous_id := tolower(previous_id)]

    pcl <- merge(pcl, old2new, by='previous_id', all.x=TRUE)
    # Average expression when two old IDs map to the same new ID.  There are also
    # cases where the old ID maps to more than 1 new IDs. The merge will
    # effectively duplicate these genes. This is not ideal but it's ok.
    pcl <- pcl[, list(gex=mean(gex)), list(array_id, gene_id)]
    pcl <- dcast(pcl, gene_id ~ array_id, value.var='gex')
    setnames(pcl, 'gene_id', gene_col)
    return(pcl)
}

read_pcl <- function(pcl) {
    dataset <- sub('\\.pcl$', '', sub('\\.gz$', '', basename(pcl)))
    dat <- fread(pcl)
    if('name' %in% names(dat)) {dat[, name := NULL]}
    if('GWeight' %in% names(dat)) {dat[, GWeight := NULL]}

    setnames(dat, 1, 'gene_id')
    dat <- melt(data=dat, id.vars= 'gene_id', value.name= 'gex', variable.name= 'array_id')
    return(dat)
}

