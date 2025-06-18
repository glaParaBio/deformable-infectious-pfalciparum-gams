library(data.table)

jaccard <- function(x, y) {
    union <- length(unique(c(x, y)))
    intx <- length(intersect(x, y))
    return(intx/union)
}

pcl <- fread('/export/projects/III-data/wcmp_bioinformatics/db291g/data/20210128_matt_infection/inputFiles/filteredInput/llinasHB3.pcl.gz')
pcl[, name := NULL]
pcl[, GWeight := NULL]
pcl <- melt(data= pcl, id.vars= 'X')
pcl <- pcl[, list(value= mean(value, na.rm= TRUE)), by= list(X, variable)]
pcl <- dcast(data= pcl, X ~ variable, var.value= 'value')
pcl <- as.matrix(pcl, rownames= 'X')
pcl <- scale(pcl)

odd <- pcl[, which(1:ncol(pcl) %% 2 == 0) - 1]
even <- pcl[, which(1:ncol(pcl) %% 2 == 0)]

rr <- rep(NA, ncol(odd))
for(i in 1:ncol(odd)) {
    r <- cor(odd[,i], even[,i], use= 'complete.obs')
    rr[i] <- r
}
timepoint_r <- data.table(odd_timepoint= colnames(odd), even_timepoint= colnames(even), corr_coef= rr)
timepoint_r 

# mat <- cor(pcl, use= 'complete.obs')
# pdf('tmp.pdf')
# corrplot(mat, diag= FALSE)
# dev.off()

odd <- odd[, which(rr > 0)]
even <- even[, which(rr > 0)]

hcor <- 0.5
clst <- list()
for(x in c('odd', 'even')){
    if(x == 'odd') {
        m <- odd
    } else if(x == 'even') {
        m <- even
    } else {
        stop()
    }
    cmat <- 1 - cor(t(m), use= 'complete.obs')
    hc <- hclust(as.dist(cmat), method= 'complete')
    cls <- cutree(hc, h= 1 - hcor)
    cls <- data.table(
        dataset= x,
        gene_id= names(cls),
        cluster= cls
    )
    clst[[length(clst) + 1]] <- cls
}
clst <- rbindlist(clst)

intx <- list()
for(k in unique(clst$cluster)) {
    genes <- clst[dataset == 'odd' & cluster == k]$gene_id
    matches <- clst[dataset != 'odd', list(odd_cluster_name= k, odd_size= length(genes), even_size= nrow(.SD), 
        in_common= sum(genes %in% .SD$gene_id), jaccard= jaccard(genes, .SD$gene_id)), by= list(even_cluster_name= cluster, dataset)]
    stopifnot(length(unique(matches$dataset)) == 1)
    matches[, dataset := NULL]
    matches <- matches[odd_size > 5 & even_size > 5]
    intx[[length(intx) + 1]] <- matches[order(-jaccard)][1]
}
intx <- rbindlist(intx)[!is.na(odd_cluster_name)][order(-jaccard)]
intx[, overlap_coef := in_common / ifelse(odd_size < even_size, odd_size, even_size)]

quit()

##################

pcl2[, name := NULL]
pcl2[, GWeight := NULL]

common <- intersect(names(pcl), names(pcl2))
pcl <- pcl[, common, with= FALSE]
pcl2 <- pcl2[, common, with= FALSE]

genes <- intersect(pcl$X, pcl2$X)
rr <- rep(NA, length(g))
for(i in 1:length(genes)) {
    g <- genes[i]
    x <- pcl[X == g, 2:ncol(pcl)]
    y <- pcl2[X == g, 2:ncol(pcl2)]
    if(nrow(x) != 1 | nrow(y) != 1) {
        next
    }
    r <- cor(unlist(x), unlist(y))
    rr[i] <- r
}

pcl[, name := NULL]
pcl[, GWeight := NULL]
setnames(pcl, 'X', 'gene_id')
sample_order <- names(pcl)[which(names(pcl) != 'gene_id')]

pcl <- melt(pcl, id.vars= 'gene_id', variable.name= 'sample_id', value.name= 'exprs')
na <- unique(pcl[is.na(exprs)]$sample_id)
stopifnot(length(na) < length(unique(pcl$sample_id))/10) # We don't drop to many samples
pcl <- pcl[!sample_id %in% na]
sample_order <- sample_order[which(!sample_order %in% na)]
pcl <- pcl[, list(exprs= mean(exprs)), by= list(gene_id, sample_id)]

mat <- as.matrix(dcast(data= pcl, gene_id ~ sample_id, value.var= 'exprs'), rownames= 'gene_id')
stopifnot(identical(colnames(mat), sample_order))
mat <- scale(mat)

reps <- list()
for(i in c(0, 1)) {
    use <- which(1:ncol(mat) %% 2 == 0) - i
    cmat <- 1 - cor(t(mat[, use]))
    hc <- hclust(as.dist(cmat))
    cls <- cutree(hc, k= 200)
    cls <- data.table(
        use= i,
        gene_id= names(cls),
        cluster= cls
    )
    reps[[length(reps) + 1]] <- cls
}
reps <- rbindlist(reps)

intx <- list()
for(k in unique(reps$cluster)) {
    genes <- reps[use == 0 & cluster == k]$gene_id
    matches <- reps[use != 0, list(first_k= k, first_size= length(genes), second_size= nrow(.SD), in_common= sum(genes %in% .SD$gene_id), jaccard= jaccard(genes, .SD$gene_id)), by= list(second_k= cluster, use)]
    stopifnot(length(unique(matches$use)) == 1)
    matches[, use := NULL]
    matches <- matches[first_size > 5 & second_size > 5]
    intx[[length(intx) + 1]] <- matches[order(-jaccard)][1]
}
intx <- rbindlist(intx)[!is.na(first_k)][order(-jaccard)]

