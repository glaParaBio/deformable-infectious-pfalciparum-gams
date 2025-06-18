library(data.table)
library(ggplot2)
library(UpSetR)

source(file.path(snakemake@scriptdir, 'utils.R'))

ss <- fread(snakemake@input[['ss']]) 

lss <- parse_sample_sheet(ss) 

stopifnot(all(!is.na(lss$value)))
lss[, value := 1]
lss <- dcast(data= lss, sample_ID ~ trait, value.var= 'value', fill= 0)
lss[, sample_ID := NULL]

pdf(snakemake@output[['upset']], width= 16/2.54, height= 12/2.54, onefile= FALSE)
upset(lss, 
      sets.x.label= 'N. samples',
      mainbar.y.label= 'Samples in intersection',
      mb.ratio = c(0.6, 0.4),
      text.scale = c(numbers_above_bars= 1.5, set_names= 2), order.by= c('freq'), decreasing= TRUE)
dev.off()
