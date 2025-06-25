library(data.table)
library(ggplot2)
library(ggbeeswarm)

source(file.path(snakemake@scriptdir, 'utils.R'))

ss <- fread(snakemake@input[['ss']]) 

lss <- parse_sample_sheet(ss) 

avg <- lss[, list(value= mean(value)), by= list(trait, Day, day_offset)]

gg <- ggplot(data= lss[day_offset == 'only_infect'], aes(x= Day, y= value)) +
    geom_line(data= avg[day_offset == 'only_infect'], aes(colour=NULL), colour= 'grey60', linewidth=1, alpha=0.5) +
    geom_quasirandom(groupOnX= TRUE, size= 1, width= 0.2) +
    scale_colour_brewer(palette='Dark2') +
    facet_wrap(~ trait, scales= 'free_y', ncol= 1) +
    scale_x_continuous(breaks= unique(lss[day_offset == 'only_infect']$Day), minor_breaks= NULL) +
    ylab(NULL) +
    theme_light() +
    theme(strip.text= element_text(colour= 'black', size= 12), legend.position='none')
ggsave(snakemake@output[['only_infect']], width= 16, height= 2 + (5 * 3), units= 'cm')

gg <- ggplot(data=lss[trait == 'Deformability'], aes(x=Day, y=dicho_value)) +
    geom_quasirandom(groupOnX=TRUE, size=1.5, width=0.4, orientation='x') +
    theme_light() +
    ylab('') +
    scale_y_discrete(expand=c(0.1, 0.1)) +
    scale_x_continuous(breaks= unique(lss[trait == 'Deformability']$Day), minor_breaks= NULL)
ggsave(snakemake@output[['deform']], width=14, height=6, units='cm')
