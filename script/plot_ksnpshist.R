#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

in_file <- 'ksnpcounts.csv'
p_title <- 'Frequency of SNPs in k-mers'
p_subtitle <- 'for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)

# Reading input file
d <- read.table(opt$file, sep=',', header=TRUE)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    # Data subsetting
    ds <- d[d$dataset == cdataset, c('k', 'count')]
    t <- tibble(k=factor(ds$k),
                count=ds$count)
    g <- ggplot(t, aes(count, colour=k, fill=k)) +
         geom_histogram(binwidth=1, alpha=0.5) +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Number of SNPs',
              y='K-mers') +
         scale_x_continuous(breaks=0:(max(t$count)+1)) +
         theme_bw()
    g <- g + 
         scale_y_continuous(limit=c(0, max(ggplot_build(g)$data[[1]]$count)*1.1), expand=c(0, 0))
    ggsave(paste('ksnpshist_', cdataset, '.pdf', sep=''))
}
