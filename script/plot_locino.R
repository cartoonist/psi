#!/usr/bin/env Rscript
library('tibble')
library('reshape2')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Starting loci'
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
    ds <- d[d$dataset == cdataset, c('pathno', 'stepsize', 'patched', 'locino', 'uniqnodes')]
    # Removing duplicates
    uds <- unique(ds[order(ds$pathno),])
    # Melting
    uds <- melt(uds, id.vars=c('pathno', 'stepsize', 'patched'))
    # Converting dataframe to tibble
    #t <- as_data_frame(uds)
    t <- tibble(npath=uds$pathno,
                stepsize=factor(uds$stepsize),
                type=factor(ifelse(uds$patched == 'yes', 'Patched', 'Full')),
                variable=factor(ifelse(uds$variable == 'locino', 'Loci', 'Nodes')),
                loci=uds$value)
    g <- ggplot(t, aes(x=npath, y=loci, colour=stepsize, linetype=variable, group=interaction(stepsize, variable))) +
         geom_point(size=2, alpha=0.7) +
         geom_line(alpha=0.7) +
         facet_grid(. ~ type) +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Paths',
              y='Counts',
              colour='Step size',
              linetype='Starting') +
         scale_x_continuous(breaks=unique(t$npath)) +
         theme_bw()
    ggsave(paste('locino_', cdataset, '.pdf', sep=''))
}
