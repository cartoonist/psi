#!/usr/bin/env Rscript
library('tibble')
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
    # Converting dataframe to tibble
    #t <- as_data_frame(uds)
    t <- tibble(npath=uds$pathno,
                stepsize=factor(uds$stepsize),
                type=factor(ifelse(uds$patched == 'yes', 'Patched', 'Full')),
                nloci=uds$locino,
                uniq=uds$uniqnodes)
    g <- ggplot(t, aes(x=npath, y=nloci, colour=stepsize, shape=type, group=interaction(stepsize, type))) +
         geom_point(size=2) +
         geom_line() +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Paths',
              y='Loci',
              colour='Step size',
              shape='Type') +
         scale_x_continuous(breaks=unique(t$npath)) +
         scale_shape_manual(values=c(2, 6)) +
         theme_bw()
    ggsave(paste('locino_', cdataset, '.pdf', sep=''))
}
