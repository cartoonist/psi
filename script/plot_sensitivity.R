#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Sensitivity'
p_subtitle <- 'for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-n', '--nofreads'), type='integer', help='Number of reads', metavar='integer'),
                    make_option(c('-l', '--readslen'), type='integer', help='Length of reads', metavar='integer'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)

if (is.null(opt$nofreads)) {
    print_help(opt_parser)
    stop('The number of reads must be specified.', call.=FALSE)
}

if (is.null(opt$readslen)) {
    print_help(opt_parser)
    stop('The length of the reads must be specified.', call.=FALSE)
}

# Reading input file
d <- read.table(opt$file, sep=',', header=TRUE)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    # Data subsetting
    ds <- d[d$dataset == cdataset, c('seedlen', 'pathno', 'chunksize', 'patched', 'stepsize', 'totalhits')]
    ds <- na.omit(ds)
    # Removing duplicates
    uds <- unique(ds[order(ds$chunksize),])
    col_readslen <- rep(opt$readslen, dim(uds)[1])
    col_nofreads <- rep(opt$nofreads, dim(uds)[1])
    # Converting dataframe to tibble
    #t <- as_data_frame(uds)
    t <- tibble(chunksize=factor(uds$chunksize),
                stepsize=uds$stepsize,
                type=factor(ifelse(uds$patched=='yes', 'Patched', 'Full')),
                npath=factor(uds$pathno),
                sensitivity=(uds$totalhits / ((col_readslen - uds$seedlen + 1) * col_nofreads)))
    g <- ggplot(t, aes(x=stepsize, y=sensitivity, colour=npath, group=npath)) +
         geom_point() +
         facet_grid(. ~ type, labeller=label_both) +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Step size',
              y='Sensitivity',
              colour='Paths') +
         scale_x_continuous(breaks=levels(t$stepsize)) +
         scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, 0.2)) +
         theme_bw()
    ggsave(paste('sensitivity_', cdataset, '.pdf', sep=''))
}
