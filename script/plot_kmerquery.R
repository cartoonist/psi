#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'K-mer query time'
p_subtitle <- 'k = %s; for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
legend_pos <- 'right'
patched_value <- 'all'  # show locino for both patched and full path index
stepsize_value <- 'all' # show all values for step size
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-n', '--nofreads'), type='integer', help='Number of reads.', metavar='integer'),
                    make_option(c('-l', '--readslen'), type='integer', help='Length of reads.', metavar='integer'),
                    make_option(c('-p', '--patched'), type='character', default=patched_value, help='Only report this `patched` value\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-e', '--stepsize'), type='character', default=patched_value, help='Only report this `stepsize` value\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-L', '--legend'), type='character', default=legend_pos, help='Legend position\n\t\t[default=\'%default\']', metavar='character'))
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
    for (cseedlen in unique(d$seedlen)) {
        # Data subsetting
        ds <- d[d$dataset == cdataset & d$seedlen == cseedlen, c('pathno', 'chunksize', 'patched', 'stepsize', 'seedfindingtime')]
        ds <- na.omit(ds)
        # Removing duplicates
        uds <- unique(ds[order(ds$chunksize),])
        if (opt$patched != 'all') {
            uds <- uds[!(uds$patched != opt$patched),]
        }
        if (opt$stepsize != 'all') {
            uds <- uds[!(uds$stepsize != opt$stepsize),]
        }
        col_readslen <- rep(opt$readslen, dim(uds)[1])
        col_nofreads <- rep(opt$nofreads, dim(uds)[1])
        # Converting dataframe to tibble
        #t <- as_data_frame(uds)
        t <- tibble(chunksize=factor(uds$chunksize),
                    stepsize=uds$stepsize,
                    type=factor(ifelse(uds$patched=='yes', 'Patched', 'Full')),
                    npath=factor(uds$pathno),
                    kmerquery=(uds$seedfindingtime / (as.integer(col_readslen / cseedlen) * col_nofreads) * 10^3))
        g <- ggplot(t, aes(x=chunksize, y=kmerquery, fill=npath)) +
             geom_col(colour='#444444', position=position_dodge()) +
             labs(title=opt$title,
                  subtitle=sprintf(opt$subtitle, cseedlen, cdataset),
                  x='Chunk size',
                  y=expression('K-mer query (ms)'),
                  fill='Paths') +
             scale_y_continuous(limit=c(0, max(t$kmerquery)*1.1), expand=c(0, 0)) +
             theme_bw() +
             theme(axis.text=element_text(size=22),
                   axis.title=element_text(size=24),
                   panel.border=element_rect(fill=NA, colour='gray20', size=3),
                   legend.position=opt$legend)
        if (opt$patched == 'all' && opt$stepsize == 'all') {
            g <- g + facet_grid(type ~ stepsize, labeller=label_both)
        }
        else if (opt$patched == 'all') {
            g <- g + facet_grid(. ~ type)
        }
        else if (opt$stepsize == 'all') {
            g <- g + facet_grid(. ~ stepsize)
        }
        ggsave(paste('kmerquery_', cdataset, '_', cseedlen, '.pdf', sep=''))
    }
}
