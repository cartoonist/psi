#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Performance per occurrence'
p_subtitle <- 'k = %s; for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
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
    for (cseedlen in unique(d$seedlen)) {
        # Data subsetting
        ds <- d[d$dataset == cdataset & d$seedlen == cseedlen, c('pathno', 'chunksize', 'patched', 'stepsize', 'totalhits', 'seedfindingtime')]
        ds <- na.omit(ds)
        # Removing duplicates
        uds <- unique(ds[order(ds$chunksize),])
        # Converting dataframe to tibble
        #t <- as_data_frame(uds)
        t <- tibble(chunksize=factor(uds$chunksize),
                    stepsize=uds$stepsize,
                    type=factor(ifelse(uds$patched=='yes', 'Patched', 'Full')),
                    npath=factor(uds$pathno),
                    perfperocc=(uds$seedfindingtime / uds$totalhits))
        g <- ggplot(t, aes(x=chunksize, y=perfperocc, fill=npath)) +
             geom_col(colour='#444444', position=position_dodge()) +
             facet_grid(type ~ stepsize, labeller=label_both) +
             labs(title=opt$title,
                  subtitle=sprintf(opt$subtitle, cseedlen, cdataset),
                  x='Chunk size',
                  y=expression('Average time per occurrence ('*mu*'s)'),
                  fill='Paths') +
             scale_y_continuous(limit=c(0, max(t$perfperocc)*1.1), expand=c(0, 0)) +
             theme_bw()
        ggsave(paste('perfperocc_', cdataset, '_', cseedlen, '.pdf', sep=''))
    }
}
