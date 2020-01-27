#!/usr/bin/env Rscript
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- '"%s" vs Parameters'
p_subtitle <- 'for "snd-%s"'
# Command-line options parser
option_list = list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                   make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title pattern\n\t\t[default=\'%default\']', metavar='character'),
                   make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'))
opt_parser = OptionParser(option_list=option_list)
# Parsing command-line arguments
opt = parse_args(opt_parser)

# Reading input file
d <- read.table(opt$file, sep=',', header=T)
# NOTE: Parameters/measures description
parameters <- c('seedlen', 'pathno', 'patched', 'chunksize', 'stepsize')
unused <- c('region', 'context', 'readsindex')
measures <- names(d)[!names(d) %in% c('dataset', parameters, unused)]
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    for (m in measures) {
        # Data subsetting
        ds <- d[d$dataset == cdataset, c(parameters, m)]
        # Renaming values under 'patched' column
        ds$patched <- ifelse(ds$patched=='yes', 1, -1)
        # Removing duplicates
        uds <- unique(ds)
        pdf(paste('pairs_', cdataset, '-', m, '.pdf', sep=''))
        pairs(uds, main=paste(sprintf(p_title, m), '--', sprintf(p_subtitle, cdataset)))
        dev.off()
    }
}
