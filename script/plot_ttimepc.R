#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Traverse time per chunk'
p_subtitle <- 'for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
# Command-line options parser
option_list = list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                   make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                   make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'))
opt_parser = OptionParser(option_list=option_list)
# Parsing command-line arguments
opt = parse_args(opt_parser)

# Reading input file
d <- read.table(opt$file, sep=',', header=T)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    # Data subsetting
    ds <- d[d$dataset == cdataset & d$chunksize < 10000, c('pathno', 'chunksize', 'stepsize', 'patched', 'traversetimeperchunk')]
    # Removing duplicates
    uds <- unique(ds[order(ds$pathno),])
    # Converting dataframe to tibble
    #t <- as_data_frame(uds)
    t <- tibble(npath=uds$pathno,
                chunk=factor(uds$chunksize),
                stepsize=factor(uds$stepsize),
                type=factor(ifelse(uds$patched=='yes', 'Patched', 'Full')),
                ttime=(uds$traversetimeperchunk/10^6))
    g <- ggplot(t, aes(x=npath, y=ttime, colour=chunk, group=chunk)) +
         geom_point() +
         geom_line() +
         facet_grid(type ~ stepsize) +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Paths',
              y='Traverse time per chunk (s)',
              colour='Chunk size') +
         scale_x_continuous(breaks=unique(t$npath)) +
         theme_bw()
    ggsave(paste('ttimepc_', cdataset, '.pdf', sep=''))
}
