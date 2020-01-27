#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Seed hits on paths'
p_subtitle <- 'for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)

# Reading input file
d <- read.table(opt$file, sep=',', header=T)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    # Data subsetting
    ds <- d[d$dataset == cdataset, c('seedlen', 'pathno', 'patched', 'stepsize', 'totalhits', 'seedsonpaths')]
    # Removing duplicates
    uds <- unique(ds[order(ds$pathno),])
    uds <- na.omit(uds)
    # Converting dataframe to tibble
    #t <- as_data_frame(uds)
    t <- tibble(seedlen=factor(uds$seedlen),
                npath=uds$pathno,
                stepsize=factor(uds$stepsize),
                type=factor(ifelse(uds$patched == 'yes', 'Patched', 'Full')),
                hits=uds$seedsonpaths/uds$totalhits)
    g <- ggplot(t, aes(x=npath, y=hits, colour=seedlen)) +
         geom_line() +
         geom_point() +
         facet_grid(stepsize ~ type, labeller=label_both) +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Paths',
              y='Ratio to total hits',
              colour='seed length') +
         scale_x_continuous(breaks=unique(t$npath)) +
         theme_bw()
    ggsave(paste('seedhits_', cdataset, '.pdf', sep=''))
}
