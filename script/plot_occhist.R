#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

in_file <- 'occs.csv'
p_title <- 'Seed occurrence count distribution'
p_subtitle <- 'of simulated reads from dataset "%s"'
legend_pos <- 'right'
default_bins <- NULL
# Command-line options parser
option_list <-
  list(make_option(c('-f', '--file'), type='character', default=in_file,
                   help='Input CSV file\n\t\t[default=\'%default\']',
                   metavar='character'),
       make_option(c('-t', '--title'), type='character', default=p_title,
                   help='Plot title\n\t\t[default=\'%default\']',
                   metavar='character'),
       make_option(c('-s', '--subtitle'), type='character', default=p_subtitle,
                   help='Plot subtitle pattern\n\t\t[default=\'%default\']',
                   metavar='character'),
       make_option(c('-L', '--legend'), type='character', default=legend_pos,
                   help='Legend position\n\t\t[default=\'%default\']',
                   metavar='character'),
       make_option(c('-b', '--bins'), type='integer', default=default_bins,
                   help='Number of bins\n\t\t[default=\'%default\']',
                   metavar='integer'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)

# Reading input file
d <- read.table(opt$file, sep=',', header=TRUE, stringsAsFactors=TRUE)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    # Data subsetting
    t <- tibble(count=d[d$dataset == cdataset, c('count')])
    g <- ggplot(t, aes(count)) +
      geom_histogram(bins=opt$bins, alpha=0.5) +
      labs(title=opt$title,
           subtitle=sprintf(opt$subtitle, cdataset),
           x='Occurrences',
           y='Frequency') +
      scale_y_continuous(trans="log10") +
      theme_bw() +
      theme(axis.text=element_text(size=22),
            axis.title=element_text(size=24),
            panel.border=element_rect(fill=NA, colour='gray20', size=3),
            legend.position=opt$legend)
    ggsave(paste('occhist_', cdataset, '.pdf', sep=''))
}
