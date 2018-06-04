#!/usr/bin/env Rscript
library('tibble')
library('reshape2')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Path indexing time'
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
    ds <- d[d$dataset == cdataset, c('pathno', 'patched', 'pathpicktime', 'pathindextime', 'pathsavetime')]
    # Removing duplicates
    uds <- unique(ds[order(ds$pathno),])
    # Melting
    uds <- melt(uds, id.vars=c('pathno', 'patched'))
    trans <- function(x) {
      if (x == 'pathpicktime') {
        'Pick'
      } else if (x == 'pathindextime') {
        'Index'
      } else {
        'Save'
      }
    }
    trans_order <- function(x) {
      if (x == 'pathpicktime') {
        '3'
      } else if (x == 'pathindextime') {
        '2'
      } else {
        '1'
      }
    }
    # Converting dataframe to tibble
    #t <- as_data_frame(uds)
    t <- tibble(npath=uds$pathno,
                type=factor(ifelse(uds$patched=='yes', 'Patched', 'Full')),
                stage=factor(sapply(uds$variable, trans)),
                stage_order=factor(sapply(uds$variable, trans_order)),
                rtime=(uds$value / 10^6))
    g <- ggplot(t, aes(x=npath, y=rtime, fill=stage_order)) +
         geom_area() +
         facet_grid(. ~ type) +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Paths',
              y=expression('Time (s)'),
              fill='Operation') +
         scale_x_continuous(breaks=unique(t$npath)) +
         scale_y_continuous(limit=c(0, max(aggregate(t$rtime, list(t$npath, t$type), FUN=sum)$x)*1.1), expand=c(0, 0)) +
         scale_fill_discrete(label=c('Save', 'Index', 'Pick')) +
         theme_bw()
    ggsave(paste('pindextime_', cdataset, '.pdf', sep=''))
}
