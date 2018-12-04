#!/usr/bin/env Rscript
library('tibble')
library('reshape2')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Path indexing time'
p_subtitle <- 'for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
legend_pos <- 'right'
size_unit <- 's'
patched_value <- 'all'  # show locino for both patched and full path index
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-p', '--patched'), type='character', default=patched_value, help='Only report this `patched` value\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-m', '--microseconds'), action="store_true", default=FALSE, help='Times are in microseconds\n\t\t[default=\'%default\']'),
                    make_option(c('-u', '--unit'), type='character', default=size_unit, help='Time unit\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-L', '--legend'), type='character', default=legend_pos, help='Legend position\n\t\t[default=\'%default\']', metavar='character'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)
opt$time_coeff <- switch(opt$unit, h=3600, m=60, s=1, ms=10^-3, us=10^-6)
if (is.null(opt$time_coeff)) {
    opt$unit = 's'
    opt$time_coeff = 1
}
opt$time_coeff <- ifelse(opt$microseconds, opt$time_coeff * 10^6, opt$time_coeff)

# Reading input file
d <- read.table(opt$file, sep=',', header=TRUE)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    # Data subsetting
    ds <- d[d$dataset == cdataset, c('pathno', 'patched', 'pathpicktime', 'pathindextime', 'pathsavetime')]
    # Removing duplicates
    uds <- unique(ds[order(ds$pathno),])
    # Melting
    uds <- melt(uds, id.vars=c('pathno', 'patched'))
    if (opt$patched != 'all') {
        uds <- uds[!(uds$patched != opt$patched),]
    }
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
                rtime=(uds$value / opt$time_coeff))
    g <- ggplot(t, aes(x=npath, y=rtime, fill=stage_order)) +
         geom_area() +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Paths',
              y=paste('Time (', opt$unit, ')', sep=''),
              fill='Operation') +
         scale_x_continuous(breaks=unique(t$npath)) +
         scale_y_continuous(limit=c(0, max(aggregate(t$rtime, list(t$npath, t$type), FUN=sum)$x)*1.1), expand=c(0, 0)) +
         scale_fill_discrete(label=c('Save', 'Index', 'Pick')) +
         theme_bw() +
         theme(axis.text=element_text(size=22),
               axis.title=element_text(size=24),
               panel.border=element_rect(fill=NA, colour='gray20', size=3),
               legend.position=opt$legend)
    if (opt$patched == 'all') {
        g <- g + facet_grid(. ~ type)
    }
    ggsave(paste('pindextime_', cdataset, '.pdf', sep=''))
}
