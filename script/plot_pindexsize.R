#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Path index size'
p_subtitle <- 'for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
legend_pos <- 'right'
size_unit <- 'MB'
patched_value <- 'all'  # show locino for both patched and full path index
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-u', '--unit'), type='character', default=size_unit, help='Size unit\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-p', '--patched'), type='character', default=patched_value, help='Only report this `patched` value\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-L', '--legend'), type='character', default=legend_pos, help='Legend position\n\t\t[default=\'%default\']', metavar='character'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)
opt$pwr <- switch(opt$unit, B=0, KB=1, MB=2, GB=3, TB=4)
if (is.null(opt$pwr)) {
    opt$unit = 'MB'
    opt$pwr = 2
}

# Reading input file
d <- read.table(opt$file, sep=',', header=TRUE)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    # Data subsetting
    ds <- d[d$dataset == cdataset, c('pathno', 'patched', 'pindexsize')]
    # Removing duplicates
    uds <- unique(ds[order(ds$pathno),])
    if (opt$patched != 'all') {
        uds <- uds[!(uds$patched != opt$patched),]
    }
    # Converting dataframe to tibble
    #t <- as_data_frame(uds)
    t <- tibble(npath=uds$pathno,
                type=factor(ifelse(uds$patched=='yes', 'Patched', 'Full')),
                size=(uds$pindexsize/1024^opt$pwr))
    g <- ggplot(t, aes(x=npath, y=size, colour=type, group=type)) +
         geom_line(size=2) +
         geom_point(size=10) +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Paths',
              y=paste('Size (', opt$unit, ')', sep=''),
              colour='Type') +
         scale_x_continuous(breaks=unique(t$npath)) +
         theme_bw() +
         theme(axis.text=element_text(size=22),
               axis.title=element_text(size=24),
               panel.border=element_rect(fill=NA, colour='gray20', size=3),
               legend.position=opt$legend)
    ggsave(paste('pindexsize_', cdataset, '.pdf', sep=''))
}
