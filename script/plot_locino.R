#!/usr/bin/env Rscript
library('tibble')
library('reshape2')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Starting loci'
p_subtitle <- 'k = %s, p = \'%s\'; for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
legend_pos <- 'right'
patched_value <- 'all'  # show locino for both patched and full path index
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-p', '--patched'), type='character', default=patched_value, help='Only report this `patched` value\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-L', '--legend'), type='character', default=legend_pos, help='Legend position\n\t\t[default=\'%default\']', metavar='character'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)

# Reading input file
d <- read.table(opt$file, sep=',', header=TRUE)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    for (cseedlen in unique(d$seedlen)) {
        # Data subsetting
        ds <- d[d$dataset == cdataset & d$seedlen == cseedlen, c('pathno', 'stepsize', 'patched', 'locino', 'uniqnodes')]
        # Removing duplicates
        uds <- unique(ds[order(ds$pathno),])
        # Melting
        uds <- melt(uds, id.vars=c('pathno', 'stepsize', 'patched'))
	if (opt$patched != 'all') {
            uds <- uds[!(uds$patched != opt$patched),]
	}
        # Converting dataframe to tibble
        #t <- as_data_frame(uds)
        t <- tibble(npath=uds$pathno,
                    stepsize=factor(uds$stepsize),
                    type=factor(ifelse(uds$patched == 'yes', 'Patched', 'Full')),
                    variable=factor(ifelse(uds$variable == 'locino', 'Loci', 'Nodes')),
                    loci=uds$value)
        g <- ggplot(t, aes(x=npath, y=loci, shape=stepsize, colour=variable, group=interaction(stepsize, variable))) +
             geom_line(size=2) +
             geom_point(size=10) +
             labs(title=opt$title,
                  subtitle=sprintf(opt$subtitle, cseedlen, opt$patched, cdataset),
                  x='Paths',
                  y='Counts',
                  colour='Step size',
                  linetype='Starting') +
             scale_x_continuous(breaks=unique(t$npath)) +
             theme_bw() +
             theme(axis.text=element_text(size=22),
                   axis.title=element_text(size=24),
                   panel.border=element_rect(fill=NA, colour='gray20', size=3),
                   legend.position=opt$legend)
        if (opt$patched == 'all') {
            g <- g + facet_grid(. ~ type)
        }
        ggsave(paste('locino_', cdataset, '_', cseedlen, '.pdf', sep=''))
    }
}
