#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Performance per occurrence'
p_subtitle <- 'k = %s; for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
legend_pos <- 'right'
patched_value <- 'all'  # show locino for both patched and full path index
stepsize_value <- 'all' # show all values for step size
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-p', '--patched'), type='character', default=patched_value, help='Only report this `patched` value\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-e', '--stepsize'), type='character', default=patched_value, help='Only report this `stepsize` value\n\t\t[default=\'%default\']', metavar='character'),
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
        ds <- d[d$dataset == cdataset & d$seedlen == cseedlen, c('pathno', 'chunksize', 'patched', 'stepsize', 'totalhits', 'seedfindingtime')]
        ds <- na.omit(ds)
        # Removing duplicates
        uds <- unique(ds[order(ds$chunksize),])
        if (opt$patched != 'all') {
            uds <- uds[!(uds$patched != opt$patched),]
        }
        if (opt$stepsize != 'all') {
            uds <- uds[!(uds$stepsize != opt$stepsize),]
        }
        # Converting dataframe to tibble
        #t <- as_data_frame(uds)
        t <- tibble(chunksize=factor(uds$chunksize),
                    stepsize=uds$stepsize,
                    type=factor(ifelse(uds$patched=='yes', 'Patched', 'Full')),
                    npath=factor(uds$pathno),
                    perfperocc=(uds$seedfindingtime / uds$totalhits * 10^6))
        g <- ggplot(t, aes(x=chunksize, y=perfperocc, group=npath, color=npath, shape=npath)) +
             geom_line(size=2) +
             geom_point(size=10) +
             labs(title=opt$title,
                  subtitle=sprintf(opt$subtitle, cseedlen, cdataset),
                  x='Chunk size',
                  y=expression('Average time per occurrence ('*mu*'s)'),
                  fill='Paths') +
             scale_y_continuous(limit=c(0, max(t$perfperocc)*1.1), expand=c(0, 0)) +
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
        ggsave(paste('perfperocc_', cdataset, '_', cseedlen, '.pdf', sep=''))
    }
}
