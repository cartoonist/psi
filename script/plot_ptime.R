#!/usr/bin/env Rscript
library('tibble')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'benchmark.csv'
p_title <- 'Path query time'
p_subtitle <- 'for simulated dataset "%s" from N. Deltacephalinicoli genome ("snd")'
# Command-line options parser
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-n', '--nofreads'), type='integer', help='Number of reads', metavar='integer'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)

if (is.null(opt$nofreads)) {
    print_help(opt_parser)
    stop('The number of reads must be specified.', call.=FALSE)
}

# Reading input file
d <- read.table(opt$file, sep=',', header=TRUE)
# Generating a plot for each dataset
for (cdataset in levels(d$dataset)) {
    # Data subsetting
    ds <- d[d$dataset == cdataset, c('pathno', 'chunksize', 'patched', 'stepsize', 'pathquerytimeperchunk')]
    ds <- na.omit(ds)
    # Removing duplicates
    uds <- unique(ds[order(ds$chunksize),])
    # Converting dataframe to tibble
    #t <- as_data_frame(uds)
    t <- tibble(chunksize=factor(uds$chunksize),
                stepsize=factor(uds$stepsize),
                type=factor(ifelse(uds$patched=='yes', 'Patched', 'Full')),
                npath=factor(uds$pathno),
                ptime=(ds$pathquerytimeperchunk * (opt$nofreads / ds$chunksize) / 10^6))
    g <- ggplot(t, aes(x=chunksize, y=ptime, fill=npath)) +
         geom_col(colour='#444444', position=position_dodge()) +
         geom_text(aes(label=round(ptime, 2)), position=position_dodge(width=0.9), vjust=0, hjust=-0.25, angle=45, size=2) +
         facet_grid(type ~ stepsize, labeller=label_both) +
         labs(title=opt$title,
              subtitle=sprintf(opt$subtitle, cdataset),
              x='Chunk size',
              y='Path query time (s)',
              fill='Paths') +
         scale_y_continuous(limits=c(0,max(t$ptime)*1.1), expand=c(0,0)) +
         theme_bw()
    ggsave(paste('ptime_', cdataset, '.pdf', sep=''))
}
