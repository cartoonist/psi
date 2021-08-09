#!/usr/bin/env Rscript
library("tibble")
library("reshape2")
library("ggplot2")
library("optparse")
library("tidyverse")
library("viridis")
library("hrbrthemes")
library("showtext")
#library("ggrepel")

## Adding Roboto Condensed and Roboto Condensed Light fonts
font_add_google(name="Roboto Condensed", family="Roboto Condensed")
font_add_google(name="Roboto Condensed", family="Roboto Condensed Light")
showtext_auto()

## Default values
in_file <- 'analyse.csv'
p_title <- 'GA/PSI Mapping Performance'
p_subtitle <- 'for "%s" dataset'
p_trim <- 0
## Command-line options parser
option_list <- list(
  make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
  make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
  make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
  make_option(c('-n', '--nofreads'), type='integer', help='Number of reads', metavar='integer'),
  make_option(c('-T', '--trim'), type='integer', default=p_trim, help='Trim this percentage of the first (high-confidence) group of mapped reads', metavar='integer')
)
opt_parser <- OptionParser(option_list=option_list)
## Parsing command-line arguments
opt <- parse_args(opt_parser)

if (is.null(opt$nofreads)) {
  print_help(opt_parser)
  stop('The number of reads must be specified.', call.=FALSE)
}
## Reading input file
data <- read.table(opt$file, header=TRUE)
data$NNAL <- rep(opt$nofreads, dim(data)[1]) - data$NRAL
dataset <- unique(data$dataset)
val_vars <- c("NVAL", "NINV")
map_vars <- c("NHUP", "NHMP", "NLUP", "NLMP", "NHUS", "NHDS", "NHMS", "NLUS", "NLDS", "NLMS", "NNAL")
tru_vars <- c("NFFM", "NFPM", "NPPM", "NFNM", "NPNM", "NNNM", "NNAL")
#map_lab <- c("NHMP", "NLMP", "NHMS", "NLMS")
d_val <- melt(data[c("run", val_vars)], id.vars="run")
d_val <- d_val %>% group_by(run) %>% mutate(cumsum = cumsum(value))
d_val$cum_percent <- d_val$cumsum / opt$nofreads
d_map <- melt(data[c("run", map_vars)], id.vars="run")
d_map <- d_map %>% group_by(run) %>% mutate(cumsum = cumsum(value))
d_map$cum_percent <- d_map$cumsum / opt$nofreads
#d_map[!(d_map$variable %in% map_lab),]$label <- NA
d_tru <- melt(data[c("run", tru_vars)], id.vars="run")
d_tru <- d_tru %>% group_by(run) %>% mutate(cumsum = cumsum(value))
d_tru$cum_percent <- d_tru$cumsum / opt$nofreads
min_hup <- min(d_map[d_map$variable == "NHUP",]$value)
min_ylim <- round(min_hup / opt$nofreads * opt$trim)
#palette <- c("#B93556FF", "#BD3853FF", "#C13A50FF", "#C63D4DFF",  # paired
#             "#3579A2FF", "#357CA3FF", "#347FA4FF", "#3483A5FF", "#3486A5FF", "#3489A6FF",  # single
#             "gray"  # unaligned
#             )
g <- ggplot(d_map, aes(fill=factor(variable, levels=rev(map_vars)), y=(value*100/opt$nofreads), x=as.factor(run))) +
  geom_bar(position="stack", stat="identity") +
  labs(title=opt$title,
       subtitle=sprintf(opt$subtitle, dataset),
       fill="Category", x="Run", y="% reads") +
  scale_fill_viridis(discrete=TRUE, direction=-1) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
#  geom_text_repel(aes(y=cumsum, label=cum_percent),
#                  vjust=1.5,
#                  nudge_x=0.22,
#                  colour="gray",
#                  direction="y") +
  coord_cartesian(ylim=c(min_ylim, 100)) +
  theme_ipsum_rc()
ggsave(paste('mapping_performance_', dataset, '.pdf', sep=''))

min_nffm <- min(d_tru[d_tru$variable == "NFFM",]$value)
min_ylim <- round(min_nffm / opt$nofreads * opt$trim)
g <- ggplot(d_tru, aes(fill=factor(variable, levels=rev(tru_vars)), y=(value*100/opt$nofreads), x=as.factor(run))) +
  geom_bar(position="stack", stat="identity") +
  labs(title=opt$title,
       subtitle=paste(sprintf(opt$subtitle, dataset), "(compared to the ground truth)"),
       fill="Category", x="Run", y="% reads") +
  scale_fill_viridis(discrete=TRUE, direction=-1) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  coord_cartesian(ylim=c(min_ylim, 100)) +
  theme_ipsum_rc()
ggsave(paste('mapping_performance_ground_truth_', dataset, '.pdf', sep=''))
