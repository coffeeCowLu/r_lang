#!/usr/bin/env Rscript
#author : luhuilong
#email  : luhl@singlerobio.com
#date   : 20200210

suppressMessages({
  library(argparser)
  library(ggpubr)
  library(ggthemes)
})

argv <- arg_parser(description = 'the volcano plot', name = 'volcano')
argv <- add_argument(argv,"--deg", help="the differentially expressed genes file")
argv <- add_argument(argv,"--labeln", help="the number of labels each", default=10)
argv <- add_argument(argv,"--name", help="the sample name")
argv <- add_argument(argv,"--outdir", help="the output directory")
argv <- parse_args(argv)

difgenes <- argv$deg
labeln <- as.numeric(argv$labeln)
name <- argv$name
outdir <- argv$outdir

# read DEG file
difgenes = read.table(difgenes, header = T, stringsAsFactors = F)
difgenes$symbol = rownames(difgenes)
difgenes$group = 'not-significant'
difgenes$group[which((difgenes$avg_logFC>0) &  (difgenes$pct.1 > difgenes$pct.2))] = 'up-regulated'
difgenes$group[which((difgenes$avg_logFC<0) &  (difgenes$pct.1 < difgenes$pct.2))] = 'down-regulated'
difgenes$dif2 = abs(difgenes$pct.1 - difgenes$pct.2)

# top n labels
difgenes$Label = ''
difgenes = difgenes[order(difgenes$avg_logFC),]
#upgenes = tail(difgenes$symbol[which(difgenes$group == 'up-regulated')], 10)
upgenes = tail(difgenes$symbol, labeln)
downgenes = head(difgenes$symbol, labeln)
#downgenes2 = head(difgenes$symbol[which(difgenes$group == 'down-regulated')], 10)
degtopn = c(upgenes, downgenes)
difgenes$Label[match(degtopn, difgenes$symbol)] = degtopn


##plot
if (length(unique(difgenes$group)) == 2){
	palette = c('#2f5688', '#CC0000')}else{
	palette = c('#2f5688', 'gray', '#CC0000')}

title = sprintf("Differentially expressed genes(%s)", name )
p = ggscatter(difgenes, x = 'avg_logFC', y = 'dif2', 
          color = 'group', palette = palette,
          ellipse =T, ellipse.level = .95,
          ellipse.type = "norm", ellipse.alpha = .1,
          ellipse.border.remove = FALSE, mean.point = FALSE,
          size = 1.7,
          label = difgenes$Label, font.label = 7, repel = T,
          title = title,
          ylab  = 'abs(pct.1-pct.2)') + 
  geom_vline(xintercept = c(-0.25, 0.25), linetype = 'dashed')+
  theme_base()+
  guides(color = guide_legend(override.aes = list(fill = NA, linetype = 0, reverse= F))) +
  theme(legend.title = element_blank(),
        legend.justification = c(.1, .97),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = 'bold', size = 12, vjust = .5, hjust = .5),
        axis.title = element_text(face = 'bold', size = 12, vjust = .5, hjust = .5),
        plot.background=element_blank())

pdf(paste(outdir,'/',name,'.deg.volcano.pdf',sep=''),title=name)
p
dev.off()

png(paste(outdir,'/',name,'.deg.volcano.png',sep=''),title=name)
p
dev.off()
