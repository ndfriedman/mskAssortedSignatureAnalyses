#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")


plotBoxPlot <- function(df, xValParam, yColParam, comps=NA, compYPositions=NA, maxY, title, yLabelText){
  
  df <- my.rename(df, xValParam, "xVal")
  df <- my.rename(df, yColParam, "yCol")
  
  p<- ggplot(df, aes(reorder(xVal, plotOrdering), y=yCol)) +
    geom_boxplot()+
    coord_cartesian(ylim=c(0, maxY))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ylab(yLabelText)
  
  #p<- p + stat_compare_means(comparisons = comps, label.y=compYPositions, tip.length=.01)#
  p<- p+ geom_jitter(shape=16, position=position_jitter(0.1))+
    ggtitle(title)
  return(p)
}


df <- read.table('~/Desktop/dataForLocalPlotting/mutburdenBoxplot.tsv',sep = '\t', header=TRUE)



numberOfRows =2
numberOfColumns =3
plt<-plot_grid(
  plotBoxPlot(df, 'label', 'nHotspots', maxY=max(df$nHotspots), title='N Hotspots', yLabelText = 'n hotspot mutations'),
  plotBoxPlot(df, 'label', 'nOncogenicMutations', maxY=max(df$nOncogenicMutation) + 5, title='N Oncogenic Drivers', yLabelText = 'n oncoKb possibly oncogenic'),
  plotBoxPlot(df, 'label', 'nOncogenicOrHotspotMutations', maxY=max(df$nOncogenicOrHotspotMutations) + 5, title='N Drivers', yLabelText = 'n oncoKb possibly oncogenic or hotspot mutations'),
  plotBoxPlot(df, 'label', 'fracHotspotsAtEnrichedMotif', maxY=1.05, title='Driver mutations at dominant sig motif', yLabelText = 'fraction driver mutations at Aging vs Hypermutator Sig Motif'),
  plotBoxPlot(df, 'label', 'Nmut', maxY=750, title='Mutation load', yLabelText = 'number of mutations in impact'),
  align='hv', nrow=numberOfRows, ncol=numberOfColumns
  )
  

ggsave('~/Desktop/testNoah.pdf', plot=plt, width = 7*numberOfColumns, height = 7*numberOfRows, units = c("in"))









