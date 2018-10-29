#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(ggpubr)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

plotBar <- function(df, title){
  plt<- ggplot(df, aes(x = variable, y=value))+
    geom_bar(stat='identity')+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=5))+
    ggtitle(title)+
    get_empty_theme()+
    ylim(0,1)
  return(plt)
}

plotBoxPlot <- function(df){
  p<- ggplot(df, aes(x=Amino_Acid_Change, y=nMutAttributedToSmoking)) + 
    geom_boxplot()+
    #ylim(0,100)
    coord_cartesian(ylim=c(0, 120))
  
  my_comparisons <- list( 
                          c("G12V", "G12D"),
                          c("G12A", "G12C"),
                          c("G12C", "G12D"),
                          c("G12C", "G12R"),
                          c("G12C", "G12S"),
                          c("G12C", "G12V"),
                          c("G13C", "G12C"), 
                          c("G12D", "G13D"))
  p<- p + stat_compare_means(comparisons = my_comparisons,label.y=c(75,80,85,90, 95,100,105,110,115), tip.length=.01)#
  p<- p+ geom_jitter(shape=16, position=position_jitter(0.1))+
  ggtitle("Number of mutations attributable to smoking by KRAS codon")
  
  return(p)
}

plot_signature_bars <- function(df){
  ggplot(df, aes(x=Amino_Acid_Change, y=smokingSpectraProb)) + 
    geom_bar(stat= "summary", fun.y = "mean")+
    get_empty_theme()
}

df <- read.table('~/Desktop/hotspotPrevalenceAnalysis/ratioData_kras_lung.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfG12D <- df[df$hotspot == 'G12D',]
dfG12C <- df[df$hotspot == 'G12C',]
dfG12V <- df[df$hotspot == 'G12V',]
dfG12R <- df[df$hotspot == 'G12R',]
dfG12A <- df[df$hotspot == 'G12A',]
dfG12S <- df[df$hotspot == 'G12S',]
dfG13D <- df[df$hotspot == 'G13D',]
dfG13C <- df[df$hotspot == 'G13C',]

fullDf <- read.table('~/Desktop/krasMutsTemp.tsv',sep = '\t', header=TRUE) 
plotBoxPlot(fullDf)
plot_signature_bars(fullDf)



plt<-plot_grid(
  plotBoxPlot(fullDf),
  plot_signature_bars(fullDf),
  align='hv', nrow=2, rel_heights = c(1,.5)
)

plt<- ggarrange(
  plotBar(dfG12C, paste("G12C", "// n cases: ")),
  plotBar(dfG12D, paste("G12D", "// n cases: ")),
  plotBar(dfG12V, paste("G12V", "// n cases: ")),
  plotBar(dfG12A, paste("G12A", "// n cases: ")),
  plotBar(dfG12S, paste("G12S", "// n cases: ")),
  plotBar(dfG12R, paste("G12R", "// n cases: ")),
  plotBar(dfG13C, paste("G13C", "// n cases: ")),
  plotBar(dfG13D, paste("G13D", "// n cases: ")),
  heights = c(.2,.2,.2,.2,.2,.2,.2,.2)
)

ggsave('~/Desktop/testNoah.pdf', plot=plt)
