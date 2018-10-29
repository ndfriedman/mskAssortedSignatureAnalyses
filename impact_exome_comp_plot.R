#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

make_correlation_plot <- function(df, corCol1, corCol2, titleString=''){
  df <- my.rename(df, corCol1, "cc1")
  df <- my.rename(df, corCol2, "cc2")
  
  Cor <- cor.test(df$cc1, df$cc2, method = "pearson", conf.level = 0.95)
  corEstimate <- Cor$"estimate"
  title <- paste(titleString, "Ncases: ",  nrow(df), " Correlation: ", round(corEstimate, digits = 2))
  
  scatter_plot <- ggplot(df, aes(cc1, cc2))+
  xlim(0,1)+
  ylim(0,1)
  scatter_plot <- scatter_plot + geom_point() + labs(x = corCol1, y = corCol2) + geom_smooth(method="lm")
  scatter_plot <- scatter_plot + ggtitle(title)
  
  return(scatter_plot)
}

dfSigsComp <- read.table('~/Desktop/sigsComp.tsv',sep = '\t', header=TRUE)
#scatterComp(dfSigsComp)

nrow(dfSigsComp)

dfLessThan4Muts <- dfSigsComp[dfSigsComp$Nmut_y < 4,]
df4To8Muts <- dfSigsComp[dfSigsComp$Nmut_y <8 & dfSigsComp$Nmut_y > 3,]
df8To12Muts <- dfSigsComp[dfSigsComp$Nmut_y <12 & dfSigsComp$Nmut_y > 7,]
df12To20Muts <- dfSigsComp[dfSigsComp$Nmut_y <20 & dfSigsComp$Nmut_y > 11,]
df20To40Muts <- dfSigsComp[dfSigsComp$Nmut_y <40 & dfSigsComp$Nmut_y > 19,]
df40PlusMuts <- dfSigsComp[dfSigsComp$Nmut_y > 39,]

finalPlot <- plot_grid(
  #row1
  make_correlation_plot(dfLessThan4Muts, "exome_mean_1", "impact_mean_1",titleString='<4muts'),
  make_correlation_plot(df4To8Muts, "exome_mean_1", "impact_mean_1",titleString='4-7muts'),
  make_correlation_plot(df8To12Muts, "exome_mean_1", "impact_mean_1",titleString='8-11muts'),
  make_correlation_plot(df12To20Muts, "exome_mean_1", "impact_mean_1",titleString='12-19muts'),
  make_correlation_plot(df20To40Muts, "exome_mean_1", "impact_mean_1",titleString='20-39muts'),
  make_correlation_plot(df40PlusMuts, "exome_mean_1", "impact_mean_1",titleString='>39muts'),
  
  make_correlation_plot(dfLessThan4Muts, "exome_mean_APOBEC", "impact_mean_APOBEC",titleString='<4muts'),
  make_correlation_plot(df4To8Muts, "exome_mean_APOBEC", "impact_mean_APOBEC",titleString='4-7muts'),
  make_correlation_plot(df8To12Muts, "exome_mean_APOBEC", "impact_mean_APOBEC",titleString='8-11muts'),
  make_correlation_plot(df12To20Muts, "exome_mean_APOBEC", "impact_mean_APOBEC",titleString='12-19muts'),
  make_correlation_plot(df20To40Muts, "exome_mean_APOBEC", "impact_mean_APOBEC",titleString='20-39muts'),
  make_correlation_plot(df40PlusMuts, "exome_mean_APOBEC", "impact_mean_APOBEC",titleString='>39muts'),
  
  make_correlation_plot(dfLessThan4Muts, "exome_mean_3", "impact_mean_3",titleString='<4muts'),
  make_correlation_plot(df4To8Muts, "exome_mean_3", "impact_mean_3",titleString='4-7muts'),
  make_correlation_plot(df8To12Muts, "exome_mean_3", "impact_mean_3",titleString='8-11muts'),
  make_correlation_plot(df12To20Muts, "exome_mean_3", "impact_mean_3",titleString='12-19muts'),
  make_correlation_plot(df20To40Muts, "exome_mean_3", "impact_mean_3",titleString='20-39muts'),
  make_correlation_plot(df40PlusMuts, "exome_mean_3", "impact_mean_3",titleString='>39muts'),
  
  make_correlation_plot(dfLessThan4Muts, "exome_mean_SMOKING", "impact_mean_SMOKING",titleString='<4muts'),
  make_correlation_plot(df4To8Muts, "exome_mean_SMOKING", "impact_mean_SMOKING",titleString='4-7muts'),
  make_correlation_plot(df8To12Muts, "exome_mean_SMOKING", "impact_mean_SMOKING",titleString='8-11muts'),
  make_correlation_plot(df12To20Muts, "exome_mean_SMOKING", "impact_mean_SMOKING",titleString='12-19muts'),
  make_correlation_plot(df20To40Muts, "exome_mean_SMOKING", "impact_mean_SMOKING",titleString='20-39muts'),
  make_correlation_plot(df40PlusMuts, "exome_mean_SMOKING", "impact_mean_SMOKING",titleString='>39muts'),
  
  make_correlation_plot(dfLessThan4Muts, "exome_mean_5", "impact_mean_5",titleString='<4muts'),
  make_correlation_plot(df4To8Muts, "exome_mean_5", "impact_mean_5",titleString='4-7muts'),
  make_correlation_plot(df8To12Muts, "exome_mean_5", "impact_mean_5",titleString='8-11muts'),
  make_correlation_plot(df12To20Muts, "exome_mean_5", "impact_mean_5",titleString='12-19muts'),
  make_correlation_plot(df20To40Muts, "exome_mean_5", "impact_mean_5",titleString='20-39muts'),
  make_correlation_plot(df40PlusMuts, "exome_mean_5", "impact_mean_5",titleString='>39muts'),
  
  make_correlation_plot(dfLessThan4Muts, "exome_mean_MMR", "impact_mean_MMR",titleString='<4muts'),
  make_correlation_plot(df4To8Muts, "exome_mean_MMR", "impact_mean_MMR",titleString='4-7muts'),
  make_correlation_plot(df8To12Muts, "exome_mean_MMR", "impact_mean_MMR",titleString='8-11muts'),
  make_correlation_plot(df12To20Muts, "exome_mean_MMR", "impact_mean_MMR",titleString='12-19muts'),
  make_correlation_plot(df20To40Muts, "exome_mean_MMR", "impact_mean_MMR",titleString='20-39muts'),
  make_correlation_plot(df40PlusMuts, "exome_mean_MMR", "impact_mean_MMR",titleString='>39muts'),
  
  make_correlation_plot(dfLessThan4Muts, "exome_mean_7", "impact_mean_7",titleString='<4muts'),
  make_correlation_plot(df4To8Muts, "exome_mean_7", "impact_mean_7",titleString='4-7muts'),
  make_correlation_plot(df8To12Muts, "exome_mean_7", "impact_mean_7",titleString='8-11muts'),
  make_correlation_plot(df12To20Muts, "exome_mean_7", "impact_mean_7",titleString='12-19muts'),
  make_correlation_plot(df20To40Muts, "exome_mean_7", "impact_mean_7",titleString='20-39muts'),
  make_correlation_plot(df40PlusMuts, "exome_mean_7", "impact_mean_7",titleString='>39muts'),
  
   
                   align='hv', ncol=6, nrow=7,
                   rel_widths = c(1,1,1,1,1),
                   rel_heights = c(1,1,1,1,1,1,1),
                   scale = 1)


ggsave('~/Desktop/testNoah.pdf', plot=finalPlot, width = 49, height = 49, units = c("in"))








