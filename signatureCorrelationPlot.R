#written by noah friedman

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
require(plyr)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/plottingScripts/corPlotUtil.R")

dfJZ = read.table('/Users/friedman/Desktop/myTestDf.tsv',sep = '\t', header=TRUE)
dfJYY = read.table('~/Desktop/myTestDf2.tsv',sep = '\t', header=TRUE)

colnames(dfJZ)

dfJZ$mean_1
dfJZ$frac_Signature.1_driver

x = c('sig.1', 
      'sig.APOBEC',
      'sig.MMR',
      'sig.7',
      'sig.4',
      'sig.10',
      'sig.11',
      'sig.14',
      'sig.17',
      'sig.18',
      'sig.22',
      'sig.23',
      'sig.27')
y = c(make_correlation_plot(dfJZ, "mean_1", "frac_Signature.1_driver"),
      make_correlation_plot(dfJZ, "mean_APOBEC", "frac_Signature.APOBEC_driver"),
      make_correlation_plot(dfJZ, "mean_MMR", "frac_Signature.MMR_driver"),
      make_correlation_plot(dfJZ, "mean_4", "frac_Signature.4_driver"),
      make_correlation_plot(dfJZ, "mean_7", "frac_Signature.7_driver"),
      make_correlation_plot(dfJZ, "mean_10", "frac_Signature.10_driver"),
      make_correlation_plot(dfJZ, "mean_11", "frac_Signature.11_driver"),
      make_correlation_plot(dfJZ, "mean_14", "frac_Signature.14_driver"),
      make_correlation_plot(dfJZ, "mean_17", "frac_Signature.17_driver"),
      make_correlation_plot(dfJZ, "mean_18", "frac_Signature.18_driver"),
      make_correlation_plot(dfJZ, "mean_22", "frac_Signature.22_driver"),
      make_correlation_plot(dfJZ, "mean_23", "frac_Signature.23_driver"),
      make_correlation_plot(dfJZ, "mean_27", "frac_Signature.27_driver"))
df <- data.frame(x,y)
p1 <- ggplot(df, aes(x, y)) + geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle("driver mut in signature context/signature magnitude correlation")+
  ylim(0,1)+
  theme( #rotate axes, change size
    axis.title.x=element_blank(), #make there be no x title
    axis.title.y=element_text(size=4),
    legend.title=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank()
  )+
  theme(plot.title = element_text(size = 5))



colnames(dfJYY)


x = c('sig.1', 
      'sig.APOBEC',
      'sig.MMR',
      'sig.7',
      'sig.4',
      'sig.10',
      'sig.11',
      'sig.14',
      'sig.17',
      'sig.18',
      'sig.22',
      'sig.23',
      'sig.27')
y = c(make_correlation_plot(dfJYY, "mean_1", "frac_Signature.1_snp"),
      make_correlation_plot(dfJYY, "mean_APOBEC", "frac_Signature.APOBEC_snp"),
      make_correlation_plot(dfJYY, "mean_MMR", "frac_Signature.MMR_snp"),
      make_correlation_plot(dfJYY, "mean_4", "frac_Signature.4_snp"),
      make_correlation_plot(dfJYY, "mean_7", "frac_Signature.7_snp"),
      make_correlation_plot(dfJYY, "mean_10", "frac_Signature.10_snp"),
      make_correlation_plot(dfJYY, "mean_11", "frac_Signature.11_snp"),
      make_correlation_plot(dfJYY, "mean_14", "frac_Signature.14_snp"),
      make_correlation_plot(dfJYY, "mean_17", "frac_Signature.17_snp"),
      make_correlation_plot(dfJYY, "mean_18", "frac_Signature.18_snp"),
      make_correlation_plot(dfJYY, "mean_22", "frac_Signature.22_snp"),
      make_correlation_plot(dfJYY, "mean_23", "frac_Signature.23_snp"),
      make_correlation_plot(dfJYY, "mean_27", "frac_Signature.27_snp"))
df <- data.frame(x,y)
p2 <- ggplot(df, aes(x, y)) + geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ggtitle("SNP in signature context/signature magnitude correlation")+
  ylim(0,1)+
  theme(plot.title = element_text(size = 5))

ggarrange(p1, p2)

dfJYY$frac_Signature.10_snp
dfJZ$frac_Signature.10_driver

















y = c(make_correlation_plot(dfJZ, "confidence_1", "frac_Signature.1_Driver"),
      make_correlation_plot(dfJZ, "confidence_APOBEC", "frac_Signature.APOBEC_Driver"),
      make_correlation_plot(dfJZ, "confidence_MMR", "frac_Signature.MMR_Driver"),
      make_correlation_plot(dfJZ, "confidence_4", "frac_Signature.4_Driver"),
      make_correlation_plot(dfJZ, "confidence_7", "frac_Signature.7_Driver"),
      make_correlation_plot(dfJZ, "confidence_10", "frac_Signature.10_Driver"),
      make_correlation_plot(dfJZ, "confidence_11", "frac_Signature.11_Driver"),
      make_correlation_plot(dfJZ, "confidence_14", "frac_Signature.14_Driver"),
      make_correlation_plot(dfJZ, "confidence_17", "frac_Signature.17_Driver"),
      make_correlation_plot(dfJZ, "confidence_18", "frac_Signature.18_Driver"),
      make_correlation_plot(dfJZ, "confidence_22", "frac_Signature.22_Driver"),
      make_correlation_plot(dfJZ, "confidence_23", "frac_Signature.23_Driver"),
      make_correlation_plot(dfJZ, "confidence_27", "frac_Signature.27_Driver"))
df <- data.frame(x,y)
ggplot(df, aes(x, y)) + geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 60, hjust = 1)) + ggtitle("driver mut in signature context/signature magnitude correlation")



