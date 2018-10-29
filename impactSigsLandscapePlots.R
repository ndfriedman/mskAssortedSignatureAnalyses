#written by Noah Friedman

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

dfOvarians <- read.table('/Users/friedman/Desktop/impactLandscapeCharts/impactDataDfOvarian.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfPancreatics <- read.table('/Users/friedman/Desktop/impactLandscapeCharts/impactDataDfPancreatic.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfBreasts <- read.table('/Users/friedman/Desktop/impactLandscapeCharts/impactDataDfBreast.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfProstates <- read.table('/Users/friedman/Desktop/impactLandscapeCharts/impactDataDfProstate.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfNonHrd <- read.table('/Users/friedman/Desktop/impactLandscapeCharts/impactDataDfNonHrd.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

plottingLevels <- c('mean_1', 'mean_APOBEC', 'mean_4', 'mean_7', 'mean_10', 'mean_11', 'mean_17', 'mean_MMR', 'other')
barColors = c(
  "#00FFFF", #age  
  "#FF0000", #brca
  "#FFA500", #smoking 
  "#FFF600", #Uv
  "#ADFF2F", #POLE
  "#2A52BE", #TMZ
  "#551A8B", #sig17
  "#267574", #mmr
  "#D3D3D3" #OTHER
)

finalPlotSansLegend <- ggarrange(
  make_bar_comparison(
    dfOvarians, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_3",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="mean_3",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="BRCA in ovarian (36.8% sufficient mutation burden)",
    mainColor="#FF1493",
    hideLegend = TRUE) + theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfOvarians, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="mean_3", 
    yAxisValParam="NmutAdj",
    yAxisFillParam="nMutClipped"),
  
  make_bar_comparison(
    dfPancreatics, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_3",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="mean_3",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    mainColor="#FF1493",
    title="BRCA in pancreatics (23.7% sufficient mutation burden)",
    hideLegend = TRUE) + theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfPancreatics, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="mean_3", 
    yAxisValParam="NmutAdj",
    yAxisFillParam="nMutClipped"),
  
  make_bar_comparison(
    dfBreasts, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_3",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="mean_3",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    mainColor="#FF1493",
    title="BRCA in breasts (31.5% sufficient mutation burden)",
    hideLegend = FALSE) + theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfBreasts, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="mean_3", 
    yAxisValParam="NmutAdj",
    yAxisFillParam="nMutClipped"),
  
  make_bar_comparison(
    dfProstates, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_3",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="mean_3",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="BRCA in prostates (22.7 % sufficient mutation burden)",
    mainColor="#FF1493",
    hideLegend = TRUE) + theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfProstates, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="mean_3", 
    yAxisValParam="NmutAdj",
    yAxisFillParam="nMutClipped"),
  
  make_bar_comparison(
    dfNonHrd, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_3",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="mean_3",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="BRCA in Non HRD cancers (44.5 % sufficient mutation burden)",
    mainColor="#FF1493",
    hideLegend = TRUE) + theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfNonHrd, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="mean_3", 
    yAxisValParam="NmutAdj",
    yAxisFillParam="nMutClipped"),
  
  heights=c(1,.2,
            1,.2,
            1,.2,
            1,.2,
            1,.2)
)

ggsave('~/Desktop/testNoah.pdf', plot=finalPlotSansLegend)


