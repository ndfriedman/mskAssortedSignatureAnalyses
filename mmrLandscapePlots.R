#written by Noah friedman

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

topMargin <- unit(c(0,-.5,-.5,-.5), "lines")
myStandardMargin <- unit(c(.3,-.5,-.5,-.5), "lines")

plottingLevels <- c('mean_1', 'mean_APOBEC', 'mean_4', 'mean_7', 'mean_10', 'mean_11', 'mean_17', 'mean_3', 'mean_14', 'other')
barColors = c(
  "#00FFFF", #age  
  "#FF0000", #APOBEC
  "#FFA500", #smoking 
  "#FFF600", #Uv
  "#ADFF2F", #POLE
  "#2A52BE", #TMZ
  "#551A8B", #sig17
  "#FF1493", #brca
  "#00DD5D", #sig14
  "#D3D3D3" #OTHER
)

plottingLevelsCancerType <- c('Colorectal Cancer', 'Endometrial Cancer', 'Prostate Cancer', 'Breast Cancer', 'Non-Small Cell Lung Cancer', 'Pancreatic Cancer', 'Esophagogastric Cancer', 'Glioma', 'Hepatobiliary Cancer', 'Other')
cancerTypeColors = c(
  "#00DFFF", #COLORECTAL  
  "#9932CC", #ENDOMETRIAL
  "#FFA500", #PROSTATE
  "#FF1493", #BREAST
  "#FF4500", #LUNG CANCER
  "#32CD32", #PANCREATIC CANCER
  "#FF0000", #ESOPH cancer
  "#003366", #GLIOMA
  "#00DD5D", #LIVER
  "#D3D3D3" #OTHER
)

dfMMRpos <- read.table('/Users/friedman/Desktop/otherLandscapePlot/exomeMmrPos.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfMMRneg <- read.table('/Users/friedman/Desktop/otherLandscapePlot/exomeMmrNeg.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

finalPlotSansLegend <- ggarrange(
  make_bar_comparison(
    dfMMRpos, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="mean_MMR",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="MMR gene mutated cancers",
    mainColor="#267574",
    hideLegend = FALSE)+ theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfMMRpos, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="mean_MMR", 
    yAxisValParam="NmutAdj",
    yAxisFillParam="nMutClipped"),
  
  make_bar_comparison(
    dfMMRneg, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="mean_MMR",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="MMR gene wildtype cancers",
    mainColor="#267574",
    hideLegend = TRUE)+ theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfMMRneg, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="mean_MMR", 
    yAxisValParam="NmutAdj",
    yAxisFillParam="nMutClipped"),
    
  heights=c(1,.2,
      1,.2)
)

legendSigs <- get_legend(make_bar_comparison(
  dfMMRpos, 
  xAxisValParam="Tumor_Sample_Barcode",
  yAxisValParam1="mean_MMR",
  yAxisValParam2="otherPredominantSigMagnitude",
  yAxisFillParam="otherPredominantSigName",
  orderingValParam="mean_MMR",
  coloringSpecCols=plottingLevels, barColorPalette=barColors,
  title="MMR gene wildtype cancers",
  mainColor="#267574",
  hideLegend = FALSE))

finalPlot <- plot_grid(finalPlotSansLegend, legendSigs,
                       align='hv', ncol=2,
                       rel_widths = c(1,.2), scale = 0.9)

ggsave('~/Desktop/testNoah.pdf', plot=finalPlot)

####################################MMR impact

dfMLH1Mut <- read.table('/Users/friedman/Desktop/otherLandscapePlot/impactMLH1Mut.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfMSH2Mut <- read.table('/Users/friedman/Desktop/otherLandscapePlot/impactMSH2Mut.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfMSH3Mut <- read.table('/Users/friedman/Desktop/otherLandscapePlot/impactMSH3Mut.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfMSH6Mut <- read.table('/Users/friedman/Desktop/otherLandscapePlot/impactMSH6Mut.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfPMS1Mut <- read.table('/Users/friedman/Desktop/otherLandscapePlot/impactPMS1Mut.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfPMS2Mut <- read.table('/Users/friedman/Desktop/otherLandscapePlot/impactPMS2Mut.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
dfMSIWildtype <- read.table('/Users/friedman/Desktop/otherLandscapePlot/impactMSIWildtype.tsv',sep = '\t', header=TRUE)

finalPlotSansLegend <- ggarrange(
  make_bar_comparison(
    dfMLH1Mut, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="orderingVal",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="MLH1 gene mutated cancers",
    mainColor="#267574",
    hideLegend = TRUE) + theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfMLH1Mut, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="orderingVal", 
    yAxisValParam="Nmut",
    yAxisFillParam="#000000")+ theme(plot.margin = myStandardMargin),
  generate_gradient_tiles(dfMLH1Mut, 
                          xAxisValParam="Tumor_Sample_Barcode", 
                          xAxisOrderingParam="orderingVal", 
                          fillArgParam="msiScore")+ theme(plot.margin = myStandardMargin),
  generate_ggplot_tiles(dfMLH1Mut, xAxisValParam="Tumor_Sample_Barcode", xAxisOrderingParam="orderingVal",
                        fillArgParam="cancer_type_adjusted", hideLegend = TRUE, mode = "cancerTypeTiles", 
                        textSizeParam=.5, tileColorPalette=cancerTypeColors, orderSpecCols=plottingLevelsCancerType),
  
  make_bar_comparison(
    dfMSH2Mut, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="orderingVal",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="MSH2 gene mutated cancers",
    mainColor="#267574",
    hideLegend = TRUE) + theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfMSH2Mut, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="orderingVal", 
    yAxisValParam="Nmut",
    yAxisFillParam="#000000")+ theme(plot.margin = myStandardMargin),
  generate_gradient_tiles(dfMSH2Mut, 
                          xAxisValParam="Tumor_Sample_Barcode", 
                          xAxisOrderingParam="orderingVal", 
                          fillArgParam="msiScore")+ theme(plot.margin = myStandardMargin),
  generate_ggplot_tiles(dfMSH2Mut, xAxisValParam="Tumor_Sample_Barcode", xAxisOrderingParam="orderingVal",
                        fillArgParam="cancer_type_adjusted", hideLegend = TRUE, mode = "cancerTypeTiles", 
                        textSizeParam=.5, tileColorPalette=cancerTypeColors, orderSpecCols=plottingLevelsCancerType),
  
  make_bar_comparison(
    dfMSH3Mut, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="orderingVal",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="MSH3 gene mutated cancers",
    mainColor="#267574",
    hideLegend = TRUE)+ theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfMSH3Mut, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="orderingVal", 
    yAxisValParam="Nmut",
    yAxisFillParam="#000000")+ theme(plot.margin = myStandardMargin),
  generate_gradient_tiles(dfMSH3Mut, 
                          xAxisValParam="Tumor_Sample_Barcode", 
                          xAxisOrderingParam="orderingVal", 
                          fillArgParam="msiScore")+ theme(plot.margin = myStandardMargin),
  generate_ggplot_tiles(dfMSH3Mut, xAxisValParam="Tumor_Sample_Barcode", xAxisOrderingParam="orderingVal",
                        fillArgParam="cancer_type_adjusted", hideLegend = TRUE, mode = "cancerTypeTiles", 
                        textSizeParam=.5, tileColorPalette=cancerTypeColors, orderSpecCols=plottingLevelsCancerType),
  
  make_bar_comparison(
    dfMSH6Mut, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="orderingVal",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="MSH6 gene mutated cancers",
    mainColor="#267574",
    hideLegend = TRUE)+ theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfMSH6Mut, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="orderingVal", 
    yAxisValParam="Nmut",
    yAxisFillParam="#000000")+ theme(plot.margin = myStandardMargin),
  generate_gradient_tiles(dfMSH6Mut, 
                          xAxisValParam="Tumor_Sample_Barcode", 
                          xAxisOrderingParam="orderingVal", 
                          fillArgParam="msiScore")+ theme(plot.margin = myStandardMargin),
  generate_ggplot_tiles(dfMSH6Mut, xAxisValParam="Tumor_Sample_Barcode", xAxisOrderingParam="orderingVal",
                        fillArgParam="cancer_type_adjusted", hideLegend = TRUE, mode = "cancerTypeTiles", 
                        textSizeParam=.5, tileColorPalette=cancerTypeColors, orderSpecCols=plottingLevelsCancerType),
  
  make_bar_comparison(
    dfPMS1Mut, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="orderingVal",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="PMS1 gene mutated cancers",
    mainColor="#267574",
    hideLegend = TRUE)+ theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfPMS1Mut, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="orderingVal", 
    yAxisValParam="Nmut",
    yAxisFillParam="#000000")+ theme(plot.margin = myStandardMargin),
  generate_gradient_tiles(dfPMS1Mut, 
                          xAxisValParam="Tumor_Sample_Barcode", 
                          xAxisOrderingParam="orderingVal", 
                          fillArgParam="msiScore")+ theme(plot.margin = myStandardMargin),
  generate_ggplot_tiles(dfPMS1Mut, xAxisValParam="Tumor_Sample_Barcode", xAxisOrderingParam="mean_MMR",
                        fillArgParam="cancer_type_adjusted", hideLegend = TRUE, mode = "cancerTypeTiles", 
                        textSizeParam=.5, tileColorPalette=cancerTypeColors, orderSpecCols=plottingLevelsCancerType),
  
  make_bar_comparison(
    dfPMS2Mut, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="orderingVal",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="PMS2 gene mutated cancers",
    mainColor="#267574",
    hideLegend = TRUE)+ theme(plot.margin = topMargin),
  generate_mut_burden_bar(
    dfPMS2Mut, 
    xAxisValParam="Tumor_Sample_Barcode", 
    xAxisOrderingParam="orderingVal", 
    yAxisValParam="Nmut",
    yAxisFillParam="#000000")+ theme(plot.margin = myStandardMargin),
  generate_gradient_tiles(dfPMS2Mut, 
                          xAxisValParam="Tumor_Sample_Barcode", 
                          xAxisOrderingParam="orderingVal", 
                          fillArgParam="msiScore")+ theme(plot.margin = myStandardMargin),
  generate_ggplot_tiles(dfPMS2Mut, xAxisValParam="Tumor_Sample_Barcode", xAxisOrderingParam="orderingVal",
                        fillArgParam="cancer_type_adjusted", hideLegend = TRUE, mode = "cancerTypeTiles", 
                        textSizeParam=.5, tileColorPalette=cancerTypeColors, orderSpecCols=plottingLevelsCancerType),
  
  #make_bar_comparison(
  #  dfMSIWildtype, 
  #  xAxisValParam="Tumor_Sample_Barcode",
  #  yAxisValParam1="mean_MMR",
  #  yAxisValParam2="otherPredominantSigMagnitude",
  #  yAxisFillParam="otherPredominantSigName",
  #  orderingValParam="orderingVal",
  #  coloringSpecCols=plottingLevels, barColorPalette=barColors,
  #  title="MSI gene wildtype",
  #  mainColor="#267574",
  #  hideLegend = TRUE)+ theme(plot.margin = topMargin),
  #generate_mut_burden_bar(
  #  dfMSIWildtype, 
  #  xAxisValParam="Tumor_Sample_Barcode", 
  #  xAxisOrderingParam="orderingVal", 
  #  yAxisValParam="Nmut",
  #  yAxisFillParam="#000000")+ theme(plot.margin = myStandardMargin),
  
  #generate_ggplot_tiles(dfMSIWildtype, xAxisValParam="Tumor_Sample_Barcode", xAxisOrderingParam="orderingVal",
  #                      fillArgParam="cancer_type_adjusted", hideLegend = TRUE, mode = "cancerTypeTiles", 
  #                      textSizeParam=.5, tileColorPalette=cancerTypeColors, orderSpecCols=plottingLevelsCancerType),
  
  heights=c(1,.35,.2,.2,
            1,.35,.2,.2,
            1,.35,.2,.2,
            1,.35,.2,.2,
            1,.35,.2,.2,
            1,.35,.2,.2)#,
            #1,.35,.2)
)

legendSigs <- get_legend(
  make_bar_comparison(
    dfMLH1Mut, 
    xAxisValParam="Tumor_Sample_Barcode",
    yAxisValParam1="mean_MMR",
    yAxisValParam2="otherPredominantSigMagnitude",
    yAxisFillParam="otherPredominantSigName",
    orderingValParam="orderingVal",
    coloringSpecCols=plottingLevels, barColorPalette=barColors,
    title="MLH1 gene mutated cancers",
    mainColor="#267574",
    hideLegend = FALSE)
)

legendSigs2 <- get_legend(generate_ggplot_tiles(dfMLH1Mut, xAxisValParam="Tumor_Sample_Barcode", xAxisOrderingParam="orderingVal",
                                                fillArgParam="cancer_type_adjusted", hideLegend = FALSE, mode = "cancerTypeTiles", 
                                                textSizeParam=.5, tileColorPalette=cancerTypeColors, orderSpecCols=plottingLevelsCancerType))

finalPlot <- plot_grid(finalPlotSansLegend, legendSigs,legendSigs2,
                       align='hv', ncol=3,
                       rel_widths = c(1,.2, .2), scale = 0.9)

ggsave('~/Desktop/testNoah.pdf', plot=finalPlot)



