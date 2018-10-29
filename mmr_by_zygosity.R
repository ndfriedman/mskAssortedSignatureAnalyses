#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)

if(!exists("foo", mode="function")) source("/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/landscape_plot_util.R")

topMargin <- unit(c(0,-.5,-.5,-.5), "lines")
myStandardMargin <- unit(c(.3,-.5,-.5,-.5), "lines")

plottingLevels <- c('1', 'APOBEC',  '3', '4', '7','10', '11', '14', '17', '18', 'other')
barColors = c(
  "#00FFFF", #age  
  "#FF0000", #APOBEC
  "#FF1493", #brca
  "#FFA500", #smoking 
  "#FFF600", #Uv
  "#ADFF2F", #POLE
  "#2A52BE", #TMZ
  "#00DD5D", #sig14
  "#551A8B", #sig17
  "#BF94E4", #sig18
  "#D3D3D3" #OTHER
)

#ALERT SORTING OF MMR IS WRING

############MLH1############

mlh1GermlineMonoallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/mlh1GermlineMonoallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
mlh1GermlineBiallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/mlh1GermlineBiallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
mlh1SomaticMonoallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/mlh1SomaticMonoallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
mlh1SomaticBiallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/mlh1SomaticBiallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

panelDataMLH1GermMono <- plot_panel(mlh1GermlineMonoallelicSigs,
                               #parameters for plotting the signature comparisson bars
                               xAxisValParam_="pid",
                               barParam1_="MMR",
                               barParam2_="otherPredominantSigMagnitude",
                               yAxisFillParam_="otherPredominantSigName",
                               orderingValParam_="orderingVal",
                               coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                               title_="MLH1 Germline Monoallelics",
                               primaryBarColor_="#267574",
                               mutBurdenBarParam_="Nmut",
                               allelicTypeColorPalette_="Spectral",
                               legendMode=FALSE,
                               gradientBarParam_="msiScore",
                               additionalBarParam_="indelSnpRatio",
                               cancerTypeColorPalette_="Set1")

panelDataMLH1GermBiallelic <- plot_panel(mlh1GermlineBiallelicSigs,
                                    #parameters for plotting the signature comparisson bars
                                    xAxisValParam_="pid",
                                    barParam1_="MMR",
                                    barParam2_="otherPredominantSigMagnitude",
                                    yAxisFillParam_="otherPredominantSigName",
                                    orderingValParam_="orderingVal",
                                    coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                    title_="MLH1 Germline Biallelics",
                                    primaryBarColor_="#267574",
                                    mutBurdenBarParam_="Nmut",
                                    allelicTypeColorPalette_="Spectral",
                                    legendMode=FALSE,
                                    gradientBarParam_="msiScore",
                                    additionalBarParam_="indelSnpRatio",
                                    cancerTypeColorPalette_="Set1")

panelDataMLH1SomaticMono <- plot_panel(mlh1SomaticMonoallelicSigs,
                                    #parameters for plotting the signature comparisson bars
                                    xAxisValParam_="pid",
                                    barParam1_="MMR",
                                    barParam2_="otherPredominantSigMagnitude",
                                    yAxisFillParam_="otherPredominantSigName",
                                    orderingValParam_="orderingVal",
                                    coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                    title_="MLH1 Somatic Monoallelics",
                                    primaryBarColor_="#267574",
                                    mutBurdenBarParam_="Nmut",
                                    allelicTypeColorPalette_="Spectral",
                                    legendMode=FALSE,
                                    gradientBarParam_="msiScore",
                                    additionalBarParam_="indelSnpRatio",
                                    cancerTypeColorPalette_="Set1")

panelDataMLH1SomaticBiallelic <- plot_panel(mlh1SomaticBiallelicSigs,
                                    #parameters for plotting the signature comparisson bars
                                    xAxisValParam_="pid",
                                    barParam1_="MMR",
                                    barParam2_="otherPredominantSigMagnitude",
                                    yAxisFillParam_="otherPredominantSigName",
                                    orderingValParam_="orderingVal",
                                    coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                    title_="MLH1 Somatic Biallelics",
                                    primaryBarColor_="#267574",
                                    mutBurdenBarParam_="Nmut",
                                    allelicTypeColorPalette_="Spectral",
                                    legendMode=FALSE,
                                    gradientBarParam_="msiScore",
                                    additionalBarParam_="indelSnpRatio",
                                    cancerTypeColorPalette_="Set1")

finalPlotSansLegend <- ggarrange(
  panelDataMLH1GermMono[[1]]+ theme(plot.margin = topMargin), panelDataMLH1GermMono[[2]]+ theme(plot.margin = myStandardMargin),panelDataMLH1GermMono[[3]]+ theme(plot.margin = myStandardMargin),panelDataMLH1GermMono[[4]],
  panelDataMLH1GermBiallelic[[1]]+ theme(plot.margin = topMargin), panelDataMLH1GermBiallelic[[2]]+ theme(plot.margin = myStandardMargin),panelDataMLH1GermBiallelic[[3]]+ theme(plot.margin = myStandardMargin),panelDataMLH1GermBiallelic[[4]],
  panelDataMLH1SomaticMono[[1]]+ theme(plot.margin = topMargin), panelDataMLH1SomaticMono[[2]]+ theme(plot.margin = myStandardMargin),panelDataMLH1SomaticMono[[3]]+ theme(plot.margin = myStandardMargin),panelDataMLH1SomaticMono[[4]],
  panelDataMLH1SomaticBiallelic[[1]]+ theme(plot.margin = topMargin), panelDataMLH1SomaticBiallelic[[2]]+ theme(plot.margin = myStandardMargin),panelDataMLH1SomaticBiallelic[[3]]+ theme(plot.margin = myStandardMargin),panelDataMLH1SomaticBiallelic[[4]],
  heights=c(1,.4,.1,.6,
            1,.4,.1,.6,
            1,.4,.1,.6,
            1,.4,.1,.6
  )
)

legendSigs <-plot_panel(mlh1SomaticBiallelicSigs,
             #parameters for plotting the signature comparisson bars
             xAxisValParam_="pid",
             barParam1_="MMR",
             barParam2_="otherPredominantSigMagnitude",
             yAxisFillParam_="otherPredominantSigName",
             orderingValParam_="orderingVal",
             coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
             title_="MLH1 Germline Monoallelics",
             primaryBarColor_="#267574",
             mutBurdenBarParam_="Nmut",
             allelicTypeColorPalette_="Spectral",
             legendMode=TRUE,
             cancerTypeColorPalette_="Set1")

finalPlot <- plot_grid(finalPlotSansLegend, legendSigs[[1]],
                       align='hv', ncol=2,
                       rel_widths = c(1,.2), scale = 0.9)

ggsave('~/Desktop/mlh1Plot.pdf', plot=finalPlot)

#############################MSH2####################################
msh2GermlineMonoallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/msh2GermlineMonoallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
msh2GermlineBiallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/msh2GermlineBiallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
msh2SomaticMonoallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/msh2SomaticMonoallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
msh2SomaticBiallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/msh2SomaticBiallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

panelDatamsh2GermMono <- plot_panel(msh2GermlineMonoallelicSigs,
                                    #parameters for plotting the signature comparisson bars
                                    xAxisValParam_="pid",
                                    barParam1_="MMR",
                                    barParam2_="otherPredominantSigMagnitude",
                                    yAxisFillParam_="otherPredominantSigName",
                                    orderingValParam_="orderingVal",
                                    coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                    title_="msh2 Germline Monoallelics",
                                    primaryBarColor_="#267574",
                                    mutBurdenBarParam_="Nmut",
                                    allelicTypeColorPalette_="Spectral",
                                    legendMode=FALSE,
                                    gradientBarParam_="msiScore",
                                    additionalBarParam_="indelSnpRatio",
                                    cancerTypeColorPalette_="Set1")

panelDatamsh2GermBiallelic <- plot_panel(msh2GermlineBiallelicSigs,
                                         #parameters for plotting the signature comparisson bars
                                         xAxisValParam_="pid",
                                         barParam1_="MMR",
                                         barParam2_="otherPredominantSigMagnitude",
                                         yAxisFillParam_="otherPredominantSigName",
                                         orderingValParam_="orderingVal",
                                         coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                         title_="msh2 Germline Biallelics",
                                         primaryBarColor_="#267574",
                                         mutBurdenBarParam_="Nmut",
                                         allelicTypeColorPalette_="Spectral",
                                         legendMode=FALSE,
                                         gradientBarParam_="msiScore",
                                         additionalBarParam_="indelSnpRatio",
                                         cancerTypeColorPalette_="Set1")

panelDatamsh2SomaticMono <- plot_panel(msh2SomaticMonoallelicSigs,
                                       #parameters for plotting the signature comparisson bars
                                       xAxisValParam_="pid",
                                       barParam1_="MMR",
                                       barParam2_="otherPredominantSigMagnitude",
                                       yAxisFillParam_="otherPredominantSigName",
                                       orderingValParam_="orderingVal",
                                       coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                       title_="msh2 Somatic Monoallelics",
                                       primaryBarColor_="#267574",
                                       mutBurdenBarParam_="Nmut",
                                       allelicTypeColorPalette_="Spectral",
                                       legendMode=FALSE,
                                       gradientBarParam_="msiScore",
                                       additionalBarParam_="indelSnpRatio",
                                       cancerTypeColorPalette_="Set1")

panelDatamsh2SomaticBiallelic <- plot_panel(msh2SomaticBiallelicSigs,
                                            #parameters for plotting the signature comparisson bars
                                            xAxisValParam_="pid",
                                            barParam1_="MMR",
                                            barParam2_="otherPredominantSigMagnitude",
                                            yAxisFillParam_="otherPredominantSigName",
                                            orderingValParam_="orderingVal",
                                            coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                            title_="msh2 Somatic Biallelics",
                                            primaryBarColor_="#267574",
                                            mutBurdenBarParam_="Nmut",
                                            allelicTypeColorPalette_="Spectral",
                                            legendMode=FALSE,
                                            gradientBarParam_="msiScore",
                                            additionalBarParam_="indelSnpRatio",
                                            cancerTypeColorPalette_="Set1")

finalPlotSansLegend <- ggarrange(
  
  panelDatamsh2GermMono[[1]]+ theme(plot.margin = topMargin), panelDatamsh2GermMono[[2]]+ theme(plot.margin = myStandardMargin),panelDatamsh2GermMono[[3]]+ theme(plot.margin = myStandardMargin),panelDatamsh2GermMono[[4]],
  panelDatamsh2GermBiallelic[[1]]+ theme(plot.margin = topMargin), panelDatamsh2GermBiallelic[[2]]+ theme(plot.margin = myStandardMargin),panelDatamsh2GermBiallelic[[3]]+ theme(plot.margin = myStandardMargin),panelDatamsh2GermBiallelic[[4]],
  panelDatamsh2SomaticMono[[1]]+ theme(plot.margin = topMargin), panelDatamsh2SomaticMono[[2]]+ theme(plot.margin = myStandardMargin),panelDatamsh2SomaticMono[[3]]+ theme(plot.margin = myStandardMargin),panelDatamsh2SomaticMono[[4]],
  panelDatamsh2SomaticBiallelic[[1]]+ theme(plot.margin = topMargin), panelDatamsh2SomaticBiallelic[[2]]+ theme(plot.margin = myStandardMargin),panelDatamsh2SomaticBiallelic[[3]]+ theme(plot.margin = myStandardMargin),panelDatamsh2SomaticBiallelic[[4]],
  heights=c(1,.4,.1,.6,
            1,.4,.1,.6,
            1,.4,.1,.6,
            1,.4,.1,.6
  )
)

legendSigs <-plot_panel(msh2SomaticBiallelicSigs,
                        #parameters for plotting the signature comparisson bars
                        xAxisValParam_="pid",
                        barParam1_="MMR",
                        barParam2_="otherPredominantSigMagnitude",
                        yAxisFillParam_="otherPredominantSigName",
                        orderingValParam_="orderingVal",
                        coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                        title_="msh2 Germline Monoallelics",
                        primaryBarColor_="#267574",
                        mutBurdenBarParam_="Nmut",
                        allelicTypeColorPalette_="Spectral",
                        legendMode=TRUE,
                        cancerTypeColorPalette_="Set1")

finalPlot <- plot_grid(finalPlotSansLegend, legendSigs[[1]],
                       align='hv', ncol=2,
                       rel_widths = c(1,.2), scale = 0.9)

ggsave('~/Desktop/msh2Plot.pdf', plot=finalPlot)

#############################MSH6##################################
msh6GermlineMonoallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/msh6GermlineMonoallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
msh6GermlineBiallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/msh6GermlineBiallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
msh6SomaticMonoallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/msh6SomaticMonoallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
msh6SomaticBiallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/msh6SomaticBiallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

panelDatamsh6GermMono <- plot_panel(msh6GermlineMonoallelicSigs,
                                    #parameters for plotting the signature comparisson bars
                                    xAxisValParam_="pid",
                                    barParam1_="MMR",
                                    barParam2_="otherPredominantSigMagnitude",
                                    yAxisFillParam_="otherPredominantSigName",
                                    orderingValParam_="orderingVal",
                                    coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                    title_="msh6 Germline Monoallelics",
                                    primaryBarColor_="#267574",
                                    mutBurdenBarParam_="Nmut",
                                    allelicTypeColorPalette_="Spectral",
                                    legendMode=FALSE,
                                    gradientBarParam_="msiScore",
                                    additionalBarParam_="indelSnpRatio",
                                    cancerTypeColorPalette_="Set1")

panelDatamsh6GermBiallelic <- plot_panel(msh6GermlineBiallelicSigs,
                                         #parameters for plotting the signature comparisson bars
                                         xAxisValParam_="pid",
                                         barParam1_="MMR",
                                         barParam2_="otherPredominantSigMagnitude",
                                         yAxisFillParam_="otherPredominantSigName",
                                         orderingValParam_="orderingVal",
                                         coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                         title_="msh6 Germline Biallelics",
                                         primaryBarColor_="#267574",
                                         mutBurdenBarParam_="Nmut",
                                         allelicTypeColorPalette_="Spectral",
                                         legendMode=FALSE,
                                         gradientBarParam_="msiScore",
                                         additionalBarParam_="indelSnpRatio",
                                         cancerTypeColorPalette_="Set1")

panelDatamsh6SomaticMono <- plot_panel(msh6SomaticMonoallelicSigs,
                                       #parameters for plotting the signature comparisson bars
                                       xAxisValParam_="pid",
                                       barParam1_="MMR",
                                       barParam2_="otherPredominantSigMagnitude",
                                       yAxisFillParam_="otherPredominantSigName",
                                       orderingValParam_="orderingVal",
                                       coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                       title_="msh6 Somatic Monoallelics",
                                       primaryBarColor_="#267574",
                                       mutBurdenBarParam_="Nmut",
                                       allelicTypeColorPalette_="Spectral",
                                       legendMode=FALSE,
                                       gradientBarParam_="msiScore",
                                       additionalBarParam_="indelSnpRatio",
                                       cancerTypeColorPalette_="Set1")

panelDatamsh6SomaticBiallelic <- plot_panel(msh6SomaticBiallelicSigs,
                                            #parameters for plotting the signature comparisson bars
                                            xAxisValParam_="pid",
                                            barParam1_="MMR",
                                            barParam2_="otherPredominantSigMagnitude",
                                            yAxisFillParam_="otherPredominantSigName",
                                            orderingValParam_="orderingVal",
                                            coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                            title_="msh6 Somatic Biallelics",
                                            primaryBarColor_="#267574",
                                            mutBurdenBarParam_="Nmut",
                                            allelicTypeColorPalette_="Spectral",
                                            legendMode=FALSE,
                                            gradientBarParam_="msiScore",
                                            additionalBarParam_="indelSnpRatio",
                                            cancerTypeColorPalette_="Set1")

finalPlotSansLegend <- ggarrange(
  
  panelDatamsh6GermMono[[1]]+ theme(plot.margin = topMargin), panelDatamsh6GermMono[[2]]+ theme(plot.margin = myStandardMargin),panelDatamsh6GermMono[[3]]+ theme(plot.margin = myStandardMargin),panelDatamsh6GermMono[[4]],
  panelDatamsh6GermBiallelic[[1]]+ theme(plot.margin = topMargin), panelDatamsh6GermBiallelic[[2]]+ theme(plot.margin = myStandardMargin),panelDatamsh6GermBiallelic[[3]]+ theme(plot.margin = myStandardMargin),panelDatamsh6GermBiallelic[[4]],
  panelDatamsh6SomaticMono[[1]]+ theme(plot.margin = topMargin), panelDatamsh6SomaticMono[[2]]+ theme(plot.margin = myStandardMargin),panelDatamsh6SomaticMono[[3]]+ theme(plot.margin = myStandardMargin),panelDatamsh6SomaticMono[[4]],
  panelDatamsh6SomaticBiallelic[[1]]+ theme(plot.margin = topMargin), panelDatamsh6SomaticBiallelic[[2]]+ theme(plot.margin = myStandardMargin),panelDatamsh6SomaticBiallelic[[3]]+ theme(plot.margin = myStandardMargin),panelDatamsh6SomaticBiallelic[[4]],
  heights=c(1,.4,.1,.6,
            1,.4,.1,.6,
            1,.4,.1,.6,
            1,.4,.1,.6
  )
)

legendSigs <-plot_panel(msh6SomaticBiallelicSigs,
                        #parameters for plotting the signature comparisson bars
                        xAxisValParam_="pid",
                        barParam1_="MMR",
                        barParam2_="otherPredominantSigMagnitude",
                        yAxisFillParam_="otherPredominantSigName",
                        orderingValParam_="orderingVal",
                        coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                        title_="msh6 Germline Monoallelics",
                        primaryBarColor_="#267574",
                        mutBurdenBarParam_="Nmut",
                        allelicTypeColorPalette_="Spectral",
                        legendMode=TRUE,
                        cancerTypeColorPalette_="Set1")

finalPlot <- plot_grid(finalPlotSansLegend, legendSigs[[1]],
                       align='hv', ncol=2,
                       rel_widths = c(1,.2), scale = 0.9)

ggsave('~/Desktop/msh6Plot.pdf', plot=finalPlot)

#############################PMS2#################################
pms2GermlineMonoallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/pms2GermlineMonoallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
pms2GermlineBiallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/pms2GermlineBiallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
pms2SomaticMonoallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/pms2SomaticMonoallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util
pms2SomaticBiallelicSigs <- read.table('/Users/friedman/Desktop/landscapePlots/pms2SomaticBiallelicSigs.tsv',sep = '\t', header=TRUE) #load a dataframe that has already been formatted properly by a python util

panelDatapms2GermMono <- plot_panel(pms2GermlineMonoallelicSigs,
                                    #parameters for plotting the signature comparisson bars
                                    xAxisValParam_="pid",
                                    barParam1_="MMR",
                                    barParam2_="otherPredominantSigMagnitude",
                                    yAxisFillParam_="otherPredominantSigName",
                                    orderingValParam_="orderingVal",
                                    coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                    title_="pms2 Germline Monoallelics",
                                    primaryBarColor_="#267574",
                                    mutBurdenBarParam_="Nmut",
                                    allelicTypeColorPalette_="Spectral",
                                    legendMode=FALSE,
                                    gradientBarParam_="msiScore",
                                    additionalBarParam_="indelSnpRatio",
                                    cancerTypeColorPalette_="Set1")

panelDatapms2GermBiallelic <- plot_panel(pms2GermlineBiallelicSigs,
                                         #parameters for plotting the signature comparisson bars
                                         xAxisValParam_="pid",
                                         barParam1_="MMR",
                                         barParam2_="otherPredominantSigMagnitude",
                                         yAxisFillParam_="otherPredominantSigName",
                                         orderingValParam_="orderingVal",
                                         coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                         title_="pms2 Germline Biallelics",
                                         primaryBarColor_="#267574",
                                         mutBurdenBarParam_="Nmut",
                                         allelicTypeColorPalette_="Spectral",
                                         legendMode=FALSE,
                                         gradientBarParam_="msiScore",
                                         additionalBarParam_="indelSnpRatio",
                                         cancerTypeColorPalette_="Set1")

panelDatapms2SomaticMono <- plot_panel(pms2SomaticMonoallelicSigs,
                                       #parameters for plotting the signature comparisson bars
                                       xAxisValParam_="pid",
                                       barParam1_="MMR",
                                       barParam2_="otherPredominantSigMagnitude",
                                       yAxisFillParam_="otherPredominantSigName",
                                       orderingValParam_="orderingVal",
                                       coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                       title_="pms2 Somatic Monoallelics",
                                       primaryBarColor_="#267574",
                                       mutBurdenBarParam_="Nmut",
                                       allelicTypeColorPalette_="Spectral",
                                       legendMode=FALSE,
                                       gradientBarParam_="msiScore",
                                       additionalBarParam_="indelSnpRatio",
                                       cancerTypeColorPalette_="Set1")

panelDatapms2SomaticBiallelic <- plot_panel(pms2SomaticBiallelicSigs,
                                            #parameters for plotting the signature comparisson bars
                                            xAxisValParam_="pid",
                                            barParam1_="MMR",
                                            barParam2_="otherPredominantSigMagnitude",
                                            yAxisFillParam_="otherPredominantSigName",
                                            orderingValParam_="orderingVal",
                                            coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                                            title_="pms2 Somatic Biallelics",
                                            primaryBarColor_="#267574",
                                            mutBurdenBarParam_="Nmut",
                                            allelicTypeColorPalette_="Spectral",
                                            legendMode=FALSE,
                                            gradientBarParam_="msiScore",
                                            additionalBarParam_="indelSnpRatio",
                                            cancerTypeColorPalette_="Set1")

finalPlotSansLegend <- ggarrange(
  
  panelDatapms2GermMono[[1]]+ theme(plot.margin = topMargin), panelDatapms2GermMono[[2]]+ theme(plot.margin = myStandardMargin),panelDatapms2GermMono[[3]]+ theme(plot.margin = myStandardMargin),panelDatapms2GermMono[[4]],
  panelDatapms2GermBiallelic[[1]]+ theme(plot.margin = topMargin), panelDatapms2GermBiallelic[[2]]+ theme(plot.margin = myStandardMargin),panelDatapms2GermBiallelic[[3]]+ theme(plot.margin = myStandardMargin),panelDatapms2GermBiallelic[[4]],
  panelDatapms2SomaticMono[[1]]+ theme(plot.margin = topMargin), panelDatapms2SomaticMono[[2]]+ theme(plot.margin = myStandardMargin),panelDatapms2SomaticMono[[3]]+ theme(plot.margin = myStandardMargin),panelDatapms2SomaticMono[[4]],
  panelDatapms2SomaticBiallelic[[1]]+ theme(plot.margin = topMargin), panelDatapms2SomaticBiallelic[[2]]+ theme(plot.margin = myStandardMargin),panelDatapms2SomaticBiallelic[[3]]+ theme(plot.margin = myStandardMargin),panelDatapms2SomaticBiallelic[[4]],
  heights=c(1,.4,.1,.6,
            1,.4,.1,.6,
            1,.4,.1,.6,
            1,.4,.1,.6
  )
)

legendSigs <-plot_panel(pms2SomaticBiallelicSigs,
                        #parameters for plotting the signature comparisson bars
                        xAxisValParam_="pid",
                        barParam1_="MMR",
                        barParam2_="otherPredominantSigMagnitude",
                        yAxisFillParam_="otherPredominantSigName",
                        orderingValParam_="orderingVal",
                        coloringSpecCols_=plottingLevels, barColorPalette_=barColors,
                        title_="pms2 Germline Monoallelics",
                        primaryBarColor_="#267574",
                        mutBurdenBarParam_="Nmut",
                        allelicTypeColorPalette_="Spectral",
                        legendMode=TRUE,
                        cancerTypeColorPalette_="Set1")

finalPlot <- plot_grid(finalPlotSansLegend, legendSigs[[1]],
                       align='hv', ncol=2,
                       rel_widths = c(1,.2), scale = 0.9)

ggsave('~/Desktop/pms2Plot.pdf', plot=finalPlot)


