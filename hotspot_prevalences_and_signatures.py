#written by Noah Friedman (a template for scripts to be excuted in the spyder environment
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

sys.path.append(pathPrefix + '/ifs/work/taylorlab/friedman/myUtils')
import analysis_utils 
import mutationSigUtils 
import maf_analysis_utils


def calculate_hotspot_data(sigsData, hotspots,
                           smokingSpectraD, #dictionaries mapping quadnucs to spectra probs 
                           observedSpectraD,
                           syntheticSpectraD,
                           minMutThesholdForSignatures = 10, #we do not consider signatures calculated in cases with fewer than this number of mutations
                           comparissonHotspotName = 'G12C' #the hotspot we compare everything to
                           ):
    listOfDicts = [] #what we will use to build up a dataframe to return
    #calculate relevant information for mutations that we are going to compare everything to
    comparissonQuadNuc = sigsData[sigsData['Amino_Acid_Change'] == comparissonHotspotName].iloc[0]['quadNuc']
    #get spectra probabilities from different distributions of spectra
    comparissonSmokingSpectraProb = smokingSpectraD[comparissonQuadNuc]
    comparissonObservedSpectraProb = observedSpectraD[comparissonQuadNuc]
    comparissonSynthethicSpectraProb = syntheticSpectraD[comparissonQuadNuc]
    comparissonHotspotMuts = sigsData[sigsData['Amino_Acid_Change'] == comparissonHotspotName]
    sigsDataEnoughMuts = sigsData[sigsData['Nmut'] > minMutThesholdForSignatures]
    cohortSmokingSigMedian = np.nanmedian(sigsDataEnoughMuts['smoking'])
    for hotspot in hotspots:
        hotspotQuadNuc = sigsData[sigsData['Amino_Acid_Change'] == hotspot].iloc[0]['quadNuc'] #get the current quad nuc
        
        #ratio calculations
        ratioToComparissonSmoking = smokingSpectraDict[hotspotQuadNuc]/comparissonSmokingSpectraProb
        ratioToComparissonObservedSpectra = observedSpectraDict[hotspotQuadNuc]/comparissonObservedSpectraProb
        ratioToComparissonSyntheticSpectra = syntheticSpectraDict[hotspotQuadNuc]/comparissonSynthethicSpectraProb
        currentHotspotMuts = sigsData[sigsData['Amino_Acid_Change'] == hotspot]
        ratioToComparissonObserved = 1.0*currentHotspotMuts.shape[0]/comparissonHotspotMuts.shape[0]
        nMutationsSmoking = np.nanmedian(currentHotspotMuts['nMutAttributedToSmoking'])
        
        #calculations about mutation signatures
        currentHotspotEnoughMuts = currentHotspotMuts[currentHotspotMuts['Nmut'] > minMutThesholdForSignatures]
        fracNotEnoughMuts = 1 - 1.0*currentHotspotEnoughMuts.shape[0]/currentHotspotMuts.shape[0]
        smokingSigAbs = np.nanmedian(currentHotspotEnoughMuts['smoking'])
        smokingSigToCohortRatio = 1.0*smokingSigAbs/cohortSmokingSigMedian
        nCases = currentHotspotMuts.shape[0]
        listOfDicts.append({
                    'hotspot': hotspot,
                    'hotspotQuadNuc': hotspotQuadNuc,
                    'ratioToCompSmoking': ratioToComparissonSmoking,
                    'ratioToCompObservedSpectra': ratioToComparissonObservedSpectra,
                    'ratioToCompSyntheticSpectra': ratioToComparissonSyntheticSpectra,
                    'ratioToCompObserved': ratioToComparissonObserved,
                    'smokingSigAbs': smokingSigAbs,
                    'nMutationsSmoking': nMutationsSmoking,
                    'smokingSigToCohortRatio': smokingSigToCohortRatio,
                    'insufficientMutFrac': fracNotEnoughMuts,
                    'nCases': nCases
                })
    df = pd.DataFrame(listOfDicts)
    return df

def get_observed_mut_spectra_kras_prob(row):
    if row['Nmut'] < 10:
        return 
    else:
        return 

hotspotMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hotspotReducedAnalysis10-19.tsv')
impactOnlyMuts = hotspotMaf[(hotspotMaf['Tumor_Sample_Barcode'].str.contains('P-')) & (hotspotMaf['Tumor_Sample_Barcode'].str.len() == 9)]
impactOnlyMuts['pid'] = impactOnlyMuts['Tumor_Sample_Barcode']
#dfSaved['Amino_Acid_Change']
impactKrasMutations = impactOnlyMuts[impactOnlyMuts['Hugo_Symbol'] == 'KRAS']
impactKrasMutations['quadNuc'] = impactKrasMutations.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

#Load impact sigs and annotate them#######################
impactSigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/pensona/dmp_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['smoking'] = impactSigs.apply(lambda row: row['mean_4'] + row['mean_18'] + row['mean_24'] + row['mean_29'], axis=1) #smoking sig is sometimes misclassified as these 
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

#Add information about hotspots to the Sigs DF
krasHotspots = ['G12D', 'G12V', 'G12C', 'G12R', 'G12A', 'G12S', 'G13D', 'G13C']
impactSigs = maf_analysis_utils.add_hotspot_maf_info_to_df(impactSigs, impactKrasMutations, krasHotspots)
impactSigsKras = impactSigs[impactSigs['Amino_Acid_Change'].notnull()]
impactSigsKras = impactSigsKras[impactSigsKras['quadNuc'].notnull()]
spectrumDicts = mutationSigUtils.convert_spectrum_file_to_dict_of_dicts(pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')

expandedDf['CCAA_prob']
expandedDf['smoking']

impactSigsKrasMerged = mutationSigUtils.merge_signature_columns(impactSigsKras, smokingMerge=True)
impactSigsKrasMerged['nMutAttributedToMMR'] = impactSigsKrasMerged.apply(lambda row: row['mean_MMR'] * row['Nmut'], axis=1)
impactSigsKrasMerged['nMutAttributedToBRCA'] = impactSigsKrasMerged.apply(lambda row: row['mean_3'] * row['Nmut'], axis=1)
impactSigsKrasLung['nMutAttributedToSmoking'] = impactSigsKrasLung.apply(lambda row: row['Nmut'] * row['smoking'], axis=1)


#mutationSigUtils.expand_df_to_include_spectrum_probs(df, spectrumDict)

#Assign some information about spectrums#################################3
spectrumFilePath='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'
spectrumDf = pd.read_table(pathPrefix + spectrumFilePath)
smokingSpectraDict = spectrumDf[spectrumDf.index == 'Signature.4'].to_dict('records')[0]

#reduce down to cancer type of interest (lung) only
impactSigsKrasLung = impactSigsKras[impactSigsKras['cancer_type'] == 'Non-Small Cell Lung Cancer']

synthetic_means_lung = [0.17, 0.06, 0.15, 0.08, 0.07, 0.04, 0.02, 0.05, 0.04, 0.01, 0.02, 0.01, 0.04, 0.01, 0.01, 0.03, 0.01, 0.02, 0.02, 0.01, 0.01, 0.0, 0.01, 0.05, 0.02, 0.01, 0.0, 0.01, 0.03, 0.03] #means in low mut burden Lung cancers
expandedDf = mutationSigUtils.expand_df_to_include_spectrum_probs(impactSigsKras, spectrumDicts)


impactSigsKrasLung['smokingSpectraProb'] = impactSigsKrasLung['quadNuc'].apply(lambda x: smokingSpectraDict[x] if x in smokingSpectraDict else None)
impactSigsKrasLung['nMutAttributedToSmoking'] = impactSigsKrasLung.apply(lambda row: row['Nmut'] * row['smoking'], axis=1)
impactSigsKrasLung['observedMutSpectraProbKras'] = impactSigsKrasLung.apply(lambda row:

impactSigsKrasLung.to_csv('~/Desktop/krasMutsTemp.tsv', sep='\t', index=False)


#impactSigsKrasLungSufficientMuts = impactSigsKrasLung[impactSigsKrasLung['Nmut'] > 10] #set 

sigNames = ['mean_' + str(i) for i in range(1,31)]
meanSigs = [np.nanmean(impactSigsKrasLungSufficientMuts[name]) for name in sigNames]
observedSpectraDict = mutationSigUtils.get_spectrum_mutation_frac_for_cohort(meanSigs, spectrumFilePath= pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
#make the synthetic sigs
meanSigsClipped = [i if i > .025 else 0 for i in meanSigs] #only take signatures present at greater thant .025
meanSigsClipped[0] = meanSigsClipped[0]*3 #up the weighting for the Aging sigs (1 and 5) because cases with fewer mutations have enriched aging/sig5 amounts
meanSigsClipped[4] = meanSigsClipped[4]*3 
normedClippedSigs = [float(i)/sum(meanSigsClipped) for i in meanSigsClipped]
synthetic_means = [0.17, 0.06, 0.15, 0.08, 0.07, 0.04, 0.02, 0.05, 0.04, 0.01, 0.02, 0.01, 0.04, 0.01, 0.01, 0.03, 0.01, 0.02, 0.02, 0.01, 0.01, 0.0, 0.01, 0.05, 0.02, 0.01, 0.0, 0.01, 0.03, 0.03] #means in low mut burden cancers
normedSyntheticSigs = [float(i)/sum(synthetic_means) for i in synthetic_means]
syntheticSpectraDict = mutationSigUtils.get_spectrum_mutation_frac_for_cohort(normedClippedSigs, spectrumFilePath= pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')

returnDf = calculate_hotspot_data(impactSigsKrasLung, krasHotspots,
                           smokingSpectraDict, #dictionaries mapping quadnucs to spectra probs 
                           observedSpectraDict,
                           syntheticSpectraDict,
                           minMutThesholdForSignatures = 10, #we do not consider signatures calculated in cases with fewer than this number of mutations
                           comparissonHotspotName = 'G12C' #the hotspot we compare everything to
                           )

impactSigsKrasLungSufficientMuts.columns.values

impactSigsKrasLungSufficientMuts.to_csv('~/Desktop/krasMutsTemp.tsv', sep='\t', index=False)

meltedDf = pd.melt(returnDf, id_vars=['hotspot'], value_vars=['ratioToCompSmoking', 'ratioToCompSyntheticSpectra', 'ratioToCompObservedSpectra'])
meltedDf.to_csv('~/Desktop/hotspotPrevalenceAnalysis/ratioData_kras_lung.tsv', sep='\t', index=False)

returnDf.to_csv('~/Desktop/noahTestHotspotAnal.tsv', sep='\t', index=False)

#print returnDf.transpose()
######################################################################3






















impactSigsKrasPancreas = impactSigsKras[impactSigsKras['cancer_type'] == 'Pancreatic Cancer']
impactSigsKrasPancreasg12R = impactSigsKrasPancreas[impactSigsKrasPancreas['Amino_Acid_Change'] == 'G12R']
impactSigsKrasPancreasNOTg12R = impactSigsKrasPancreas[impactSigsKrasPancreas['Amino_Acid_Change'] != 'G12R']
krasG12Rids = set(impactSigsKrasPancreasg12R['Tumor_Sample_Barcode'])
krasNotG12Rids = set(impactSigsKrasPancreasNOTg12R['Tumor_Sample_Barcode'])


np.nanmedian(impactSigsKrasPancreasg12R['mean_3']), np.nanmedian(impactSigsKrasPancreasNOTg12R['mean_3'])

F = open('/Users/friedman/Desktop/file.txt','w')
for entry in set(impactSigsKrasPancreasNOTg12R['Tumor_Sample_Barcode']):
    F.write(entry + '\n')
    
brcaDataDf = pd.read_table('/Users/friedman/Desktop/gray-93017-exome-annotation.tsv')
pancreatics = brcaDataDf[brcaDataDf['cancer_type'] == 'Pancreatic Cancer']
 = pancreatics[pancreatics['gene'].notnull()]

pancreatics[pancreatics['dmp_sample'].isin(krasNotG12Rids)]

np.nanmean(pancreatics[pancreatics['dmp_sample'].isin(krasG12Rids)]['brca_signature'])


np.nanmean(pancreatics['brca_signature'])

print max(impactSigsKrasLung['smokingSpectraProb'])
print set(impactSigsKrasLung[impactSigsKrasLung['Amino_Acid_Change'] == 'G13C']['smokingSpectraProb'])







sigNames = ['mean_' + str(i) for i in range(1,31)]

#isolate out the lung cancers with sufficient mutations
lungCancerKrasSigs = impactSigsKrasMuts[impactSigsKrasMuts['cancer_type'] == 'Non-Small Cell Lung Cancer']
lungCancerKrasSufficientMutSigs = lungCancerKrasSigs[lungCancerKrasSigs['Nmut'] >= 10]
lungCancerKrasSufficientMutIds = set(lungCancerKrasSufficientMutSigs['pid'])
meanSigs = [np.nanmean(lungCancerKrasSufficientMutSigs[name]) for name in sigNames]






lungCancerKrasSigs.columns.values


mutSigFracs = mutationSigUtils.get_spectrum_mutation_frac_for_cohort(meanSigs, spectrumFilePath= pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')

krasHotspotMotifDict = {'ACTC':'G12D', 'ACAC':'G12V', 'CCAA':'G12C', 'CCGA':'G12R', 'ACGC':'G12A', 'CCTA':'G12S',
                        'GCTC':'G13D', 'CCAA':'G13C'} #kras mutations in g12/g13
krasg12_13Hotspots = set(krasHotspotMotifDict.keys())
trinucsHotspots = mutSigFracs.to_dict()
trinucsHotspots = {k: v for k, v in trinucsHotspots.items() if k in krasg12_13Hotspots}

#normalize the dictionary
factor=1.0/sum(trinucsHotspots.itervalues())
for k in trinucsHotspots:
  trinucsHotspots[k] = trinucsHotspots[k]*factor

krasg12_13muts = impactKrasMutations[impactKrasMutations['Amino_Acid_Change'].isin(set(krasHotspotMotifDict.values()))]
cntrObj = Counter(krasg12_13muts['Amino_Acid_Change'])

#normalize the counter 
total = sum(cntrObj.values(), 0.0)
for key in cntrObj:
    cntrObj[key] /= total


