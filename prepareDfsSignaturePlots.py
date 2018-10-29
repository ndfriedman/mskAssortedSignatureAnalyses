#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'
    
sys.path.append(pathPrefix + '/ifs/work/taylorlab/friedman/myUtils')
import analysis_utils 
import mutationSigUtils   
    
def get_msi_score_dict():
    msiList = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/IM_MSIscores_092818_.txt')
    msiList['tid'] = msiList['DMP'].apply(lambda x: x[:13])
    msiScoreDict = dict(zip(list(msiList['tid']), list(msiList['%'])))
    return msiScoreDict  

def find_second_most_common(row, primarySig, returnMode, 
                            sigNamesToSpecify = set(['mean_1', 'mean_3', 'mean_4', 'mean_7', 'mean_10', 'mean_11','mean_14', 'mean_17', 'mean_MMR', 'mean_APOBEC']) #a set of signatures we actually mark on the chart
                            ):
    colNames = row.to_dict().keys()
    signatureColumns = [i for i in list(row.keys()) if 'mean' in i]
    rowSigsOnly = row[signatureColumns]
    rowAsDict = rowSigsOnly.to_dict()
    items = rowAsDict.items()
    sortedItems = sorted(items, key=lambda x: x[1], reverse=True)
    if sortedItems[0][0] == primarySig:
        if returnMode == 'name':
            sigName = sortedItems[1][0]
            if sigName in sigNamesToSpecify:
                return sigName
            else:
                return 'other'
        else:
            return sortedItems[1][1]
    else:
        if returnMode == 'name':
            sigName = sortedItems[0][0]
            if sigName in sigNamesToSpecify:
                return sigName
            else:
                return 'other' 
        else:
            return sortedItems[0][1]
        
#define an ordering for the ggplot plot based on the first signature to hit a clipping thresh
#add functionality for low mut burden
def ordering_function_clip_mode(row, clipThresh = .15, sigsToOrderBy = ['mean_MMR', 'mean_1', 'mean_APOBEC', 'mean_3', 'mean_4', 'mean_7', 'mean_10']):
    orderingNum = len(sigsToOrderBy)
    for i in range(len(sigsToOrderBy)):
        curSigToConsider = sigsToOrderBy[i]    
        if row[curSigToConsider] > clipThresh and (curSigToConsider == row['otherPredominantSigName'] or curSigToConsider == 'mean_MMR'):
            return orderingNum + row[curSigToConsider]
        else:
            orderingNum -=1
    return 0
    
def ordering_function_dom_sig_mode(row, sigsToOrderBy = ['mean_APOBEC', 'mean_3', 'mean_4', 'mean_7', 'mean_10', 'mean_14', 'mean_11', 'mean_17']):
    orderingNum = len(sigsToOrderBy) + 2
    if row['Nmut'] < 10: return -1 #cases with less than 10 mutations go at the far side
    
    if row['mean_MMR'] > row['otherPredominantSigMagnitude']:
        return orderingNum + row['mean_MMR']
    orderingNum -= 1
    if row['otherPredominantSigName'] == 'mean_1': #sort age dominant cases by mean MMR
        return orderingNum + row['mean_MMR']
    orderingNum -= 1
    for i in range(len(sigsToOrderBy)):
        curSigToConsider = sigsToOrderBy[i]    
        if curSigToConsider == row['otherPredominantSigName']:
            return orderingNum + row[curSigToConsider]
        else:
            orderingNum -=1
    return 0

exomeDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/pensona/dmp_sigs/mixed_impact_wes_2017_data_mutations_extended.sigs.tab.txt')
df = pd.read_table(pathPrefix + '/ifs/work/taylorlab/pensona/dmp_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')

df['pid'] = df['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
df['tid'] = df['Tumor_Sample_Barcode'].apply(lambda x: x[:13])
df['cancer_type'] = df['pid'].apply(lambda x: cDict[x] if x in cDict else None)
df = mutationSigUtils.merge_signature_columns(df, mode='Stratton')
df['otherPredominantSigName'] = df.apply(lambda row: find_second_most_common(row, 'mean_APOBEC', 'name'), axis=1)
df['otherPredominantSigMagnitude'] = df.apply(lambda row: find_second_most_common(row, 'mean_APOBEC', 'magnitude'), axis=1)

exomeDf['pid'] = exomeDf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
exomeDf['tid'] = exomeDf['Tumor_Sample_Barcode'].apply(lambda x: x[:13])
exomeDf['cancer_type'] = exomeDf['pid'].apply(lambda x: cDict[x] if x in cDict else None)
exomeDf = mutationSigUtils.merge_signature_columns(exomeDf, mode='Stratton')
exomeDf['otherPredominantSigName'] = exomeDf.apply(lambda row: find_second_most_common(row, 'mean_APOBEC', 'name'), axis=1)
exomeDf['otherPredominantSigMagnitude'] = exomeDf.apply(lambda row: find_second_most_common(row, 'mean_APOBEC', 'magnitude'), axis=1)

clippingThresh = 750 
exomeDf['NmutAdj'] = exomeDf['Nmut'].apply(lambda x: clippingThresh if x >= clippingThresh else x)
exomeDf['nMutClipped'] = exomeDf['NmutAdj'].apply(lambda x: False if x < clippingThresh else True)

impactInExomeDf = df[df['tid'].isin(set(exomeDf['tid']))]
apobecImpactDict = dict(zip(list(impactInExomeDf['tid']), list(impactInExomeDf['confidence_APOBEC'])))
exomeDf['apobecImpactMag'] = exomeDf['tid'].apply(lambda x: -.00001 if x not in apobecImpactDict else apobecImpactDict[x])


exomeTids = set(exomeDf['tid'])

maf = pd.read_table(pathPrefix + '/ifs/res/taylorlab/ang46/ext/dmp/mixedpact_ang46/data_mutations_unfiltered_reviewed_oncokb.txt')
maf['tid'] = maf['Tumor_Sample_Barcode'].apply(lambda x: x[:13])

mafExomeCasesOnly = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedMafImpactVarsFromExome_trinuc') #LOAD the version I generated
mafExomeCasesOnly['Ref_Tri']
#mafExomeCasesOnly = maf[maf['tid'].isin(exomeTids)]
mafTids = set(mafExomeCasesOnly['tid'])
mafExomeCasesOnly['quadNuc'] = mafExomeCasesOnly.apply(lambda row: 
    mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
    if row['Variant_Type'] == 'SNP'
    else None
    , axis=1)

for tid in mafTids:
    if tid in apobecFracDict and tid in agingFracDict:
        if apobecFracDict[tid] != None and agingFracDict[tid] != None:
            print apobecFracDict[tid] + agingFracDict[tid]


exomeDf['apobecSignificantMutPortion'] = exomeDf['tid'].apply(lambda x: apobecFracDict[x] if x in apobecFracDict else None)
exomeDf['agingSignificantMutPortion'] = exomeDf['tid'].apply(lambda x: agingFracDict[x] if x in agingFracDict else None)
exomeDf['smokingSignificantMutPortion'] = exomeDf['tid'].apply(lambda x: smokingFracDict[x] if x in smokingFracDict else None)


apobecRelatedCancers = set(['Breast Cancer', 'Bladder Cancer', '','Head and Neck Cancer', 'Cervical Cancer', 'Uterine Sarcoma', 'Esophagogastric Cancer', 'Small Cell Lung Cancer', 'Non-Small Cell Lung Cancer'])
otherApobecCancerTypes = set(['Head and Neck Cancer', 'Cervical Cancer', 'Uterine Sarcoma', 'Esophagogastric Cancer'])
bladders = exomeDf[exomeDf['cancer_type'] == 'Bladder Cancer']
breasts = exomeDf[exomeDf['cancer_type'] == 'Breast Cancer']
lungCancers = exomeDf[exomeDf['cancer_type'].str.contains('Lung')]
otherApobecRelatedCancers = exomeDf[exomeDf['cancer_type'].isin(otherApobecCancerTypes)]
nonApobecCancers = exomeDf[~exomeDf['cancer_type'].isin(apobecRelatedCancers)]

#dfInferrableSignatures = df[df['Nmut'] > 10]
#bladders = dfInferrableSignatures[dfInferrableSignatures['cancer_type'] == 'Bladder Cancer']
#breasts = dfInferrableSignatures[dfInferrableSignatures['cancer_type'] == 'Breast Cancer']

bladders.to_csv('/Users/friedman/Desktop/otherLandscapePlot/bladders.tsv', sep='\t', index=False)
breasts.to_csv('/Users/friedman/Desktop/otherLandscapePlot/breasts.tsv', sep='\t', index=False)
lungCancers.to_csv('/Users/friedman/Desktop/otherLandscapePlot/lungs.tsv', sep='\t', index=False)
otherApobecRelatedCancers.to_csv('/Users/friedman/Desktop/otherLandscapePlot/otherApobecs.tsv', sep='\t', index=False)
nonApobecCancers.to_csv('/Users/friedman/Desktop/otherLandscapePlot/nonApobecCancers.tsv', sep='\t', index=False)

#############MMR analyses

#-----------------------##########MMR exome
msiList = pd.read_table(pathPrefix + '/ifs/work/taylorlab/pensona/MSI/cohort_list.txt')
msiList['tid'] = msiList['SAMPLE_ID'].apply(lambda x: x[:13])

msiScoreDict = dict(zip(list(msiList['tid']), list(msiList['msisensor'])))
msiSensorDict = dict(zip(list(msiList['tid']), list(msiList['msisensor'])))

mlh1List = set(msiList[(msiList['MLH1'] == 'missense') | (msiList['MLH1'] == 'truncating')]['tid'])
msh2List = set(msiList[(msiList['MSH2'] == 'missense') | (msiList['MSH2'] == 'truncating')]['tid'])
msh3List = set(msiList[(msiList['MSH3'] == 'missense') | (msiList['MSH3'] == 'truncating')]['tid'])
msh6List = set(msiList[(msiList['MSH6'] == 'missense') | (msiList['MSH6'] == 'truncating')]['tid'])
pms1List = set(msiList[(msiList['PMS1'] == 'missense') | (msiList['PMS1'] == 'truncating')]['tid'])
pms2List = set(msiList[(msiList['PMS2'] == 'missense') | (msiList['PMS2'] == 'truncating')]['tid'])

casesWithMMRGeneMutations = mlh1List| msh2List | msh3List | msh6List | pms1List | pms2List

mmrClippingThresh = 5000 # the average mmr case has waaaaay more mutation
exomeDf['NmutAdj'] = exomeDf['Nmut'].apply(lambda x: mmrClippingThresh if x >= mmrClippingThresh else x)
exomeDf['nMutClipped'] = exomeDf['NmutAdj'].apply(lambda x: False if x < mmrClippingThresh else True)

exomeMMRDf = exomeDf[exomeDf['tid'].isin(casesWithMMRGeneMutations)]
exomeNoMMRDf = exomeDf[~exomeDf['tid'].isin(casesWithMMRGeneMutations)]

exomeMMRDf.to_csv('/Users/friedman/Desktop/otherLandscapePlot/exomeMmrPos.tsv', sep='\t', index=False)
exomeNoMMRDf.to_csv('/Users/friedman/Desktop/otherLandscapePlot/exomeMmrNeg.tsv', sep='\t', index=False)

#--------------------------#########MMR impact

impactDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/pensona/dmp_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')

impactDf['pid'] = impactDf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
impactDf['tid'] = impactDf['Tumor_Sample_Barcode'].apply(lambda x: x[:13])
impactDf['cancer_type'] = impactDf['pid'].apply(lambda x: cDict[x] if x in cDict else None)
impactDf = mutationSigUtils.merge_signature_columns(impactDf, mode='Stratton')
impactDf['otherPredominantSigName'] = impactDf.apply(lambda row: find_second_most_common(row, 'mean_MMR', 'name'), axis=1)
impactDf['otherPredominantSigMagnitude'] = impactDf.apply(lambda row: find_second_most_common(row, 'mean_MMR', 'magnitude'), axis=1)
impactDf['msiScore'] = impactDf['tid'].apply(lambda x: msiScoreDict[x] if x in msiScoreDict else None)
impactDf['msiSensor'] = impactDf['tid'].apply(lambda x: msiSensorDict[x] if x in msiSensorDict else None)


impactDfMutsAdj = impactDf.copy()
mutThresh = 10

impactDfMutsAdj['mean_MMR'] = impactDfMutsAdj.apply(lambda row: 0 if row['Nmut'] < mutThresh else row['mean_MMR'], axis=1)
impactDfMutsAdj['otherPredominantSigMagnitude'] = impactDfMutsAdj.apply(lambda row: 0 if row['Nmut'] < mutThresh else row['otherPredominantSigMagnitude'], axis=1)

cancerTypesToMark = set(['Colorectal Cancer', 'Endometrial Cancer', 'Prostate Cancer', 'Breast Cancer', 'Non-Small Cell Lung Cancer', 'Pancreatic Cancer', 'Esophagogastric Cancer', 'Glioma', 'Hepatobiliary Cancer'])
impactDfMutsAdj['cancer_type_adjusted'] = impactDfMutsAdj['cancer_type'].apply(lambda x: x if x in cancerTypesToMark else 'Other')
impactDfMutsAdj['orderingVal'] = impactDfMutsAdj.apply(lambda row: ordering_function_dom_sig_mode(row), axis=1)
msiScoreDict = get_msi_score_dict()
impactDfMutsAdj['msiScore'] = impactDfMutsAdj['tid'].apply(lambda x: msiScoreDict[x] if x in msiScoreDict else None)

print max(list(impactMLH1Mut['msiScore']))

impactMLH1Mut =  impactDfMutsAdj[impactDfMutsAdj['tid'].isin(mlh1List)] 
impactMSH2Mut =  impactDfMutsAdj[impactDfMutsAdj['tid'].isin(msh2List)] 
impactMSH3Mut =  impactDfMutsAdj[impactDfMutsAdj['tid'].isin(msh3List)] 
impactMSH6Mut =  impactDfMutsAdj[impactDfMutsAdj['tid'].isin(msh6List)] 
impactPMS1Mut =  impactDfMutsAdj[impactDfMutsAdj['tid'].isin(pms1List)] 
impactPMS2Mut =  impactDfMutsAdj[impactDfMutsAdj['tid'].isin(pms2List)] 
impactMSIWildtype = impactDfMutsAdj[~impactDfMutsAdj['tid'].isin(casesWithMMRGeneMutations)]

impactMSINotWT = impactDfMutsAdj[impactDfMutsAdj['tid'].isin(casesWithMMRGeneMutations)]


impactMLH1Mut.to_csv('/Users/friedman/Desktop/otherLandscapePlot/impactMLH1Mut.tsv', sep='\t', index=False)
impactMSH2Mut.to_csv('/Users/friedman/Desktop/otherLandscapePlot/impactMSH2Mut.tsv', sep='\t', index=False)
impactMSH3Mut.to_csv('/Users/friedman/Desktop/otherLandscapePlot/impactMSH3Mut.tsv', sep='\t', index=False)
impactMSH6Mut.to_csv('/Users/friedman/Desktop/otherLandscapePlot/impactMSH6Mut.tsv', sep='\t', index=False)
impactPMS1Mut.to_csv('/Users/friedman/Desktop/otherLandscapePlot/impactPMS1Mut.tsv', sep='\t', index=False)
impactPMS2Mut.to_csv('/Users/friedman/Desktop/otherLandscapePlot/impactPMS2Mut.tsv', sep='\t', index=False)
impactMSIWildtype.to_csv('/Users/friedman/Desktop/otherLandscapePlot/impactMSIWildtype.tsv', sep='\t', index=False)

############BRCA somatic status analyses


