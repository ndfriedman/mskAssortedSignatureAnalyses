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

exomeSigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/pensona/dmp_sigs/mixed_impact_wes_2017_data_mutations_extended.sigs.tab.txt')
impactSigs = pd.read_table(pathPrefix + '/ifs/work/taylorlab/pensona/dmp_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
exomeSigs['tid'] = exomeSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:13])
impactSigs['tid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:13])
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)

impactSigsComp = impactSigs[impactSigs['tid'].isin(set(exomeSigs['tid']))]
exomeSigsComp = exomeSigs[exomeSigs['tid'].isin(set(impactSigsComp['tid']))]

impactSigsComp = mutationSigUtils.merge_signature_columns(impactSigsComp, smokingMerge=True)
exomeSigsComp = mutationSigUtils.merge_signature_columns(exomeSigsComp, smokingMerge=True)

lungCancers = impactSigsComp[impactSigsComp['cancer_type'] == 'Non-Small Cell Lung Cancer']
lowMutBurdenLung = set(lungCancers[lungCancers['Nmut'] < 10]['tid'])
exomeLowMutBurdenLung = exomeSigsComp[exomeSigsComp['tid'].isin(lowMutBurdenLung)]
exomeLowMutBurdenLung = exomeLowMutBurdenLung[exomeLowMutBurdenLung['Nmut'] > 10]

exomeLowMutBurdenLung['Nmut']

means= []
sigColumnsNoMerge = []
for i in range(1,31):
    print i 
    means.append(round(np.nanmean(exomeLowMutBurdenLung['mean_' + str(i)]),2))
print means

sigColumns = [i for i in exomeSigsComp.columns.values if 'mean' in i]
renameDict = {}
for col in sigColumns:
    renameDict[col] = 'exome_' + col
exomeSigsComp = exomeSigsComp.rename(columns=renameDict)

renameDict = {}
for col in sigColumns:
    renameDict[col] = 'impact_' + col
impactSigsComp = impactSigsComp.rename(columns=renameDict)


sigColumns = [i for i in impactSigsComp.columns.values if 'mean' in i]
tupleList = []
for col in sigColumns:
    tupleList.append((col, np.nanmean(impactSigsComp[col])))
    
mergeDf = exomeSigsComp.merge(impactSigsComp, on='tid')
mergeDf.to_csv('~/Desktop/sigsComp.tsv', sep='\t', index=False)

mergeDf.columns.values