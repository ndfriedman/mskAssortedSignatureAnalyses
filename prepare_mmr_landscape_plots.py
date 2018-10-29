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

#the bi-allelic = 0 is heterozygous, 1 is bi-allelic

otherMMRdata = pd.read_table('~/Desktop/MMR_MSI_merge_uniq_pts_signature.txt')
otherMMRdata[otherMMRdata['nmut'].notnull()].shape

sigDfAnonymized = pd.read_table(pathPrefix + '/ifs/res/taylorlab/jonssonp/msk_impact_gml_som/data/all-somatic-mutations-concatenated.mutsig.txt')
sigDfAnonymized['pid'] = sigDfAnonymized['Sample Name'].apply(lambda x: x[:-8])
#keep track of duplicate ids for later
casesWithDuplicates = set(sigDfAnonymized['pid']) - set(sigDfAnonymized.drop_duplicates(subset=['pid'], keep=False)['pid'])

sigDfAnonymized2 = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/All.dmp_somatic_data_mutations_unfiltered.mafAnno.oncokb.hotspots.sigs_summ.txt')
pivotMonster = sigDfAnonymized2.pivot_table(values='mean', index='Tumor_Sample_Barcode', columns='Signature')
pivotMonster['Tumor_Sample_Barcode'] = pivotMonster.index
sigsDf = pivotMonster
sigsDf['pid'] = sigsDf['Tumor_Sample_Barcode'].apply(lambda x: x[:-8])

mmrDataDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/unique_MMR-4genes.txt')
mmrGermlineMonoallelic = set(mmrDataDf[(mmrDataDf['bi-allelic'] == 0) & (mmrDataDf['Mutation_type'] == 'Germline')]['Patient_ID'])
mmrGermlineBiallelic = set(mmrDataDf[(mmrDataDf['bi-allelic'] > 0) & (mmrDataDf['Mutation_type'] == 'Germline')]['Patient_ID'])
mmrSomaticMonoallelic = set(mmrDataDf[(mmrDataDf['bi-allelic'] == 0) & (mmrDataDf['Mutation_type'] == 'Somatic')]['Patient_ID'])
mmrSomaticBiallelic = set(mmrDataDf[(mmrDataDf['bi-allelic'] > 0) & (mmrDataDf['Mutation_type'] == 'Somatic')]['Patient_ID'])

mmrDataDf['pid'] = mmrDataDf['Patient_ID']

mmrRatioDict = dict(zip(mmrDataDf['Patient_ID'], mmrDataDf['ratio_indel_SNP']))

mmrDataDf[mmrDataDf['nmut'].notnull()].shape, mmrDataDf.shape

nmutDict = dict(zip(otherMMRdata['Patient_ID'], otherMMRdata['nmut']))
sigsDf['Nmut'] = sigsDf['pid'].apply(lambda x: nmutDict[x] if x in nmutDict else None)

mlh1 = set(mmrDataDf[mmrDataDf['Hugo_Symbol'] == 'MLH1']['Patient_ID'])
msh2 = set(mmrDataDf[mmrDataDf['Hugo_Symbol'] == 'MSH2']['Patient_ID'])
msh3 = set(mmrDataDf[mmrDataDf['Hugo_Symbol'] == 'MSH3']['Patient_ID'])
msh6 = set(mmrDataDf[mmrDataDf['Hugo_Symbol'] == 'MSH6']['Patient_ID'])
pms1 = set(mmrDataDf[mmrDataDf['Hugo_Symbol'] == 'PMS1']['Patient_ID'])
pms2 = set(mmrDataDf[mmrDataDf['Hugo_Symbol'] == 'PMS2']['Patient_ID'])


sigsDf['indelSnpRatio'] = sigsDf['pid'].apply(lambda x: mmrRatioDict[x] if x in mmrRatioDict else None)

signatureColumns = ['1', '10', '11', '12', '14', '16', '17', '18', '19', '22', '23', '24', '25', '27', '28', '29',
       '3', '30', '4', '5', '7', '8', '9', 'APOBEC', 'MMR']
sigsToLookAt = ['1', 'APOBEC',  '3', '4', '7','10', '11', '14', '17', '18']

sigsDfReducedCols = sigsDf.drop(columns= ['6', '15', '20', '21', '26', '2', '13', 'MMRa', 'MMRb'], axis=1)

sigsDfReducedCols.columns.values

sigsDf['otherPredominantSigName'] = sigsDf.apply(lambda row:
    mutationSigUtils.find_second_most_common_signature(row, 'MMR', 'name', sigNamesToSpecify = sigsToLookAt, signatureColumns=signatureColumns), axis=1)
sigsDf['otherPredominantSigMagnitude'] = sigsDf.apply(lambda row:
    mutationSigUtils.find_second_most_common_signature(row, 'mean_3', 'magnitude', sigNamesToSpecify = sigsToLookAt, signatureColumns=signatureColumns), axis=1)

sigsDf['orderingVal'] = sigsDf.apply(lambda row:
   signature_landscape_plot_prep_util.ordering_function_dom_sig_mode(row, 'MMR', ageSigReOrderMode=True, sigsToOrderBy = sigsToLookAt), axis=1)

sigsDf['MMR']
sigsDf['otherPredominantSigMagnitude']     
    
mmrGermlineMonoallelicSigs['orderingVal']

mlh1SomaticBiallelicSigs['6']
    
mmrGermlineMonoallelicSigs = sigsDf[sigsDf['pid'].isin(mmrGermlineMonoallelic)]
mmrGermlineBiallelicSigs = sigsDf[sigsDf['pid'].isin(mmrGermlineBiallelic)]
mmrSomaticMonoallelicSigs = sigsDf[sigsDf['pid'].isin(mmrSomaticMonoallelic)] 
mmrSomaticBiallelicSigs = sigsDf[sigsDf['pid'].isin(mmrSomaticBiallelic)]

mlh1GermlineMonoallelicSigs = mmrGermlineMonoallelicSigs[mmrGermlineMonoallelicSigs['pid'].isin(mlh1)]
mlh1GermlineBiallelicSigs = mmrGermlineBiallelicSigs[mmrGermlineBiallelicSigs['pid'].isin(mlh1)]
mlh1SomaticMonoallelicSigs = mmrSomaticMonoallelicSigs[mmrSomaticMonoallelicSigs['pid'].isin(mlh1)]
mlh1SomaticBiallelicSigs = mmrSomaticBiallelicSigs[mmrSomaticBiallelicSigs['pid'].isin(mlh1)]

msh2GermlineMonoallelicSigs = mmrGermlineMonoallelicSigs[mmrGermlineMonoallelicSigs['pid'].isin(msh2)]
msh2GermlineBiallelicSigs = mmrGermlineBiallelicSigs[mmrGermlineBiallelicSigs['pid'].isin(msh2)]
msh2SomaticMonoallelicSigs = mmrSomaticMonoallelicSigs[mmrSomaticMonoallelicSigs['pid'].isin(msh2)]
msh2SomaticBiallelicSigs = mmrSomaticBiallelicSigs[mmrSomaticBiallelicSigs['pid'].isin(msh2)]

msh6GermlineMonoallelicSigs = mmrGermlineMonoallelicSigs[mmrGermlineMonoallelicSigs['pid'].isin(msh6)]
msh6GermlineBiallelicSigs = mmrGermlineBiallelicSigs[mmrGermlineBiallelicSigs['pid'].isin(msh6)]
msh6SomaticMonoallelicSigs = mmrSomaticMonoallelicSigs[mmrSomaticMonoallelicSigs['pid'].isin(msh6)]
msh6SomaticBiallelicSigs = mmrSomaticBiallelicSigs[mmrSomaticBiallelicSigs['pid'].isin(msh6)]

pms2GermlineMonoallelicSigs = mmrGermlineMonoallelicSigs[mmrGermlineMonoallelicSigs['pid'].isin(pms2)]
pms2GermlineBiallelicSigs = mmrGermlineBiallelicSigs[mmrGermlineBiallelicSigs['pid'].isin(pms2)]
pms2SomaticMonoallelicSigs = mmrSomaticMonoallelicSigs[mmrSomaticMonoallelicSigs['pid'].isin(pms2)]
pms2SomaticBiallelicSigs = mmrSomaticBiallelicSigs[mmrSomaticBiallelicSigs['pid'].isin(pms2)]

#_____________________________________________________________________________________________________

mlh1GermlineMonoallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/mlh1GermlineMonoallelicSigs.tsv', sep='\t', index=False)
mlh1GermlineBiallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/mlh1GermlineBiallelicSigs.tsv', sep='\t', index=False)
mlh1SomaticMonoallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/mlh1SomaticMonoallelicSigs.tsv', sep='\t', index=False)
mlh1SomaticBiallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/mlh1SomaticBiallelicSigs.tsv', sep='\t', index=False)

msh2GermlineMonoallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/msh2GermlineMonoallelicSigs.tsv', sep='\t', index=False)
msh2GermlineBiallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/msh2GermlineBiallelicSigs.tsv', sep='\t', index=False)
msh2SomaticMonoallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/msh2SomaticMonoallelicSigs.tsv', sep='\t', index=False)
msh2SomaticBiallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/msh2SomaticBiallelicSigs.tsv', sep='\t', index=False)

msh6GermlineMonoallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/msh6GermlineMonoallelicSigs.tsv', sep='\t', index=False)
msh6GermlineBiallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/msh6GermlineBiallelicSigs.tsv', sep='\t', index=False)
msh6SomaticMonoallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/msh6SomaticMonoallelicSigs.tsv', sep='\t', index=False)
msh6SomaticBiallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/msh6SomaticBiallelicSigs.tsv', sep='\t', index=False)

pms2GermlineMonoallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/pms2GermlineMonoallelicSigs.tsv', sep='\t', index=False)
pms2GermlineBiallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/pms2GermlineBiallelicSigs.tsv', sep='\t', index=False)
pms2SomaticMonoallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/pms2SomaticMonoallelicSigs.tsv', sep='\t', index=False)
pms2SomaticBiallelicSigs.to_csv('/Users/friedman/Desktop/landscapePlots/pms2SomaticBiallelicSigs.tsv', sep='\t', index=False)


mlh1GermlineMonoallelicSigs.columns.values



