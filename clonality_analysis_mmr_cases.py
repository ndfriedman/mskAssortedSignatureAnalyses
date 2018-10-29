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

import matplotlib.pyplot as plt
    
mafWithClonalityInfo = pd.read_table(pathPrefix + '/ifs/res/taylorlab/bandlamc/somatic_germline/somatic_germline_analysis/cohort_mafs/May_2018/All.dmp_somatic_data_mutations_extended.mafAnno.oncokb.hotspots.zygosity.maf')
mafWithClonalityInfo['ccf_Mcopies']
mafWithClonalityInfo['Tumor_Sample_Barcode']
mafWithClonalityInfo['pid'] = mafWithClonalityInfo['Tumor_Sample_Barcode'].apply(lambda x: x[:-8])

mmrDataDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/unique_MMR-4genes.txt')

otherMMRdata = pd.read_table('~/Desktop/MMR_MSI_merge_uniq_pts_signature.txt')
otherMMRdata['Top_Signature'] = otherMMRdata['Top_Signature'].apply(lambda x: 'MMR' if x in set(['Signature.6', 'Signature.15', 'Signature.20', 'Signature.21', 'Signature.26']) else x)
biallelicNoMMRsig = set(otherMMRdata[otherMMRdata['Top_Signature'] != 'MMR']['Patient_ID']) & set(otherMMRdata[otherMMRdata['bi-allelic'] >= 1]['Patient_ID']) & set(otherMMRdata[otherMMRdata['Hugo_Symbol'].isin(set(['MLH1', 'MSH2', 'MSH6', 'PMS1']))]['Patient_ID'])

mmrDataOtherSigs = otherMMRdata[otherMMRdata['Patient_ID'].isin(biallelicNoMMRsig)]

uvs = mmrDataOtherSigs[mmrDataOtherSigs['Top_Signature'] == 'Signature.7']
uvs['MMRa_contribution']
uvs['Patient_ID']
uvIds = set(uvs['Patient_ID'])

apobecs = mmrDataOtherSigs[(mmrDataOtherSigs['Top_Signature'] == 'Signature.2') | (mmrDataOtherSigs['Top_Signature'] == 'Signature.13')]
apobecs['MMRa_contribution']
apobecIds = set(apobecs['Patient_ID'])

mafBiallelicNoMMRCases = mafWithClonalityInfo[mafWithClonalityInfo['pid'].isin(biallelicNoMMRsig)]

caseIDs = set(mafBiallelicNoMMRCases['pid'])
for caseId in apobecIds:
    localDf = mafWithClonalityInfo[mafWithClonalityInfo['pid'] == caseId]
    fig = plt.hist(localDf['ccf_Mcopies'], range=[0,1])
    filename = caseId + '_clonality.png'
    plt.title(filename + ' Nmut: ' + str(localDf.shape[0]))
    plt.savefig('/Users/friedman/Desktop/mmrClonality/' + filename)
    plt.clf()
    

vals = [1,1,2,3,1,2,3,7]
plt.hist(vals, num_bins=7, range=[0,7])