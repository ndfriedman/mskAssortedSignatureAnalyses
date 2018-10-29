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

hotspotMaf = pd.read_table(pathPrefix + '/ifs/res/taylorlab/hotspot_mafs/supermaf_2018/maf_42k_FINAL0802.txt')

hotspotMaf['Reference_Amino_Acid']
hotspotMaf[hotspotMaf['Hugo_Symbol'] == 'KRAS'].shape

hotspotMaf.columns.values

KRAS = hotspotMaf[hotspotMaf['Hugo_Symbol'] == 'KRAS']

g13s = KRAS[KRAS['Amino_Acid_Position'] == '61']
for variant in set(g13s['Variant_Amino_Acid']):
    krasTypeMuts = g13s[g13s['Variant_Amino_Acid'] == variant]
    nCases = krasTypeMuts.shape[0]
    if nCases > 50:
        for refTri in set(krasTypeMuts['Ref_Tri']):
            mutsAtRefTri = krasTypeMuts[krasTypeMuts['Ref_Tri']  == refTri]
            if mutsAtRefTri.shape[0] > 50:
                print variant, refTri, mutsAtRefTri.shape[0], set(mutsAtRefTri['Tumor_Seq_Allele2']) 

tempMaf = hotspotMaf[hotspotMaf['Hugo_Symbol'].isin(set(['TP53', 'KRAS', 'PIK3CA', 'PTEN', 'AKT1', 'ERBB2', 'CTNNB1']))]

tempMaf.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hotspotReducedAnalysis10-19.tsv', sep='\t', index=False)



#dfSaved = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/hotspotReducedAnalysis10-19.tsv')
#dfSaved['Amino_Acid_Change']
Counter(dfSaved['Tumor_Sample_Barcode'])











#a146s--->MMR??

a146s = KRAS[KRAS['Amino_Acid_Position'] == '146']
a146s[a146s['Tumor_Seq_Allele2'] == ''].shape
Counter(a146s['Ref_Tri'])
print set(a146s['Ref_Tri']), set(a146s['Tumor_Seq_Allele2'])