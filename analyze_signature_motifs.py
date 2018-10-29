#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 12:22:36 2018

@author: friedman
"""

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
import signature_attribution_util 

def calculate_fraction_cols(df, signatureColumnsToCalculate, spectraDict, suffix='_driver', accessCol='driverTrinucs'):
    for col in signatureColumnsToCalculate:
        driverColName = 'frac_' + col + suffix
        df[driverColName] = df[accessCol].apply(lambda x:
            signature_attribution_util.calculate_driver_fraction(x, col, spectraDict))
    return df

def enumerate_snp_dicts(df):
    snpDict = {}
    tids = set(df['tid'])
    for tid in tids:
        localMaf = df[df['tid'] == tid]
        snpDict[tid] = list(localMaf['quadNuc'])
    return snpDict

signaturesToIgnore= set(['Signature.3','Signature.5','Signature.8','Signature.9','Signature.12','Signature.16','Signature.19','Signature.24','Signature.25','Signature.28','Signature.29','Signature.30'])
spectraD = mutationSigUtils.get_enriched_spectra_for_signatures(pathPrefix=pathPrefix, signaturesToIgnore=signaturesToIgnore)
signaturesToPayAttentionTo = spectraD.keys()

impactMutationsMaf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')
impactMutationsMaf['tid'] = impactMutationsMaf['Tumor_Sample_Barcode'].apply(lambda x: x[:13])
impactMutationsMaf['pid'] = impactMutationsMaf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
impactMutationsMaf['cancer_type'] = impactMutationsMaf['pid'].apply(lambda x: cDict[x] if x in cDict else None)


impactMutationsMaf['quadNuc'] = impactMutationsMaf.apply(lambda row: 
    mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
    if row['Variant_Type'] == 'SNP'
    else None
    , axis=1)
    
impactMutationsMaf['isSNPdriver'] = impactMutationsMaf.apply(lambda row:
    signature_attribution_util.annotate_possible_snp_driver(row), axis=1)
impactMutationsSNPDrivers = impactMutationsMaf[impactMutationsMaf['isSNPdriver'] == True]
snpDriverDict = enumerate_snp_dicts(impactMutationsSNPDrivers)

impactSNPs = impactMutationsMaf[impactMutationsMaf['Variant_Type'] == 'SNP']
impactSNPsDict = enumerate_snp_dicts(impactSNPs)
    
impactSignaturesDf = pd.read_table(pathPrefix + '/ifs/work/taylorlab/pensona/dmp_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSignaturesDf['tid'] = impactSignaturesDf['Tumor_Sample_Barcode'].apply(lambda x: x[:13])

#mnake a column called driver trinucs from the dictionary we did above. clean up so all entries are non except those where there are drivers
impactSignaturesDf['driverTrinucs'] = impactSignaturesDf['tid'].apply(lambda x: 
    snpDriverDict[x] if x in snpDriverDict else None)
impactSignaturesDf['snps'] = impactSignaturesDf['tid'].apply(lambda x: 
    impactSNPsDict[x] if x in impactSNPsDict else None)
#Remove empty entries
impactSignaturesDf['driverTrinucs'] = impactSignaturesDf['driverTrinucs'].apply(lambda x: None if x == None else x if len(x) > 0 else None)
impactSignaturesDf['snps'] = impactSignaturesDf['snps'].apply(lambda x: None if x == None else x if len(x) > 0 else None)

impactSignaturesDfDriverTrinucsOnly = impactSignaturesDf[impactSignaturesDf['driverTrinucs'].notnull()]

smokingSpectra = [] #make the smoking spectra by enumerating all C->A mutations
for letter1 in ['A', 'C', 'T', 'G']:
    for letter2 in ['A', 'C', 'T', 'G']:
        smokingSpectra.append(letter1 + 'CA' + letter2)
spectraD['Signature.4'] = set(smokingSpectra)

impactSignaturesDfDriverTrinucsOnly = calculate_fraction_cols(impactSignaturesDfDriverTrinucsOnly, signaturesToPayAttentionTo, spectraD)
impactSignaturesDfDriverTrinucsOnly = mutationSigUtils.merge_signature_columns(impactSignaturesDfDriverTrinucsOnly, mode='Stratton')
impactSignaturesDfDriverTrinucsOnly.to_csv('/Users/friedman/Desktop/myTestDf.tsv', sep='\t', index=False)


impactSignaturesDfSNPsOnly = impactSignaturesDf[impactSignaturesDf['snps'].notnull()]
impactSignaturesDfSNPsOnly = calculate_fraction_cols(impactSignaturesDfSNPsOnly, signaturesToPayAttentionTo, spectraD, suffix='_snp', accessCol='snps')
impactSignaturesDfSNPsOnly = mutationSigUtils.merge_signature_columns(impactSignaturesDfSNPsOnly, mode='Stratton')
impactSignaturesDfSNPsOnly.to_csv('~/Desktop/myTestDf2.tsv', sep='\t', index=False)

impactMutationsSNPDrivers.to_csv('~/Desktop/mafHardProcess.bak', sep='\t', index=False)
    
impactSignaturesDfSNPsOnly.columns.values
    