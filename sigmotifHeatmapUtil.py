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

def sig_ordering_function(x):
    sigNumber = x.strip('Signature ')
    return sigNumber

def trinuc_ordering_function(x):
    d = {'CA': 1, 'CG': 2, 'CT': 3, 'TA': 4, 'TC':5, 'TG':6}
    valD = {'A':0, 'C':1, 'G':2, 'T':3}
    change = x[1] + x[5]
    return d[change] + .1*valD[x[0]] + .01*valD[x[2]]

dfLongFormat = pd.read_table(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/sigMotifsMelted.tsv')
dfLongFormat['sigOrderingCol'] = dfLongFormat['X1'].apply(lambda x: sig_ordering_function(x))
dfLongFormat['trinucOrderingCol'] = dfLongFormat['X2'].apply(lambda x: trinuc_ordering_function(x))

dfLongFormat.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/sigMotifsMelted_adj.tsv', sep='\t', index=False)

dfLongFormat['X2']

dfLongFormat['value']

dfLongFormatAdj = dfLongFormat.copy()
#dfLongFormatAdj['value'] = dfLongFormatAdj['value'].apply(lambda x: x if x > .005 and x < .05 else 0)

krasHotspotMotifs = set(['ACC>ATC', 'ACC>AAC', 'CCA>CAA', 'CCA>CGA', 'ACC>AGC', 'CCA>CTA'])

dfLongFormatAdj['value'] = dfLongFormatAdj.apply(lambda row: row['value'] if row['X2'] in krasHotspotMotifs else 0, axis=1)
dfLongFormatAdj.to_csv(pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/sigMotifsMelted_adjSignif.tsv', sep='\t', index=False)




