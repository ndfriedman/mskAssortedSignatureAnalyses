#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 12:06:48 2018

@author: friedman
"""

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

hotspotDf = pd.read_table(pathPrefix + '/home/gavrilae/noahmaf.txt')
hotspotDf['Ref_Tri']
hotspotDf['quadNuc'] = hotspotDf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2']), axis=1)

hotspotDf.columns.values
hotspotDf['Reference_Allele']
hotspotDf['Tumor_Seq_Allele2']

hotspotDf['Amino_Acid_Change']

Counter(hotspotDf['quadNuc'])

print Counter(hotspotDf[hotspotDf['quadNuc'] == 'ACTC']['Amino_Acid_Change'])

print Counter(hotspotDf['Amino_Acid_Change'])

print Counter(hotspotDf[hotspotDf['Amino_Acid_Change'] == 'G12R']['oncotree_detailed'])