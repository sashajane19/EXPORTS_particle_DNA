#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sasha Kramer
MBARI

"""
###Running permanova on EXPORTS 18S data
#Using Aitchison distance

#Install packages needed for analysis
import pandas as pd
import numpy as np
import os
from skbio.stats.distance import permanova, DistanceMatrix
from sklearn.neighbors import DistanceMetric
from statsmodels.sandbox.stats.multicomp import multipletests

dist = DistanceMetric.get_metric('euclidean')

#Import data
os.chdir('~/18S/aitchison/')

NP_feat = pd.read_csv('f1_np_featP.csv',header=None)
NP_feat = np.asmatrix(NP_feat)

NP_ty=pd.read_csv('f1_np_ty.csv',header=None)
NP_ty = np.asarray(NP_ty)

#NP samples type
ind12 = NP_ty[np.where(NP_ty==[1,2])[0]]
print(permanova(DistanceMatrix(dist.pairwise(NP_feat[:,np.where(NP_ty==[1,2])[0]].T)),ind12,permutations=9999))

ind13 = NP_ty[np.where(NP_ty==[1,3])[0]]
print(permanova(DistanceMatrix(dist.pairwise(NP_feat[:,np.where(NP_ty==[1,3])[0]].T)),ind13,permutations=9999))

ind23 = NP_ty[np.where(NP_ty==[2,3])[0]]
print(permanova(DistanceMatrix(dist.pairwise(NP_feat[:,np.where(NP_ty==[2,3])[0]].T)),ind23,permutations=9999))

##Bonferroni correcting pvals
NPty_pvals=[]
NPty_pvals.append(permanova(DistanceMatrix(dist.pairwise(NP_feat[:,np.where(NP_ty==[1,2])[0]].T)),ind12,permutations=9999)['p-value'])
NPty_pvals.append(permanova(DistanceMatrix(dist.pairwise(NP_feat[:,np.where(NP_ty==[1,3])[0]].T)),ind13,permutations=9999)['p-value'])
NPty_pvals.append(permanova(DistanceMatrix(dist.pairwise(NP_feat[:,np.where(NP_ty==[2,3])[0]].T)),ind23,permutations=9999)['p-value'])

NPty_padjusted = multipletests(NPty_pvals, alpha=0.001, method='bonferroni')

print(NPty_padjusted)

#Repeat as needed for all replicate samples and types (e.g., NP particle types, NA sample types, etc.)