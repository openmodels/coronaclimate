#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

import numpy as np
import scipy

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

resultfiles = [\
	'stats_firstdiff_L10_pol1_nomob_10.csv', # dlog as dep. variable, on weather
	'stats_firstdiff_mobility_L10_pol1_mobpca_99.csv', # dlog as dep. variable, on mobility
	'stats_firstdiff_mortality_L10_pol1_nomob_10.csv', # mortality as dep. variable, on weather
	'stats_firstdiff_contemp_mobility_pol1_mobility_pca1_10.csv' # mobility as dep. variable, on weather
]

datapath = './results'
for resultfile in resultfiles:
	df = pd.read_csv(os.path.join(datapath, resultfile))

	# model M is full model; compare with empty model after subtraction of group means
	SSM = df['SSTwithin'] - df['SSE'] # sample within-group variance minus sum of squared errors
	DFM = df['ncoeffswithin'] # number of variables of model
	MSM = SSM / DFM

	SSE = df['SSE'] # unexplained variance
	DFE = df['nobs'] - df['ncoeffs']# + df['ncoeffswithin'] # dof of error
	MSE = SSE / DFE

	F = MSM / MSE

	# H0 : all beta = 0; reject if p < 0.05
	p_value = 1. - scipy.stats.f.cdf(F, DFM, DFE)

	print(p_value)

## =======================
## effect of adding mobility while controlling for weather

datapath = './results'
resultfile = 'stats_firstdiff_L10_pol1_mobpca_10.csv'
df_1 = pd.read_csv(os.path.join(datapath, resultfile)) # full model
resultfile = 'stats_firstdiff_L10_pol1_mobcntr_10.csv'
df_2 = pd.read_csv(os.path.join(datapath, resultfile)) # restricted model

RSS_F = df_1['SSE'] # unexplained variance
RSS_R = df_2['SSE'] # unexplained variance

DF_F = df_1['nobs'] - df_1['ncoeffs']# + df['ncoeffswithin'] # dof of error
DF_R = df_2['nobs'] - df_2['ncoeffs']# + df['ncoeffswithin'] # dof of error

F = ((RSS_R - RSS_F) / (DF_R - DF_F)) / (RSS_F / DF_F)

# H0 : all beta = 0; reject if p < 0.05
p_value = 1. - scipy.stats.f.cdf(F, DF_R - DF_F, DF_F)

## =======================
## effect of adding weather while controlling for mobility

datapath = './results'
resultfile = 'stats_firstdiff_L10_pol1_mobpca_10.csv'
df_1 = pd.read_csv(os.path.join(datapath, resultfile)) # full model
resultfile = 'stats_firstdiff_mobility_L10_pol1_mobpca_99.csv'
df_2 = pd.read_csv(os.path.join(datapath, resultfile)) # restricted model

RSS_F = df_1['SSE'] # unexplained variance
RSS_R = df_2['SSE'] # unexplained variance

DF_F = df_1['nobs'] - df_1['ncoeffs']# + df['ncoeffswithin'] # dof of error
DF_R = df_2['nobs'] - df_2['ncoeffs']# + df['ncoeffswithin'] # dof of error

F = ((RSS_R - RSS_F) / (DF_R - DF_F)) / (RSS_F / DF_F)

# H0 : all beta = 0; reject if p < 0.05
p_value = 1. - scipy.stats.f.cdf(F, DF_R - DF_F, DF_F)
