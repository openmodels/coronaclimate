#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## ======================================================= ##

RESULTPATH = './results/'

## ======================================================= ##
## with cumulative coefficients
## ======================================================= ##

df_models = pd.read_csv('modelselection_allvariables.csv')

df_models = df_models.loc[~df_models['variables'].str.contains('wbgt'), :]
df_models = df_models.loc[~df_models['variables'].str.contains('q'), :]
df_models = df_models.loc[~df_models['variables'].str.contains('absh'), :]
df_models = df_models.loc[~df_models['variables'].apply(lambda x: 'r' in x.split(',')), :]

variable_sets = df_models['n'].unique()
df_all_coeffs = pd.DataFrame()
for i, variable_set in enumerate(variable_sets):
	variables_str = df_models.loc[df_models['n'] == int(variable_set), 'variables'].values
	if ',' in variables_str[0]:
		variables = variables_str[0].split(',')
	else:
		variables = list(variables_str)
	for j, variable in enumerate(variables):
		if (i == 0) & (j == 0):
			df_all_coeffs = pd.read_csv(os.path.join(RESULTPATH, 'coeffscum_firstdiff_modelselection_L10_pol1_L-cum_nomob_{0:d}_{1:s}.csv'.format(variable_set, variable)))
			df_all_coeffs['variable_set'] = variable_set
			df_all_coeffs['variables'] = variables_str
			df_all_coeffs['variable'] = variable
		else:
			df = pd.read_csv(os.path.join(RESULTPATH, 'coeffscum_firstdiff_modelselection_L10_pol1_L-cum_nomob_{0:d}_{1:s}.csv'.format(variable_set, variable)))
			df['variable_set'] = variable_set
			df['variables'] = variables_str
			df['variable'] = variable
			df_all_coeffs = df_all_coeffs.append(df, ignore_index=True)

colors = ['r', 'm', 'y', 'b']
for v, var in enumerate(['t2m', 'utci', 'ssrd', 'tp']):

	fig, ax = plt.subplots(figsize=(4,3))
	for i, variable_set in enumerate(df_all_coeffs['variable_set'].unique()):
		dfs = df_all_coeffs.loc[df_all_coeffs['variable_set'] == variable_set, :]
		if var in dfs['variable'].values:
			dfss = dfs.loc[dfs['variable'] == var, :]
			y = dfss['coeffs']
			yerr = dfss['se'] * 1.96
		else:
			y = np.nan
			yerr = np.nan

		ax.plot(i, y, 'o', color=colors[v], label=var)
		ax.errorbar(i, y, yerr, color=colors[v], markersize=5.)

	xlims = ax.get_xlim()
	ax.plot(xlims, [0., 0.], 'k--', lw=0.5)
	ax.set_xlim(xlims)

	ax.set_xticks(np.arange(0, np.size(df_all_coeffs['variable_set'].unique()), 1))
	ax.set_xticklabels(df_all_coeffs['variables'].unique(), rotation=90.)
	ax.set_xlabel('Combination of variables')
	ax.set_ylabel('Cumulative effect (10 lags)')
	#sns.despine(ax=ax, offset=1., right=True, top=True)
	fig.savefig(os.path.join('./figures', 'multicollinearity_L10_{0:s}.pdf'.format(var)), bbox_inches='tight')
