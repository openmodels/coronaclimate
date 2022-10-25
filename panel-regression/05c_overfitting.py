#!/usr/bin/env python

import os
import numpy as np 
import pandas as pd
import datetime

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

df_variables = pd.read_csv('modelselection_allvariables.csv')
df_variables['variables_list'] = df_variables['variables'].apply(lambda x: x.split(','))
df_variables['variables_n'] = df_variables['variables_list'].apply(lambda x: len(x))

df_variables = df_variables.loc[~df_variables['variables'].str.contains('wbgt'), :]
df_variables = df_variables.loc[~df_variables['variables'].str.contains('q'), :]
df_variables = df_variables.loc[~df_variables['variables'].str.contains('absh'), :]
df_variables = df_variables.loc[~df_variables['variables'].apply(lambda x: 'r' in x.split(',')), :]
df_variables = df_variables.reset_index(drop=True)

df_variables['r2within'] = np.nan
df_variables['delta2'] = np.nan
for i, n in enumerate(df_variables['n'].values):
	df = pd.read_csv('./results/stats_firstdiff_modelselection_L10_pol1_nomob_{0:d}.csv'.format(n))
	df_variables.loc[i, 'r2within'] = df['r2within'].values[0]
	df_variables.loc[i, 'delta2'] = df['delta2'].values[0]

y1 = df_variables['r2within'] * 1000
y2 = df_variables['delta2'] * 1000

corr = np.corrcoef(y1, y2)[0][1]

fig, ax = plt.subplots(figsize=(5, 4))
for i, value in enumerate(y1):
	x = value
	y = y2[i]
	ax.plot(x, y, 'o', markersize=14, alpha=0.7, label='{0:d}: {1:s}'.format(i, df_variables.loc[i, 'variables']))
	ax.annotate(s=str(i), xy=(x, y), xycoords='data', va='center', ha='center')
ax.set_xlabel('R2 within * 1000')
ax.set_ylabel('Crossvalidation prediction error * 1000')
ax.annotate(s='Pearson corr.: {0:3.2f}'.format(corr), xy=(0.95, 0.99), xycoords='axes fraction', va='top', ha='right')
ax.legend(loc='upper right', bbox_to_anchor=(1.5, 1.05))
sns.despine(ax=ax, offset=0.05, right=True, left=False, top=True)
fig.savefig('./figures/scatter_modelfit.pdf', bbox_inches='tight')
fig.savefig('./figures/scatter_modelfit.png', bbox_inches='tight')
