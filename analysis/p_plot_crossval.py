#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

[some description]

"""

__author__ = "Manuel Linsenmeier"
__email__ = "m.linsenmeier@lse.ac.uk"
__version__ = "0.0.1"
__status__ = "Draft"

import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


datapath = '~/Dropbox/Coronavirus and Climate/results'
datafile = 'crossval_results_update.csv'
df = pd.read_csv(os.path.join(datapath, datafile))

#df = pd.concat([df, df2])
#df.to_csv('../results/crossval_results.csv', index=False)

for country in df['country'].unique():
	for growth in df['growth'].unique():
		df_red = df.loc[(df['country'] == country) & (df['growth'] == growth), :]

		Z = df_red.pivot(index='last', columns='first', values='rsqr').values

		fig, ax = plt.subplots(figsize=(6,6))
		ax.set_title('result of cross-validation: $R^{2}$ \n country: ' + country + ' \n time period: ' + growth)
		cax = ax.imshow(Z, origin='lower')
		sns.despine(ax=ax, offset=1., right=True, top=True)
		ax.set_xlabel("First day for averaging (days until recorded)")
		ax.set_ylabel("Last day for averaging (days until recorded)")
		ax.set_xticks(range(0, Z.shape[0], 2))
		ax.set_yticks(range(0, Z.shape[1], 2))
		ax.set_xticklabels([str(i+1) for i in ax.get_xticks()])
		ax.set_yticklabels([str(i+1) for i in ax.get_yticks()])
		cbar = plt.colorbar(cax)
		for i in range(0, 10):
			x, y = df_red.sort_values(by='rsqr', ascending=False).loc[:, ['first', 'last']].values[i]
			#print(country, i, x, y)
			ax.annotate(s=str(i+1), xy=(x-1, y-1), xycoords='data', ha='center', va='center')
		fig.savefig('../../figures/crossval_{0:s}_{1:s}_update.png'.format(country, growth).replace(' ', '-'), bbox_inches='tight')
"""

datapath = '~/Dropbox/Coronavirus and Climate'
datafile = 'crossval_results_lags.csv'
df = pd.read_csv(os.path.join(datapath, datafile))

#df = pd.concat([df, df2])
#df.to_csv('../results/crossval_results.csv', index=False)

for lagmode in df['lag.mode'].unique():
	for nobs in df['nobs.min'].unique():

		df_red = df.loc[(df['lag.mode'] == lagmode) & (df['nobs.min'] == nobs), :]

		Z = df_red.pivot(index='last', columns='first', values='rsqr').values

		fig, ax = plt.subplots(figsize=(6,6))
		fig.suptitle('result of cross-validation: $R^{2}$, ' + '{0:s} lags, min. obs.: {1:d}'.format(lagmode, nobs))
		cax = ax.imshow(Z, origin='lower', vmin=0., vmax=0.0015)
		sns.despine(ax=ax, offset=1., right=True, top=True)
		ax.set_xlabel("First day of averaging period (days until onset of symptoms)")
		ax.set_ylabel("Last day of averaging period (days until onset of symptoms)")
		ax.set_xticks(range(0, Z.shape[0], 2))
		ax.set_yticks(range(0, Z.shape[1], 2))
		ax.set_xticklabels([str(i+1) for i in ax.get_xticks()])
		ax.set_yticklabels([str(i+1) for i in ax.get_yticks()])
		cbar = plt.colorbar(cax)
		for i in range(0, 10):
			x, y = df_red.sort_values(by='rsqr', ascending=False).loc[:, ['first', 'last']].values[i]
			#print(country, i, x, y)
			ax.annotate(s=str(i+1), xy=(x-1, y-1), xycoords='data', ha='center', va='center')
		fig.savefig('../../figures/crossval_{0:s}-lags_nobs{1:d}.png'.format(lagmode, nobs), bbox_inches='tight')
"""
