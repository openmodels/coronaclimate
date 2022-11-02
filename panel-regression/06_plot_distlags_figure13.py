#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

###

NLAGS = 20

###

## plot single plots
ylimdict = {\
	'absh': [-0.1, 0.07],
	't2m': [-0.08, 0.15],
	'utci': [-0.06, 0.02],
	'tp': [-0.02, 0.02],
	'ssrd': [-0.03, 0.03],
	'r': [-0.1, 0.1]
}

variablesdict = {
	'absh': 'Absolute humidity [g m-3]',
	'r': 'Relative humidity [\%]',
	't2m': 'Temperature [degree K]',
	'utci': 'UTCI [-]',
	'tp': 'Precipitation [mm day-1]',
	'ssrd': 'Solar radiation [W m-2]'
}

variables = ['t2m', 'utci', 'tp', 'ssrd'		]

for i, variable in enumerate(variables):

	fig, ax = plt.subplots(figsize=(6,4))

	xi = 0.5

	## plot coeffs of lags
	for lag in range(1, NLAGS+1, 1):
		resultfile = 'coeffs_firstdiff_L20_pol1_nomob_10.csv'
		try:
			df = pd.read_csv(os.path.join('./results', resultfile))
			if df.iloc[:, 0].apply(lambda x: x[:2] == variable[:2]).sum() < 1:
				continue
		except:
			print('File not found: ', os.path.join('./results', resultfile))
			continue
		x = np.array([xi, xi+3-0.1])
		y = np.array([float(df.loc[df.iloc[:, 0] == (variable+'_L'+str(lag)), 'Estimate'].values)]*2)
		yerr = np.array([float(df.loc[df.iloc[:, 0] == variable+'_L'+str(lag), 'Cluster.s.e.'].values)]*2)
		ax.plot(x, y, color='b', lw=1.)
		ax.fill_between(x, y1=y-1.96*yerr, y2=y+1.96*yerr, color='b', alpha=0.3, lw=0.)	
		xi += 3.

	## plot cumulative coeffs
	resultfile = 'coeffscum_firstdiff_L20_pol1_L-cum5_nomob_10_{0:s}.csv'.format(variable)
	try:
		df = pd.read_csv(os.path.join('./results', resultfile))
	except:
		print('File not found: ', os.path.join('./results', resultfile))
		continue
	x = -4.
	y = df['coeffs'].values
	yerr = df['se'].values
	ax.plot(x, y, 'o', color='b', lw=1., markersize=5.)
	ax.errorbar(x, y, yerr=1.96*yerr, fmt='o', capthick=2, markersize=0., lw=1., color='b')

	resultfile = 'coeffscum_firstdiff_L20_pol1_L-cum_nomob_10_{0:s}.csv'.format(variable)
	df = pd.read_csv(os.path.join('./results', resultfile))
	x = -9.
	y = df['coeffs'].values
	yerr = df['se'].values
	ax.plot(x, y, 'o', color='b', lw=1., markersize=5.)
	ax.errorbar(x, y, yerr=1.96*yerr, fmt='o', capthick=2, markersize=0., lw=1., color='b')

	ax.set_xlim(-14., NLAGS * 3 + 1)
	ax.set_xticks([-9., -4.] + list(np.arange(1, NLAGS * 3 + 1, 6)))
	ax.set_xticklabels(['c', 'c5'] + [str(i) for i in np.arange(1, NLAGS * 3 + 1, 6)])
	ax.set_xlabel('Days before reporting')

	ax.annotate(s=variablesdict[variable], xy=(0.95, 0.99), xycoords='axes fraction', ha='right', va='top')
	ax.plot([-14., NLAGS * 3.], [0., 0.], 'k-', lw=0.5)
	ax.set_ylabel('Coefficient of lag')
	sns.despine(ax=ax, offset=1., right=True, top=True)
	fig.savefig(os.path.join('./figures/', 'coeffs_distlags_{0:s}.png'.format(variable)), bbox_inches='tight')

