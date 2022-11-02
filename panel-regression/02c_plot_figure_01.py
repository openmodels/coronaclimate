#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

variables = ['t2m', 'tp', 'ssrd', 'utci']
RESULTPATH = './results'
FIGUREPATH = './figures'

## ======================================================= ##
## ...
## ======================================================= ##

unitsdict = {
	'absh': 'g m-3',
	'r': '\%',
	't2m': 'degree K',
	'utci': 'index',
	'tp': 'mm day-1',
	'ssrd': 'W m-2'
}

## ======================================================= ##

means_sds = pd.read_csv(os.path.join(RESULTPATH, 'means-sds_firstdiff_L10.csv'))
for var in variables:
	means_sds[var] = [0., means_sds.loc[1:, [col for col in means_sds.columns if var in col]].mean(axis=1)[1]]
sds = means_sds.iloc[1, -4:]
sds['tp'] = sds['tp'] * 24 * 1000

## ======================================================= ##

variables = ['t2m', 'utci', 'ssrd', 'tp']
variable_ids = {'t2m': 4,
				'utci': 1,
				'ssrd': 2,
				'tp': 3}
colors = ['r', 'm', 'y', 'b']

## ======================================================= ##

mobility = 'nomob'
variable_set = 10
nlags = 10

## ======================================================= ##
## ...
## ======================================================= ##

fig = plt.figure(figsize=(10,4))

gs_left = fig.add_gridspec(1, 1)
gs_right = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[1, 1], height_ratios=[1, 1], wspace=0.2)

ax = fig.add_subplot(gs_left[:, 0])
for ivar, variable in enumerate(variables):
	df = pd.read_csv(os.path.join(RESULTPATH, 'coeffscum_firstdiff_L05_pol1_L-cum_{0:s}_{1:d}_{2:s}.csv'.format(mobility, variable_set, variable)))
	coef = df.iloc[0, 1]
	se = df.iloc[0, 2]
	ax.plot(1. * ivar - 0.2, coef, 'o', color=colors[ivar], markersize=5.)
	ax.errorbar(1. * ivar - 0.2, coef, 1.96 * se, color=colors[ivar], markersize=5.)
	print('L05, {0:s}: {1:5.4f} {2:5.4f}; [{3:5.4f}, {4:5.4f}]'.format(variable, coef / sds.loc[variable], 1.96 * se / sds.loc[variable],
		coef / sds.loc[variable] - 1.96 * se / sds.loc[variable], coef / sds.loc[variable] + 1.96 * se / sds.loc[variable]))
for ivar, variable in enumerate(variables):
	df = pd.read_csv(os.path.join(RESULTPATH, 'coeffscum_firstdiff_L10_pol1_L-cum_{0:s}_{1:d}_{2:s}.csv'.format(mobility, variable_set, variable)))
	coef = df.iloc[0, 1]
	se = df.iloc[0, 2]
	ax.plot(1. * ivar, coef, 'o', color=colors[ivar], markersize=5.)
	ax.errorbar(1. * ivar, coef, 1.96 * se, color=colors[ivar], markersize=5.)
	#print('L10, {0:s}: {1:5.4f} {2:5.4f}'.format(variable, coef, se))
for ivar, variable in enumerate(variables):
	df = pd.read_csv(os.path.join(RESULTPATH, 'coeffscum_firstdiff_L10_pol1_L-cum_{0:s}_{1:d}_{2:s}.csv'.format(mobility, variable_set, variable)))
	coef = df.iloc[0, 1]
	se = df.iloc[0, 2]
	ax.plot(1. * ivar + 0.2, coef, 'o', color=colors[ivar], markersize=5.)
	ax.errorbar(1. * ivar + 0.2, coef, 1.96 * se, color=colors[ivar], markersize=5.)
	#print('L20, {0:s}: {1:5.4f} {2:5.4f}'.format(variable, coef, se))
ax.plot(ax.get_xlim(), [0., 0.], 'k-', lw=0.5)
ax.set_xlim(-0.5, len(variables)-0.5)
ax.set_ylabel('Cumulative change in growth rate')
sns.despine(ax=ax, offset=1., right=True, top=True)
ax.set_xticks(np.arange(len(variables)))
ax.set_xticklabels(variables)


colors = ['r', 'm', 'y', 'b']
NLAGS = 5
for ivar, variable in enumerate(variables[:2]):
	ax = fig.add_subplot(gs_right[0, ivar])

	## plot coeffs of lags
	xi = 0.5

	for lag in range(1, NLAGS+1, 1):#, 'pol2']):
		resultfile = 'coeffs_firstdiff_L{0:02d}_pol1_nomob_10.csv'.format(NLAGS)
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
		ax.plot(x, y, color=colors[ivar], lw=1.)
		ax.fill_between(x, y1=y-1.96*yerr, y2=y+1.96*yerr, color=colors[ivar], alpha=0.3, lw=0.)	
		xi += 3.

	if ivar == 1:
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")
		sns.despine(ax=ax, offset=1., right=False, left=True, top=True)
	else:
		sns.despine(ax=ax, offset=1., right=True, top=True)

	ax.set_xlim(0., NLAGS * 3 + 1)
	ax.set_xticks(list(np.arange(2, NLAGS * 3 + 1, 6)))
	ax.set_xticklabels([str(i) for i in np.arange(2, NLAGS * 3 + 1, 6)])

	ax.plot([0., NLAGS * 3.], [0., 0.], 'k-', lw=0.5)

	# add a bit more space for labels
	ylims = ax.get_ylim()
	ax.set_ylim(np.array(ylims) * 1.2)

	t = ax.annotate(s='{0:s} ({1:s})'.format(variable, unitsdict[variable]), xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top')
	t.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

for ivar, variable in enumerate(variables[2:]):
	ax = fig.add_subplot(gs_right[1, ivar])

	## plot coeffs of lags
	xi = 0.5

	for lag in range(1, NLAGS+1, 1):#, 'pol2']):
		resultfile = 'coeffs_firstdiff_L{0:02d}_pol1_nomob_10.csv'.format(NLAGS)
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
		ax.plot(x, y, color=colors[2+ivar], lw=1.)
		ax.fill_between(x, y1=y-1.96*yerr, y2=y+1.96*yerr, color=colors[2+ivar], alpha=0.3, lw=0.)	
		xi += 3.

	if ivar == 1:
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")
		sns.despine(ax=ax, offset=1., right=False, left=True, top=True)
	else:
		sns.despine(ax=ax, offset=1., right=True, top=True)

	ax.set_xlim(0., NLAGS * 3 + 1)
	ax.set_xticks(list(np.arange(2, NLAGS * 3 + 1, 6)))
	ax.set_xticklabels([str(i) for i in np.arange(2, NLAGS * 3 + 1, 6)])

	# add a bit more space for labels
	ylims = ax.get_ylim()
	ax.set_ylim(np.array(ylims) * 1.2)

	ax.plot([0., NLAGS * 3.], [0., 0.], 'k-', lw=0.5)

	t = ax.annotate(s='{0:s} ({1:s})'.format(variable, unitsdict[variable]), xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top')#set_ylabel(variable)
	t.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

fig.text(0.7, 0.03, "Days before reporting (central day of time interval)", rotation="horizontal", va="top", ha="center")
fig.text(0.42, 0.5, "Change in growth rate", rotation="vertical", va="center")

# now the plots are on top of each other, we'll have to adjust their edges so that they won't overlap
gs_left.update(right=0.4)
gs_right.update(left=0.5)

# also, we want to get rid of the horizontal spacing in the left gridspec
gs_left.update(wspace=0)
gs_right.update(wspace=0.2)
gs_right.update(hspace=0.3)

fig.savefig(os.path.join(FIGUREPATH, 'figure_01.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(FIGUREPATH, 'figure_01.png'), bbox_inches='tight')
