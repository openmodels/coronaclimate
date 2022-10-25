#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RESULTPATH = './results/'
FIGUREPATH = './figures/'

variables = ['t2m', 'utci', 'ssrd', 'tp']
variable_ids = {'t2m': 4,
				'utci': 1,
				'ssrd': 2,
				'tp': 3}
colors = ['r', 'm', 'y', 'b']
nbins = 5

bin_bounds = pd.read_csv('variables_bins_L10.csv')
bin_bounds.loc[bin_bounds['variable'] == 't2m', ['bin_max', 'bin_min']] = bin_bounds.loc[bin_bounds['variable'] == 't2m', ['bin_max', 'bin_min']].values - 273.15
bin_bounds.loc[bin_bounds['variable'] == 'tp', ['bin_max', 'bin_min']] = np.abs(bin_bounds.loc[bin_bounds['variable'] == 'tp', ['bin_max', 'bin_min']].values * 1000. * 24.)

fig = plt.figure(figsize=(8,4))

gs_right = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[1, 1], height_ratios=[1, 1], wspace=0.2)

missing_bin = int(np.median(np.arange(1, nbins+1, 1)))
for ivar, variable in enumerate(variables[:2]):
	ax = fig.add_subplot(gs_right[0, ivar])
	x = []
	y = []
	yerr = []
	for ibin in range(1, nbins+1, 1):
		x_i = ibin
		if ibin == missing_bin:
			y_i = 0.
			yerr_i = 0.
		else:
			df = pd.read_csv(os.path.join(RESULTPATH, 'coeffscum_distlag_bins_L10_bins_L-cum_nomob_10.csv'))
			df = df.loc[df['variable'] == (variable + '_' + str(ibin)), :]
			y_i = float(df['coef'].values)
			yerr_i = float(df['se'].values)
		x.append(x_i)
		y.append(y_i)
		yerr.append(yerr_i)
	y = np.array(y)
	yerr = np.array(yerr)
	ax.plot(x, y, 'o-', color=colors[ivar], lw=1., markersize=5.)
	ax.fill_between(x, y1=y-1.96*yerr, y2=y+1.96*yerr, color=colors[ivar], alpha=0.3, lw=0.)
	ax.plot(ax.get_xlim(), [0., 0.], 'k-', lw=0.5)
	if ivar == 1:
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")
		sns.despine(ax=ax, offset=1., right=False, left=True, top=True)
	else:
		sns.despine(ax=ax, offset=1., right=True, top=True)
	ax.set_xticks(np.arange(0.5, nbins + 1.5, 1))
	bin_indices = bin_bounds['variable'] == variable
	ax.set_xticklabels(['{0:3.0f}'.format(s) for s in list(np.hstack((bin_bounds.loc[bin_indices, 'bin_min'].values, bin_bounds.loc[bin_indices, 'bin_max'].values[-1])))], rotation=0., size='x-small')
	t = ax.annotate(s=variable, xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top')#set_ylabel(variable)
	t.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='none'))

for ivar, variable in enumerate(variables[2:]):
	ax = fig.add_subplot(gs_right[1, ivar])
	x = []
	y = []
	yerr = []
	for ibin in range(1, nbins+1, 1):
		x_i = ibin
		if ibin == missing_bin:
			y_i = 0.
			yerr_i = 0.
		else:
			df = pd.read_csv(os.path.join(RESULTPATH, 'coeffscum_distlag_bins_L10_bins_L-cum_nomob_10.csv'))
			df = df.loc[df['variable'] == (variable + '_' + str(ibin)), :]
			y_i = float(df['coef'].values)
			yerr_i = float(df['se'].values)
		x.append(x_i)
		y.append(y_i)
		yerr.append(yerr_i)
	y = np.array(y)
	yerr = np.array(yerr)
	ax.plot(x, y, 'o-', color=colors[2+ivar], lw=1., markersize=5.)
	ax.fill_between(x, y1=y-1.96*yerr, y2=y+1.96*yerr, color=colors[2+ivar], alpha=0.3, lw=0.)
	ax.plot(ax.get_xlim(), [0., 0.], 'k-', lw=0.5)
	if ivar == 1:
		ax.yaxis.tick_right()
		ax.yaxis.set_label_position("right")
		sns.despine(ax=ax, offset=1., right=False, left=True, top=True)
	else:
		sns.despine(ax=ax, offset=1., right=True, top=True)
	ax.set_xticks(np.arange(0.5, nbins + 1.5, 1))
	bin_indices = bin_bounds['variable'] == variable
	if variable == 'tp':
		ax.set_xticklabels(['{0:3.0f}'.format(s) for s in list(np.hstack((bin_bounds.loc[bin_indices, 'bin_min'].values, bin_bounds.loc[bin_indices, 'bin_max'].values[-1])))], rotation=0., size='x-small')
	else:
		ax.set_xticklabels(['{0:3.0f}'.format(s) for s in list(np.hstack((bin_bounds.loc[bin_indices, 'bin_min'].values, bin_bounds.loc[bin_indices, 'bin_max'].values[-1])))], rotation=0., size='x-small')
	t = ax.annotate(s=variable, xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top')#set_ylabel(variable)
	t.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='none'))
	print(variable, y)

fig.text(0.97, 0.5, "Cumulative effect of bin\n (relative to missing bin)", rotation="vertical", va="center")

# now the plots are on top of each other, we'll have to adjust their edges so that they won't overlap
gs_right.update(left=0.5)

fig.savefig(os.path.join(FIGUREPATH, 'figure_nonlinear_L10_bins.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(FIGUREPATH, 'figure_nonlinear_L10_bins.png'), bbox_inches='tight')