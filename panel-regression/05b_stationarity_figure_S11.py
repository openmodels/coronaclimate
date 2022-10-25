#!/usr/bin/env python

import os
import numpy as np 
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.stattools import kpss

## ===============================
## convenience functions
## ===============================

def adfuller_y(y, mode, minobs=10):
	valid = np.isfinite(y) & pd.notnull(y)
	nobs = np.sum(valid)
	result = [np.nan, np.nan]
	if nobs > minobs:
		try:
			result = adfuller(y[valid], regression=mode, autolag="AIC")
		except:
			pass
	return result 

def kpss_y(y, mode, lags='auto', minobs=10):
	valid = np.isfinite(y) & pd.notnull(y)
	nobs = np.sum(valid)
	result = [np.nan, np.nan]
	if lags != 'auto':
		minobs = np.max([(lags + 2 + 1)*2. + 1., minobs])
	if nobs > minobs:
		try:
			result = kpss(y[valid], regression=mode, nlags=lags)
		except:
			pass
	return result 

## ===============================
## prepare data
## ===============================

df = pd.read_csv("./data/panel_prepared_distlags_L05_dependent.csv")
df = df.loc[df['lowest_level'] == 1, :]
df['dow'] = df['days'] % 7
df['dlog'] = df['dlog'].replace([-np.inf, np.inf], np.nan)
df = df.sort_values(['regid', 'days'], ascending=True)

columns = ['dlog']

# create a matrix
df['regid'] = pd.factorize(df['regid'])[0]
dfm = df.loc[:, ['regid', 'days', 'dlog']].pivot(index='days', values='dlog', columns='regid')

# drop units and time periods with only missing values
dfm = dfm.loc[dfm.isnull().sum(axis=1) != dfm.shape[1], dfm.isnull().sum(axis=0) != dfm.shape[0]]

# set all values to missing values that do not belong to longest uninterrupted period of that unit
for col in dfm.columns:
	dfm['t'] = dfm[col].notnull().groupby(dfm[col].isnull().cumsum()).transform('count')
	dfm.loc[dfm['t'] != dfm['t'].max(), col] = np.nan 

dfm = dfm.drop(columns=['t'])
dfm = dfm.loc[:, dfm.count() > 14.]

## ===============================
## analyse
## ===============================

for i, transformation in enumerate(['nodiff', 'firstdiff', 'seconddiff']):
	dft = [dfm, dfm.diff(), dfm.diff().diff()][i]
	#dft = dft.dropna()
	
	mode = 'nc'
	df_adf_stats_p_1 = dft.apply(lambda x: adfuller_y(x, mode)[1], axis=0)
	mode = 'c'
	df_adf_stats_p_2 = dft.apply(lambda x: adfuller_y(x, mode)[1], axis=0)
	mode = 'ct'
	df_adf_stats_p_3 = dft.apply(lambda x: adfuller_y(x, mode)[1], axis=0)

	mode = 'c'
	df_kpss_stats_p_2 = dft.apply(lambda x: kpss_y(x, mode)[1], axis=0)
	mode = 'ct'
	df_kpss_stats_p_3 = dft.apply(lambda x: kpss_y(x, mode)[1], axis=0)

	df_all = pd.concat([df_adf_stats_p_1,
						df_adf_stats_p_2,
						df_adf_stats_p_3,
						df_kpss_stats_p_2,
						df_kpss_stats_p_3], axis=1, ignore_index=True)
	df_all.columns = ['ADF nc', 'ADF c', 'ADF ct', 'KPSS c', 'KPSS ct']
	df_all.to_csv('./results/stationarity_{0:s}.csv'.format(transformation))

for i, transformation in enumerate(['nodiff', 'firstdiff', 'seconddiff']):
	df_all = pd.read_csv('./results/stationarity_{0:s}.csv'.format(transformation))
	df_all = df_all.drop(columns=['regid'])
	columns = df_all.columns
	columns = 'ADF ct'
	for column in columns:
		fig, ax = plt.subplots()
		ax.hist(df_all[column].values, bins=[0., 0.001, 0.01, 0.05, 0.1, 1.])
		ax.set_xscale('log')
		ax.set_xlabel('p-value')
		ax.set_ylabel('relative frequency')
		ax.annotate(s=column, xy=(0.01, 0.99), xycoords='axes fraction', ha='left', va='top')
		ax.annotate(s='number of units: {0:d}'.format(df_all.shape[0]), xy=(0.99, 0.99), xycoords='axes fraction', ha='right', va='top')
		yticks = ax.get_yticks()
		ax.set_yticklabels(['{0:3.2f}'.format(tick / df_all.shape[0]) for tick in yticks])
		ylim = ax.get_ylim()
		ax.plot([0.05, 0.05], ylim, 'k--', lw=0.5)
		ax.set_ylim(ylim)
		sns.despine(ax=ax, offset=0.05, right=True, left=False, top=True)
		fig.savefig('./figures/hist_{0:s}_{1:s}.pdf'.format(transformation, column))


