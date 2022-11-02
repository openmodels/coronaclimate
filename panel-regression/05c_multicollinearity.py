#!/usr/bin/env python

import os
import numpy as np 
import pandas as pd
import datetime

from statsmodels.stats.outliers_influence import variance_inflation_factor

for transformation in ['', '_firstdiff', '_secdiff']:

	df = pd.read_csv("./data_raw/panel_all.csv")
	df['regid'] = df['Country'].astype(str) +'-'+ df['Region'].astype(str) +'-'+ df['Locality'].astype(str)
	df = df.loc[df['lowest_level'] == 1, :]
	df['superset'] = df['Country']
	df.loc[(pd.isnull(df['Region']) & pd.isnull(df['Locality'])), 'superset'] = 'global'
	df['Date'] = df['Date'].apply(lambda x: datetime.datetime.strptime(x, "%Y-%m-%d"))
	df['days'] = (df['Date'] - datetime.datetime(2020, 1, 1)).apply(lambda x: x.days)
	df['week'] = np.floor(df['days'] / 7.)
	df['dow'] = df['days'] % 7

	models = pd.read_csv('modelselection.csv')

	df_all = pd.DataFrame()
	for i, model in enumerate(models['variables']):
		variables = model.split(',')
		#print(variables)

		if len(variables) < 2:
			continue

		dft = df.sort_values(['regid', 'days'], ascending=True)
		if transformation == '_firstdiff':
			for col in variables:
				dft.loc[:, col] = dft.groupby('regid')[col].diff()
		if transformation == '_secdiff':
			for col in variables:
				dft.loc[:, col] = dft.groupby('regid')[col].diff().diff()

		# data for this model
		X = dft.loc[:, variables].replace([-np.inf, np.inf], np.nan)
		X = X.dropna()
		print(X.shape)

		# VIF dataframe
		df_vif = pd.DataFrame()
		df_vif["feature"] = X.columns

		# calculating VIF for each feature
		df_vif["VIF"] = [variance_inflation_factor(X.values, i)
							for i in range(len(X.columns))]
		df_vif['model'] = models.loc[models['variables'] == model, 'n'].values[0]
		df_all = pd.concat([df_all, df_vif], axis=0, ignore_index=True)
	df_all.to_csv('./results/vif_contemp{0:s}.csv'.format(transformation))


	df_all = pd.DataFrame()
	for i, model in enumerate(models['variables']):
		variables = model.split(',')
		#print(variables)

		if len(variables) < 2:
			continue
		
		# subtract group means (fixed effects)
		dft = df.replace([-np.inf, np.inf], np.nan)
		groups = [['regid', 'dow'], ['regid', 'week'], ['superset', 'Date']]
		for col in variables:
			for group in groups:
				colmean = dft[col].mean()
				dft[col] = dft[col] - dft.groupby(group)[col].transform('mean') + colmean

		dft = dft.sort_values(['regid', 'days'], ascending=True)
		if transformation == '_firstdiff':
			for col in variables:
				dft.loc[:, col] = dft.groupby('regid')[col].diff()
		if transformation == '_secdiff':
			for col in variables:
				dft.loc[:, col] = dft.groupby('regid')[col].diff().diff()

		# data for this model
		X = dft.loc[:, variables].replace([-np.inf, np.inf], np.nan)
		X = X.dropna()

		# VIF dataframe
		df_vif = pd.DataFrame()
		df_vif["feature"] = X.columns

		# calculating VIF for each feature
		df_vif["VIF"] = [variance_inflation_factor(X.values, i)
							for i in range(len(X.columns))]
		df_vif['model'] = models.loc[models['variables'] == model, 'n'].values[0]
		df_all = pd.concat([df_all, df_vif], axis=0, ignore_index=True)
	df_all.to_csv('./results/vif_dm_contemp{0:s}.csv'.format(transformation))

