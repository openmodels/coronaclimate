#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import datetime
import socket

import numpy as np
import pandas as pd
from itertools import accumulate

##

NBINS = 5
DATA_RAW = './data_raw/'
DATA = './data/'
GAMMA = 0.1234567901234568

##

datafile = 'panel_all.csv'
df = pd.read_csv(os.path.join(DATA_RAW, datafile))

df['regid'] = df['Country'].astype(str) +'-'+ df['Region'].astype(str) +'-'+ df['Locality'].astype(str)
df = df.sort_values(by=["regid", "Date"], ascending=True)

variables_sum = []

variables_weather = ['absh', 'de', 'q', 'r', 'ssrd', 't2m', 'tp', 'wbgt', 'utci']
variables_mobility = ['mobility_pca1', 'mobility_pca2', 'mobility_pca3', 'mobility_retail_and_recreation', 'mobility_grocery_and_pharmacy', 'mobility_parks', 'mobility_transit_stations', 'mobility_workplaces', 'mobility_residential']
variables_mean = variables_weather + variables_mobility

## create bins
df_bins = pd.DataFrame()

variables_weather_bins = []
for column in variables_weather:

	if column == 'tp':
		bounds = np.array([0., 0.01, 1., 5., 10., 20.]) / (24. * 1000.)
	else:
		bounds = np.percentile(df[column].values, np.linspace(0., 100., NBINS+1))

	for ibound, bound in enumerate(bounds[:-1]):

		df_bins = df_bins.append({'variable': column, 'bin_id': ibound, 'bin_min': bounds[ibound], 'bin_max': bounds[ibound+1]}, ignore_index=True)

		newcolname = column+'_{0:d}'.format(ibound+1)
		df.loc[:, newcolname] = ((df.loc[:, column] >= bound) & (df.loc[:, column] < bounds[ibound+1])).astype(int)
		variables_weather_bins.append(newcolname)

print(variables_weather_bins)
df_bins.to_csv('variables_contemp_bins.csv', index=False)

df = df.sort_values(by=["regid", "Date"], ascending=True)

## first differences
"""
for ilag in range(1, int(NDAYS / LAGS_LENGTH)+1, 1):

	for column in variables_sum + variables_mean:
		df.loc[:, column+'_L'+str(ilag)] = df.loc[:, column+'_L'+str(ilag)] - df.loc[:, column+'_L'+str(ilag+1)]
"""

#df = df.drop(columns=variables_sum+variables_mean)

df["d.Confirmed"] = df["Confirmed"] - df.groupby('regid')["Confirmed"].shift(1)
df.loc[df["d.Confirmed"].isnull(), "d.Confirmed"] = 0.
accumulated = df.groupby('regid')['d.Confirmed'].apply(lambda x: list(accumulate(x, lambda v0, dc: v0 * (1. - float(GAMMA)) + dc)))
df["Infected"] = np.nan
for regid in df['regid'].unique():
	df.loc[df['regid'] == regid, "Infected"] = accumulated[regid]

df["Confirmed.4w"] = df["Confirmed"] - df.groupby('regid')["Confirmed"].shift(28)
df['dlog.4w'] = np.log(df['Confirmed.4w']) - np.log(df.groupby('regid')["Confirmed.4w"].shift(1))
df = df.drop(columns=['Confirmed.4w'])
df['dlog.4w'] = df['dlog.4w'].replace([np.inf, -np.inf], np.nan)

df['dlog'] = np.log(df["Infected"]) - np.log(df.groupby('regid')["Infected"].shift(1))
df = df.drop(columns=['d.Confirmed'])
df['dlog'] = df['dlog'].replace([np.inf, -np.inf], np.nan)

df['Confirmed.lag'] = df.groupby('regid')['Confirmed'].shift(1)
df['dlog.confirmed'] = np.log(df['Confirmed']) - np.log(df['Confirmed.lag'])
df = df.drop(columns=['Confirmed.lag'])
df['dlog.confirmed'] = df['dlog.confirmed'].replace([np.inf, -np.inf], np.nan)

df['Deaths.lag'] = df.groupby('regid')['Deaths'].shift(1)
df['dlog.deaths'] = np.log(df['Deaths']) - np.log(df['Deaths.lag'])
df['mort.pc'] = (df['Deaths'] - df['Deaths.lag']) / df['population']
df = df.drop(columns=['Deaths.lag'])
df['dlog.deaths'] = df['dlog.deaths'].replace([np.inf, -np.inf], np.nan)
df['mort.pc'] = df['mort.pc'].replace([np.inf, -np.inf], np.nan)

df["Confirmed.4w"] = df["Confirmed"] - df.groupby('regid')["Confirmed"].shift(28)
df["Confirmed.2w"] = df["Confirmed"] - df.groupby('regid')["Confirmed"].shift(14)
df["mort.4w"] = (df['Deaths'] - df.groupby('regid')['Deaths'].shift(1)) / df["Confirmed.4w"]
df["mort.2w"] = (df['Deaths'] - df.groupby('regid')['Deaths'].shift(1)) / df["Confirmed.2w"]
df["log.mort.4w"] = np.log(df["mort.4w"])
df = df.drop(columns=["Confirmed.4w", "Confirmed.2w"])
df['mort.4w'] = df['mort.4w'].replace([np.inf, -np.inf], np.nan)
df['mort.2w'] = df['mort.2w'].replace([np.inf, -np.inf], np.nan)
df['log.mort.4w'] = df['log.mort.4w'].replace([np.inf, -np.inf], np.nan)

df['superset'] = df['Country']
df.loc[(pd.isnull(df['Region']) & pd.isnull(df['Locality'])), 'superset'] = 'global'
df['Date'] = df['Date'].apply(lambda x: datetime.datetime.strptime(x, "%Y-%m-%d"))
df['days'] = (df['Date'] - datetime.datetime(2020, 1, 1)).apply(lambda x: x.days)
df['week'] = np.floor(df['days'] / 7.)

## first differences
for column in ['dlog', 'dlog.confirmed', 'mort.4w', 'mort.2w', 'log.mort.4w'] + variables_weather_bins + variables_weather + variables_mobility:
	df[column] = df[column] - df.groupby('regid')[column].shift(1)

df = df.loc[df['Country'].notnull(), :]

## population density
dft = df.loc[df['lowest_level'] == 1, :].groupby('Country')['population_density'].quantile(0.90)
for country in df['Country'].unique():
	if pd.isnull(country):
		continue
	index = df['Country'] == country
	df.loc[index, 'urban'] = (df.loc[index, 'population_density'] >= dft.loc[country]) & (df.loc[index, 'Region'].notnull() | df.loc[index, 'Locality'].notnull())
	df.loc[index, 'rural'] = (df.loc[index, 'population_density'] < dft.loc[country]) & (df.loc[index, 'Region'].notnull() | df.loc[index, 'Locality'].notnull())
df['urban'] = df['urban'].astype(int)
df['rural'] = df['rural'].astype(int)

datafile = 'panel_prepared_firstdiff_contemp.csv'
df.to_csv(os.path.join(DATA, datafile), index=False)
