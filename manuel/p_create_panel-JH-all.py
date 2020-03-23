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
import socket

import numpy as np
import pandas as pd
import datetime


MACHINE = socket.gethostname()

if 'manuel' in MACHINE:
	DATAPATH = '/home/manuel/ResearchData/covid-19/john-hopkins'

# read in JH data
datapath = os.path.join(DATAPATH, "COVID-19/csse_covid_19_data/csse_covid_19_time_series")
datafile = "time_series_19-covid-Confirmed.csv"
df_JH_confirmed = pd.read_csv(os.path.join(datapath, datafile))

# read in JH data
datapath = os.path.join(DATAPATH, "COVID-19/csse_covid_19_data/csse_covid_19_time_series")
datafile = "time_series_19-covid-Deaths.csv"
df_JH_deaths = pd.read_csv(os.path.join(datapath, datafile))

# read in JH data
datapath = os.path.join(DATAPATH, "COVID-19/csse_covid_19_data/csse_covid_19_time_series")
datafile = "time_series_19-covid-Recovered.csv"
df_JH_recovered = pd.read_csv(os.path.join(datapath, datafile))


df_all = pd.DataFrame(columns=['Country', 'Region', 'Locality', 'Date', 'Confirmed', 'Deaths', 'Recovered', 'Source'])

dates = [datetime.datetime.strptime(t, "%m/%d/%y") for t in df_JH_confirmed.columns[4:]]
n = np.size(dates)

for i in np.arange(df_JH_confirmed.shape[0]):
	country = df_JH_confirmed.iloc[i, :].loc['Country/Region']
	region = df_JH_confirmed.iloc[i, :].loc['Province/State']
	confirmed = df_JH_confirmed.iloc[i, 4:].values

	indices_deaths = ((df_JH_deaths['Country/Region'] == country) & (df_JH_deaths['Province/State'] == region))
	if np.size(df_JH_deaths.loc[indices_deaths, :]) > 0:
		deaths = df_JH_deaths.loc[indices_deaths, :].iloc[:, 4:].values[0]
	else:
		deaths = [np.nan]*n

	indices_recovered = ((df_JH_recovered['Country/Region'] == country) & (df_JH_recovered['Province/State'] == region))
	if np.size(df_JH_recovered.loc[indices_recovered, :]) > 0:
		recovered = df_JH_recovered.loc[indices_recovered, :].iloc[:, 4:].values[0]
	else:
		recovered = [np.nan]*n

	df_all = pd.concat([df_all, pd.DataFrame({'Country': [country]*n, 'Region': [region]*n, 'Locality': [np.nan]*n, 'Date': dates, 'Confirmed': confirmed, 'Deaths': deaths, 'Recovered': recovered})], axis=0)

df_all['Source'] = 'John Hopkins'

# write to file
datapath = DATAPATH
datafile = "panel_john-hopkins.csv"
df_all.to_csv(os.path.join(datapath, datafile), index=False)

