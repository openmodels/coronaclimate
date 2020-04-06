#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

[ some description ]

"""

__author__ = "Manuel Linsenmeier"
__email__ = "m.linsenmeier@lse.ac.uk"
__version__ = "0.0.1"
__status__ = "Draft"

import sys
import os
import datetime
import socket

import numpy as np
import pandas as pd
import pickle

import geopandas as gpd
import xarray as xr

LABEL_WEIGHTED = 'w2020'

MACHINE = socket.gethostname()

if 'manuel' in MACHINE:
	DATAPATH = '/home/manuel/ResearchData/'
	SCRIPTPATH = os.getcwd()
else:
	DATAPATH = os.path.join(os.path.expanduser("~"), 'data')
	SCRIPTPATH = os.getcwd()

## load the shapefile as a GeoDataFrame using geopandas
datapath = os.path.join(DATAPATH, 'gadm/gadm36_levels_shp')
datafile = "gadm36_0.shp"
gdf_regions = gpd.read_file(os.path.join(datapath, datafile))
ID_COLUMN = 'GID_0'

# load the climate data
datapath = os.path.join(DATAPATH, 'ecmwf-covid')
ifile = 'era5_2020-01-01_2020-03-17_t2m-tp-d2m-sp-ssrd_daymean.nc'

# do not read in data for the two poles, because population data there not defined and therefore cell not included in weights
ds_clim = xr.open_dataset(os.path.join(datapath, ifile)).sel(latitude=slice(89.9, -89.9))
ds_clim = ds_clim.assign_coords(longitude=(((ds_clim.longitude + 180) % 360) - 180)).sortby('longitude')

# check whether panel file already exists
datapath = os.path.join(DATAPATH, 'ecmwf-covid')
ofile = 'panel_climate.csv'
if os.path.exists(os.path.join(datapath, ofile)):
	df_panel = pd.read_csv(os.path.join(datapath, ofile))
	first_date_string = df_panel['TIME'].iloc[-1]
	first_date = datetime.strptime(first_date_string, "%Y-%m-%d")
else:
	# create an empty dataset to store results in form of a panel
	df_panel = pd.DataFrame(columns=['UNIT_ID', 'TIME'])
	first_date = datetime.datetime(2020, 1, 1)

# get first and last date
t = datetime.datetime.now()
last_date = datetime.datetime(t.year, t.month, t.day) - datetime.timedelta(days=5)

# loop over days, compute the statistic for all grid cells, then loop over regions and apply region masks
# first, compute statistic
for date in [first_date + delta for delta in [datetime.timedelta(days=d) for d in range(0, (last_date-first_date).days, 1)]]:

	date_string = date.strftime("%Y-%m-%d")
	statistic_01 = ds_clim.sel(time=date_string)['t2m'].mean()

	for i, geo_ID in enumerate(gdf_regions[ID_COLUMN]):

		# read in the masks as xarray datasets
		datapath = os.path.join(DATAPATH, 'ecmwf-covid/region-masks/countries')
		ifile = 'weights_{0:s}.p'.format(geo_ID)
		with open(os.path.join(datapath, ifile), 'rb') as ifp:
			mask_unweighted = pickle.load(ifp)

		ifile = 'weights_{0:s}_{1:s}.p'.format(geo_ID, LABEL_WEIGHTED)
		with open(os.path.join(datapath, ifile), 'rb') as ifp:
			mask_weighted = pickle.load(ifp)
			
		df_panel = df_panel.append({
			'TIME': date_string, 'UNIT_ID': geo_ID, 
			'T': (statistic_01 * mask_weighted).sum().values, 'T_w': (statistic_01 * mask_unweighted).sum().values,
			}, ignore_index=True)

# write results to output
datapath = os.path.join(DATAPATH, 'ecmwf-covid')
ofile = 'panel_climate.csv'
df_panel.to_csv(os.path.join(datapath, ofile), index=False)
