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

import pickle
import socket

import xarray as xr
import geopandas as gpd
import regionmask

## == set parameters for new grid

MAKE_REGIONMASKS = True

FILE_LABEL = 'global_agg15min_2020_ERA'
LABEL_WEIGHTED = 'w2020'

MACHINE = socket.gethostname()

if 'manuel' in MACHINE:
	DATAPATH = '/home/manuel/ResearchData/'
else:
	DATAPATH = os.path.join(os.path.expanduser("~"), 'data')

## == create masks with population weights for each geometry object (shape)

if MAKE_REGIONMASKS == True:

	# 1. create an xarray-dataset with population that has the same shape as the climate data 
	# read in population data and climate data and 
	datapath = os.path.join(DATAPATH, 'gpw/')
	ifile = 'gpw_v4_{0:s}.nc'.format(FILE_LABEL)
	ds_pop = xr.open_dataset(os.path.join(datapath, ifile)).sel(latitude=slice(89.9, -89.9))

	# 2. load the shapefiles and create two masks for each region - population weighted and not-weighted; write them as Python pickle file
	## load the shapefile as a GeoDataFrame using geopandas
	datapath = os.path.join(DATAPATH, 'gadm/gadm36_levels_shp')
	datafile = "gadm36_0.shp"
	gdf_shapes = gpd.read_file(os.path.join(datapath, datafile))

	# set id-column
	ID_COLUMN = 'GID_0'

	# create the mask (a dataset with one multidimensional array, whose values are the "numbers" of the geometry objects to which a grid point belongs)
	nuts_mask_polygons = regionmask.Regions(list(gdf_shapes.geometry.values), name='regions_mask', names=list(gdf_shapes[ID_COLUMN]), abbrevs=list(gdf_shapes[ID_COLUMN]))
	mask = nuts_mask_polygons.mask(ds_pop, lat_name='latitude', lon_name='longitude')

	datapath = os.path.join(DATAPATH, 'ecmwf-covid/region-masks/countries')
	for i, region in enumerate(gdf_shapes[ID_COLUMN]):

		mask_unweighted = ((mask == float(i)) * 1.) / np.sum((mask == float(i)) * 1.)
		mask_weighted = ((mask == float(i)) * ds_pop['population']) / np.sum((mask == float(i)) * ds_pop['population'])
		
		# for some regions, there is no grid point with centroid inside their geometry
		# for these regions, take the "nearest neighbour"
		if mask_unweighted.sum().values < 0.8:
			clat = np.array(gdf_shapes.loc[i, 'geometry'].centroid.coords.xy)[1, 0]
			clon = np.array(gdf_shapes.loc[i, 'geometry'].centroid.coords.xy)[0, 0]
			ilat = np.argmin(np.abs(ds_pop['latitude'].values - clat))
			ilon = np.argmin(np.abs(ds_pop['longitude'].values - clon))
			mask_unweighted[ilat, ilon] = 1.
			mask_weighted[ilat, ilon] = 1.
			"""
			print('Using nearest neighbour for geometry: {region:s}; \n \
				coordinates of centroid: {true_lat:4.3f}, {true_lon:4.3f}; coordinates of nearest neighbour: {nn_lat:4.3f}, {nn_lon:4.3f}'.format(
					**{'region': region, 'true_lat': clat, 'true_lon': clon,
					'nn_lat': ds_pop['latitude'].values[ilat], 'nn_lon': ds_pop['longitude'].values[ilon]
					})
				)
			"""
		ofile = 'weights_{0:s}.p'.format(region)
		with open(os.path.join(datapath, ofile), 'wb') as ofp:
			pickle.dump(mask_unweighted, ofp)

		ofile = 'weights_{0:s}_{1:s}.p'.format(region, LABEL_WEIGHTED)
		with open(os.path.join(datapath, ofile), 'wb') as ofp:
			pickle.dump(mask_weighted, ofp)
