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
import copy
import socket

import numpy as np
import subprocess
import datetime

import xarray as xr
import cdsapi

# === 

DOWNLOAD = True
PROCESS = True

# === 

variables = ['2m_temperature', 'total_precipitation', '2m_dewpoint_temperature', 'surface_pressure', 'surface_solar_radiation_downwards']
variables_short = ['t2m', 'tp', 'd2m', 'sp', 'ssrd'] # see https://apps.ecmwf.int/codes/grib/param-db

first_date = datetime.datetime(2020, 3, 14)
t = datetime.datetime.now()
last_date = datetime.datetime(t.year, t.month, t.day) - datetime.timedelta(days=6)
#datetime.datetime(2020, 3, 10)

# === 

MACHINE = socket.gethostname()

if 'manuel' in MACHINE:
	DATAPATH = '/home/manuel/ResearchData/ecmwf-covid'
	SCRIPTPATH = os.getcwd()
else:
	DATAPATH = '/users/linsenme/data/ecmwf-covid'
	SCRIPTPATH = os.getcwd()

# === 

def correct_expver_era5(ifile):
	ds = xr.open_dataset(ifile)
	if 'expver' not in ds.coords:
		return
	else:
		ds_new = copy.deepcopy(ds)
		for var_short in ds.variables:
			if (var_short != 'expver') and ('expver' in ds[var_short].coords):

				array_expver1 = ds[var_short].sel(expver=1)
				array_expver5 = ds[var_short].sel(expver=5)

				# merge the two arrays based on missing values (assuming they are perfect complements)
				merged_array = np.where(array_expver1 == np.nan, array_expver5, array_expver1)

				# debug
				#merged_array_two = np.where(array_expver5 == np.nan, array_expver1, array_expver5)
				#print(np.sum(merged_array - merged_array_two))

				ds_new = ds_new.drop_vars(var_short)
				ds_new = ds_new.assign({
					var_short: (
			    		tuple(['time', 'latitude', 'longitude']),
						merged_array,
						ds[var_short].attrs
						)
					})

		# remove the dimension no longer required
		ds_new = ds_new.drop_vars('expver')
		ds_new.to_netcdf(ifile)
		return

# === 

# check whether dates specifiec correctly
if (first_date.hour + last_date.hour) != 0:
	raise Exception("Please set the hour of both the first and the last date to 0.")

# get datetime object for every 3 hour in time period; hours are 1, 4, 7, ..., 22, corresponding to model intervals 0-1, 3-4, ..., 21-22
delta = last_date - first_date
dates = [first_date + datetime.timedelta(hours=h) for h in range(1, ((delta.days + 1) * 24) + 1, 3)]

# get lists of years and months (to loop over them when downloading the netcdf files)
years = list(set([date.year for date in dates]))
months = list(set([date.month for date in dates]))

# create empty list of time periods
time_periods = []

# loop over months
for month in months:

	# get dates of the month
	dates_this_month = [date for date in dates if (date.month == month)]
	first_date_this_month = dates_this_month[0]
	last_date_this_month = dates_this_month[-1]

	# append time period to list of time periods (used below for processing the files)
	time_periods.append('{first_date:s}_{last_date:s}'.format(
		**{'first_date': first_date_this_month.strftime("%Y-%m-%d"), 'last_date': last_date_this_month.strftime("%Y-%m-%d")}))

	# create lists of days and hours to pass to API
	days = list(set([date.day for date in dates_this_month]))
	hours = list(set([date.hour for date in dates_this_month]))

	# loop over variables and download monthly netcdf files
	if DOWNLOAD == True:

		for i, variable in enumerate(variables):

			variable_short = variables_short[i]

			newfilename = 'era5_{first_date:s}_{last_date:s}_{variable_short:s}.nc'.format(
				**{'first_date': first_date_this_month.strftime("%Y-%m-%d"), 'last_date': last_date_this_month.strftime("%Y-%m-%d"),
				'variable_short': variable_short})
			newfile = os.path.join(DATAPATH, newfilename)

			if not os.path.exists(newfile):

				c = cdsapi.Client()

				print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
				print('Downloading now: month ', month, ' , variable: ', variable)

				tempfile = "temp.nc"
				
				r = c.retrieve(
				    "reanalysis-era5-single-levels",
				    {
				       "variable": variable,
				       "product_type": "reanalysis",
				       "year": years,
				       "month": month,
				       "day": days,
				       "time": hours,
				       "format": "netcdf",
				    },
				)
				r.download(tempfile)
				correct_expver_era5(tempfile)
				os.rename(tempfile, newfile)

# now do some post-processing
time_periods = list(set(time_periods))
merged_files = []

if PROCESS == True:

	# create one file with all variables and daily means for the whole time period
	# merge the files for the different time periods
	merged_file_allperiods = os.path.join(
		DATAPATH,
		'era5_{first_date:s}_{last_date:s}_{variables_string:s}_daymean.nc'.format(
		**{
			'first_date': first_date.strftime("%Y-%m-%d"),
			'last_date': last_date.strftime("%Y-%m-%d"),
		 	'variables_string': '-'.join(variables_short)
		 }))

	if not os.path.exists(merged_file_allperiods):

		# create one cdo command to merge all monthly files for every variable
		cdo_commands_mergetime = ''
		tempfiles = []
		for variable_short in variables_short:

			nc_files_variable = ' '.join([os.path.join(DATAPATH, 'era5_' + time_period + '_' + variable_short + '.nc') for time_period in time_periods])
			tempfile = variable_short + '.nc'
			cdo_commands_mergetime += 'cdo -b F32 mergetime ' + nc_files_variable + ' ' + os.path.join(DATAPATH, tempfile) + '; '
			tempfiles.append(os.path.join(DATAPATH, tempfile))

		nc_files_allperiods = ' '.join(tempfiles)

		# define the equations to compute humidity variables
		humidity_equations_file = os.path.join(SCRIPTPATH, 'eq_humidity.txt')
		lines = ['r = 10 ^ (7.591386 * ((d2m/(d2m + 240.7263)) - (t2m/(t2m + 240.7263)))) ; \n',
				'absh = ((6.1120 * 2.71828 ^ ((17.670 * (t2m - 273.15))/(t2m + 243.5 - 273.15)) * r * 2.1674))/(t2m) ; \n']
		with open(humidity_equations_file, 'w') as ofp:
			ofp.writelines(lines)

		# merge different variables in one file using these cdo commands, compute relative and absolute humidity for every time step, then compute daily means
		cdo_command = '{cdo_commands_mergetime:s} cdo -b F32 daymean -shifttime,-1hour -merge {nc_files_allperiods:s} -exprf,{humidity_equations_file:s} -merge {nc_files_allperiods:s} {merged_file_allperiods:s}'.format(
			**{'cdo_commands_mergetime': cdo_commands_mergetime, 'nc_files_allperiods': nc_files_allperiods, 'humidity_equations_file': humidity_equations_file, 'merged_file_allperiods': merged_file_allperiods})
		result = subprocess.check_output(cdo_command, shell=True)
		print(result)
		for tempfile in tempfiles:
			os.remove(tempfile)
