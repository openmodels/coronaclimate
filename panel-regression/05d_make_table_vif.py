#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

import numpy as np
import pandas as pd

## ================================== ##

var2name = {\
		'absh': 'x',
		'q': 'x',
		'de': 'x',
		'r': 'x',
		'ssrd': 'x',
		'tp': 'x',
		'wbgt': 'x',
		'utci': 'x',
			}

variables = var2name.keys()
variables = ['t2m', 'utci', 'ssrd', 'tp']

## ================================== ##

methods = ['contemp', 'dm_contemp', 'dm_contemp_firstdiff', 'dm_contemp_secdiff']

# read in data
df = pd.DataFrame()

for method in methods:
	ifile = 'vif_{0:s}.csv'.format(method)
	datapath = './results'	
	df_this = pd.read_csv(os.path.join(datapath, ifile))
	df_this['method'] = method
	df_this = pd.pivot_table(df_this, index=['model', 'method'], columns='feature', values='VIF').reset_index()
	#df_this.columns[:2] = df_this.columns[0][:2]
	df = pd.concat([df, df_this], ignore_index=True)

df.index = df['method']

## ================================== ##

labels = ['levels', 'levels', 'first-differences', 'second-differences']

df = df.loc[df['model'] == 12, :]

# make table with vif

line0 = '\\begin{tabular}{ L{3.0cm} '
for x in methods:
	line0 = line0 + r' R{2.2cm} '
line0 = line0 + ' }\\\\ \n'

line1 = '\\toprule \\midrule \n'
line2 = 'Variable'
for i, x in enumerate(methods):
	line2 = line2 + ' & ' + '{0:s}'.format(labels[i])
line2 = line2 + '\\\\ \n'
line3 = '\\midrule \n'
lines_header = [line0, line1, line2, line3]

lines_core = []
for i, x1 in enumerate(variables):
	line = '{0:s}'.format(x1)
	for x2 in methods:
		line = line + ' & ' + '{0:7.2f}'.format(df.loc[x2, x1]).replace('nan', '')
	line = line + '\\\\ \n'
	lines_core.append(line)

line5 = '\\midrule\n'
line6 = '\\bottomrule\n'
line7 = '\\end{tabular}\n'
lines_footer = [line5, line6, line7]

lines = lines_header + lines_core + lines_footer

with open('./tables/table_vif.tex', 'w') as ofp:
	ofp.writelines(lines)
