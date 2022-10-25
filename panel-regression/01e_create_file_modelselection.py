#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import itertools

import numpy as np
import pandas as pd

variables = ['tp', 'ssrd', 'utci', 't2m']

all_sets = []
for r in range(1, 5, 1):
	this_set = list(itertools.combinations(variables, r))
	all_sets = all_sets + this_set

all_sets = [','.join(this_set) for this_set in all_sets]
df = pd.DataFrame(all_sets)
df.columns = ['variables']
df['n'] = df.index.values
df.to_csv('./modelselection_allvariables.csv', index=False)
