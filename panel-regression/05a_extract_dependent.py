#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

import numpy as np
import pandas as pd

df = pd.read_csv('./data/panel_prepared_distlags_L05.csv')
df = df.loc[:, ['dlog', 'dlog.confirmed', 'regid', 'Date', 'superset', 'lowest_level', 'week', 'days']]
df.to_csv('./data/panel_prepared_distlags_L05_dependent.csv')

df = pd.read_csv('./data/panel_prepared_distlags_firstdiff_L05.csv')
df = df.loc[:, ['dlog', 'dlog.confirmed', 'regid', 'Date', 'superset', 'lowest_level', 'week', 'days']]
df.to_csv('./data/panel_prepared_distlags_firstdiff_L05_dependent.csv')
