#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 15:36:20 2020

@author: anademenezes
"""

import requests

r = requests.get('https://covidtrackerapi.bsg.ox.ac.uk/api/stringency/date-range/2020-01-02/2020-04-20')
print(r.content)
open('temp.txt', 'wb').write(r.content)