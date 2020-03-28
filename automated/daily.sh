#!/bin/bash

cd /home/jrising/research/coronavirus/coronaclimate/cases/
python p_create_panel-JH-all.py
/home/jrising/added/Dropbox-Uploader/dropbox_uploader.sh upload panel_john-hopkins.csv "Coronavirus and Climate/cases/john-hopkins/"
