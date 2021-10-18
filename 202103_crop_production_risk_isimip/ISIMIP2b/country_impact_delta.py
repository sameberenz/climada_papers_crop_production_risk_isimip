#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 14:03:10 2021

@author: eberenzs
"""



import os
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join
from pathlib import Path

import matplotlib.pyplot as plt

import crop_config as co

data_path = '/Users/eberenzs/Documents/Projects/crop_analysis/data/ISIMIP2b'
output_dir = Path(os.path.join(data_path, 'output/Impact/'))

irr_strings  = {
#    'noirr': "noirr_",
#    'firr': "firr_",
    'noirr-fullexp': "noirr-everywhere_",
    'firr-fullexp': "firr-everywhere_",
    'combi': "combi_",
    }
irr_colors = {
    'combi': "k",
    'noirr': "r",
    'firr': "g",
    'noirr-fullexp': '#FF00FF',
    'firr-fullexp': 'b',
    }
var = {
    'combi': "Production",
    'noirr': "Production_noirr",
    'firr': "Production_firr",
    'noirr-fullexp': 'Production_noirr-everywhere',
    'firr-fullexp': 'Production_firr-everywhere',
    }
labels = {
    'combi': "combi",
    'noirr': "Production_noirr",
    'firr': "Production_firr",
    'noirr-fullexp': 'no irrigation',
    'firr-fullexp': 'full irrigation',
    }


ggcm = co.ggcms[co.input_version][2]
gcm = co.gcms[co.input_version][2].lower()
scenario = "rcp60"
crop = "soy"
country = "ITA"
baseline = pd.read_csv(Path(data_path) / "output" / "baseline_exposure_irr.csv",
                       encoding="ISO-8859-1", header=0)

impact_dict = dict()
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(2,1,1)
for irr_key in irr_strings:
    irr_str = irr_strings[irr_key]
    filename = 'impact_' + irr_str + ggcm +'_'+ gcm +'_'+ scenario + '_'+ crop +'.csv'

    impact_dict[irr_key] = pd.read_csv(output_dir / filename, encoding="ISO-8859-1", header=0)
    values = impact_dict[irr_key][country] + baseline.loc[(baseline.crop==crop) & (baseline.variable==var[irr_key])][country].values[0]

    ax.plot(impact_dict[irr_key]['Year'], values, linestyle='-', color=irr_colors[irr_key], label=labels[irr_key], alpha=.7, zorder=5)

values =  impact_dict['noirr-fullexp'][country] - impact_dict['firr-fullexp'][country]

ax.plot(impact_dict[irr_key]['Year'], values, linestyle=':', color='r', label='Delta', alpha=.9, zorder=5)

irr_ratio = baseline.loc[(baseline.crop==crop) & (baseline.variable=="Irrigation_ratio")][country].values[0]
values = values * (1-irr_ratio)
irr_ratio = baseline.loc[(baseline.crop==crop) & (baseline.variable=="Irrigation_ratio")][country].values[0]
ax.plot(impact_dict[irr_key]['Year'], values, linestyle='-', color='r', label='Delta scaled', alpha=.9, zorder=5)
ax.hlines(0, impact_dict[irr_key]['Year'].min(),  impact_dict[irr_key]['Year'].max())
ax.legend(bbox_to_anchor=(1,1), loc=2)
ax.set_ylabel('Crop production [t/y]')
plt.title("%s %s %s %s %s" %(crop, country, scenario, ggcm, gcm))
plt.show()

### CALC DELTA and DELTA SCALED
"""
filenames_firr = [f.name for f in output_dir.iterdir()
                 if f.is_file() and f.name.startswith('impact_%s' %(irr_strings['firr-fullexp']))]
filenames_noirr = [f.name for f in output_dir.iterdir()
                 if f.is_file() and f.name.startswith('impact_%s' %(irr_strings['noirr-fullexp']))]
filenames_firr.sort()
filenames_noirr.sort()
for idx, fn_firr in enumerate(filenames_firr):
    fn_noirr = fn_firr.replace("_firr-","_noirr-")
    fn_delta = fn_firr.replace("firr-everywhere","delta")
    fn_delta_scaled = fn_firr.replace("firr-everywhere","delta-scaled")
    
    data_firr = pd.read_csv(output_dir / fn_firr, encoding="ISO-8859-1", header=0)
    data_noirr = pd.read_csv(output_dir / fn_noirr, encoding="ISO-8859-1", header=0)
    delta = data_noirr-data_firr
    delta['Year'] = data_firr['Year']
    try:
        delta = delta.drop(columns=['0'])
    except:
        print('')
    try:
        delta = delta.drop(columns=['Unnamed: 0'])
    except:
        print('')
    delta.to_csv(output_dir/ fn_delta, index=0)
    
    for col in delta.columns:
        if len(col) != 3:
            continue
        else:
            irr_ratio = baseline.loc[(baseline.crop==fn_noirr[-7:-4]) & (baseline.variable=="Irrigation_ratio")][col].values[0]
            delta[col] = delta[col] * (1-irr_ratio)
    delta.to_csv(output_dir/ fn_delta_scaled, index=0)
"""
