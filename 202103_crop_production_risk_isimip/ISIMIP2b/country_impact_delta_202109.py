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

"""
ggcm = co.ggcms[co.input_version][2]
gcm = co.gcms[co.input_version][2].lower()
scenario = "rcp60"
crop = "soy"
country = "ITA"
"""
ggcm = 'gepic'
gcm = co.gcms[co.input_version][0].lower()
scenario = "historical"
crop = "mai"
country = "BGD"
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

# values =  impact_dict['noirr-fullexp'][country] - impact_dict['firr-fullexp'][country]
delta_rel =  impact_dict['noirr-fullexp'][country] / baseline.loc[(baseline.crop==crop) & (baseline.variable==var['noirr-fullexp'])][country].values[0] \
            - impact_dict['firr-fullexp'][country] / baseline.loc[(baseline.crop==crop) & (baseline.variable==var['firr-fullexp'])][country].values[0]
delta_combi = delta_rel * baseline.loc[(baseline.crop==crop) & (baseline.variable==var['combi'])][country].values[0]


ax.plot(impact_dict[irr_key]['Year'], delta_combi, linestyle=':', color='r', label='Delta', alpha=.9, zorder=5)

irr_ratio = baseline.loc[(baseline.crop==crop) & (baseline.variable=="Irrigation_ratio")][country].values[0]
delta_combi_scaled = delta_combi * (1-irr_ratio)
irr_ratio = baseline.loc[(baseline.crop==crop) & (baseline.variable=="Irrigation_ratio")][country].values[0]
ax.plot(impact_dict[irr_key]['Year'], delta_combi_scaled, linestyle='-', color='r', label='Delta scaled', alpha=.9, zorder=5)
ax.hlines(0, impact_dict[irr_key]['Year'].min(),  impact_dict[irr_key]['Year'].max())
ax.legend(bbox_to_anchor=(1,1), loc=2)
ax.set_ylabel('Crop production [t/y]')
plt.title("%s %s %s %s %s" %(crop, country, scenario, ggcm, gcm))
plt.show()

### CALC DELTA and DELTA SCALED

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
    fn_delta_rel = fn_firr.replace("firr-everywhere","delta-rel")
    fn_delta_rel_scaled = fn_firr.replace("firr-everywhere","delta-rel-scaled")
    crop = fn_firr.split(sep='_')[-1].split(sep='.')[0]

    data_firr = pd.read_csv(output_dir / fn_firr, encoding="ISO-8859-1", header=0)
    data_noirr = pd.read_csv(output_dir / fn_noirr, encoding="ISO-8859-1", header=0)
    delta_rel_combi = pd.read_csv(output_dir / fn_noirr, encoding="ISO-8859-1", header=0)
    for cntry in delta_rel_combi.columns:
        if len(cntry) != 3:
            print()
            continue
        delta_rel_combi[cntry] = (data_noirr[cntry] / \
                                  baseline.loc[(baseline.crop==crop) & (baseline.variable==var['noirr-fullexp'])][cntry].values[0] \
                                 - data_firr[cntry] 
                                 baseline.loc[(baseline.crop==crop) & (baseline.variable==var['firr-fullexp'])][cntry].values[0]) \
                                * baseline.loc[(baseline.crop==crop) & (baseline.variable==var['combi'])][cntry].values[0]
    
    delta = data_noirr-data_firr
    delta['Year'] = data_firr['Year']
    try:
        delta = delta.drop(columns=['0'])
        delta_rel_combi = delta_rel_combi.drop(columns=['0'])
    except:
        print('')
    try:
        delta = delta.drop(columns=['Unnamed: 0'])
        delta_rel_combi = delta_rel_combi.drop(columns=['Unnamed: 0'])
    except:
        print('')
    delta.to_csv(output_dir/ fn_delta, index=0)
    delta_rel_combi.to_csv(output_dir/ fn_delta_rel, index=0)
    
    for col in delta.columns:
        if len(col) != 3:
            continue
        else:
            irr_ratio = baseline.loc[(baseline.crop==crop) & (baseline.variable=="Irrigation_ratio")][col].values[0]
            delta_rel_combi[col] = delta_rel_combi[col] * (1-irr_ratio)
            delta[col] = delta[col] * (1-irr_ratio)
    delta.to_csv(output_dir/ fn_delta_scaled, index=0)
    delta_rel_combi.to_csv(output_dir/ fn_delta_rel_scaled, index=0)


