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


ggcm = co.ggcms[co.input_version][1]
gcm = co.gcms[co.input_version][2].lower()
scenario = "rcp60"
crop = "ric"
country = "CHN"
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

values = impact_dict['firr-fullexp'][country] - impact_dict['noirr-fullexp'][country]

ax.plot(impact_dict[irr_key]['Year'], values, linestyle=':', color='r', label='Delta', alpha=.9, zorder=5)
ax.hlines(0, impact_dict[irr_key]['Year'].min(),  impact_dict[irr_key]['Year'].max())
ax.legend(bbox_to_anchor=(1,1), loc=2)
ax.set_ylabel('Crop production [t/y]')
plt.title("%s %s %s %s %s" %(crop, country, scenario, ggcm, gcm))
plt.show()
