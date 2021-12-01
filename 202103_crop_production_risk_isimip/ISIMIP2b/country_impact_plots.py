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

data_path = Path('/Users/eberenzs/Documents/Projects/crop_analysis/data/ISIMIP2b')
output_dir = data_path / 'output'/ 'Impact'
plot_dir = data_path / 'plots'
irr_strings  = {
#    'noirr': "noirr_",
#    'firr': "firr_",
    'noirr-fullexp': "noirr-everywhere_",
    'firr-fullexp': "firr-everywhere_",
    'combi': "combi_",
    'delta_rel_scaled': "delta-rel-scaled_"
    }
irr_colors = {
    'combi': "k",
    'noirr': 'r',
    'firr': "g",
    'noirr-fullexp': 'r',
    'firr-fullexp': 'b',
    'delta_rel_scaled': '#FF00FF'
    }
var = {
    'combi': "Production",
    'noirr': "Production_noirr",
    'firr': "Production_firr",
    'noirr-fullexp': 'Production_noirr-everywhere',
    'firr-fullexp': 'Production_firr-everywhere',
    'delta_rel_scaled': 'Production'
    }
labels = {
    'combi': "combi",
    'noirr': "Production_noirr",
    'firr': "Production_firr",
    'noirr-fullexp': 'no irrigation',
    'firr-fullexp': 'full irrigation',
    'delta_rel_scaled': 'water scarcity signal'
    }
styles = {
    'combi': "-",
    'noirr': '--',
    'firr': "--",
    'noirr-fullexp': ':',
    'firr-fullexp': ':',
    'delta_rel_scaled': '-'
    }

#ggcm = co.ggcms[co.input_version][1]
#gcm = co.gcms[co.input_version][2].lower()
makeplots = True
savefigures= True
scenario = "rcp60"
#crops = ["mai", "ric", "soy", "whe"]
crops = ["mai" "soy"]
#countries = ["CHN", "IND", "RUS", "USA", "FRA", 'TUR', 'ESP', 'TUN', 'ITA', 'KEN', 'ZAF']
countries = ["FRA", 'TUR']
baseline = pd.read_csv(Path(data_path) / "output" / "baseline_exposure_irr_202112.csv",
                       encoding="ISO-8859-1", header=0)

impact_dict = dict()
crop_column = list()
gcm_column = list()
ggcm_column=list()
scen_column= list()


for crop in crops:
    crop_column.append(crop)
    scen_column.append(scenario)
df = pd.DataFrame(columns =['crop', 'scenario'])
df['crop'] = crop_column
df['scenario'] = scen_column

for country in countries:
    print(country)
    country_mean = list()
    country_std = list()
    for crop in crops:
        change_rel = []
        for gcm in co.gcms[co.input_version]:
            for ggcm in co.ggcms[co.input_version]:
                try:
                    filename = 'impact_' + "combi_" + ggcm +'_'+ gcm +'_'+ scenario + '_'+ crop +'.csv'
                    xxx = pd.read_csv(output_dir / filename, encoding="ISO-8859-1", header=0)
                except FileNotFoundError:
                    continue
                if makeplots:
                    fig = plt.figure(figsize=(10,8))
                    ax = fig.add_subplot(2,1,1)
                minimum = 1e30
                for irr_key in irr_strings:
                    irr_str = irr_strings[irr_key]
                    filename = 'impact_' + irr_str + ggcm +'_'+ gcm +'_'+ scenario + '_'+ crop +'.csv'
                    impact_dict[irr_key] = pd.read_csv(output_dir / filename, encoding="ISO-8859-1", header=0)
        
                    values = impact_dict[irr_key][country] + baseline.loc[(baseline.crop==crop) & (baseline.variable==var[irr_key])][country].values[0]
                    minimum = np.min([minimum, np.min(values)])
                    if makeplots:
                        ax.plot(impact_dict[irr_key]['Year'], values, color=irr_colors[irr_key], \
                                linestyle = styles[irr_key], label=labels[irr_key], alpha=.7, zorder=5)
                    if irr_key == 'delta_rel_scaled':
                        end_of_century_mean = np.mean(values[-20:]) # 2080-2099
                        baseline_combi = baseline.loc[(baseline.crop==crop) & (baseline.variable==var['combi'])][country].values[0]
                        change_rel.append((end_of_century_mean-baseline_combi)/baseline_combi)
                        std = np.std(values[-20:]/baseline_combi)
                        if makeplots:
                            ax.hlines(baseline_combi, impact_dict['delta_rel_scaled']['Year'].min()-2,  impact_dict['delta_rel_scaled']['Year'].max()+2, color='k', alpha=.25, linestyle='-', label='baseline')
        
        #values = impact_dict['firr-fullexp'][country] - impact_dict['noirr-fullexp'][country]
        #ax.plot(impact_dict[irr_key]['Year'], values, linestyle=':', color='r', label='Delta', alpha=.9, zorder=5)
        #ax.hlines(0, impact_dict[irr_key]['Year'].min(),  impact_dict[irr_key]['Year'].max())
                if makeplots:
                    s = f' end of century signal:\n mean = {change_rel[-1]*100:+.2f}%, \u03C3 = {std*100:.2f}%'
                    ax.legend(bbox_to_anchor=(1,1), loc=2)
                    ax.set_ylabel('Crop production [t/y]')
                    ax.text(2075, 1.05*minimum, s, color=irr_colors["delta_rel_scaled"])
                    plt.title("%s %s %s %s %s" %(crop, country, scenario, ggcm, gcm))
                    plt.show()
                    if savefigures:
                        fig.savefig(plot_dir / f'signal_eoc_{crop}_{country}_{scenario}_{ggcm}_{gcm}.png', \
                                dpi=300, facecolor='w', edgecolor='w', \
                                orientation='landscape', papertype=None, format='png', \
                                transparent=True, bbox_inches=None, pad_inches=0., \
                                frameon=None, metadata=None)
        country_mean.append(np.mean(change_rel))
        country_std.append(np.std(change_rel))
    df[f'{country}'] = country_mean
    df[f'{country} std'] = country_std
#    df.to_csv(data_path / 'output' / 'signal_eoc.csv')
        
        
        

