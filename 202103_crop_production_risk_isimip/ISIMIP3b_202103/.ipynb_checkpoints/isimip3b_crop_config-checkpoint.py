#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 13:44:10 2020

@author: eberenzs

config.py
"""

import os
import numpy as np

from pathlib import Path

# Input data source version, e.g. ISIMIP generation:
input_version = 'ISIMIP3b'

# print info on mutliple steps for debugging?
debug = True

# Path: parent directory everything
parent_dir = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#data_dir = parent_dir / 'data'

if os.path.dirname(os.path.dirname(os.path.abspath(__file__))).startswith('/cluster/'):
    on_cluster = True
    subfolder = False
    parent_dir = Path('/cluster/work/climate/eberenzs/crop_analysis')
    make_plots = False
    save_plots = False

    init_haz = False # True
    if init_haz: save_haz = True
    init_exp = False # initiate and save exposure set?

    run_calc_impact_sets = False
    load_impact_mats = False
    if load_impact_mats: combine_crop_mats = True
    run_calc_country_impacts = True # required once for next steps
    calc_country_statistics = True
    calc_country_statistics_rel2bin = True
    if calc_country_statistics or calc_country_statistics_rel2bin:
        overlapping = True
    gmt_bins_stat = np.arange(0, 7.5, 1) # for impact matrics only, not used for countries
    gmt_bins_stat_i_ref_bin = 1 # for impact matrics only, not used for countries

else: # local
    on_cluster = False
    subfolder = '20201121_06'
    
    make_plots = True 

    init_haz = False # True
    if init_haz: save_haz = True
    init_exp = False # initiate exposure set?

    run_calc_impact_sets = False
    load_impact_mats = False #TODO :
    if load_impact_mats: combine_crop_mats = False #TODO :
    run_calc_country_impacts = False # required once for next steps
    calc_country_statistics = True #TODO :
    calc_country_statistics_rel2bin = True
    if calc_country_statistics or calc_country_statistics_rel2bin:
        overlapping = True

    gmt_bins_stat = np.arange(0, 6.5, 1) # for impact matrics only, not used for countries
    gmt_bins_stat_i_ref_bin = 1 # for impact matrics only, not used for countries
print(f'on_cluster: {on_cluster}')
#----------------------------------


work_dir = parent_dir / input_version
in_dir = parent_dir / 'data' / input_version / 'input' # Path: input data directory
out_dir = parent_dir / 'data' / input_version / 'output' # Path: input data directory
res_dir = parent_dir / 'results' / input_version # Path: results directory

haz_in_dir = in_dir / 'Hazard' #'Hazard_test'
exp_in_dir = in_dir / 'Exposure'
histmean_dir = out_dir / 'Hist_mean'
impact_dir = out_dir / 'Impact' #/ '20201104_crop_full_04/' : special folder

haz_dir = out_dir / 'Hazard'
haz_dir_comb = out_dir / 'Hazard_combined'
exp_dir = out_dir / 'Exposure'
exp_norm_dir = out_dir / 'Exposure_norm'
imp_mat_dir = out_dir / 'Impact_mat'
stats_dir = out_dir / 'Stats'
stats_cntry_dir = out_dir / 'Stats_countries'
imp_mat_crop_dir = {'mai': imp_mat_dir / 'mai',
                'ric': imp_mat_dir / 'ric',
                'soy': imp_mat_dir / 'soy',
                'whe': imp_mat_dir / 'whe',
                'combi': imp_mat_dir / 'combi'
                }

haz_in_sub_dir = {
                'mai-firr': haz_in_dir / 'mai-firr',
                'ric-firr': haz_in_dir / 'ric-firr',
                'soy-firr': haz_in_dir / 'soy-firr',
                'whe-firr': haz_in_dir / 'whe-firr',
                'mai-noirr': haz_in_dir / 'mai-noirr',
                'ric-noirr': haz_in_dir / 'ric-noirr',
                'soy-noirr': haz_in_dir / 'soy-noirr',
                'whe-noirr': haz_in_dir / 'whe-noirr',
                }

plot_dir = out_dir / 'Plots'

if subfolder:
    impact_dir = impact_dir / subfolder
    stats_dir = stats_dir / subfolder
    stats_cntry_dir = stats_cntry_dir / subfolder
    baseline_exp_file_path = out_dir / subfolder
else:
    baseline_exp_file_path = None

for path in [work_dir, in_dir, haz_in_dir, exp_in_dir, out_dir, histmean_dir,
             impact_dir, res_dir, haz_dir, haz_dir_comb, exp_dir, exp_norm_dir,
             stats_dir, stats_cntry_dir, imp_mat_dir, plot_dir]:
    path.mkdir(parents=True, exist_ok=True) # make directories if they don't exist
for key in imp_mat_crop_dir:
    imp_mat_crop_dir[key].mkdir(parents=True, exist_ok=True)
for key in haz_in_sub_dir:
    haz_in_sub_dir[key].mkdir(parents=True, exist_ok=True)
# --------------------------------------------------------------------------- #

"""Hazard and exposure generation"""
# booleans:
# init_haz = True # True
# if init_haz: save_haz = True
# init_exp = True # initiate exposure set?# Warning: exposures look different depending on which hist_mean sets are available in the histmean_dir

# year ranges:
yearrange_his = (1850, 2014) # historical run year range (if None, set to default in hazard creation)
yearrange_mean = (1983, 2013) # reference period for hist_mean, hazard, and exposure creation
yearrange_fao_normalization = (2008, 2018)  # The mean of FAO crop production quantity over
                                            # these years is used for normalization of the
                                            # exposure value (normalization per country and crop)
combine_subcrops = True # combine ri1+ri2=ric and swh+wwh=whe?

# default bounding box:
bbox = (-180, -85, 180, 85)
cp_unit = 't' # Crop potential (exposure) unit, either 't' or 'USD'

"""Impact calculation"""
impact_mats_filename_list = None # ['imp_full_gepic_2005soc_co2_whe_firr.npz', 'imp_full_pepic_2005soc_co2_whe_firr.npz'] # None

normalize_exposure = True
rename_event_names = True
save_combined_haz = False
save_combined_full_imp = True
if input_version == 'ISIMIP3b':
    co2_fert = 'default'
    soc = '2015soc'
elif 'ISIMIP2' in input_version:
    co2_fert = 'co2' # '2005co2' or 'co2' for ISIMIP2b
    soc = '2005soc'
else:
    co2_fert = None
    soc = None
event_names_file_path = imp_mat_dir / 'event_names_all.csv'

""" GMT binning """

gmt_dir = parent_dir / 'util' / 'gmt' / input_version
gcms = {'ISIMIP2b':['GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC5'],
        'ISIMIP3b': ['GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL']}
experiments = {'ISIMIP2b': ['historical', 'rcp26', 'rcp60', 'rcp85'], # rcp85 and rcp45 can be added if needed
               'ISIMIP3b': ['historical', 'ssp126', 'ssp585']
               }
pi_control_yr = {'ISIMIP2b': {'GFDL-ESM2M': (1661, 2099),
                              'HadGEM2-ES': (1661, 2299),
                              'IPSL-CM5A-LR': (1661, 2299),
                              'MIROC5': (1661, 2299)},
                 'ISIMIP3b': (1601, 2100)}
exp_yr1 = {'ISIMIP2b': {'historical': (1861, 2005),
                       'rcp26': (2006, 2099),
                       'rcp45': (2006, 2099),
                       'rcp60': (2006, 2099),
                       'rcp85': (2006, 2099)
                       },
           'ISIMIP3b': {'historical': (1850, 2014),
                       'ssp126': (2015, 2100),
                       'ssp585': (2015, 2100),
                       },
           }
exp_yr2 = {'ISIMIP2b': {'historical': (1861, 2005),
                       'rcp26': (2006, 2299),
                       'rcp45': (2006, 2299),
                       'rcp60': (2006, 2299),
                       'rcp85': (2006, 2299)}}
gmt_bins = np.arange(-1.5, 10., 1.) # np.arange(-2, 13, .5)
gmt_bins_overlapping = [
    np.arange(-.5, 8.5, 1),
    np.arange(-.25, 8.75, 1),
    np.arange(0.0, 9., 1),
    np.arange(0.25, 9.25, 1),
    ] # TODO: implement overlaping bins?


""" Statistics """
crop_types = ['mai', 'ric', 'soy', 'whe']
ggcms = {'ISIMIP2b': ['gepic', 'pepic', 'lpjml', 'clm45'],
         'ISIMIP3b': ['acea', 'crover', 'cygma1p74', 'lpjml', 'epic-iiasa', 'pepic', 'promet', 'simplace-lintul5'],
         # 'ISIMIP3b': ['AquaCrop-ACEA', 'CROVER', 'CYGMA1p74', 'LPJmL', 'EPIC-IIASA', 'PEPIC', 'PROMET', 'SIMPLACE-LINTUL5+'],
         }
"""
    statistics_list:
        mean, std and quantiles, IQR can be compared in absolute termns between bins (e.g. mean going up)
        count_rel... are always relative to the same bin only, e.g. irrespective of what the mean is, how often 10% below or lower (for ('count_rel_mean', -.1))
"""
reference_bin = .5
statistics_list = [ # for statistics per COUNTRY
    'mean',
    'std',
    ('count_rel_mean', -.1),
    ('count_rel_mean', -.2),
    1/100, 1/80, 1/60, 1/50, 1/40, 1/30, 1/20, 1/10, 1/5, 1/3, 1/2, 2/3, 4/5, 9/10, 19/20, 29/30, 39/40, 49/50, 59/60, 79/80, 99/100,
    'IQR'] # ['mean', 'std', ('count_rel', -.1), 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975]
statistics_list_ref = ['mean', 0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.5]

statistics_list_mat =  [ # for statistics per GRID CELL
    ('count_rel_mean', -.2),
    0.025,
    0.05,
    0.1,
    0.5,
    'mean',
    'std',
    'IQR',
    # 'max==0',
    ]
statistics_list_mat_ref = [
    ('count_rel_mean', -.2),
    0.025,
    0.05,
    0.1,
    0.5,
    'mean'
    ]

""" Plotting """
colors4 = ['#e41a1c','#377eb8','#4daf4a','#984ea3']
colors8 = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666']
colors8_heating = ['#377eb8', # blue
                   '#4daf4a', # green
                   '#984ea3', # purple
                   '#e41a1c', # red
                   '#ff7f00', # orange
                   '#f781bf', # rose
                   '#ffff33', # yellow
                   '#a65628', # brown
                   ]
colors8_heating_bins = dict()
for i, col in enumerate(colors8_heating):
    colors8_heating_bins[i] = col


colors_ggcm = {'ISIMIP2b': {'gepic': colors4[0],
                            'pepic': colors4[1],
                            'lpjml': colors4[2],
                            'clm45': colors4[3]},
               }
colors_ggcm['ISIMIP3b'] = dict()
for i, ggcm in enumerate(ggcms['ISIMIP3b']):
    colors_ggcm['ISIMIP3b'][ggcm] = colors8[i]
