#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 13:44:10 2020

@author: Samuel Eberenz

Config script setting parameters and variables for all main_isimip3b and all
functions called there.

Users can manipulate parameters here, for example directory paths, etc.
"""

import os
import numpy as np
from pathlib import Path

# Input data source version, i.e., ISIMIP generation:
input_version = 'ISIMIP3b' # Example: 'ISIMIP3b'

# print info on mutliple steps for debugging?
debug = True # boolean

# Path: parent directory for everything
parent_dir = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Default: Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Kcal per ton, converison factors:
KCAL_PER_TON = {'mai': 3.56e6,
                'ric': 2.80e6,
                'soy': 3.35e6,
                'whe': 3.34e6,
                }
# "caloric contents of the “as purchased” biomass (i.e. including the water content)
#of 3.56kcal/g for maize, 2.8kcal/g for rice, 3.35kcal/g for soybean and of 3.34kcal/g
# for spring and winter wheat, following FAO (2001)." ( Mueller et al., 2021)

# set parameters specific for computation on unix cluster:
if os.path.dirname(os.path.dirname(os.path.abspath(__file__))).startswith('/cluster/'):
    on_cluster = True # run on cluster?
    subfolder = False # save output in customized subfolder? if yes, provide folder name as str, else: False

    # might require manual parent_dir definition:
    #parent_dir = Path('/cluster/work/climate/eberenzs/crop_analysis')
    parent_dir = Path('/cluster/work/climate/eberenzs/climada_papers_crop_production_risk_isimip/202103_crop_production_risk_isimip/')
    ## booleans (True or False):
    make_plots = False
    save_plots = False

    init_haz = False # True # Initiate hazard sets from NetCDF?
    if init_haz: save_haz = True # save hazard sets?
    init_exp = False # initiate and save exposure set?

    run_calc_impact_sets = False # True # calculate full gridded impacts sets?
    load_impact_mats = False # load full gridded impact sets from files?
    if load_impact_mats: combine_crop_mats = False # combine crops across full gridded impact sets?
    run_calc_country_impacts = False # compute impacts per country from hazard and exposures? (required once for next steps)
    calc_country_statistics = True # compute country-level statistics?
    calc_country_statistics_rel2bin = True # compute country-level statistics relative to reference bin?
    overlapping = True # use overallping GMT bins?
    save_country_impact_binned = True # save country-level impacts binned per GMT level?
    calc_country_statistics_kcal = True # compute country-level statistics in kcal?
    save_results_table = True

else: # if run on local machine:
    on_cluster = False
    subfolder = False

    make_plots = True

    init_haz = False # True
    if init_haz: save_haz = True
    init_exp = False # initiate exposure set?

    run_calc_impact_sets = False
    load_impact_mats = False
    if load_impact_mats: combine_crop_mats = False
    run_calc_country_impacts = False # required once for next steps
    calc_country_statistics = False
    calc_country_statistics_rel2bin = False
    overlapping = True
    save_country_impact_binned = False
    calc_country_statistics_kcal = False
    save_results_table = True

print(f'on_cluster: {on_cluster}')
#----------------------------------

work_dir = parent_dir / input_version
in_dir = parent_dir / 'data' / input_version / 'input' # Path: input data directory
out_dir = parent_dir / 'data' / input_version / 'output' # Path: input data directory
res_dir = parent_dir / 'results' / input_version # Path: results directory

haz_in_dir = in_dir / 'Hazard'
# haz_in_dir: contains gridded yearly crop yield data (netCDF files) per crop-irrigation combi
# e.g., whe-firr/lpjml_gfdl-esm4_w5e5_historical_2015soc_default_yield-wwh-firr_global_annual_1850_2014.nc
exp_in_dir = in_dir / 'Exposure'
# exp_in_dir: contains gridded landuse data from ISIMIP and FAO stats
# FAOSTAT_data_production_quantity.csv; FAOSTAT_data_country_codes.csv;
# histsoc_landuse-15crops_annual_1850_2014.nc; 2015soc_landuse-15crops_annual_1850_2014.nc;
gmt_dir = in_dir / 'gmt'
# gmt_dir: contains global mean temperature (GMT) info per year sorted by climate scenario and climate model (GCM)

# Output directory for intermediate and final outputs:
histmean_dir = out_dir / 'Hist_mean'
impact_dir = out_dir / 'Impact' #/ '20201104_crop_full_04/' : customize folder
impact_dir_noirr = impact_dir / 'noirr'
impact_dir_firr = impact_dir / 'firr'
impact_dir_noirr_fullexp = impact_dir / 'noirr_fullexp'
impact_dir_firr_fullexp = impact_dir / 'firr_fullexp'
impact_dir_detr = impact_dir / 'Detrended'

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

# result sub-directories:
plot_dir = res_dir / 'plots'

if subfolder: # if customized subfolder provided above
    impact_dir = impact_dir / subfolder
    stats_dir = stats_dir / subfolder
    stats_cntry_dir = stats_cntry_dir / subfolder
    baseline_exp_file_path = out_dir / subfolder
else:
    baseline_exp_file_path = out_dir

# create directories if not yet existing
for path in [work_dir, in_dir, haz_in_dir, exp_in_dir, out_dir, histmean_dir,
             impact_dir, impact_dir_noirr, impact_dir_firr, impact_dir_noirr_fullexp,
             impact_dir_firr_fullexp, impact_dir_detr, res_dir, haz_dir,
             haz_dir_comb, exp_dir, exp_norm_dir,
             stats_dir, stats_cntry_dir, imp_mat_dir, plot_dir]:
    path.mkdir(parents=True, exist_ok=True) # make directories if they don't exist
for key in imp_mat_crop_dir:
    imp_mat_crop_dir[key].mkdir(parents=True, exist_ok=True)
for key in haz_in_sub_dir:
    haz_in_sub_dir[key].mkdir(parents=True, exist_ok=True)
# --------------------------------------------------------------------------- #

""" Hazard and exposure generation: Parameters & Variables """

# year ranges:
yearrange_his = (1850, 2014) # historical run year range (if None, set to default in hazard creation)
yearrange_mean = (1983, 2013) # reference period used for hist_mean, hazard, and exposure creation
yearrange_fao_normalization = (2008, 2018)  # The mean of FAO crop production quantity over
                                            # these years is used for normalization of the
                                            # exposure value (normalization per country and crop)
combine_subcrops = True # combine ri1+ri2=ric and swh+wwh=whe?

# default bounding box:
bbox = (-180, -85, 180, 85)
cp_unit = 't' # Crop potential (exposure) unit, either 't' or 'USD'


exp_from_isimip = False
# if exp_from_isimip==False: exposure created from MIRCA area and SPAM/Ray yield:
filename_area = 'cultivated_area_MIRCA_GGCMI.nc4'
filename_yield = 'spam_ray_yields.nc4'
# which crop index in the netcdf corresponds to which crop type:
crop_idx_area = {'mai': 1,
                 'whe': 2,
                 'soy': 3,
                 'ric': 4}
crop_idx_yield = {'mai': 1,
                  'whe': 2,
                  'soy': 4,
                  'ric': 3}
varnames_area = {'noirr': 'cultivated area rainfed',
                 'firr': 'cultivated area irrigated',
                 'combi': 'cultivated area all'}
varnames_yield = {'noirr': 'yield.rf',
                 'firr': 'yield.ir',
                 'combi': 'yield.tot'}

""" Impact calculation: Parameters & Variables """
impact_mats_filename_list = None # ['imp_full_gepic_2005soc_co2_whe_firr.npz', 'imp_full_pepic_2005soc_co2_whe_firr.npz'] # None

normalize_exposure = False
rename_event_names = True
save_combined_haz = False
save_combined_full_imp = True
if input_version == 'ISIMIP3b': # set specific parameters depending on gridded input files
    co2_fert = 'default'
    soc = '2015soc'
elif 'ISIMIP2' in input_version:
    co2_fert = 'co2' # '2005co2' or 'co2' for ISIMIP2b
    soc = '2005soc'
else:
    co2_fert = None
    soc = None
event_names_file_path = imp_mat_dir / 'event_names_all.csv'

irr_type_list = [None] # [None, 'noirr', 'firr', 'noirr_fullexp', 'firr_fullexp']

detrend_country_impacts = True # if True, detrending is done; but both are kept, detrended and not detrended imapct data
detrend_country_impacts_mean_corr = True # reset/correct the mean value per GMT bin after detrending? Default: True

""" GMT binning: Parameters & Variables """
# names of input climate models (GCMs):
gcms = {'ISIMIP2b':['GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC5'],
        'ISIMIP3b': ['GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL']}
# names of input climate scenarios / experiments:
experiments = {'ISIMIP2b': ['historical', 'rcp26', 'rcp60', 'rcp85'], # rcp85 and rcp45 can be added if needed
               'ISIMIP3b': ['historical', 'ssp126', 'ssp585']
               }
# year range for pi_control:
pi_control_yr = {'ISIMIP2b': {'GFDL-ESM2M': (1661, 2099),
                              'HadGEM2-ES': (1661, 2299),
                              'IPSL-CM5A-LR': (1661, 2299),
                              'MIROC5': (1661, 2299)},
                 'ISIMIP3b': (1601, 2100)}
# year range for pi_control:
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

## for country-level impact:
# gmt_bins = np.array([0, 1, 1.5, 2.5, 3.5, 4.5, 5.5]) # bin limits, not center points(!) np.arange(-1.5, 10., 1.)
gmt_bins = np.array([-1, 0, 1, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
gmt_bins_overlapping = [
    np.arange(-.5, 8.5, 1),
    np.arange(-.25, 8.75, 1),
    np.arange(0.0, 9., 1),
    np.arange(0.25, 9.25, 1),
    ]

## for gridded impact:
# gmt_bins_stat: provide bin center points(!) for impact matrics only, not used for countries:
gmt_bins_stat = np.array([.5, 2, 3, 4, 5]) #np.arange(0, 7.5, 1)
gmt_bins_stat_i_ref_bin = 0 # for impact matrics only, not used for countries


""" Statistics calculation: Parameters & Variables """
crop_types = ['mai', 'ric', 'soy', 'whe']
crop_types_tons = ['mai', 'ric', 'soy', 'whe', 'combi-tons'] # ['combi-tons'] #
crop_types_kcal = ['mai-kcal', 'ric-kcal', 'soy-kcal', 'whe-kcal', 'combi-kcal']
ggcms = {'ISIMIP2b': ['gepic', 'pepic', 'lpjml', 'clm45'],
         'ISIMIP3b': ['crover', 'cygma1p74', 'epic-iiasa', 'isam', 'ldndc', 'lpjml', 'pepic', 'pdssat', 'promet', 'simplace-lintul5'], # May 2021: 10 models, excluding ACEA, DSSAT-Pythia, LPJ GUESS
         'ISIMIP3b_allcrops': ['crover', 'cygma1p74', 'epic-iiasa', 'isam', 'ldndc', 'lpjml', 'pepic', 'pdssat', 'promet'], # May 2021: 9 models, excluding simplace (no rice), check others...
         #'ISIMIP3b': ['acea', 'crover', 'cygma1p74', 'epic-iiasa', 'lpjml', 'pepic', 'promet', 'simplace-lintul5'], # 8 models
         #'ISIMIP3b_allcrops': ['crover', 'cygma1p74', 'epic-iiasa', 'lpjml', 'pepic', 'promet'], # 6 models
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
#statistics_list_ref = ['mean', 0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.5]
statistics_list_ref = ['mean', ('count_rel_mean', -.1), ('count_rel_mean', -.2),
                       0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.5]

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
statistics_list_mat_ref = [ # for statistics per GRID CELL, compared to reference bin
    ('count_rel_mean', -.2),
    0.025,
    0.05,
    0.1,
    0.5,
    'mean'
    ]


""" Plotting: Parameters and color lists """
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
    colors_ggcm['ISIMIP3b'][ggcm] = '#377eb8'
#    colors_ggcm['ISIMIP3b'][ggcm] = colors8[i] # TODO, more colors?
