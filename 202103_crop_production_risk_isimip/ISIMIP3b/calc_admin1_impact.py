#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 10:15:43 2021

@author: eberenzs
"""


import os
import numpy as np
import pandas as pd
from pathlib import Path
import pycountry

# from climada_petals.hazard.relative_cropyield import RelativeCropyield
#from climada_petals.entity.exposures.crop_production import CropProduction
from climada.hazard import Hazard
from climada.entity import ImpactFuncSet
from climada.entity.exposures import Exposures
from climada_petals.entity import ImpfRelativeCropyield
from climada.engine import Impact



import crop_helper_code as chc


### TODO: set paths and constants:
# set your input and output paths manually:
## INPUT   
# folder where hazard hdf5 files are stored:
haz_path = Path("/Users/eberenzs/Documents/Projects/climada_papers_crop_production_risk_isimip/202103_crop_production_risk_isimip/data/ISIMIP3b/output/Hazard")
# folder where exposure hdf5 files are stored:
exp_path = Path("/Users/eberenzs/Documents/Projects/climada_papers_crop_production_risk_isimip/202103_crop_production_risk_isimip/data/ISIMIP3b/output/Exposure")
## OUTPUT:
timeseries_path = Path("/Users/eberenzs/Documents/Projects/climada_papers_crop_production_risk_isimip/202103_crop_production_risk_isimip/data/ISIMIP3b/output/Timeseries")
exposure_baseline_path = Path("/Users/eberenzs/Documents/Projects/climada_papers_crop_production_risk_isimip/202103_crop_production_risk_isimip/data/ISIMIP3b/output/Exporsure_baseline_region")


### TODO: define your model run (this ca go in a loop later with each variable as list - but for now manually):
crop_type = 'mai' # crop type, either 'mai' or 'whe' for maize or wheat
ggcm = 'pepic' # global gridded crop model GGCM
gcm = 'gfdl-esm4' # global climate model GCM
scenario = 'historical' # climate scenario, either 'historical', 'ssp126', or 'ssp585'

### the following code defines filenames based on input above:
# year range is defined by scenario:
if scenario == 'historical':
    yearrange = '1850-2014'
else:
    yearrange = '2015-2100'

# filename hazard full irrigation:
fn_haz_firr = f'haz_{ggcm}_{gcm}_{scenario}_2015soc_default_{crop_type}-firr_{yearrange}.hdf5'
# example filename: haz_pepic_gfdl-esm4_historical_2015soc_default_whe-firr_1850-2014.hdf5
print(fn_haz_firr)
# filename hazard no irrigation:
fn_haz_noirr = f'haz_{ggcm}_{gcm}_{scenario}_2015soc_default_{crop_type}-noirr_{yearrange}.hdf5'

# filenames exposure:
fn_exp_firr = f'crop_production_{crop_type}-firr_spamray-mirca.hdf5'
fn_exp_noirr = f'crop_production_{crop_type}-noirr_spamray-mirca.hdf5'
###########

# define country:
# Country code as ISO3alpha, e.g., 'FRA' for France, 'DEU' for Germany, 'USA', ...
# Find all codes online at https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3
country_iso3a = 'FRA' # iso3alpha or name of country
country_iso3a = pycountry.countries.lookup(country_iso3a).alpha_3
country_iso3n = int(pycountry.countries.lookup(country_iso3a).numeric)

# get list of admin1 regions:
admin1_names = chc.get_admin1_names(country_iso3a)
print(admin1_names)

# select admin1 regions to calculate impact for, either whole list or select single one(s), e.g.:
# admin1_names_sel = admin1_names # or:
admin1_names_sel = ['Ain', 'Ardèche', 'Drôme', 'Isère', 'Loire', 'Rhône', 'Savoie', 'Haute-Savoie'] # this is an example for region Rhone-Alpes, see https://en.wikipedia.org/wiki/Rh%C3%B4ne-Alpes

#Import impact function
if_cp = ImpactFuncSet()
if_def = ImpfRelativeCropyield()
if_def.set_relativeyield()
if_cp.append(if_def)
if_cp.check()

# All calculation needs to be done once for full irrigaton (firr) and once for no irrigation (noirr) and added up at the end:
# load hazard and exposure files:
haz_firr = Hazard.from_hdf5(haz_path / fn_haz_firr)
haz_noirr = Hazard.from_hdf5(haz_path / fn_haz_noirr)

exp_firr = Exposures.from_hdf5(exp_path / fn_exp_firr)
exp_noirr = Exposures.from_hdf5(exp_path / fn_exp_noirr)
# assign centroids for impact calc:
exp_firr.assign_centroids(haz_firr, threshold=20)
exp_noirr.assign_centroids(haz_noirr, threshold=20)
exp_firr.gdf['impf_'] = 1
exp_noirr.gdf['impf_'] = 1
# crop exposures to country (region_id coresponds to country_iso3n, i.e. numerical code):
exp_firr_cntry = Exposures(exp_firr.gdf.loc[exp_firr.gdf.region_id == country_iso3n])
exp_noirr_cntry = Exposures(exp_noirr.gdf.loc[exp_noirr.gdf.region_id == country_iso3n])
exp_firr_cntry.set_geometry_points()
exp_noirr_cntry.set_geometry_points()

# make plot of country exposure:
exp_firr_cntry.plot_scatter(pop_name=False)


column_names = ['Year']
impact_per_region = np.zeros([len(admin1_names_sel)+len(column_names), len(haz_firr.event_id)])
impact_per_region[0, :] = haz_firr.event_name
exp_per_region = np.zeros([len(admin1_names_sel), 1])

for admin1_idx, admin1_name in enumerate(admin1_names_sel):
    # get shape of admin1 region:
    shape_admin1 = chc.get_admin1_shape(country_iso3a, admin1_name)
    if shape_admin1 is None:
        raise ValueError(f'invalid shape name : {admin1_name} for country {country_iso3a}')
    column_names.append(str(admin1_name))
    # crop exposures to region:
    exp_firr_admin1 = Exposures(exp_firr_cntry.gdf.loc[exp_firr_cntry.gdf.geometry.within(shape_admin1)])
    exp_noirr_admin1 = Exposures(exp_firr_cntry.gdf.loc[exp_firr_cntry.gdf.geometry.within(shape_admin1)])
    # sum baseline production (exposure) per region:
    exp_per_region[admin1_idx, 0] = sum(exp_firr_admin1.gdf.value) + sum(exp_noirr_admin1.gdf.value)

    # calculate impact:
    impact_firr = Impact()
    impact_firr.calc(exp_firr_admin1, if_cp, haz_firr, save_mat=True)
    impact_noirr = Impact()
    impact_noirr.calc(exp_noirr_admin1, if_cp, haz_noirr, save_mat=True)
    # write to output matrix:
    impact_per_region[admin1_idx+1, :] = np.nan_to_num(impact_firr.at_event)+np.nan_to_num(impact_noirr.at_event)
# transform to dataframe
dataframe_impact = pd.DataFrame(impact_per_region.T, columns=column_names)
dataframe_impact = dataframe_impact.astype({"Year": int})
dataframe_exposure = pd.DataFrame(exp_per_region.T, columns=admin1_names_sel)
# save results as from dataframe to CSV:
impact_filename = 'impact_combi_'+ggcm+'_'+gcm+'_'+scenario+'_'+crop_type+'.csv'
exp_per_region_filename = 'exposure_baseline_'+crop_type+'.csv'
dataframe_impact.to_csv(timeseries_path / impact_filename, index=0)
dataframe_exposure.to_csv(exposure_baseline_path / exp_per_region_filename, index=0)
