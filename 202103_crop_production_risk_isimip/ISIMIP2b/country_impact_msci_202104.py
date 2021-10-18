#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 14:55:48 2020
@author: carmensteinmann

Calculate the impact per country in t/y and save to a csv file. The data_path has to be adapted to  
the path containing the hazard and exposure files created by the script "Main_create_haz_exp". 
"""

import os
import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join

from climada.hazard.relative_cropyield import RelativeCropyield
import pycountry
from climada.hazard.relative_cropyield import set_multiple_rc_from_isimip
from climada.entity.exposures import Exposures
from climada.entity.exposures.crop_production import CropProduction, normalize_several_exp, init_full_exp_set_isimip
from climada.entity import ImpactFuncSet, ImpfRelativeCropyield
from climada.engine import Impact


"""Set data path and desired computations"""
data_path = '/Users/eberenzs/Documents/Projects/crop_analysis/data/ISIMIP2b'

create_haz = False
create_exp = False
compute_impact = False
calc_firr_noirr_impact = False
normalize = False   #normalize the exposure data with FAO crop production per country
baseline = False     #generate excel sheet with exposure values per country


"""Computations"""
input_haz = os.path.join(data_path, 'input', 'Hazard')
input_exp = os.path.join(data_path, 'input', 'Exposure')

path_hist_mean = data_path + 'output/Hist_mean/'
output_dir = os.path.join(data_path, 'output')

haz_path = os.path.join(output_dir, 'Hazard')

if not os.path.exists(os.path.join(data_path,'Output', 'Impact')):
    os.mkdir(os.path.join(data_path,'Output', 'Impact'))

"""Create hazard and exposure files"""    
if create_haz:
    filelist_haz = set_multiple_rc_from_isimip(input_dir=input_haz, output_dir=output_dir)
#, combine_subcrops=False
if create_exp:
    
    filename_yield = 'spam_ray_yields.nc4'
    filename_area = 'cultivated_area_MIRCA_GGCMI.nc4'

    # crop layers and variable names in default input files:
    layers_yield = {'mai': 1, 'whe': 2, 'soy': 4, 'ric': 3}
    layers_area = {'mai': 1, 'whe': 2, 'soy': 3, 'ric': 4}

    varnames_yield = {'noirr': 'yield.rf',
                      'firr': 'yield.ir',
                     'all': 'yield.tot'}
    varnames_area = {'noirr': 'cultivated area rainfed',
                     'firr': 'cultivated area irrigated',
                     'all': 'cultivated area all'}

    ### TODO: init exposure from spam/mirca

    """Calculate impact"""

if compute_impact:
    #list climate scenarios
    all_haz_files = [f for f in listdir(haz_path) if (isfile(join(haz_path, f))) if not \
                 f.startswith('.')]
    list_scenarios = list()
    for file, filename in enumerate(all_haz_files): 
        scenario = filename.split('_')[3]
        if scenario not in list_scenarios:
            list_scenarios.append(scenario)

    #initialize first hazard to assign exposures to
    haz = RelativeCropyield()
    haz.read_hdf5(os.path.join(haz_path, all_haz_files[0]))
    
    #normalize exposures
    crop_list, countries_list, ratio_list, \
    exp_firr_norm, exp_noirr_norm, fao_cp_list, exp_tot_cp_list = normalize_several_exp(input_dir=input_exp, output_dir=output_dir)
    
    # filenames exposured (not normalized, only firr):
    filenames_exp = [f for f in listdir(os.path.join(output_dir, 'Exposure')) if
                     (isfile(join(os.path.join(output_dir, 'Exposure'), f))) if not
                     f.startswith('.') if 'firr' in f]

    #Import impact function
    if_cp = ImpactFuncSet()
    if_def = ImpfRelativeCropyield()
    if_def.set_relativeyield()
    if_cp.append(if_def)
    if_cp.check()
    

    #loop over crop and climate scenario
    
    """####
    # OVERWRITE crop_list hardcode:
    crop_list = ['soy'] # TODO remove
    countries_list[0] = np.array([32, 50, 380]) # TODO remove
    countries_list[1] = np.array([32, 50, 380]) # TODO remove
    countries_list[2] = np.array([32, 50, 380]) # TODO remove
    countries_list[3] = np.array([32, 50, 380]) # TODO remove
    ####"""
    exp_per_country = np.empty([len(crop_list), len(countries_list[0])])
    

    for crop_idx, cropname in enumerate(crop_list):
        items_exp = filenames_exp[crop_idx].split('_')
        filename_noirr = items_exp[0] + '_' + items_exp[1] + '_' + items_exp[2].split('-')[0] +\
        '-' + 'noirr' + '_' + items_exp[3]

        # load exposures
        exp_firr = CropProduction()
        exp_noirr = CropProduction()
        exp_combi = CropProduction()
        if normalize: #use normalized exposure
            exp_firr = exp_firr_norm[crop_idx]
            exp_noirr = exp_noirr_norm[crop_idx]
            exp_combi.read_hdf5(os.path.join(output_dir, 'Exposure', filename_noirr))
            exp_combi.gdf['value'] = exp_firr.gdf['value'] + exp_noirr.gdf['value']
        else: #use original exposure before normalization
            exp_firr.read_hdf5(os.path.join(output_dir, 'Exposure', filenames_exp[crop_idx]))
            exp_noirr.read_hdf5(os.path.join(output_dir, 'Exposure', filename_noirr))
            exp_combi.read_hdf5(os.path.join(output_dir, 'Exposure', filename_noirr))
            exp_combi.gdf['value'] += exp_firr.gdf['value']

        exp_firr.assign_centroids(haz, threshold=20)
        exp_noirr.gdf['centr_RC'] = exp_firr.gdf.centr_RC
        exp_combi.gdf['centr_RC'] = exp_firr.gdf.centr_RC
        
        exp_noirr.check()
        exp_firr.check()
        exp_combi.check()

        for _, scenario_name in enumerate(list_scenarios):
            #generate list of all full irrigation hazards
            filenames_haz = [f for f in listdir(haz_path) if (isfile(join(haz_path, f))) if not \
                 f.startswith('.') if '-firr_' in f if "_%s-" %(cropname) in f if scenario_name in f]

            for hazfile, hazfile_name in enumerate(filenames_haz):
                print(hazfile_name)
                # read full irrigation hazard
                haz_firr = RelativeCropyield()
                haz_firr.read_hdf5(os.path.join(haz_path, hazfile_name))
                
                # read corresponding non irrigation hazard
                items_haz = hazfile_name.split('_')
                haz_noirr = RelativeCropyield()
                filename_noirr = items_haz[0]+'_'+items_haz[1]+'_'+items_haz[2]+'_'+\
                items_haz[3]+'_'+items_haz[4]+'_'+items_haz[5]+'_'+(items_haz[6]).split('-')[0]+ \
                '-noirr_'+items_haz[7]
                haz_noirr.read_hdf5(os.path.join(haz_path, filename_noirr))

                #calculate impact per country
                column_names = ['Year']
                impact_final = np.zeros([len(countries_list[crop_idx])+len(column_names), len(haz_firr.event_id)])
                impact_final[0, :] = haz_firr.event_name
                impact_final_firr = np.zeros([len(countries_list[crop_idx])+len(column_names), len(haz_firr.event_id)])
                impact_final_firr[0, :] = haz_firr.event_name
                impact_final_noirr = np.zeros([len(countries_list[crop_idx])+len(column_names), len(haz_firr.event_id)])
                impact_final_noirr[0, :] = haz_firr.event_name

                if calc_firr_noirr_impact:
                    impact_final_firr_fullexp = np.zeros([len(countries_list[crop_idx])+len(column_names), len(haz_firr.event_id)])
                    impact_final_firr_fullexp[0, :] = haz_firr.event_name
                    impact_final_noirr_fullexp = np.zeros([len(countries_list[crop_idx])+len(column_names), len(haz_firr.event_id)])
                    impact_final_noirr_fullexp[0, :] = haz_firr.event_name

                #for country_idx, country_name in enumerate((countries_list[crop_idx])[1:],1): #enumerate([152, 156])
                for country_idx, country_name in enumerate(countries_list[crop_idx]): #enumerate([152, 156])
                    try:
                        #column_names.append(iso_cntry.get(country_name).alpha3) # old climada
                        column_names.append(pycountry.countries.lookup(f'{country_name:03}').alpha_3)
                    except LookupError:
                        column_names.append(str(country_name))
                    print(column_names[-1])

                    exp_per_country[crop_idx,country_idx] = sum(exp_noirr.gdf.loc[exp_noirr.gdf.region_id == country_name].value)+\
                        sum(exp_firr.gdf.loc[exp_noirr.gdf.region_id == country_name].value)
                    
                    impact_firr = Impact()
                    exp_cntry_firr = Exposures(exp_firr.gdf.loc[exp_firr.gdf.region_id == country_name])
                    impact_firr.calc(exp_cntry_firr, if_cp, haz_firr, save_mat=True)  
    
                    impact_noirr = Impact()
                    exp_cntry_noirr = Exposures(exp_noirr.gdf.loc[exp_noirr.gdf.region_id == country_name])
                    impact_noirr.calc(exp_cntry_noirr, if_cp, haz_noirr, save_mat=True) 

                    impact_final_firr[country_idx+1, :] = np.nan_to_num(impact_firr.at_event)
                    impact_final_noirr[country_idx+1, :] = np.nan_to_num(impact_noirr.at_event)
                    impact_final[country_idx+1, :] = np.nan_to_num(impact_firr.at_event)+np.nan_to_num(impact_noirr.at_event)

                    if calc_firr_noirr_impact:
                        exp_cntry = Exposures(exp_combi.gdf.loc[exp_combi.gdf.region_id == country_name])
                        impact_firr_combiexp = Impact()
                        impact_firr_combiexp.calc(exp_cntry, if_cp, haz_firr, save_mat=True)
                        impact_final_firr_fullexp[country_idx+1, :] = np.nan_to_num(impact_firr_combiexp.at_event)

                        impact_noirr_combiexp = Impact()
                        impact_noirr_combiexp.calc(exp_cntry, if_cp, haz_noirr, save_mat=True)
                        impact_final_noirr_fullexp[country_idx+1, :] = np.nan_to_num(impact_noirr_combiexp.at_event)
                #save impact to .csv file
                impact_filename = 'impact_combi_'+items_haz[1]+'_'+items_haz[2]+'_'+scenario_name+'_'+cropname+'.csv'
                dataframe = pd.DataFrame(impact_final.T, columns=column_names)
                dataframe = dataframe.astype({"Year": int})
                dataframe.to_csv(os.path.join(output_dir, 'Impact', impact_filename))

                impact_filename = 'impact_firr_'+items_haz[1]+'_'+items_haz[2]+'_'+scenario_name+'_'+cropname+'.csv'
                dataframe = pd.DataFrame(impact_final_firr.T, columns=column_names)
                dataframe = dataframe.astype({"Year": int})
                dataframe.to_csv(os.path.join(output_dir, 'Impact', impact_filename))

                impact_filename = 'impact_noirr_'+items_haz[1]+'_'+items_haz[2]+'_'+scenario_name+'_'+cropname+'.csv'
                dataframe = pd.DataFrame(impact_final_noirr.T, columns=column_names)
                dataframe = dataframe.astype({"Year": int})
                dataframe.to_csv(os.path.join(output_dir, 'Impact', impact_filename))

                if calc_firr_noirr_impact:
                    impact_filename_firr = 'impact_firr-fullexp_'+items_haz[1]+'_'+items_haz[2]+'_'+scenario_name+'_'+cropname+'.csv'
                    dataframe = pd.DataFrame(impact_final_firr_fullexp.T, columns=column_names)
                    dataframe = dataframe.astype({"Year": int})
                    dataframe.to_csv(os.path.join(output_dir, 'Impact_irr_diff', impact_filename_firr))
                    impact_filename_noirr = 'impact_noirr-fullexp_'+items_haz[1]+'_'+items_haz[2]+'_'+scenario_name+'_'+cropname+'.csv'
                    dataframe = pd.DataFrame(impact_final_noirr_fullexp.T, columns=column_names)
                    dataframe = dataframe.astype({"Year": int})
                    dataframe.to_csv(os.path.join(output_dir, 'Impact_irr_diff', impact_filename_noirr))

"""Generate excel sheet containing exposure baseline"""
if baseline:
    croplist = list()
    rows2 = list()
    rows_per_crop = ['Production', 'Production_firr', 'Production_noirr', 'Irrigation_ratio']
    n_rows_per_crop = len(rows_per_crop)
    data = np.empty([len(crop_list)*n_rows_per_crop, len(column_names)-1])
    for crop_idx, cropname in enumerate(crop_list):
        
        
        items_exp = filenames_exp[crop_idx].split('_')
        filename_noirr = items_exp[0] + '_' + items_exp[1] + '_' + items_exp[2].split('-')[0] +\
        '-' + 'noirr' + '_' + items_exp[3]

        # load exposures
        exp_firr = CropProduction()
        exp_noirr = CropProduction()
        exp_combi = CropProduction()
        if normalize: #use normalized exposure
            exp_firr = exp_firr_norm[crop_idx]
            exp_noirr = exp_noirr_norm[crop_idx]
            exp_combi.read_hdf5(os.path.join(output_dir, 'Exposure', filename_noirr))
            exp_combi.gdf['value'] = exp_firr.gdf['value'] + exp_noirr.gdf['value']
        else: #use original exposure before normalization
            exp_firr.read_hdf5(os.path.join(output_dir, 'Exposure', filenames_exp[crop_idx]))
            exp_noirr.read_hdf5(os.path.join(output_dir, 'Exposure', filename_noirr))
            exp_combi.read_hdf5(os.path.join(output_dir, 'Exposure', filename_noirr))
            exp_combi.gdf['value'] += exp_firr.gdf['value']

        exp_firr.assign_centroids(haz, threshold=20)
        exp_noirr.gdf['centr_RC'] = exp_firr.gdf.centr_RC
        exp_combi.gdf['centr_RC'] = exp_firr.gdf.centr_RC

        exp_noirr.check()
        exp_firr.check()
        exp_combi.check()

        for country_idx, country_name in enumerate(countries_list[crop_idx]):
            exp_cntry_firr = Exposures(exp_firr.gdf.loc[exp_firr.gdf.region_id == country_name])
            exp_cntry_noirr = Exposures(exp_noirr.gdf.loc[exp_noirr.gdf.region_id == country_name])
            exp_cntry = Exposures(exp_combi.gdf.loc[exp_combi.gdf.region_id == country_name])
            data[n_rows_per_crop*crop_idx, country_idx] = exp_cntry.gdf.value.sum()
            data[n_rows_per_crop*crop_idx+1, country_idx] = exp_cntry_firr.gdf.value.sum()
            data[n_rows_per_crop*crop_idx+2, country_idx] = exp_cntry_noirr.gdf.value.sum()
            data[n_rows_per_crop*crop_idx+3, country_idx] = exp_cntry_firr.gdf.value.sum() / exp_cntry.gdf.value.sum()
        rows2.extend(rows_per_crop)
        croplist.extend([cropname]*n_rows_per_crop)
    dataframe1 = pd.DataFrame(croplist)
    dataframe2 = pd.DataFrame(rows2)
    dataframe3 = pd.DataFrame(data, columns=column_names[1:])
    result = pd.concat([dataframe1, dataframe2, dataframe3], axis=1, sort=False)
    result.to_csv(os.path.join(output_dir, 'baseline_exposure_irr.csv'))
"""
if baseline:
    cp = list()
    croplist = list()
    data = np.empty([len(crop_list), len(column_names)-1])
    for crop_idx, cropname in enumerate(crop_list):
        cp.extend(['ISIMIP'])
        croplist.extend([cropname])
        data[crop_idx,:] = (exp_tot_cp_list[crop_idx])[1:]
        # data[3*crop_idx+1,:] = (fao_cp_list[crop_idx])[1:] #print fao crop production into the same excel sheet
        # data[3*crop_idx+2,:] = (exp_per_country[crop_idx])[1:] #print the used exposure into the same excel sheet (to cross check)
        
    saveto = os.path.join(output_dir, 'baseline_exposure.csv')
    data1 = {'Crop': croplist}
    data2 = {'Crop Production': cp}
    dataframe1 = pd.DataFrame(data1)
    dataframe2 = pd.DataFrame(data2)
    dataframe3 = pd.DataFrame(data, columns=column_names[1:])
    result = pd.concat([dataframe1, dataframe2, dataframe3], axis=1, sort=False)
    result.to_csv(saveto)
"""
