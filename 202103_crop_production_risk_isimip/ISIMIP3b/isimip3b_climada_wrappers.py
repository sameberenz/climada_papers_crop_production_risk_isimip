#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 14:11:31 2020

wrapper script using CLIMADA classes and methods to calculate hazard, exposures,
and impact from gridded crop yield data. Required by main script.

@author: Samuel Eberenz
"""
import os
import numpy as np
import pandas as pd
import copy
from scipy import sparse

from climada import util
import climada.util.coordinates as u_coord
from climada.engine import Impact
from climada.hazard import Hazard
from climada.entity import ImpactFuncSet, IFRelativeCropyield
import climada.hazard.relative_cropyield as rc
import climada.entity.exposures.crop_production as cp

import isimip3b_crop_config as co

# set logging level: ERROR --> warnings and info not returned to output
util.config.setup_logging('ERROR')

def init_hazs_and_exps(input_dir=co.haz_in_dir, return_data=False,
                       init_exp=co.init_exp, init_haz=co.init_haz,
                       exp_from_isimip=co.exp_from_isimip):
    """
    init hazard and exposure set from all files provided in input folders

    Parameters
    ----------
    input_dir : Path, optional
        DESCRIPTION. The default is co.haz_in_dir.
    return_data : boolean, optional
        Return data? The default is False. (heavy on memory)
    init_exp : boolean, optional
        initiate exposures? The default is co.init_exp.
    init_haz : boolean, optional
        Initiate hazard set? The default is co.init_haz.

    Returns
    -------
    haz_filenames (list of str): list of filenames
    haz_sets (list of Hazard): list of generated output data (hazards and historical mean)
    exp_filenames(list of str): all filenames of saved exposure files
    exp_sets (list Exposures): list containing all Exposure instances

    """
    """

    Returns:
        haz_filenames (list of str): list of filenames
        haz_sets (list of Hazard): list of generated output data (hazards and historical mean)
        exp_filenames(list of str): all filenames of saved exposure files
        exp_sets (list Exposures): list containing all Exposure instances"""

    if init_haz:
        haz_filenames, haz_sets = rc.set_multiple_rc_from_isimip(input_dir=input_dir,
                                               output_dir=co.out_dir,
                                               bbox=co.bbox, isimip_run=co.input_version,
                                               yearrange_his=co.yearrange_his,
                                               yearrange_mean=co.yearrange_mean,
                                               return_data=return_data, save=co.save_haz,
                                               combine_subcrops=co.combine_subcrops)
    else:
        haz_filenames, haz_sets = (None, None)
    if init_exp and exp_from_isimip:
        print('WARNING: Exposure values depend on choice of hist_mean sets are available in %s' %(co.histmean_dir))
        exp_filenames, exp_sets = cp.init_full_exp_set_isimip(input_dir=co.exp_in_dir,
                                                     hist_mean_dir=co.histmean_dir,
                                                     output_dir=co.out_dir, bbox=co.bbox,
                                                     isimip_version=co.input_version,
                                                     unit=co.cp_unit, return_data=return_data)
    elif init_exp:
        exp_filenames, exp_sets = init_full_exp_set_from_nc(return_data=return_data)
    else:
        exp_filenames, exp_sets = (None, None)
    
    return haz_filenames, haz_sets, exp_filenames, exp_sets

def init_full_exp_set_from_nc(crop_types=co.crop_types,
                              input_dir=co.exp_in_dir,
                              output_dir=co.exp_dir,
                              return_data=False):
    """wrapper around cp.set_from_area_and_yield_nc4(), loop over all crops"""
    filename_list = list()
    output_list = list()
    for crop_type in crop_types:
        for irr in ('noirr', 'firr'):
            exp = cp.CropProduction()
            exp.set_from_area_and_yield_nc4(crop_type,
                                            co.crop_idx_yield[crop_type],
                                            co.crop_idx_area[crop_type],
                                            bbox=co.bbox,
                                            input_dir=input_dir,
                                            filename_yield=co.filename_yield,
                                            filename_area=co.filename_area,
                                            yield_var=co.varnames_yield[irr],
                                            area_var=co.varnames_area[irr])
            exp.gdf['value'] = np.nan_to_num(exp.gdf.value) # replace NaN by 0.0
            filename = ('crop_production_' + crop_type + '-'+ irr + 
                        '_spamray-mirca.hdf5')
            filename_list.append(filename)
            exp.write_hdf5(str(output_dir / filename))
            if return_data:
                output_list.append(exp)
    return filename_list, output_list

def init_full_exp_set_from_nc_USD(crop_types=co.crop_types,
                              input_dir=co.exp_in_dir,
                              output_dir=co.exp_dir,
                              return_data=False):
    """wrapper around cp.set_from_area_and_yield_nc4(), loop over all crops"""
    filename_list = list()
    output_list = list()
    for crop_type in crop_types:
        for irr in ('noirr', 'firr'):
            exp = cp.CropProduction()
            exp.set_from_area_and_yield_nc4(crop_type,
                                            co.crop_idx_yield[crop_type],
                                            co.crop_idx_area[crop_type],
                                            bbox=co.bbox,
                                            input_dir=input_dir,
                                            filename_yield=co.filename_yield,
                                            filename_area=co.filename_area,
                                            yield_var=co.varnames_yield[irr],
                                            area_var=co.varnames_area[irr])
            exp.gdf['value'] = np.nan_to_num(exp.gdf.value) # replace NaN by 0.0
            filename = ('crop_production_' + crop_type + '-'+ irr + 
                        '_spamray-mirca_USD.hdf5')
            filename_list.append(filename)
            exp.set_value_to_usd(yearrange=(2000, 2009), input_dir=co.exp_in_dir)
            exp.gdf = exp.gdf.drop(columns=['tonnes_per_year'])
            exp.write_hdf5(str(output_dir / filename))
            if return_data:
                output_list.append(exp)
    return filename_list, output_list

def import_imp_func_set():
    """imports and returns impact function set
    
    Returns:
        impf_cp (ImpactFuncSet): instance comntaining 1 impact function for crop"""
    impf_cp = ImpactFuncSet()
    impf_def = IFRelativeCropyield()
    impf_def.set_relativeyield()
    impf_cp.append(impf_def)
    impf_cp.check()
    return impf_cp

def item_list_from_haz_filenames(haz_filenames, idx_part, keep_duplicates=True,
                                 input_version=co.input_version):
    """
    
    

    Parameters
    ----------
    haz_filenames : list
        DESCRIPTION.
    idx_part : int
         for ISIMIP2b/3b: 0: type, 1: ggcm, 2: forcing, 3: scenario, 4: soc, 5: co2, 6: crop/irr.
    keep_duplicates : bool, optional
        DESCRIPTION. The default is True.
    input_version : str, optional
        DESCRIPTION. The default is co.input_version.

    Returns
    -------
    list_items : TYPE
        DESCRIPTION.
    pure_haz_files : TYPE
        DESCRIPTION.

    """
    if not haz_filenames:
        haz_filenames = os.listdir(co.haz_dir)
    
    pure_haz_files = [f for f in haz_filenames if (os.path.isfile(co.haz_dir / f)) if not \
                     f.startswith('.') if not f.startswith('hist_mean')]

    list_items = list()
    if not input_version in ('ISIMIP2b', 'ISIMIP3b'):
        print('Warning! input not ISIMIP2b/3b? check file name splitting!')
    for file, filename in enumerate(pure_haz_files): 
        item = filename.split('_')[idx_part]
        if keep_duplicates or (item not in list_items):
            list_items.append(item)
    return list_items, pure_haz_files

def calc_country_impacts(normalize=True, haz_filenames=None, irr=None, combine_exposure_irr=False,
                         out_dir=None):
    """Compute impact (= crop production) deviation by country and export to CSV.

    Optional Parameters:
        normalize (boolean): Normalize exposure values with FAO production data per country?
        haz_filenames (list of str): file names of hazard files to be considered
        irr (str): either None (default) or 'firr' or 'noirr', if only 1 irrigation type should be considered
        combine_exposure_irr (boolean): combine firr+noirr exposure for impact calculation? only relevant if irr is not None.

    Returns:
        baseline_exposure_df (DataFrame with): table of baseline exposure values per country and crop 
    """

    if out_dir is None:
        if irr == 'noirr' and combine_exposure_irr:
            out_dir = co.impact_dir_noirr_fullexp
        elif irr == 'noirr':
            out_dir = co.impact_dir_noirr
        elif irr == 'firr' and combine_exposure_irr:
            out_dir = co.impact_dir_firr_fullexp
        elif irr == 'firr':
            out_dir = co.impact_dir_firr
        else:
            out_dir = co.impact_dir

    if not haz_filenames:
        haz_filenames = [f for f in os.listdir(co.haz_dir) if (os.path.isfile(os.path.join(co.haz_dir, f))) 
                     if f.endswith('.hdf5') if f.startswith('haz_')]

    if combine_exposure_irr and irr is None:
        print("warning: combine_exposure_irr should not be True for irr==None. Double counting of exposure values.")

    #list climate scenarios
    list_scenarios, haz_filenames = item_list_from_haz_filenames(haz_filenames, 3, keep_duplicates=False)

    #initialize first hazard to assign exposures to
    haz = rc.RelativeCropyield()
    haz.read_hdf5(co.haz_dir / haz_filenames[0])

    #normalize exposures
    if normalize:
        crop_list, countries_list, ratio_list, \
        exp_firr_norm, exp_noirr_norm, fao_cp_list, exp_tot_cp_list = cp.normalize_several_exp(input_dir=co.exp_in_dir, output_dir=co.out_dir)
    else:
        crop_list, countries_list, ratio_list, \
        _, _, fao_cp_list, exp_tot_cp_list = cp.normalize_several_exp(input_dir=co.exp_in_dir, output_dir=co.out_dir)

    filenames_exp = [f for f in os.listdir(co.exp_dir) if (os.path.isfile(co.exp_dir / f)) if not
                     f.startswith('.') if 'firr' in f]

    impf_cp = import_imp_func_set()

    exp_per_country = np.empty([len(crop_list), len(countries_list[0])])
    exp_per_country_firr = np.empty([len(crop_list), len(countries_list[0])])
    exp_per_country_noirr = np.empty([len(crop_list), len(countries_list[0])])
    #loop over crop and climate scenario
    column_names = list()
    for crop_idx, cropname in enumerate(crop_list):
        if normalize:
            exp_firr = exp_firr_norm[crop_idx]
            exp_noirr = exp_noirr_norm[crop_idx]
        else:
            items_exp = filenames_exp[crop_idx].split('_')
            exp_firr = cp.CropProduction()
            exp_firr.read_hdf5(co.exp_dir / filenames_exp[crop_idx])
        
            filename_noirr = items_exp[0] + '_' + items_exp[1] + '_' + items_exp[2].split('-')[0] +\
            '-' + 'noirr' + '_' + items_exp[3]
            exp_noirr = cp.CropProduction()
            exp_noirr.read_hdf5(co.exp_dir / filename_noirr)

        exp_firr.assign_centroids(haz, threshold=20)
        exp_noirr.gdf['centr_RC'] = exp_firr.gdf.centr_RC
        if combine_exposure_irr: # sum irrigated + non irrigated exposure
            exp_combi = cp.CropProduction()
            exp_combi.read_hdf5(co.exp_dir / filenames_exp[crop_idx])
            exp_combi.gdf.value = exp_combi.gdf.value.values + exp_noirr.gdf.value.values
            exp_combi['centr_RC'] = exp_firr.gdf.centr_RC
        print(f'{cropname}: loop over climate scenarios ...')
        for scenario_name in list_scenarios:
            """loop over climate scenarios"""
            print(scenario_name)
            #generate list of all full irrigation hazards
            haz_filenames_use = [f for f in haz_filenames if (os.path.isfile(co.haz_dir / f)) if not \
                 f.startswith('.') if '-firr_' in f if not f.startswith('hist_mean') if \
                 f'_{cropname}-' in f if scenario_name in f]

            print('... loop over hazard files ...')
            for hazfile, hazfile_name in enumerate(haz_filenames_use):
                """loop over hazard files"""
                print(hazfile_name)
                #initialize full irrigation hazard
                haz_firr = rc.RelativeCropyield()
                haz_firr.read_hdf5(co.haz_dir / hazfile_name)
                
                #initialize corresponding non irrigation hazard
                items_haz = hazfile_name.split('_')
                haz_noirr = rc.RelativeCropyield()

                if co.input_version in ('ISIMIP2b', 'ISIMIP3b'):
                    filename_noirr = items_haz[0]+'_'+items_haz[1]+'_'+items_haz[2]+'_'+\
                    items_haz[3]+'_'+items_haz[4]+'_'+items_haz[5]+'_'+(items_haz[6]).split('-')[0]+ \
                    '-noirr_'+items_haz[7]

                elif co.input_version in ('ISIMIP2a'):
                    filename_noirr = items_haz[0]+'_'+items_haz[1]+'_'+items_haz[2]+'_'+\
                    items_haz[3]+'_'+items_haz[4]+'_'+(items_haz[5]).split('-')[0]+ \
                    '-noirr_'+items_haz[6]

                haz_noirr.read_hdf5(co.haz_dir / filename_noirr)

                """calculate impact:"""
                #calculate impact per country
                impact_final = np.zeros([len(countries_list[crop_idx]), len(haz_firr.event_id)])
                impact_final[0, :] = haz_firr.event_name
                column_names = ['Year']

                for country_idx, country_name in enumerate((countries_list[crop_idx])[1:],1): #enumerate([152, 156])
                    """loop over countries"""
                    try:
                        column_names.append(u_coord.country_to_iso(country_name,
                                                                   representation="alpha3"))
                    except KeyError:
                        column_names.append(str(country_name))
                    print(column_names[-1])

                    exp_per_country_noirr[crop_idx,country_idx] = np.sum(exp_noirr.gdf.loc[exp_noirr.gdf.region_id == country_name].value)
                    exp_per_country_firr[crop_idx,country_idx] = np.sum(exp_firr.gdf.loc[exp_noirr.gdf.region_id == country_name].value)
                    exp_per_country[crop_idx,country_idx] = np.sum(exp_noirr.gdf.loc[exp_noirr.gdf.region_id == country_name].value)+\
                            np.sum(exp_firr.gdf.loc[exp_noirr.gdf.region_id == country_name].value)

                    impact_firr = Impact()
                    impact_noirr = Impact()
                    if combine_exposure_irr:
                        exp_tmp = cp.CropProduction()
                        exp_tmp.gdf = exp_combi.gdf.loc[exp_firr.gdf.region_id == country_name]
                        exp_tmp.check()
                        impact_firr.calc(exp_tmp, impf_cp, haz_firr, save_mat=True)
                        exp_tmp = cp.CropProduction()
                        exp_tmp.gdf = exp_combi.gdf.loc[exp_noirr.gdf.region_id == country_name]
                        exp_tmp.check()
                        impact_noirr.calc(exp_tmp, impf_cp, haz_noirr, save_mat=True)
                    else:
                        exp_tmp = cp.CropProduction()
                        exp_tmp.gdf = exp_firr.gdf.loc[exp_firr.gdf.region_id == country_name]
                        exp_tmp.check()
                        impact_firr.calc(exp_tmp, impf_cp, haz_firr, save_mat=True)
                        exp_tmp = cp.CropProduction()
                        exp_tmp.gdf = exp_noirr.gdf.loc[exp_noirr.gdf.region_id == country_name]
                        exp_tmp.check()
                        impact_noirr.calc(exp_tmp, impf_cp, haz_noirr, save_mat=True)
                    del exp_tmp
                    if irr == 'noirr':
                        impact_final[country_idx, :] = np.nan_to_num(impact_noirr.at_event)
                    elif irr == 'firr':
                        impact_final[country_idx, :] = np.nan_to_num(impact_firr.at_event)
                    else:
                        impact_final[country_idx, :] = np.nan_to_num(impact_firr.at_event) + np.nan_to_num(impact_noirr.at_event)

                #save impact to .csv file
                if co.input_version in ['ISIMIP2b']:
                    impact_filename = 'impact_'+items_haz[1]+'_'+items_haz[2]+'_'+items_haz[5]+'_'+scenario_name+'_'+cropname+'.csv'
                else:
                    impact_filename = 'impact_'+items_haz[1]+'_'+items_haz[2]+'_'+scenario_name+'_'+cropname+'.csv'
                dataframe = pd.DataFrame(impact_final.T, columns=column_names)
                dataframe = dataframe.astype({"Year": int})
                print(f'... saving impact per country to {impact_filename} .................')
                dataframe.to_csv(out_dir / impact_filename, index=False)

    """ generate spreadsheet to compare isimip crop production, fao crop production, and normalized exposure """
    if column_names:
        num_ = 5
        cp_list = list()
        croplist = list()
        data = np.empty([num_*len(crop_list), len(column_names)-1])
        for crop_idx, cropname in enumerate(crop_list):
            cp_list.extend(['ISIMIP', 'FAO', 'Exposure', 'Exposure_firr', 'Exposure_noirr'])
            croplist.extend(num_ * [cropname])
            data[num_*crop_idx,:] = (exp_tot_cp_list[crop_idx])[1:]
            data[num_*crop_idx+1,:] = (fao_cp_list[crop_idx])[1:]
            if irr == 'noirr' and not combine_exposure_irr:
                data[num_*crop_idx+2,:] = (exp_per_country_noirr[crop_idx])[1:]
            elif irr == 'firr' and not combine_exposure_irr:
                data[num_*crop_idx+2,:] = (exp_per_country_firr[crop_idx])[1:]
            else:
                data[num_*crop_idx+2,:] = (exp_per_country[crop_idx])[1:]
            data[num_*crop_idx+3,:] = (exp_per_country_firr[crop_idx])[1:]
            data[num_*crop_idx+4,:] = (exp_per_country_noirr[crop_idx])[1:]
    
        data1 = {'Crop': croplist}
        data2 = {'Crop Production': cp_list}
        dataframe1 = pd.DataFrame(data1)
        dataframe2 = pd.DataFrame(data2)
        dataframe3 = pd.DataFrame(data, columns=column_names[1:])
        baseline_exposure_df = pd.concat([dataframe1, dataframe2, dataframe3], axis=1, sort=False)
        print('... saving baseline_exposure.csv.................\n')
        baseline_exposure_df.to_csv(co.baseline_exp_file_path / 'baseline_exposure.csv', index=False)
        baseline_exposure_df.to_csv(out_dir / 'baseline_exposure.csv', index=False)
        return out_dir, baseline_exposure_df
    return out_dir, None

def init_haz_files_df(haz_filenames, yr_exp=co.yearrange_mean):
    """
    Initiate hazard files and organize in dataframe

    Parameters
    ----------
    haz_filenames : list
        list of hazard filenames.
    yr_exp : tuple, optional
        year range of mean. The default is co.yearrange_mean.

    Returns
    -------
    haz_df : DataFrame
    """

    haz_filenames = item_list_from_haz_filenames(haz_filenames, 0, keep_duplicates=True)[1]
    haz_df = pd.DataFrame(data=haz_filenames, columns=['haz_fn'])
    haz_fn_parts = ['type', 'ggcm', 'forcing', 'scenario', 'soc', 'co2', 'crop-irr', 'end']
    part_indices = [1, 2, 3, 4, 5, 6, 7]
    for idx_part in part_indices:
        list_tmp = item_list_from_haz_filenames(haz_filenames, idx_part,
                                                keep_duplicates=True)[0]
        if haz_fn_parts[idx_part]=='crop-irr':
            haz_df['crop'] = [item.split('-')[0] for item in list_tmp]
            haz_df['irr'] = [item.split('-')[1] for item in list_tmp]
        elif haz_fn_parts[idx_part]=='end':
            haz_df['yearrange'] = [(int(item[0:4]), int(item[5:9])) for item in list_tmp]
        else:
            haz_df[haz_fn_parts[idx_part]] = list_tmp
    haz_df['exp_fn'] = None
    for idx in haz_df.index:
        if co.exp_from_isimip:
            haz_df['exp_fn'].loc[idx] = f'crop_production_{haz_df.crop[idx]}-{haz_df.irr[idx]}_spamray-mirca.hdf5'
        else:
            haz_df['exp_fn'].loc[idx] = f'crop_production_{haz_df.crop[idx]}-{haz_df.irr[idx]}_spamray-mirca.hdf5'
    return haz_df

def rename_events(event_sets, filenames, ignore_irr=True):
    """rename_events (boolean): add run specifications to event_names?
    ---> renaming required before combining hazards of same time period
    
    Parameters:
        event_sets (list of Hazard or Impact instances)
        filenames (list of str)
    
    Optional Parameters:
        ingore_irr (boolean): Default = True

    Reurns:
        filenames (list of str): truncated filenames as added to the event names
        changed (list of boolean): indicates wheather the event set has been valid and changed by position in list
    """
    if not isinstance(event_sets, list):
        event_sets = [event_sets]
    if not isinstance(filenames, list):
        filenames = [filenames]

    changed = list()
    for idx, event_set in enumerate(event_sets):
        if isinstance(event_set, Hazard) or isinstance(event_set, Impact):
            if filenames[idx][-5:]=='.hdf5':
                filenames[idx] = filenames[idx][:-14]
            if ignore_irr and filenames[idx][-5:]=='firr_':
                filenames[idx] = filenames[idx][:-6] + '_'
            elif ignore_irr and filenames[idx][-6:]=='noirr_':
                filenames[idx] = filenames[idx][:-7] + '_'
            if filenames[idx][0:3]=='haz':
                filenames[idx]  = filenames[idx][4:]
            new_names = [filenames[idx] + name[-4:] for name in event_set.event_name]
            event_set.event_name = new_names
            changed.append(True)
        else:
            changed.append(False)
    return filenames, changed

def calc_impact_sets(normalize=True, haz_filenames=None, ignore_irr=True, debug=True,
                     return_all_mats=False):
    """
    calculate full crop impact sets

    Parameters
    ----------
    normalize : bool, optional
        Normalize exposure values with FAO production data per country? The default is True.
    haz_filenames : list of str, optional
        file names of hazard files to be considered. The default is None.
    ignore_irr : bool, optional
        DESCRIPTION. The default is True.
    debug : bool, optional
        DESCRIPTION. The default is True.
    return_all_mats : bool, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    imp_filenames : list of str
        all filenames of saved impact files.
    impact_sets : list
        list containing all Impact instances.
    haz_df : DataFrame
        DESCRIPTION.
    combi_df : DataFrame
        DESCRIPTION.
    coord_exp : TYPE
        DESCRIPTION.
    event_name_df_all : DataFrame
        event names with extra data, per imp_filename, order corresponds to order of events in imp_set.
    """

    #df all hazard sets:
    haz_df = init_haz_files_df(haz_filenames, yr_exp=co.yearrange_mean)
    # list_scenarios = list(np.unique(haz_df.scenario))
    
    group_by = ['co2', 'soc', 'crop', 'irr', 'ggcm']
    # sceanrio + gcm are pooled... change here and in for loops below to analyse GCMs seperately
    haz_df['group_id'] = None
    haz_df['both_irr_id'] = None
    group_cnt = 0
    both_irr_id_cnt = 0
    for co2 in np.unique(haz_df.co2):
        for soc in np.unique(haz_df.soc):
            for crop in np.unique(haz_df.crop):
                for ggcm in np.unique(haz_df.ggcm):
                    if haz_df.loc[(haz_df.co2==co2) & (haz_df.soc==soc) & (haz_df.crop==crop) & (haz_df.ggcm==ggcm)].shape[0]:
                        haz_df.loc[(haz_df.co2==co2) & (haz_df.soc==soc) & (haz_df.crop==crop) & (haz_df.ggcm==ggcm), 'both_irr_id'] = both_irr_id_cnt
                        both_irr_id_cnt +=1
                    for irr in np.unique(haz_df.irr):
                        if haz_df.loc[(haz_df.co2==co2) & (haz_df.soc==soc) & (haz_df.crop==crop) & (haz_df.irr==irr) & (haz_df.ggcm==ggcm)].shape[0]:
                            haz_df.loc[(haz_df.co2==co2) & (haz_df.soc==soc) & (haz_df.crop==crop) & (haz_df.irr==irr) & (haz_df.ggcm==ggcm), 'group_id'] = group_cnt
                            group_cnt +=1
    # test:
    for var in group_by:
        for idg in np.sort(haz_df.group_id.unique()):
            if haz_df.loc[haz_df.group_id==idg, var].unique().size != 1:
                print(f'grouping of hazard produced error for group {idg} in category {var}')
                raise ValueError
    haz_df.loc[haz_df.group_id==0, 'ggcm'].unique()
    # normalize exposures
    if normalize:
        crop_list, countries_list, ratio_list, \
        exp_firr_norm, exp_noirr_norm, fao_cp_list, exp_tot_cp_list = cp.normalize_several_exp(input_dir=co.exp_in_dir, output_dir=co.out_dir)
    else:
        crop_list, countries_list, ratio_list, \
        _, _, fao_cp_list, exp_tot_cp_list = cp.normalize_several_exp(input_dir=co.exp_in_dir, output_dir=co.out_dir)

    impf_cp = import_imp_func_set()
    impact_sets = list()
    imp_filenames = list()

    combi_df = pd.DataFrame(index=np.sort(haz_df.group_id.unique()), columns=['ggcm', 'soc', 'co2', 'crop', 'irr', 'group_id', 'both_irr_id'])
    combi_df.group_id = combi_df.index
    # loop over hazard groups
    for idg in np.sort(haz_df.group_id.unique()):
        for col in combi_df.columns:
            if col not in ['group_id']:
                combi_df.loc[combi_df.group_id==idg, col] = haz_df.loc[haz_df.group_id==idg, col].unique()[0]
    event_name_df_all = pd.DataFrame(columns=['ggcm', 'soc', 'co2', 'crop', 'irr', 'event_id', 'event_name', 'fn_imp'])
    for idg in np.sort(haz_df.group_id.unique()):
        print('\n ')
        haz_list = list()
        # read & combine hazard sets
        for idf, filename in enumerate(list(haz_df.loc[haz_df.group_id==idg, 'haz_fn'])):
            haz_list.append(rc.RelativeCropyield())
            haz_list[idf].read_hdf5(co.haz_dir / filename)
            if co.rename_event_names:
                rename_events(haz_list[idf], filename, ignore_irr=ignore_irr) # rename event names
            if idf > 0:
                haz_list[0].append(haz_list[idf])
            if debug: print(f'>> single hazard intensity range: {haz_list[idf].intensity.min()} to {haz_list[idf].intensity.max()}')
        haz_combined = copy.deepcopy(haz_list[0])
        if debug: print(f'>> haz_combined intensity range: {haz_combined.intensity.min()} to {haz_combined.intensity.max()}')
        del haz_list
        if co.save_combined_haz:
            fn_out = 'haz_combined_%s_%s_%s_%s_%s.hdf5' %(haz_df.loc[haz_df.group_id==idg, group_by[4]].unique()[0],
                                                          haz_df.loc[haz_df.group_id==idg, group_by[1]].unique()[0],
                                                          haz_df.loc[haz_df.group_id==idg, group_by[0]].unique()[0],
                                                          haz_df.loc[haz_df.group_id==idg, group_by[2]].unique()[0],
                                                          haz_df.loc[haz_df.group_id==idg, group_by[3]].unique()[0])
                                                                                                                 
            haz_combined.write_hdf5(co.haz_dir_comb / fn_out)
        # select fitting exposures and assign to centroids:
        crop_type = haz_df.loc[haz_df.group_id==idg, 'crop'].unique()[0]
        if normalize:
            if haz_df.loc[haz_df.group_id==idg, 'irr'].unique()[0] == 'firr':
                exp = exp_firr_norm[crop_list.index(crop_type)]
            else:
                exp = exp_noirr_norm[crop_list.index(crop_type)]
        else:
            exp = cp.CropProduction()
            exp.read_hdf5(co.exp_dir / haz_df.loc[haz_df.group_id==idg, 'exp_fn'].unique()[0])
        exp.assign_centroids(haz_combined, threshold=20)
        exp.check()

        # calc impact:
        if idg==0 or return_all_mats:
            impact_sets.append(Impact())
        else:
            impact_sets[-1] = Impact()
        impact_sets[-1].calc(exp, impf_cp, haz_combined, save_mat=True)

        fn_ = '%s_%s_%s_%s_%s' %(haz_df.loc[haz_df.group_id==idg, group_by[4]].unique()[0],
                                 haz_df.loc[haz_df.group_id==idg, group_by[1]].unique()[0],
                                 haz_df.loc[haz_df.group_id==idg, group_by[0]].unique()[0],
                                 haz_df.loc[haz_df.group_id==idg, group_by[2]].unique()[0],
                                 haz_df.loc[haz_df.group_id==idg, group_by[3]].unique()[0])
        fn_imp_out = f'imp_full_{fn_}.npz'
        imp_filenames.append(fn_imp_out)

        event_name_df = pd.DataFrame(data=impact_sets[-1].event_id, columns=['event_id'])
        event_name_df['event_name'] = impact_sets[-1].event_name
        event_name_df['fn_imp'] = fn_imp_out
        for item in group_by:
            event_name_df[item] = haz_df.loc[haz_df.group_id==idg, item].unique()[0]
        event_name_df_all = event_name_df_all.append(event_name_df)
        if co.save_combined_full_imp:
            print(f'saving impact sparse matrix: {fn_imp_out} \n')
            sparse.save_npz(co.imp_mat_crop_dir[crop_type] / fn_imp_out, impact_sets[-1].imp_mat)
            
    coord_exp = pd.DataFrame(data=impact_sets[-1].coord_exp, columns=['lat', 'lon'])
    event_name_df_all = event_name_df_all.reset_index(drop=True)
    if co.save_combined_full_imp:
         event_name_df_all.to_csv(co.event_names_file_path, index=False)
         haz_df.to_csv(co.imp_mat_dir / 'haz_df.csv', index=False)
         combi_df.to_csv(co.imp_mat_dir / 'combi_df.csv', index=False)
         coord_exp.to_csv(co.imp_mat_dir / 'coord_exp.csv', index=False)

    if not return_all_mats: impact_sets=None
    return imp_filenames, impact_sets, haz_df, combi_df, coord_exp, event_name_df_all

def import_imp_mat(input_dir=co.imp_mat_dir, mat_dir=co.imp_mat_dir, filenames=None):
    """
    read impact sparse matrix from input_dir and metadata as DFs, return lists

    Parameters
    ----------
    input_dir : TYPE, optional
        DESCRIPTION. The default is co.imp_mat_dir.
    mat_dir : TYPE, optional
        DESCRIPTION. The default is co.imp_mat_dir.
    filenames : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    imp_mat_list = list()
    if not filenames:
        filenames = [f for f in os.listdir(mat_dir) if (os.path.isfile(os.path.join(mat_dir, f))) 
                     if f.endswith('.npz') if f.startswith('imp_full_')]
    elif isinstance(filenames, str):
        filenames = [filenames]
    for f in filenames:
        imp_mat_list.append(sparse.load_npz(os.path.join(mat_dir, f)))

    if os.path.isfile(os.path.join(input_dir, 'haz_df.csv')):
        haz_df = pd.read_csv(os.path.join(input_dir, 'haz_df.csv'), encoding="ISO-8859-1", header=0)
    else:
        haz_df = None

    if os.path.isfile(os.path.join(input_dir, 'combi_df.csv')):
        combi_df = pd.read_csv(os.path.join(input_dir, 'combi_df.csv'), encoding="ISO-8859-1", header=0)
    else:
        combi_df = None

    if os.path.isfile(os.path.join(input_dir, 'coord_exp.csv')):
        coord_exp = pd.read_csv(os.path.join(input_dir, 'coord_exp.csv'), encoding="ISO-8859-1", header=0)
    else:
        coord_exp = None

    event_name_df = pd.read_csv(co.event_names_file_path, encoding="ISO-8859-1", header=0)

    return filenames, imp_mat_list, haz_df, combi_df, coord_exp, event_name_df.reset_index(drop=True)

def imp_mat_add_mean_cp(imp_mat_filenames, impact_sets_or_mats, debug=co.debug):
    """
    loop over imp_sets, load CP exposure and add to imp_mat 

    Parameters
    ----------
    imp_mat_filenames : TYPE
        DESCRIPTION.
    impact_sets_or_mats : TYPE
        DESCRIPTION.
    debug : TYPE, optional
        DESCRIPTION. The default is co.debug.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    impact_sets_or_mats : TYPE
        DESCRIPTION.

    """
    if len(imp_mat_filenames) != len(impact_sets_or_mats):
        print('error: imp_mat_filenames and impact_sets_or_mats need to correspond')
        raise ValueError
    if co.normalize_exposure:
        crop_list, _, _, exp_firr_norm, exp_noirr_norm, _, _ = cp.normalize_several_exp(input_dir=co.exp_in_dir, output_dir=co.out_dir)
    for idm, imp_fn in enumerate(imp_mat_filenames):
        _, _, ggcm, soc, co2, crop, irr = imp_fn.split('_')
        irr = irr.split('.')[0]
        if co.normalize_exposure and irr=='firr':
            exp = np.nan_to_num(exp_firr_norm[crop_list.index(crop)].value)

        elif co.normalize_exposure and irr=='noirr':
            exp = np.nan_to_num(exp_noirr_norm[crop_list.index(crop)].value)
        elif co.exp_from_isimip:
            exp = cp.CropProduction()
            exp.read_hdf5(co.exp_dir / f'crop_production_{crop}-{irr}_{co.yearrange_mean[0]}-{co.yearrange_mean[1]}.hdf5')
            exp = np.nan_to_num(exp.gdf.value)
        else:
            exp = cp.CropProduction()
            exp.read_hdf5(co.exp_dir / f'crop_production_{crop}-{irr}_spamray-mirca.hdf5')
            exp = np.nan_to_num(exp.gdf.value)
        if debug: print(f'\n >> exposure range: {exp.min()} to {exp.max()}')
        if isinstance(impact_sets_or_mats[idm], Impact):
            impact_sets_or_mats[idm].imp_mat = sparse.csr_matrix(np.nan_to_num(impact_sets_or_mats[idm].imp_mat))
            if debug: print(f' >> Impact range BEFORE addition: {impact_sets_or_mats[idm].imp_mat.min()} to {impact_sets_or_mats[idm].imp_mat.max()}')
            impact_sets_or_mats[idm].imp_mat = sparse.csr_matrix(impact_sets_or_mats[idm].imp_mat + exp)
            if debug: print(f' >> Impact range AFTER addition: {impact_sets_or_mats[idm].imp_mat.min()} to {impact_sets_or_mats[idm].imp_mat.max()}\n')
        elif isinstance(impact_sets_or_mats[idm], sparse.csr_matrix):
            # add exposure array N to sparse MxN
            impact_sets_or_mats[idm] = sparse.csr_matrix(np.nan_to_num(impact_sets_or_mats[idm]))
            if debug: print(f' >> Impact range BEFORE addition: {impact_sets_or_mats[idm].min()} to {impact_sets_or_mats[idm].max()}')
            impact_sets_or_mats[idm] = sparse.csr_matrix(impact_sets_or_mats[idm] + exp)
            if debug: print(f' >> Impact range AFTER addition: {impact_sets_or_mats[idm].min()} to {impact_sets_or_mats[idm].max()}\n')
    return impact_sets_or_mats

def imp_mat_combine_irr(imp_mat_filenames, impact_sets_or_mats, combi_df, event_name_df, debug=co.debug):
    """
    combine each two impact sets with same both_irr_id by adding fitting events (i.e., same event_name)

    Parameters
    ----------
    imp_mat_filenames : TYPE
        DESCRIPTION.
    impact_sets_or_mats : TYPE
        DESCRIPTION.
    combi_df : TYPE
        DESCRIPTION.
    event_name_df : TYPE
        DESCRIPTION.
    debug : TYPE, optional
        DESCRIPTION. The default is co.debug.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    combined_imp_mats : TYPE
        DESCRIPTION.
    combi_id : TYPE
        DESCRIPTION.
    event_name_df : TYPE
        DESCRIPTION.

    """

    both_irr_id = list()
    combi_id = list()
    combined_imp_mats = list()
    if not 'both_irr_id' in event_name_df.columns: event_name_df['both_irr_id'] = np.nan
    for idm, imp_fn in enumerate(imp_mat_filenames):
        _, _, ggcm, soc, co2, crop, irr = imp_fn.split('_')
        # irr = irr.split('.')[0]
        id_max = combi_df.loc[(combi_df.ggcm==ggcm) & (combi_df.soc==soc) & (combi_df.co2==co2) & (combi_df.crop==crop), 'both_irr_id'].values.max()
        id_min = combi_df.loc[(combi_df.ggcm==ggcm) & (combi_df.soc==soc) & (combi_df.co2==co2) & (combi_df.crop==crop), 'both_irr_id'].values.min()
        if not id_max==id_min:
            print('Error: both_irr_id not unique')
            raise ValueError
        event_name_df.loc[(event_name_df.ggcm==ggcm) & (event_name_df.soc==soc) & (event_name_df.co2==co2) & (event_name_df.crop==crop), 'both_irr_id'] = id_min
        both_irr_id.append(id_min)
        if isinstance(impact_sets_or_mats[idm], Impact):
            impact_sets_or_mats[idm] = impact_sets_or_mats[idm].imp_mat
    print(both_irr_id)
    for idb in np.unique(both_irr_id):
        idx_list = list(filter(lambda x: both_irr_id[x] == idb, range(len(both_irr_id))))
        if len(idx_list) != 2:
            print('Error: idx_list should have 2 elements, not %i' %(len(idx_list)))
            raise ValueError
        names_noirr = np.array(event_name_df.loc[(event_name_df.both_irr_id==idb) & (event_name_df.irr=='noirr'), 'event_name'])
        names_firr = np.array(event_name_df.loc[(event_name_df.both_irr_id==idb) & (event_name_df.irr=='firr'), 'event_name'])
        if names_firr.size != names_noirr.size:
            print('Error: firr and noirr impact sets need to have same number of events')
            raise ValueError
        elif not np.array([x in names_noirr for x in names_firr]).all():
            print('Error: firr and noirr impact sets need to have same event_names (independent of order)')
            raise ValueError
        elif not np.array_equal(names_noirr, names_firr): # reorder one of the impact sets if not same order already
            print('\n Reordering events in firr matrix according to noirr matrix...')
            names_map = [np.where(names_firr == name)[0][0] for i, name in enumerate(names_noirr)]
            if '_firr.' in imp_mat_filenames[idx_list[0]] and not '_firr.' in imp_mat_filenames[idx_list[1]]:
                firr_i = 0
            elif '_firr.' in imp_mat_filenames[idx_list[1]] and not '_firr.' in imp_mat_filenames[idx_list[0]]:
                firr_i = 1
            else:
                print('Error: matrices with same both_irr_id should be 1 firr and 1 noirr impact set')
                raise ValueError
            impact_sets_or_mats[idx_list[firr_i]] = sparse.csr_matrix(impact_sets_or_mats[idx_list[firr_i]][names_map, :]) # Sort along ax=0
            event_name_df.loc[(event_name_df.both_irr_id==idb) & (event_name_df.irr=='firr'), 'event_name'] = names_noirr
        # make sum:
        combined_imp_mats.append(sparse.csr_matrix(impact_sets_or_mats[idx_list[0]] + impact_sets_or_mats[idx_list[1]]))
        combi_id.append(idb)
        if debug: print(f' >> Impact 1 range: {impact_sets_or_mats[idx_list[0]].min()} to {impact_sets_or_mats[idx_list[0]].max()}')
        if debug: print(f' >> Impact 2 range: {impact_sets_or_mats[idx_list[1]].min()} to {impact_sets_or_mats[idx_list[1]].max()}')
        if debug: print(f' >> Impact 1+2 range noirr+firr: {combined_imp_mats[-1].min()} to {combined_imp_mats[-1].max()}\n')

    return combined_imp_mats, combi_id, event_name_df

def event_name_df_add_columns(event_name_df):
    """
    add column to dataframe containing event names

    Parameters
    ----------
    event_name_df : DataFrame

    Returns
    -------
    event_name_df : DataFrame

    """
    event_name_df = event_name_df.reset_index(drop=True)
    for idt, item in enumerate(['ggcm', 'gcm', 'experiment', 'soc', 'co2', 'crop', 'year']):
        if not item in event_name_df.columns:
            if item=='year':
                event_name_df[item] = [int(event_name.split('_')[idt]) for event_name in event_name_df.event_name]
            else:
                event_name_df[item] = [event_name.split('_')[idt] for event_name in event_name_df.event_name]
    return event_name_df

def import_imp_mat_combine_crops(combi_df, input_dir=co.imp_mat_crop_dir['combi']):
    """
    read impact sparse matrix from input_dir and metadata as DFs, return lists

    Parameters
    ----------
    combi_df : DataFrame
        dataframe containing combined gridded impacts.
    input_dir : Path, optional
        The default is co.imp_mat_crop_dir['combi'].

    Raises
    ------
    ValueError

    Returns
    -------
    combined_imp_mats : matrix
    cc_id_list : list
    combi_df : DataFrame
    event_name_df : DataFrame
    """

    combi_df['combine_crop_id'] = int(-999)
    crop_combi_cnt = 0
    # attribute unique combine_crop_id to each combi of co2, soc, ggcm, accross crop types
    for co2 in np.unique(combi_df.co2):
        for soc in np.unique(combi_df.soc):
                for ggcm in np.unique(combi_df.ggcm):
                    if combi_df.loc[(combi_df.co2==co2) & (combi_df.soc==soc) & (combi_df.ggcm==ggcm)].shape[0]:
                        combi_df.loc[(combi_df.co2==co2) & (combi_df.soc==soc) & (combi_df.ggcm==ggcm), 'combine_crop_id'] = int(crop_combi_cnt)
                        crop_combi_cnt +=1
    count = 0
    combined_imp_mats = list()
    cc_id_list = list()
    # loop over new id (combine_crop_id), to get 1 combined imp_mat per id:
    for cc_id_i, cc_id in enumerate(combi_df.combine_crop_id.unique()):
        print(cc_id)
        imp_mat_list = list()
        combi_id = list()
        # load impact matrices combined per crop and noirr/irr, and load event_name_df
        for idb in combi_df.loc[combi_df.combine_crop_id==cc_id].both_irr_id.unique():
            crop = combi_df.loc[combi_df.both_irr_id==idb].crop.unique()
            if  len(crop)==1:
                crop = crop[0]
            else:
                print(f'import_imp_mat_combine_crops(): crop should be unique: {crop}')
                raise ValueError
            df_tmp =  pd.read_csv(input_dir / f'event_names_combi_{crop}.csv', encoding="ISO-8859-1", header=0)
            df_tmp['crop_combi_id'] = cc_id
            if count==0:
                event_name_df = df_tmp[df_tmp.both_irr_id == idb]
            else:
                event_name_df = event_name_df.append(df_tmp[df_tmp.both_irr_id == idb])
            imp_mat_list.append(sparse.load_npz(input_dir / f'imp_full_combi_id_{int(idb):02.0f}.npz'))
            count += 1
            combi_id.append(idb)
        event_name_df = event_name_df.reset_index(drop=True)
        event_name_df['event_name_no_crop'] = ''
        for idx in event_name_df.index:
            event_name_df.loc[idx, 'event_name_no_crop'] = event_name_df.event_name[idx][:-8] + event_name_df.event_name[idx][-4:]
        event_names_per_crop = list()
        for idb in combi_id:
            event_names_per_crop.append(np.array(event_name_df.loc[(event_name_df.both_irr_id==idb) & (event_name_df.irr=='noirr'), 'event_name_no_crop']))
            print(f'cc_id {cc_id}, both_irr_id {idb}: N events = {len(event_names_per_crop[-1])}')
        # for idb_i, idb in enumerate(combi_id):
        # combine matrices:
        # check whether all have the same number of events
        if np.unique([len(x) for x in event_names_per_crop]).size != 1:
            print(f'Skipping cc_id {cc_id}: impact sets of each crop within combi need to have same number of events')
            continue
        elif len(event_names_per_crop) < 2:
            print(f'Skipping cc_id {cc_id}: only one crop type available, nothing to combine')
            continue
        elif not np.array([x in event_names_per_crop[0] for x in event_names_per_crop[1]]).all():
            print('Skipping cc_id {cc_id}: impact sets need to have same event_names (independent of crop & order)')
            continue
        elif len(imp_mat_list) != len(event_names_per_crop):
            print('Error: length of imp_mat_list and event_names_per_crop does not agree')
            print(imp_mat_list)
            print(event_names_per_crop)
            raise ValueError
        # else, if all checks are ok:
        cc_id_list.append(cc_id)
        print('\n Reordering events and summing up...')
        combined_imp_mats.append(imp_mat_list[0])
        for id_nl, name_list in enumerate(event_names_per_crop):
            if id_nl==0:
                print(imp_mat_list[id_nl].shape)
                names_map = None
            else: # sort
                names_map = [np.where(name_list == name)[0][0] for i, name in enumerate(event_names_per_crop[0])]

                imp_mat_list[id_nl] = sparse.csr_matrix(imp_mat_list[id_nl][names_map, :]) # Sort along ax=0
                # sum sorted matrices:
                combined_imp_mats[-1] += imp_mat_list[id_nl]
        event_name_df = event_name_df.reset_index(drop=True)
        event_name_df[(event_name_df.irr != 'firr') & (event_name_df.both_irr_id == combi_id[0])].to_csv(input_dir / f'event_names_per_crop_cc_id_{cc_id:02.0f}.csv', index=False)
        for idb_i, idb in enumerate(combi_id):
            if idb_i > 0: # drop all event_name rows except for the first in the combination
                event_name_df = event_name_df[event_name_df.both_irr_id != idb]
    combi_df.to_csv(co.imp_mat_dir / 'combine_crop_combi_df.csv', index=False)
    return combined_imp_mats, cc_id_list, combi_df, event_name_df

def country_impact_combine_crops(imp_dir=None, crops=None, unit=None, baseline_exp_file_path=None):
    """
    Combined country-level imapcts per across crop types

    Parameters
    ----------
    imp_dir : Path, optional
        Impact directory. The default is None.
    crops : list, optional
        list of crop types. The default is None.
    unit : str, optional
        unit of crop production impact. The default is None.
    baseline_exp_file_path : Path, optional
        Path to baseline exposures file. The default is None.

    Raises
    ------
    FileExistsError

    Returns
    -------
    baseline_path : Path
        Path to baseline exposures file containibng combined exposure.
    """
    if imp_dir is None:
        imp_dir = co.impact_dir
    if unit is None:
        unit='kcal'
    if crops is None:
        crops = co.crop_types
    if not baseline_exp_file_path: baseline_exp_file_path = co.impact_dir / 'baseline_exposure.csv'
    if not os.path.isfile(baseline_exp_file_path):
        if os.path.isfile(baseline_exp_file_path / 'baseline_exposure.csv'):
            baseline_exp_file_path = baseline_exp_file_path / 'baseline_exposure.csv'
        else:
            baseline_exp_file_path = co.out_dir / 'baseline_exposure.csv'
    if not os.path.isfile(baseline_exp_file_path):
        raise FileExistsError('baseline_exposure.csv: File not found.')
    production_baseline = pd.read_csv(baseline_exp_file_path, encoding="ISO-8859-1", header=0)
    production_baseline = production_baseline.fillna(0)
    production_baseline_combi = pd.DataFrame(columns=production_baseline.columns, index=[0])
    if len(crops) == 1:
        production_baseline_combi.loc[0, 'Crop'] = f'{crops[0]}-{unit}'
    else:
        production_baseline_combi.loc[0, 'Crop'] = f'combi-{unit}'
    production_baseline_combi.loc[0, 'Crop Production'] = 'Exposure'

    for i_crop, crop in enumerate(crops):
        for cntry in production_baseline.columns[2:]:
            if i_crop == 0:
                production_baseline_combi.loc[0, cntry] = 0
            if unit=='kcal':
                factor = co.KCAL_PER_TON[crop] # cp.KCAL_PER_TON[crop]
            else:
                factor = 1

            production_baseline_combi.loc[0, cntry] += factor * production_baseline.loc[(production_baseline['Crop']==crop) & (production_baseline['Crop Production'] == 'Exposure'), cntry].values[0]
            
    if len(crops) == 1:
        baseline_path = co.out_dir / f'baseline_exposure_{crop}-{unit}.csv'
    else:
        baseline_path = co.out_dir / f'baseline_exposure_combi-{unit}.csv'
    production_baseline_combi.to_csv(baseline_path, index=False)

    for ggcm in co.ggcms[co.input_version]:
        for gcm in co.gcms[co.input_version]:
            for scenario in co.experiments[co.input_version]:
                imp_df_dict = dict()
                i_crop = 0
                for crop in crops:
                    if unit=='kcal':
                        factor = co.KCAL_PER_TON[crop] # cp.KCAL_PER_TON[crop]
                    else:
                        factor = 1
                    # impact_cygma1p74_mri-esm2-0_historical_soy.csv
                    fn = f'impact_{ggcm}_{gcm.lower()}_{scenario}_{crop}.csv'
                    if os.path.isfile(imp_dir / fn):
                        imp_df_dict[crop] = pd.read_csv(imp_dir / fn, encoding="ISO-8859-1", header=0)
                        if i_crop == 0:
                            combi_df = imp_df_dict[crop]
                            for cntry in combi_df.columns[1:]:
                                combi_df[cntry] = combi_df[cntry]*factor
                        else:
                            for cntry in combi_df.columns[1:]:
                                combi_df[cntry] = combi_df[cntry] + imp_df_dict[crop][cntry]*factor
                        i_crop += 1
                if i_crop > 0:
                    if len(crops) == 1:
                        combi_df.to_csv(imp_dir / f'impact_{ggcm}_{gcm.lower()}_{scenario}_{crop}-{unit}.csv', index=False)
                    elif i_crop == len(crops): # only if data of model combination available for all crop types:
                        combi_df.to_csv(imp_dir / f'impact_{ggcm}_{gcm.lower()}_{scenario}_combi-{unit}.csv', index=False)
                    elif i_crop < len(crops): # only if data of model combination available for all crop types:
                        combi_df.to_csv(imp_dir / f'impact_{ggcm}_{gcm.lower()}_{scenario}_combiincomplete-{unit}.csv', index=False)
    return baseline_path