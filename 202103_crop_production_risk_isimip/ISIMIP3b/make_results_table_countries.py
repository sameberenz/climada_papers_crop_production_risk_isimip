#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 16:37:07 2020

Functions for laoding and postprocessing of statistics on country level produced in main script

@author: eberenzs
"""

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

import isimip3b_crop_config as co

# from climada.entity.exposures.crop_production import KCAL_PER_TON
from isimip3b_crop_config import KCAL_PER_TON


def load_country_stat_single(fn, directory=co.stats_cntry_dir):
    if directory:
        if isinstance(directory, str): directory = Path(directory)
        stats_df = pd.read_csv(directory / fn, encoding="ISO-8859-1", header=0)
    else:
        stats_df = pd.read_csv(fn, encoding="ISO-8859-1", header=0)
    return stats_df

def load_country_stat_combi(fn_str_end='.csv',
                        fn_str_start='stats_country_binned',
                        crop=co.crop_types[-1],
                        ggcm=None,
                        stat=None,
                        ref_bin=0.5,
                        directory=None,
                       ):
    fn_list = list()
    fn_list.append(f'{fn_str_start}_self_{crop}{fn_str_end}') # self
    fn_list.append(f'{fn_str_start}_rel2bin_{10*ref_bin:02.0f}_{crop}{fn_str_end}') # rel2bin

    try:
        df_self = load_country_stat_single(fn_list[0], directory=directory)
    except FileNotFoundError:
        print(fn_list[0])
        df_self = None
    try:
        df_rel2bin = load_country_stat_single(fn_list[1], directory=directory)
    except FileNotFoundError:
        print(fn_list[1])
        df_rel2bin = None
    if stat:
        df_self = df_self.loc[df_self.stat==stat]
        df_rel2bin = df_rel2bin.loc[df_rel2bin.stat==stat]
    if ggcm:
        df_self = df_self.loc[df_self.ggcm==ggcm]
        df_rel2bin = df_rel2bin.loc[df_rel2bin.ggcm==ggcm]
    return df_self, df_rel2bin, fn_list

def load_mat_stat_single(fn, directory=co.stats_dir, coord_exp=None, replace_zeros=False, calc_inverse=True):

    if coord_exp is None: coord_exp = co.imp_mat_dir / 'coord_exp.csv'
    if not isinstance(coord_exp, pd.DataFrame):
        coord_exp = pd.read_csv(co.imp_mat_dir / 'coord_exp.csv', encoding="ISO-8859-1", header=0)
        print(coord_exp.shape)
    if directory:
        if isinstance(directory, str): directory = Path(directory)
        data = pd.read_csv(directory / fn, encoding="ISO-8859-1", header=None)
    else:
        data = pd.read_csv(fn, encoding="ISO-8859-1", header=None)
    data = np.array(data[0].values).reshape(len(coord_exp.lat.unique()), len(coord_exp.lon.unique()))
    dataset = xr.Dataset(
        {
            "statistic": (("lat", "lon"), data),
        },
        {"lat": coord_exp.lat.unique(), "lon": coord_exp.lon.unique()},
    )
    if replace_zeros:
        # replace all values equal to 0 with np.nan
        dataset = dataset.where(dataset['statistic'] != 0)
    if calc_inverse:
        dataset['inverse'] = dataset.statistic**(-1)
    return dataset

def save_results_table(ref_bin=0.5, gmt_bins=None, stats_dir=None, stats=None, baseline_exp=None,
                       res_dir=co.res_dir, res_fn='results_table_countries.csv'):
    """
    Load country level statistics and combine into a results table saved as CSV

    Parameters:
    ----------
    ref_bin : float, optional
        historical reference GMT bin. The default is 0.5.
    gmt_bins : list, optional
        list of GMT bins to extract. The default is [0.5, 2, 4].
    stats_dir : Path, optional
        input directory. The default is None.
    stats : list, optional
        list of statistics to extract. The default is
        ['mean', 'std', 0.025, "('count_rel_mean', -0.1)"].
    baseline_exp : DataFrame, optional
        baseline exposures dataframe. The default is None.
    res_dir : Path, optional
        results directory. The default is co.res_dir.
    res_fn : str, optional
        results filename. The default is 'results_table_countries.csv'.
    Returns
    -------
    DataFrame

    """
    if gmt_bins is None:
        gmt_bins = [0.5, 2, 4]# np.arange(0,7,.5) # [0.5, 2, 4] #
    if stats_dir is None:
        stats_dir = co.stats_cntry_dir / 'overlapping' / 'Detrended' #'overlapping_20210124_detr_mean'
    if baseline_exp is None:
        try:
            baseline_exp = pd.read_csv(co.impact_dir / 'baseline_exposure_irr_level.csv',
                               encoding="ISO-8859-1", header=0)
        except FileNotFoundError:
            baseline_exp = pd.read_csv(co.impact_dir / 'baseline_exposure.csv',
                               encoding="ISO-8859-1", header=0)
    if stats is None: 
        # stats = ['mean', 'std', 0.025, "('count_rel_mean', -0.1)", "('count_rel_mean', -0.2)"]
        stats = ['mean', 'std', 0.025, "('count_hist_mean', -0.1)", "('count_rel_mean', -0.1)"]
    
    pop_data_path = co.in_dir / 'worldbank_population_v2020.csv'
    pop_year = '2018'
    # loop over crops, get stat DataFrames in dicts
    df_self_dict = dict()
    df_rel_dict = dict()
    for i, crop in enumerate(co.crop_types_kcal):
        df_self_dict[crop], df_rel_dict[crop], fn_list = load_country_stat_combi(crop=crop, ggcm=None, stat=None,
                                                                ref_bin=ref_bin, directory=stats_dir,
                                                                                    )
        print(fn_list)
    try:
        df_self_dict['combi-tons'],_ , _ = load_country_stat_combi(crop='combi-tons', ggcm=None, stat=None,
                                                                ref_bin=ref_bin, directory=stats_dir,
                                                                                    )
        combi_tons = True
    except:
        combi_tons = False
    print(f'combi-tons: {combi_tons}')
    print(crop)
    print(df_self_dict)
    countries = list(df_self_dict[crop].columns[5:])
    
    print(countries)
    
    # only for statistics in abs. values (t/y)
    stats_abs = [stat for stat in list(df_self_dict[crop].stat.unique()) if not 'count_rel_mean' in stat]
    pop_df = pd.read_csv(pop_data_path, encoding="ISO-8859-1", header=0)
    print(pop_df[pop_year].sum())
    
    pop_list = list()
    for cntry in countries:
        if pop_df.loc[pop_df['Country Code']==cntry, pop_year].size == 1:
            pop_list.append(pop_df.loc[pop_df['Country Code']==cntry, pop_year].values[0])
        else:
            pop_list.append(np.nan)
    # set global GLB population:
    pop_list[-1] = pop_df[pop_year].sum()
    
    combi_kcal_glb_exp = 0
    combi_kcal_glb_firr_exp = 0
    for i_cr, crop in enumerate(co.crop_types_kcal[:]):
        print(crop)
        column_keys = ['crop', 'ref_bin', 'country', f'population {pop_year}', 'irrigation_ratio']
        column_arrays = dict()
        
        column_arrays['crop'] = crop
        column_arrays['ref_bin'] = ref_bin
        column_arrays['country'] = countries
        column_arrays[f'population {pop_year}'] = pop_list
        if not 'combi' in crop:
            combi = False
            cr = crop.split('-')[0]
            irr_glb = np.nansum(baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'Exposure_firr'), 'AFG':].values.squeeze()) / \
                np.nansum(baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'Exposure'), 'AFG':].values.squeeze())
            combi_kcal_glb_exp += KCAL_PER_TON[cr] * np.nansum(baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'Exposure'), 'AFG':].values.squeeze())
            combi_kcal_glb_firr_exp += KCAL_PER_TON[cr] * np.nansum(baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'Exposure_firr'), 'AFG':].values.squeeze())
            if 'irr_ratio' in baseline_exp['Crop Production'].values:
                column_arrays['irrigation_ratio'] = list(baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'irr_ratio'), 'AFG':].values.squeeze()) 
            else:
                column_arrays['irrigation_ratio'] = list(baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'Exposure_firr'), 'AFG':].values.squeeze() / \
                    baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'Exposure'), 'AFG':].values.squeeze())
            column_arrays['irrigation_ratio'].append(irr_glb)
        else:
            combi = True
            cr = crop
            if 'irr_ratio' in baseline_exp['Crop Production'].values:
                column_arrays['irrigation_ratio'] = list(baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'irr_ratio'), 'AFG':].values.squeeze())
            else:
                column_arrays['irrigation_ratio'] = list(baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'Exposure_firr'), 'AFG':].values.squeeze() / \
                    baseline_exp.loc[(baseline_exp['Crop'] == cr) & (baseline_exp['Crop Production'] == 'Exposure'), 'AFG':].values.squeeze())
            column_arrays['irrigation_ratio'].append(combi_kcal_glb_firr_exp / combi_kcal_glb_exp)
    
        if combi and (co.input_version=='ISIMIP3b'):
            ggcms_sel = np.sort(co.ggcms['ISIMIP3b_allcrops'])
            for iggcm, ggcm in enumerate(ggcms_sel):
                if iggcm == 0:
                    df_self_cr = df_self_dict[crop].loc[df_self_dict[crop].ggcm==ggcm]
                    df_rel_cr = df_rel_dict[crop].loc[df_rel_dict[crop].ggcm==ggcm]
                else:
                    df_self_cr = df_self_cr.append(df_self_dict[crop].loc[df_self_dict[crop].ggcm==ggcm])
                    df_rel_cr = df_rel_cr.append(df_rel_dict[crop].loc[df_rel_dict[crop].ggcm==ggcm])
        else:
            df_self_cr = df_self_dict[crop]
            df_rel_cr = df_rel_dict[crop]
        print(crop)
        print(np.sort(df_self_cr.ggcm.unique()))
        print(df_self_cr.shape)

        # general statistics:
        for stat in stats:
            stat = str(stat)
            if 'count_hist_mean' in stat:
                continue
            i_cntry_start = 3
            column_keys.append(f'ref. {stat} {ref_bin:1.1f}C [kcal p.c./y]')
            df_ref_ = df_self_cr.loc[(df_self_cr.stat==stat) & (df_self_cr.bin==ref_bin)].reset_index(drop=True).select_dtypes(np.number)
            column_arrays[column_keys[-1]] = np.nan_to_num(df_ref_.median()[i_cntry_start:].values) / np.array(pop_list)
            
            if combi_tons and combi:
                column_keys.append(f'ref. {stat} {ref_bin:1.1f}C [t/y]')
                df_ref_combi_tons = df_self_dict['combi-tons'].loc[(df_self_dict['combi-tons'].stat==stat) & (df_self_dict['combi-tons'].bin==ref_bin)].reset_index(drop=True).select_dtypes(np.number)
                column_arrays[column_keys[-1]] = np.nan_to_num(df_ref_combi_tons.median()[i_cntry_start:].values)
            else:
                column_keys.append(f'ref. {stat} {ref_bin:1.1f}C [t/y]')
                column_arrays[column_keys[-1]] = np.nan_to_num(df_ref_.median()[i_cntry_start:].values) / KCAL_PER_TON[cr]
        # PP and PR:
        for gmt_bin in gmt_bins:
            for stat in stats:
                if isinstance(stat, str) and 'count_rel_mean' in stat:
                    dev_ = float(stat.split('-')[-1][:-1])
                    column_keys.append(f'PPr (IQR) {1-dev_:1.1f}*mean {gmt_bin:1.1f}C') # -3
                    column_keys.append(f'PPr {1-dev_:1.1f}*mean {gmt_bin:1.1f}C') #-2
                    column_keys.append(f'PRr {1-dev_:1.1f}*mean {gmt_bin:1.1f}C') # -1
                    #column_keys.append(f'PR* {1-dev_:1.1f}*mean {gmt_bin:1.1f}C') #
                    i_cntry_start = 3
                    df_ref_ = df_self_cr.loc[(df_self_cr.stat==stat) & (df_self_cr.bin==ref_bin)].reset_index(drop=True).select_dtypes(np.number)
                    df_ = df_self_cr.loc[(df_self_cr.stat==stat) & (df_self_cr.bin==gmt_bin)].reset_index(drop=True).select_dtypes(np.number)
                    key_ref_pp = f'PPr {1-dev_:1.1f}*mean {ref_bin:1.1f}C'
                    key_ref_pr = f'PRr {1-dev_:1.1f}*mean {ref_bin:1.1f}C'
                elif isinstance(stat, str) and 'count_hist_mean' in stat:
                    dev_ = float(stat.split('-')[-1][:-1])
                    column_keys.append(f'PPh (IQR) {1-dev_:1.1f}*mean {gmt_bin:1.1f}C') # -3
                    column_keys.append(f'PPh {1-dev_:1.1f}*mean {gmt_bin:1.1f}C') #-2
                    column_keys.append(f'PRh {1-dev_:1.1f}*mean {gmt_bin:1.1f}C') # -1
                    #column_keys.append(f'PR* {1-dev_:1.1f}*mean {gmt_bin:1.1f}C') #
                    i_cntry_start = 3
                    df_ref_ = df_rel_cr.loc[(df_rel_cr.stat==stat.replace('_hist_', '_rel_')) & (df_rel_cr.bin==ref_bin)].reset_index(drop=True).select_dtypes(np.number)
                    df_ = df_rel_cr.loc[(df_rel_cr.stat==stat.replace('_hist_', '_rel_')) & (df_rel_cr.bin==gmt_bin)].reset_index(drop=True).select_dtypes(np.number)
                    key_ref_pp = f'PPh {1-dev_:1.1f}*mean {ref_bin:1.1f}C'
                    key_ref_pr = f'PRh {1-dev_:1.1f}*mean {ref_bin:1.1f}C'
                elif isinstance(stat, float) or isinstance(stat, int):
                    column_keys.append(f'PPh (IQR) {100*stat}th p {gmt_bin:1.1f}C')
                    column_keys.append(f'PPh {100*stat}th p {gmt_bin:1.1f}C')
                    column_keys.append(f'PRh {100*stat}th p {gmt_bin:1.1f}C')
                    #column_keys.append(f'PR* {100*stat}th p {gmt_bin:1.1f}C')
                    key_ref_pp = f'PPh {100*stat}th p {ref_bin:1.1f}C'
                    key_ref_pr = f'PRh {100*stat}th p {ref_bin:1.1f}C'
                    stat = str(stat)
                    i_cntry_start = 3
                    df_ref_ = df_rel_cr.loc[(df_rel_cr.stat==stat) & (df_rel_cr.bin==ref_bin)].reset_index(drop=True).select_dtypes(np.number)
                    df_ = df_rel_cr.loc[(df_rel_cr.stat==stat) & (df_rel_cr.bin==gmt_bin)].reset_index(drop=True).select_dtypes(np.number)
                else:
                    continue
    
                # Probability ratio: df_  / df_ref_:
                df_pr = df_.div(df_ref_)
                data_array = np.nan_to_num(df_pr.median()[i_cntry_start:].values)
                # probability of reference GMT bin:
                data_array_ref = np.nan_to_num(df_ref_.median()[i_cntry_start:].values)
                # probability of GMT bin in [%]:
                data_array_bin = np.nan_to_num(df_.median()[i_cntry_start:].values)
                # inetrquartile range of probability of GMT bin:
                data_array_bin_IQR = np.nan_to_num((df_.quantile(.75)-df_.quantile(.25))[i_cntry_start:].values)
    
    
                # handle data where PR at ref is 0 and bin is not 0:
                for i, _ in enumerate(data_array):
                    if (data_array_ref[i] == 0) and (data_array_bin[i] == 0):
                        data_array[i] = np.nan
                    elif (data_array_ref[i] == 0) and (data_array_bin[i] > 0):
                        data_array[i] = np.inf
                column_arrays[column_keys[-2]] = np.round(data_array_bin * 100, 1)
                column_arrays[column_keys[-1]] = np.round(data_array, 1)
                #column_arrays[column_keys[-1]] = np.round(data_array, 1)
                seperator = 1.
                ggcm_consistency_array = ((df_pr > seperator).sum() > .7*df_pr.shape[0])[i_cntry_start:].values + ((df_pr <= seperator).sum() > .7*df_pr.shape[0])[i_cntry_start:].values
                negative_aggreement = ((df_pr <= seperator).sum() / df_pr.shape[0])[i_cntry_start:].values
                positive_aggreement = ((df_pr > seperator).sum() / df_pr.shape[0])[i_cntry_start:].values
                column_keys.append('MA ' + column_keys[-1]) # model agreement
                column_arrays[column_keys[-1]] = np.round(100*np.maximum.reduce([negative_aggreement,positive_aggreement]), 1)
    
                data_array_bin_IQR_str = list()
                for i_cntry, (pp, pr, ma) in  enumerate(zip(column_arrays[column_keys[-3]], column_arrays[column_keys[-2]], column_arrays[column_keys[-1]])):
                    # criterion: model agreement > 70%; 
                    if ma >= 70 and gmt_bin != ref_bin:
                        # criterion: PR change corresponding to PP change in median
                        if (pr > 1) and (pp > column_arrays[key_ref_pp][i_cntry]):
                            model_agg = '*'
                        elif (pr < 1) and (pp < column_arrays[key_ref_pp][i_cntry]):
                            model_agg = '*'
                        else:
                            model_agg = ''
                    else:
                        model_agg = ''
                    data_array_bin_IQR_str.append(f'{pp:1.1f}{model_agg} ({100*data_array_bin_IQR[i_cntry]:1.1f})')
                column_arrays[column_keys[-4]] = data_array_bin_IQR_str
    # =============================================================================
    #             column_keys.append('* ' + column_keys[-1])
    #             column_arrays[column_keys[-1]] = ggcm_consistency_array
    # 
    #             column_keys.append('PA ' + column_keys[-1])
    #             column_arrays[column_keys[-1]] = positive_aggreement
    # 
    #             column_keys.append('NA ' + column_keys[-2])
    #             column_arrays[column_keys[-1]] = negative_aggreement
    # =============================================================================
                
    
        df_crop = pd.DataFrame(index=np.arange(len(countries)), columns=column_keys)

        for col in column_keys:
            try:
                df_crop[col] = column_arrays[col]
            except ValueError as err:
                print(col)
                print(err)
        if i_cr==0:
            df = df_crop
        else:
            df = df.append(df_crop)
            df = df.reset_index(drop=True)
    df.to_csv(res_dir / res_fn, index=False)
    return df, baseline_exp
