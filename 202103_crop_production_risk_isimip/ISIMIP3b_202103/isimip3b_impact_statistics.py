#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:44:37 2020

Functions for calculating and crop production statistics and detrending impacts.
Requires gmt_binning and crop_config; called by main script.

@author: eberenzs
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

import isimip3b_crop_config as co
import isimip3b_gmt_binning as gmt_binning

def single_stat_country(df, stat):
    """
    calculate single statistic per country from dataframe
    
    Parameters
    ----------
    df : DataFrame
        Impacts per country and event.
    stat : str
        statistic to calculate.

    Returns
    -------
    results : float or None

    """
    if stat=='mean':
        return df.mean()
    elif stat=='std':
        return df.std()
    elif isinstance(stat, float) and stat>=0 and stat<=1:
        return df.quantile(q=stat)
    elif stat=='IQR':
        return df.quantile(q=.75)-df.quantile(q=.25)
    elif isinstance(stat, tuple):
        if 'count_rel' in stat[0]:
            if isinstance(stat[1], float):
                if 'mean' in stat[0]: reference = df.mean()
                elif 'median' in stat[0]: reference = df.median()
                results = reference*0
                for idx in results.index:
                    if idx=='Year':
                        results.loc[idx] = np.nan
                    else:
                        results.loc[idx] = ((df[idx]/reference[idx])<(1+stat[1])).sum() / df[idx].size
    else: results = None
    return results

def stats_country_binned_self(impact_comb, statistics_list=None, crop=None, save=True,
                              out_dir=None):
    """Compute statistics per country and globally for a given combined country-impact dataframe
    
    Parameters:
        impact_comb (df): DataFrame with country impacts combined by GMT bin
        Statistics are relative to the same bin and GGCM

    Optional Parameters:
        crop (str): crop type, e.g., 'mai'
        save (boolean): save output?
        statistics_list(list): c.f. crop_config.statistics_list (default)
        out_dir (Path): output directory for saving statistics as CSV

    Returns:
        stats_df :  DataFrame
            containing statistics per country
    """
    if out_dir is None: out_dir = co.stats_cntry_dir
    if not statistics_list: statistics_list = co.statistics_list

    stats_df = pd.DataFrame(columns=['stat']+impact_comb.columns.tolist())
    print('GMT bins:')
    print(np.sort(impact_comb.bin.unique()))
    for gmt_bin in np.sort(impact_comb.bin.unique()): # loop over bins
        for ggcm in impact_comb.ggcm.unique(): # loop over GGCMs
            stats_df_tmp = pd.DataFrame(index=np.arange(len(statistics_list)),
                                        columns=stats_df.columns)

            # loop over statistics:
            for idx, stat in enumerate(statistics_list):
                # filter by bin, ggcm, and stat and calculate statistic
                stats_df_tmp.loc[idx] = single_stat_country(impact_comb.loc[(impact_comb.bin==gmt_bin) & (impact_comb.ggcm==ggcm)].drop(columns=['bin']), stat)
                stats_df_tmp.loc[idx, 'stat'] = stat
            stats_df_tmp['N_years'] = impact_comb.loc[(impact_comb.bin==gmt_bin) & (impact_comb.ggcm==ggcm)].shape[0]
            stats_df_tmp['bin'] = gmt_bin
            stats_df_tmp['ggcm'] = ggcm
            stats_df = stats_df.append(stats_df_tmp)
    cols = stats_df.columns.tolist()[1:4] + stats_df.columns.tolist()[-1:] + stats_df.columns.tolist()[0:1] + stats_df.columns.tolist()[6:-1]
    if save:
        stats_df[cols].reset_index(drop=True).to_csv(out_dir / f'stats_country_binned_self_{crop}.csv',
                                                     index=False)
    return stats_df[cols].reset_index(drop=True)

def stats_country_binned_rel2bin(impact_comb, reference_bin=None, statistics_list=None,
                                 save=True, crop=None, out_dir=None):
    """Compute statistics per country and globally for a given combined country-impact dataframe
    Statistics are relative to the statistics of a reference_bin of the same GGCM.
    To answer question: what is the percentile of years below a certain percentile value from reference

    Parameters:
        impact_comb (df): DataFrame with country impacts combined by GMT bin

    Optional Parameters:
        crop (str): crop type, e.g., 'mai'
        save (boolean): save output?
        reference_bin (float): reference GMT bin for statistics relative to bin
        statistics_list(list): c.f. crop_config.statistics_list (default)
        out_dir (Path): output directory for saving statistics as CSV

    Returns:
        stats_df :  DataFrame
            containing statistics per country
    """
    if out_dir is None: out_dir = co.stats_cntry_dir
    if not reference_bin: reference_bin = co.reference_bin
    if not statistics_list: statistics_list = co.statistics_list_ref
    # 1. use stats_country_binned_self to get reference statistics
    stats_self = stats_country_binned_self(impact_comb, statistics_list=statistics_list, save=False)
    stats_ref = stats_self.loc[stats_self.bin==reference_bin]
    # 2. get statistics relative to reference statistics of reference_bin and same ggcm
    stats_df = pd.DataFrame(columns=stats_self.columns)

    for gmt_bin in np.sort(impact_comb.bin.unique()): # loop over bins
        print(f'{gmt_bin} ...')
        for ggcm in impact_comb.ggcm.unique(): # loop over GGCMs
            stats_df_tmp = pd.DataFrame(index=np.arange(len(statistics_list)),
                                        columns=stats_df.columns)
        
            # loop over statistics:
            for idx, stat in enumerate(statistics_list):
                for cntry in stats_df.columns.tolist()[5:]: # loop over countries
                    stat_ref_val = stats_ref.loc[(stats_ref.ggcm==ggcm) & (stats_ref.stat==stat), cntry].values[0]
                    # filter impact df by bin, ggcm, and stat and calculate statistic:
                    df_ = impact_comb.loc[(impact_comb.bin==gmt_bin) & (impact_comb.ggcm==ggcm), cntry]
                    stats_df_tmp.loc[idx, cntry] = (df_<stat_ref_val).sum() / df_.size
                stats_df_tmp.loc[idx, 'stat'] = stat
            stats_df_tmp['N_years'] = impact_comb.loc[(impact_comb.bin==gmt_bin) & (impact_comb.ggcm==ggcm)].shape[0]
            stats_df_tmp['bin'] = gmt_bin
            stats_df_tmp['ggcm'] = ggcm
            stats_df = stats_df.append(stats_df_tmp)
    if save:
        stats_df.reset_index(drop=True).to_csv(out_dir / f'stats_country_binned_rel2bin_{10*reference_bin:02.0f}_{crop}.csv',
                                               index=False)
    return stats_df.reset_index(drop=True)


def single_stat_mat(mat, stat):
    """
    calculate single statistics per grid cell for gridded impact matrix

    Parameters
    ----------
    mat : matrix
        containing impact values.
    stat : str
        statistic to calculate, e.g., 'mean'

    Returns
    -------
    results : matrix


    """
    if stat=='mean':
        return np.array(mat.mean(axis=0))[0]
    elif stat=='std':
        return np.array(mat.todense().std(axis=0))[0]
    elif isinstance(stat, float) and stat>=0 and stat<=1:
        return np.quantile(mat.todense(), stat, axis=0)
    elif stat=='IQR':
        return np.quantile(mat.todense(), .75, axis=0)-np.quantile(mat.todense(), .25, axis=0)
    elif stat=='max==0':
        return np.array(mat.todense().max(axis=0))[0] == 0
    elif isinstance(stat, tuple):
        if 'count_rel' in stat[0]:
            if isinstance(stat[1], float):
                if 'mean' in stat[0]: reference = np.array(mat.mean(axis=0))[0]
                elif 'median' in stat[0]: reference = np.quantile(mat.todense(), .5, axis=0)
                results = np.array((mat < (reference * (1+stat[1]))).sum(axis=0)).squeeze()/mat.shape[0]
    else: results = None
    return results


def stats_mat_binned_self(imp_mat, statistics_list=None, mask_zeros=True):
    """
    wrapper calling single_stat_mat()

    Parameters
    ----------
    imp_mat : TYPE
        DESCRIPTION.
    statistics_list : TYPE, optional
        DESCRIPTION. The default is None.
    mask_zeros : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    stat_array_out_list : TYPE
        DESCRIPTION.
    statistics_list : TYPE
        DESCRIPTION.

    """
    if not statistics_list: statistics_list = co.statistics_list_mat
    if mask_zeros: mask = single_stat_mat(imp_mat, 'max==0')
    stat_array_out_list = list()
    for idx, stat in enumerate(statistics_list):
        # filter by bin, ggcm, and stat and calculate statistic
        stat_array_out_list.append(single_stat_mat(imp_mat, stat))
        if mask_zeros: stat_array_out_list[-1][mask] = np.nan
    return stat_array_out_list, statistics_list

def stats_mat_binned_rel2bin(imp_mat, imp_mat_ref, statistics_list=None, mask_zeros=True):
    """
    wrapper calling single_stat_mat() for statistics relative to reference bin

    Parameters
    ----------
    imp_mat : TYPE
        DESCRIPTION.
    imp_mat_ref : TYPE
        DESCRIPTION.
    statistics_list : TYPE, optional
        DESCRIPTION. The default is None.
    mask_zeros : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    stat_array_out_list : TYPE
        DESCRIPTION.
    statistics_list : TYPE
        DESCRIPTION.

    """
    if not statistics_list: statistics_list = co.statistics_list_mat_ref
    if mask_zeros: mask = single_stat_mat(imp_mat, 'max==0')
    stat_array_out_list = list()
    for idx, stat in enumerate(statistics_list):
        # filter by bin, ggcm, and stat and calculate statistic
        if isinstance(stat, tuple):
            if 'count_rel' in stat[0]:
                if isinstance(stat[1], float):
                    if 'mean' in stat[0]: stat_ref = 'mean'
                    elif 'median' in stat[0]: stat_ref = 0.5
            stat_array_ref = single_stat_mat(imp_mat_ref, stat_ref)*(1+stat[1])
        else:
            stat_array_ref = single_stat_mat(imp_mat_ref, stat)
        stat_array_out_list.append(
            np.array((imp_mat < stat_array_ref).sum(axis=0)).squeeze()/imp_mat.shape[0]
            )
        if mask_zeros: stat_array_out_list[-1][mask] = np.nan
    return stat_array_out_list, statistics_list

def stats_mat_binned_wrapper(event_name_df, imp_mats, combi_id, gmt_bins=None, i_ref_bin=0, crops_combi=False):
    """
    wrapper calling stats_mat_binned_rel2bin() and stats_mat_binned_self()
    

    Parameters
    ----------
    event_name_df : TYPE
        DESCRIPTION.
    imp_mats : TYPE
        DESCRIPTION.
    combi_id : TYPE
        DESCRIPTION.
    gmt_bins : TYPE, optional
        DESCRIPTION. The default is None.
    i_ref_bin : TYPE, optional
        DESCRIPTION. The default is 0.
    crops_combi : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    stats_self : TYPE
        DESCRIPTION.

    """
    if crops_combi:
        id_var = 'crop_combi_id'
    else:
        id_var = 'both_irr_id'
    
    if gmt_bins is None: gmt_bins = [0.5, 1.5]
    stats_event_count_df = pd.DataFrame(columns=['combi_id', 'ggcm', 'soc', 'co2', 'crop', 'gmt_bin', 'N'],
                                     index=np.arange(len(combi_id) * len(gmt_bins)))
    stats_self = dict()
    stats_rel = dict()
    counter = 0
    crops = list()
    for ii_, id_ in enumerate(combi_id):
        stats_self[id_] = dict()
        stats_rel[id_] = dict()
        ggcm = event_name_df.loc[event_name_df[id_var]==id_, 'ggcm'].unique()[0]
        soc = event_name_df.loc[event_name_df[id_var]==id_, 'soc'].unique()[0]
        co2 = event_name_df.loc[event_name_df[id_var]==id_, 'co2'].unique()[0]
        if crops_combi:
            crop = 'combi'
        else:
            crop = event_name_df.loc[event_name_df[id_var]==id_, 'crop'].unique()[0]
        crops.append(crop)
        imp_mats_binned_ref = gmt_binning.bin_imp_mat(imp_mats[ii_], event_name_df, id_, gmt_bins[i_ref_bin])[0]
        for gmt_bin in gmt_bins:
            stats_event_count_df.loc[counter, 'combi_id'] = id_
            stats_event_count_df.loc[counter, 'ggcm'] = ggcm
            stats_event_count_df.loc[counter, 'soc'] = soc
            stats_event_count_df.loc[counter, 'co2'] = co2
            stats_event_count_df.loc[counter, 'crop'] = crop
            stats_event_count_df.loc[counter, 'gmt_bin'] = gmt_bin

            imp_mats_binned_tmp, _ = gmt_binning.bin_imp_mat(imp_mats[ii_], event_name_df, id_, gmt_bin)
            stats_event_count_df.loc[counter, 'N'] = imp_mats_binned_tmp.shape[0]
            counter += 1
            print(f'gmt_bin {gmt_bin} +/- 0.5: N={imp_mats_binned_tmp.shape[0]}')
            if imp_mats_binned_tmp.shape[0]==0: continue # skip stats if bin is empty
            stats_self[id_][gmt_bin], statistics_list = stats_mat_binned_self(imp_mats_binned_tmp)
            stats_rel[id_][gmt_bin], statistics_list_rel = stats_mat_binned_rel2bin(imp_mats_binned_tmp, imp_mats_binned_ref)
            for ids, stats_self_arr in enumerate(stats_self[id_][gmt_bin]):
                np.savetxt(co.stats_dir / f'stats_self_{ggcm}_{soc}_{co2}_{crop}_{statistics_list[ids]}_bin_{gmt_bin:1.1f}.csv', stats_self_arr, delimiter=",")
            for ids, stats_rel_arr in enumerate(stats_rel[id_][gmt_bin]):
                np.savetxt(co.stats_dir / f'stats_rel_{ggcm}_{soc}_{co2}_{crop}_{statistics_list_rel[ids]}_bin_{gmt_bin:1.1f}_ref_{gmt_bins[i_ref_bin]:1.1f}.csv', stats_rel_arr, delimiter=",")
            
    if np.unique(crops).size==1:
        stats_event_count_df.to_csv(co.stats_dir / f'stats_event_count_df_crop_{crops[0]}.csv', index=False)
    else:
        stats_event_count_df.to_csv(co.stats_dir / f'stats_event_count_df_combi_id_{combi_id}.csv', index=False)
    return stats_self

def polynomial_detrending_impact_csv(input_dir, filename, output_dir=None, order=2,
                                     save=True, save_prefix=None):
    """Detrends the impact per country using a polynomial of a specified order (default 2) and saves
    the impact hazard to the output path.
    Method following this web article:
    https://towardsdatascience.com/removing-non-linear-trends-from-timeseries-data-b21f7567ed51
    (accessed 2020/12).
    Brute force version looping over centroids in hazard set.

    Parameters:
        input_dir (string): input directory
        filename (string): hazard file (in .hdf5 format)

    Optional Parameters:
        output_dir (string): output directory, only required if save==True
        order (int): order of the polynomial to be used; default: 2
        save (boolean): if True (default) the detrended hazard is saved to the output_dir
        save_prefix (boolean or str): define prefix for output filename.
            If None, prefix is set to '' for no prefix.

    Returns:
        impact_df (DataFrame): Impacts table detrended per country over time
    """
    if isinstance(input_dir, str):
        input_dir = Path(input_dir)
    if save:
        if output_dir is None:
            output_dir = co.impact_dir_detr
        if save_prefix is None:
            save_prefix = ''
        elif isinstance(save_prefix, str) and not (save_prefix=='') and not save_prefix[-1]=='-':
            save_prefix = save_prefix + '-'

    impact_df = pd.read_csv(input_dir / filename, encoding="ISO-8859-1", header=0)


    # use the years as the x-axis for fitting the polynom to the hazard intensity
    x_ax = list(impact_df['Year'])
    x_ax = np.reshape(x_ax, (len(x_ax), 1))

#    for centroid in range(intensity.shape[1]):
    # loop over columns (countries):
    for cntry in impact_df.columns[1:]:
        data_array = impact_df[cntry].values

        #Detrend with x-order polynomial
        poly_features = PolynomialFeatures(degree=order)
        x_poly = poly_features.fit_transform(x_ax)
        model = LinearRegression()
        model.fit(x_poly, data_array)
        trend = model.predict(x_poly)

        impact_df[cntry] = [data_array[i] - trend[i] for i in range(
            0, len(data_array[:]))]

    if save:
        impact_df.to_csv(output_dir / (save_prefix + filename),
                                              index=False)
    return impact_df

def polynomial_detrending_impact_multiple(input_dir=None, filename_prefix=None,
                                          output_dir=None, order=2,
                                          save=True, save_prefix=None):
    """
    wrapper calling polynomial_detrending_impact_csv multiple times for all files in input_dir

    Parameters
    ----------
    input_dir : Path, optional
        DESCRIPTION. The default is None.
    filename_prefix : str, optional
        DESCRIPTION. The default is None.
    output_dir : Path, optional
        DESCRIPTION. The default is None.
    order : int, optional
        DESCRIPTION. The default is 2.
    save : boolean, optional
        DESCRIPTION. The default is True.
    save_prefix : str, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    if input_dir is None:
        input_dir = co.impact_dir
    if filename_prefix is None:
        filename_prefix = 'impact_'
    filenames = [f for f in os.listdir(input_dir) if (os.path.isfile(os.path.join(input_dir, f))) 
                     if f.endswith('.csv') if f.startswith('impact_')]
    for filename in filenames:
        polynomial_detrending_impact_csv(input_dir, filename, output_dir=output_dir,
                                         order=order,
                                         save=save, save_prefix=save_prefix)
