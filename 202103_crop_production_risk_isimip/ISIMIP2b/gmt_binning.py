#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 11:05:53 2020

@author: eberenzs
"""
import os
import numpy as np
import pandas as pd
import copy
from scipy import sparse

import matplotlib.pyplot as plt

import crop_config as co


def make_delta_gmt_df(input_version=co.input_version, lower_bound_only=False):
    """make dataframe with delta GMT and GMT bin per:
        global climate model (gcm) climate scenario, and year
    """

    for g, gcm in enumerate(co.gcms[input_version]):
        # get piControl baseline temperature
        exp = 'piControl'
        yr = co.pi_control_yr[input_version][gcm]
        gmt_path = co.gmt_dir / gcm / f'tas_day_{gcm}_{exp}_r1i1p1_EWEMBI_{yr[0]}-{yr[1]}.fldmean.yearmean.txt'
        gmtpiControlmean = np.mean(np.loadtxt(gmt_path, skiprows=1, usecols=(1,)))

        # load and store temperature_change data
        for e, exp in enumerate(co.experiments[input_version]):

            try:
                yr = co.exp_yr1[input_version][exp]
                gmt_path = co.gmt_dir / gcm / f'tas_day_{gcm}_{exp}_r1i1p1_EWEMBI_{yr[0]}-{yr[1]}.fldmean.yearmean.txt'
                gmtallperiods = np.loadtxt(gmt_path, skiprows=1, usecols=(1,))
            except OSError:
                yr = co.exp_yr2[input_version][exp]
                gmt_path = co.gmt_dir / gcm / f'tas_day_{gcm}_{exp}_r1i1p1_EWEMBI_{yr[0]}-{yr[1]}.fldmean.yearmean.txt'
                gmtallperiods = np.loadtxt(gmt_path, skiprows=1, usecols=(1,))

            df_tmp = pd.DataFrame(index=np.arange(np.diff(yr)+1),
                                              columns=['gcm', 'experiment', 'year', 'dGMT'])
            df_tmp.gcm = gcm
            df_tmp.experiment = exp
            df_tmp.year = np.arange(yr[0], yr[1]+1)
            df_tmp.dGMT = gmtallperiods - gmtpiControlmean

            if not g and not e:
                temperature_change = copy.deepcopy(df_tmp)
            else:
                temperature_change = temperature_change.append(df_tmp)
    temperature_change['bin'] = pd.cut(temperature_change['dGMT'], co.gmt_bins)
    bin_count = pd.DataFrame(index=np.arange(len(temperature_change.bin.unique())+1),
                             columns=['bin', 'lower_bound', 'count'])
    temperature_change = temperature_change.reset_index(drop=True)
    temperature_change['bin_left'] = np.nan
    for idx in temperature_change.index:
        try:
            temperature_change.loc[idx, 'bin_left'] = temperature_change.loc[idx, 'bin'].left
        except AttributeError:
            temperature_change.loc[idx, 'bin_left'] = -999

    for i, gmt_bin in enumerate(np.sort(temperature_change.bin_left.unique())):
        bin_count.loc[i, 'bin'] = temperature_change.loc[temperature_change.bin_left==gmt_bin, 'bin'].unique()[0]
        bin_count.loc[i, 'lower_bound'] = gmt_bin
        bin_count.loc[i, 'count'] = temperature_change.loc[temperature_change.bin_left==gmt_bin].shape[0]
    temperature_change.gcm = temperature_change.gcm.str.lower() # make GCM names lower case
    if lower_bound_only:
        temperature_change = temperature_change.reset_index(drop=True)
        temperature_change['bin'] = copy.deepcopy(temperature_change['bin_left'])
        temperature_change = temperature_change.drop(columns='bin_left')
    return temperature_change.reset_index(drop=True), bin_count

def map_run_to_gmt_bin(gcm, scenario, year, gmt_map_df, lower_bound_only=True):
    gmt_bin = gmt_map_df.loc[(gmt_map_df.gcm==gcm) & (gmt_map_df.experiment==scenario) &
                             (gmt_map_df.year==year), 'bin'].unique()
    if gmt_bin.size > 1:
        print('Error: gmt mapping not exact, got multiple results for map_run_to_gmt_bin for:')
        print(f' gcm={gcm}, scenario={scenario}, year={year}.')
        raise ValueError
    elif gmt_bin.size == 0:
        print('Error: gmt mapping not succesfull, got no results for map_run_to_gmt_bin for:')
        print(f' gcm={gcm}, scenario={scenario}, year={year}.')
        raise ValueError
    if lower_bound_only:
        return gmt_bin[0].left
    return gmt_bin[0]

def bin_country_impact(crop_type, co2='co2', absolute_production=True,
                       baseline_exp_file_path=None):
    """
    load country impact per model run as computed with climada_wrappers.calc_country_impacts()
    and rearrange by: GMT bin, return 1 DataFrame per 'crop' type and 'co2' parameter,
    differentiating 'bin' and 'ggcm' in extra column.

    Parameters:
        crop_type (str): 'whe', 'ric', 'soy', or 'mai'

    Optional Parameters:
        co2 (str): CO2 fertilization parameter from ISIMIP, e.g. '2005co2' or 'co2' (default)
    
    Returns:
        impact_comb = DataFrame
    """
    if not baseline_exp_file_path: baseline_exp_file_path = co.impact_dir / 'baseline_exposure.csv'
    impact_comb = pd.DataFrame()
    impact_in_filenames = [f for f in os.listdir(co.impact_dir) if (os.path.isfile(co.impact_dir / f)) if 
                     f.startswith('impact_') if f'_{crop_type}.csv' in f if f'_{co2}_' in f]
    if not impact_in_filenames and co2=='co2':
        impact_in_filenames = [f for f in os.listdir(co.impact_dir) if (os.path.isfile(co.impact_dir / f)) if 
                     f.startswith('impact_') if f'_{crop_type}.csv' in f]

    if absolute_production:
        production_baseline = pd.read_csv(baseline_exp_file_path)
        production_baseline = production_baseline.fillna(0)
    # loop over impact files: load and combine in one DataFrame:
    for idf, fn in enumerate(impact_in_filenames):

        df_tmp = pd.read_csv(co.impact_dir / fn)
        if absolute_production:
        # add reference mean value since impact is deviation from this value
            for cntry in df_tmp.columns[1:]:
                df_tmp[cntry] += production_baseline.loc[(production_baseline.Crop==crop_type) &
                                                         (production_baseline['Crop Production']=='Exposure'),
                                                         cntry].values[0]
        # add model names and scenario, rearrange df:
        df_tmp['ggcm'] = fn.split('_')[1]
        df_tmp['gcm'] = fn.split('_')[2]
        df_tmp['scenario'] = fn.split('_')[-2]
        cols = df_tmp.columns.tolist()[-3:] + df_tmp.columns.tolist()[:-3]
        df_tmp = df_tmp[cols]
        if not idf:
            impact_comb = df_tmp
        else:
            impact_comb = impact_comb.append(df_tmp,
                               ignore_index=False, verify_integrity=False)
    impact_comb['GLB'] = impact_comb.iloc[:, 1:].sum(axis=1)
    gmt_map_df = make_delta_gmt_df()[0]

    impact_comb['bin'] = impact_comb.apply(lambda x: 
                            map_run_to_gmt_bin(x['gcm'], x['scenario'],x['Year'], gmt_map_df),
                            axis=1)
    cols = impact_comb.columns.tolist()[-1:] + impact_comb.columns.tolist()[:-1]
    return impact_comb[cols].reset_index(drop=True), impact_in_filenames

def plot_GMT_bin_hist(gmt_bins=None, by='year'):
    """plot histogram of bins per by. by is either 'year' or 'dGMT'"""
    gmt_map, _ = make_delta_gmt_df(lower_bound_only=True)
    if not gmt_bins and not gmt_bins == 0:
        gmt_bins = list(np.sort(gmt_map.bin.unique()))
    elif not isinstance(gmt_bins,list):
        gmt_bins = [gmt_bins]
    f_h = list()
    ax_h = list()
    for idx, _ in enumerate(gmt_bins):
        if not (isinstance(gmt_bins[idx], float) or isinstance(gmt_bins[idx], int)):
            gmt_bins[idx] = gmt_bins[idx].left
        f_h.append(plt.figure(facecolor='w'))
        ax_h.append(f_h[idx].add_subplot(1,1,1))
        gmt_map.loc[gmt_map.bin==gmt_bins[idx], by].plot.hist(ax=ax_h[idx])
        plt.xlabel(by)
        plt.title(f'Histogram GMT bin ({gmt_bins[idx]}, {gmt_bins[idx]+.5}]')
    return f_h, ax_h

def bin_imp_events(event_df, check_plot=False):
    """ cut up imp_mat and recombine by bins
    return sparse matrices and meta data per event etc."""

    gmt, bin_count = make_delta_gmt_df(lower_bound_only=True)
    ### add bin to event_name_df:
    event_df['gmt_lower'] = np.nan
    event_df = event_df.reset_index(drop=True)
    for idx, _ in enumerate(event_df.index):
        bin_lower = gmt.loc[(gmt.gcm==event_df.gcm[idx]) &
                            (gmt.experiment==event_df.experiment[idx]) &
                            (gmt.year==event_df.year[idx]), 'bin'].values
        if bin_lower.size > 1:
            print('Error bin identification is not unique.')
            raise ValueError
        event_df.loc[idx, 'gmt_lower'] = bin_lower[0]
    if check_plot:
        plt.scatter(event_df.year, event_df.gmt_lower)
    return event_df #, event_names, etc


def bin_imp_mat(impact_mat_combi, event_name_df, both_irr_id, gmt_bin):
    """ cut up imp_mat and recombine by bins
    return sparse matrices and meta data per event etc."""
# TODO: debug
    if not (isinstance(gmt_bin, float) or isinstance(gmt_bin, int)):
        gmt_bin = gmt_bin.left
    events_df = event_name_df.loc[(event_name_df.both_irr_id==both_irr_id)]
    if np.array_equal(events_df.loc[events_df.irr=='firr'].event_name, events_df.loc[events_df.irr=='noirr'].event_name):
        events_df = events_df.loc[events_df.irr=='noirr']

    # cut matrix_
    impact_mat_binned = sparse.csr_matrix(impact_mat_combi[(events_df.gmt_lower == gmt_bin).values, :])
    
    # a bit like in cw: impact_sets_or_mats[idx_list[firr_i]] = sparse.csr_matrix(impact_sets_or_mats[idx_list[firr_i]]
    #                                                                               [names_map, :]) # Sort along ax=0
    return impact_mat_binned, events_df.loc[events_df.gmt_lower == gmt_bin]
