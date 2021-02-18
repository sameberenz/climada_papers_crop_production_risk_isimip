#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 11:05:53 2020

funtions for binning impacts (gridded and per country) by global mean temperature (GMT) bin

@author: eberenzs
"""
import os
import numpy as np
import pandas as pd
import copy
from scipy import sparse

import matplotlib.pyplot as plt

import isimip3b_crop_config as co


def make_delta_gmt_df(input_version=co.input_version, lower_bound_only=False, bin_center_only=False,
                      gmt_bins=None):
    """
    make dataframe with delta GMT and GMT bin per global climate model (gcm) climate scenario, and year

    Parameters
    ----------
    input_version : str, optional
        DESCRIPTION. The default is co.input_version.
    lower_bound_only : bool, optional
        DESCRIPTION. The default is False.
    bin_center_only : bool, optional
        DESCRIPTION. The default is False.
    gmt_bins : list or array, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    temperature_change : DataFrame
    bin_count : DataFrame
    """

    if gmt_bins is None: gmt_bins = co.gmt_bins
    for g, gcm in enumerate(co.gcms[input_version]):
        # get piControl baseline mean temperature
        exp = 'piControl'
        if input_version == 'ISIMIP2b':
            yr = co.pi_control_yr[input_version][gcm]
            gmt_path = co.gmt_dir / gcm / f'tas_day_{gcm}_{exp}_r1i1p1_EWEMBI_{yr[0]}-{yr[1]}.fldmean.yearmean.txt'
        elif input_version == 'ISIMIP3b':
            yr = co.pi_control_yr[input_version]
            gmt_path = co.gmt_dir / exp / gcm
            gmt_path = gmt_path / os.listdir(gmt_path)[0]
        gmtpiControlmean = np.mean(np.loadtxt(gmt_path, skiprows=1, usecols=(1,)))

        # load and store temperature_change data
        for e, exp in enumerate(co.experiments[input_version]):
            if input_version == 'ISIMIP2b':
                try:
                    yr = co.exp_yr1[input_version][exp]
                    gmt_path = co.gmt_dir / gcm / f'tas_day_{gcm}_{exp}_r1i1p1_EWEMBI_{yr[0]}-{yr[1]}.fldmean.yearmean.txt'
                    gmtallperiods = np.loadtxt(gmt_path, skiprows=1, usecols=(1,))
                except OSError:
                    yr = co.exp_yr2[input_version][exp]
                    gmt_path = co.gmt_dir / gcm / f'tas_day_{gcm}_{exp}_r1i1p1_EWEMBI_{yr[0]}-{yr[1]}.fldmean.yearmean.txt'
                    gmtallperiods = np.loadtxt(gmt_path, skiprows=1, usecols=(1,))
            elif input_version == 'ISIMIP3b':
                yr = co.exp_yr1[input_version][exp]
                gmt_path = co.gmt_dir / exp / gcm
                gmt_path = gmt_path / os.listdir(gmt_path)[0]
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
    temperature_change['bin'] = pd.cut(temperature_change['dGMT'], gmt_bins)
    bin_count = pd.DataFrame(index=np.arange(len(temperature_change.bin.unique())+1),
                             columns=['bin', 'lower_bound', 'bin_center', 'count'])
    temperature_change = temperature_change.reset_index(drop=True)
    temperature_change['bin_left'] = np.nan
    temperature_change['bin_center'] = np.nan
    for idx in temperature_change.index: # bin GMT changes in temperature_change:
        try:
            temperature_change.loc[idx, 'bin_left'] = temperature_change.loc[idx, 'bin'].left
            temperature_change.loc[idx, 'bin_center'] = (temperature_change.loc[idx, 'bin'].left + temperature_change.loc[idx, 'bin'].right)/2
        except AttributeError:
            temperature_change.loc[idx, 'bin_left'] = -999
            temperature_change.loc[idx, 'bin_center'] = -999

    for i, gmt_bin in enumerate(np.sort(temperature_change.bin_center.unique())):
        bin_count.loc[i, 'bin'] = temperature_change.loc[temperature_change.bin_center==gmt_bin, 'bin'].unique()[0]
        bin_count.loc[i, 'lower_bound'] = temperature_change.loc[temperature_change.bin_center==gmt_bin, 'bin_left'].unique()[0]
        bin_count.loc[i, 'bin_center'] = gmt_bin
        bin_count.loc[i, 'count'] = temperature_change.loc[temperature_change.bin_center==gmt_bin].shape[0]
    temperature_change.gcm = temperature_change.gcm.str.lower() # make GCM names lower case
    if lower_bound_only:
        temperature_change = temperature_change.reset_index(drop=True)
        temperature_change['bin'] = copy.deepcopy(temperature_change['bin_left'])
        temperature_change = temperature_change.drop(columns='bin_left')
    elif bin_center_only:
        temperature_change = temperature_change.reset_index(drop=True)
        temperature_change['bin'] = copy.deepcopy(temperature_change['bin_center'])
        temperature_change = temperature_change.drop(columns='bin_center')
    bin_count.to_csv(co.out_dir / 'bin_count.csv', index=False)
    return temperature_change.reset_index(drop=True), bin_count

def map_run_to_gmt_bin(gcm, scenario, year, gmt_map_df, lower_bound_only=False, bin_center_only=True):
    """
    maps model run for GCM and scenario to GMT bin.

    Parameters
    ----------
    gcm : str
        climate model.
    scenario : str
        climate scenario / experiment.
    year : int
        year.
    gmt_map_df : DataFrame
        as returned by make_delta_gmt_df().
    lower_bound_only : boolean, optional
        The default is False.
    bin_center_only : boolean, optional
        The default is True.

    Raises
    ------
    ValueError

    Returns
    -------
    gmt_bin : 
        bin of level of global warming
    """
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
    if not (isinstance(gmt_bin[0], int) or isinstance(gmt_bin[0], float)):
        if lower_bound_only:
            return gmt_bin[0].left
        if bin_center_only:
            return (gmt_bin[0].left + gmt_bin[0].right)/2
    return gmt_bin[0]

def bin_country_impact(crop_type, co2='co2', absolute_production=True,
                       baseline_exp_file_path=None, gmt_bins=None, drop_nobin=False,
                       impact_dir=None, country_impact_combi_mean_ref=None):

    """
    load country impact per model run as computed with climada_wrappers.calc_country_impacts()
    and rearrange by: GMT bin, return 1 DataFrame per 'crop' type and 'co2' parameter,
    differentiating 'bin' and 'ggcm' in extra column.
    Also, baseline exposure value are added to impact values.

    Parameters:
        crop_type (str): 'whe', 'ric', 'soy', or 'mai'; or 'combi-tons', also '-kcal'

    Optional Parameters:
        absolute_production : TYPE, optional
            DESCRIPTION. The default is True.
        baseline_exp_file_path : TYPE, optional
            DESCRIPTION. The default is None.
        gmt_bins : TYPE, optional
            DESCRIPTION. The default is None.
        drop_nobin : TYPE, optional
            DESCRIPTION. The default is False.
        impact_dir : TYPE, optional
            DESCRIPTION. The default is None.
        co2 (str): CO2 fertilization parameter from ISIMIP, e.g. '2005co2' or 'co2' (default)
        country_impact_combi_mean_ref: For mean correction of detrended data: If DataFrame with binned country impacts
            (without detrending) are provided, the mean per GMT bin, gcm and ggcm is reset to mean over bin and ggcm (not gcm!) of reference
            X_corrected = X_detrended - mean(X_detrended) + mean(X)

    Returns:
        impact_comb = DataFrame with columns: 'bin', 'ggcm', 'gcm', 'scenario', 'Year', ...country ISO3..., 'GLB'
        impact_in_filenames: list
    """
    if impact_dir is None:
        impact_dir = co.impact_dir
    if gmt_bins is None: gmt_bins = co.gmt_bins
    if not baseline_exp_file_path: baseline_exp_file_path = impact_dir / 'baseline_exposure.csv'
    if not os.path.isfile(baseline_exp_file_path):
        if os.path.isfile(baseline_exp_file_path / 'baseline_exposure.csv'):
            baseline_exp_file_path = baseline_exp_file_path / 'baseline_exposure.csv'
        else:
            baseline_exp_file_path = co.out_dir / 'baseline_exposure.csv'
    if not os.path.isfile(baseline_exp_file_path): raise FileExistsError('baseline_exposure.csv: File not found.')
    impact_comb = pd.DataFrame()
    impact_in_filenames = [f for f in os.listdir(impact_dir) if (os.path.isfile(impact_dir / f)) if 
                     f.startswith('impact_') if f'_{crop_type}.csv' in f if f'_{co2}_' in f]
    if not impact_in_filenames and co2==co.co2_fert:
        impact_in_filenames = [f for f in os.listdir(impact_dir) if (os.path.isfile(impact_dir / f)) if 
                     f.startswith('impact_') if f'_{crop_type}.csv' in f]

    if absolute_production:
        production_baseline = pd.read_csv(baseline_exp_file_path, encoding="ISO-8859-1", header=0)
        production_baseline = production_baseline.fillna(0)

    # loop over impact files: load and combine in one DataFrame:
    for idf, fn in enumerate(impact_in_filenames):

        df_tmp = pd.read_csv(impact_dir / fn, encoding="ISO-8859-1", header=0)
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
    gmt_map_df = make_delta_gmt_df(bin_center_only=True, gmt_bins=gmt_bins)[0]

    impact_comb['bin'] = impact_comb.apply(lambda x: 
                            map_run_to_gmt_bin(x['gcm'], x['scenario'],x['Year'], gmt_map_df),
                            axis=1)

    cols = impact_comb.columns.tolist()[-1:] + impact_comb.columns.tolist()[:-1]
    impact_comb = impact_comb[cols].reset_index(drop=True)

    if drop_nobin:
        impact_comb = impact_comb.loc[impact_comb.bin>-99]

    impact_comb = impact_comb.reset_index(drop=True)

    if country_impact_combi_mean_ref is not None:
        print(f'Detrended country impacts in bins: correct mean value based on given stats ({crop_type}) \n') 
        # correct mean (mean_orig) value based on given dataframe
        # (with original, un-detrended binned country impacts)
        # requires non-detrended impacts

        for gmt_bin in impact_comb.bin.unique():
            # print(gmt_bin)
            for ggcm in impact_comb.ggcm.unique():
                # print(ggcm)
                for gcm in impact_comb.gcm.unique():
                    for col in impact_comb.columns:
                        if len(col)==3 and (col not in ['bin', 'Year', 'year', 'ggcm', 'gcm', 'N_years', 'scenario']):
                            # for mean orig (no detrending) GCMs are pooled, as later for stats
                            mean_orig = country_impact_combi_mean_ref.loc[(country_impact_combi_mean_ref.bin==gmt_bin) &
                                                                          (country_impact_combi_mean_ref.ggcm==ggcm), col].mean()
                            # for detrended mean, each gcm is corrected individually 
                            mean_detr = impact_comb.loc[(impact_comb.bin==gmt_bin) &
                                                        (impact_comb.ggcm==ggcm) &
                                                        (impact_comb.gcm==gcm), col].mean()
    
                            # std_detr = impact_comb.loc[(impact_comb.bin==gmt_bin) &
                            #                             (impact_comb.ggcm==ggcm), col].std()
                            # std_orig = country_impact_combi_mean_ref.loc[(country_impact_combi_mean_ref.bin==gmt_bin) &
                            #                                               (country_impact_combi_mean_ref.ggcm==ggcm), col].std()
                            # correct impact values by mean:
                            impact_comb.loc[(impact_comb.bin==gmt_bin) & (impact_comb.ggcm==ggcm) & (impact_comb.gcm==gcm), col] = \
                                impact_comb.loc[(impact_comb.bin==gmt_bin) & (impact_comb.ggcm==ggcm) & (impact_comb.gcm==gcm), col] - \
                                    mean_detr + mean_orig
                            # std_corr = impact_comb.loc[(impact_comb.bin==gmt_bin) &
                            #                             (impact_comb.ggcm==ggcm), col].std()
                            # mean_corr = impact_comb.loc[(impact_comb.bin==gmt_bin) &
                            #                             (impact_comb.ggcm==ggcm), col].mean()
                        # if col=='USA': # check
                        #     print(f'{col}: mean (detrend): {int(mean_detr):1.4e}')
                        #     print(f'{col}: mean (original): {int(mean_orig):1.4e}')
                        #     print(f'{col}: mean (corrected): {int(mean_corr):1.4e}')
                        #     print(f'{col}: std (detrended): {int(std_detr):1.4e}')
                        #     print(f'{col}: std (original): {int(std_orig):1.4e}')
                        #     print(f'{col}: std (corrected): {int(std_corr):1.4e}')
    return impact_comb.reset_index(drop=True), impact_in_filenames

def plot_GMT_bin_hist(gmt_bins=None, by='year'):
    """
    plot histogram of bins per by. by is either 'year' or 'dGMT

    Parameters
    ----------
    gmt_bins : list, optional
        list of GMT bins. The default is None.
    by : str, optional
        x-axis. The default is 'year'.

    Returns
    -------
    f_h : figure handle
    ax_h : axes handle
    """

    gmt_map, _ = make_delta_gmt_df(bin_center_only=True)
    if gmt_bins is None:
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
        plt.title(f'Histogram GMT bin {int(gmt_bins[idx])} +/- 0.5 K]')
    return f_h, ax_h

def bin_imp_events(event_df, check_plot=False):
    """
    cut up imp_mat and recombine by bins
    return sparse matrices and meta data per event etc.

    Parameters
    ----------
    event_df : DataFrame
    check_plot : bool, optional
        make plot? The default is False.

    Raises
    ------
    ValueError

    Returns
    -------
    None.
    """

    gmt, bin_count = make_delta_gmt_df(bin_center_only=True)
    ### add bin to event_name_df:
    event_df['gmt_center'] = np.nan
    event_df = event_df.reset_index(drop=True)
    for idx, _ in enumerate(event_df.index):
        bin_center = gmt.loc[(gmt.gcm==event_df.gcm[idx]) &
                            (gmt.experiment==event_df.experiment[idx]) &
                            (gmt.year==event_df.year[idx]), 'bin'].values
        if bin_center.size > 1:
            print('Error bin identification is not unique.')
            raise ValueError
        event_df.loc[idx, 'gmt_center'] = bin_center[0]
    if check_plot:
        plt.scatter(event_df.year, event_df.gmt_center)
    return event_df #, event_names, etc


def bin_imp_mat(impact_mat_combi, event_name_df, both_irr_id, gmt_bin):
    """
    cut up imp_mat and recombine by bins
    return sparse matrices and meta data per event etc.

    Parameters
    ----------
    impact_mat_combi : matrix, either sparse or not
        DESCRIPTION.
    event_name_df : DataFrame
        DESCRIPTION.
    both_irr_id : int
        DESCRIPTION.
    gmt_bin : float
        DESCRIPTION.

    Returns
    -------
    None.

    """

    if not (isinstance(gmt_bin, float) or isinstance(gmt_bin, int)):
        gmt_bin = (gmt_bin.left+gmt_bin.right)/2
    events_df = event_name_df.loc[(event_name_df.both_irr_id==both_irr_id)]
    if np.array_equal(events_df.loc[events_df.irr=='firr'].event_name, events_df.loc[events_df.irr=='noirr'].event_name):
        events_df = events_df.loc[events_df.irr=='noirr']

    # cut matrix_
    impact_mat_binned = sparse.csr_matrix(impact_mat_combi[(events_df.gmt_center == gmt_bin).values, :])
    
    # a bit like in cw: impact_sets_or_mats[idx_list[firr_i]] = sparse.csr_matrix(impact_sets_or_mats[idx_list[firr_i]]
    #                                                                               [names_map, :]) # Sort along ax=0
    return impact_mat_binned, events_df.loc[events_df.gmt_center == gmt_bin]
