#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:44:37 2020

@author: eberenzs
"""

import numpy as np
import pandas as pd
# import copy
import matplotlib.pyplot as plt

import crop_config as co
import gmt_binning

def single_stat_country(df, stat):
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

def stats_country_binned_self(impact_comb, statistics_list=None, crop=None):
    """Compute statistics per country and globally for a given combined country-impact dataframe
    
    Parameters:
        impact_comb (df): DataFrame with country impacts combined by GMT bin
        Statistics are relative to the same bin and GGCM

    Optional Parameters:
        statistics_list(list): c.f. crop_config.statistics_list (default)
    """

    if not statistics_list: statistics_list = co.statistics_list

    stats_df = pd.DataFrame(columns=['stat']+impact_comb.columns.tolist())

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
    cols = stats_df.columns.tolist()[1:4] + stats_df.columns.tolist()[-1:] + stats_df.columns.tolist()[0:1] + stats_df.columns.tolist()[7:-1]
    stats_df[cols].reset_index(drop=True).to_csv(co.stats_cntry_dir / f'stats_country_binned_self_{crop}.csv', index=False)
    return stats_df[cols].reset_index(drop=True)


def plot_stats_country_bin_EFC(input_dir=co.stats_cntry_dir, crop=None, country=None,
                               frequencies=None, rp=True, log_x=False,
                               gmt_bin=2.0, ref_bin=0.5, ggcms=None, deviation_from=None,
                               save_dir = None):

    if not frequencies:
        frequencies=[1/60, 1/40, 1/20, 1/10, 1/5, 1/2] # [1/100, 1/80, 1/60, 1/40, 1/20, 1/10, 1/5, 1/2]
    if country is None: country='GLB'
    if ggcms is None: ggcms = co.ggcms[co.input_version]
    if rp:
        frequencies.reverse()
        x_array = np.array(frequencies)**(-1)
    else: x_array = frequencies
    if isinstance(crop, str):
        crop = [crop]
    elif crop is None:
        crop = co.crop_types

    # plotting option:
    ls_list = {'gmt_bin': '-',
               'ref_bin': '--'}
    alpha_ggcm = .4
    if deviation_from is not None and '_rel' in deviation_from:
        ylabel = 'Rel.deviation from mean production [%]'
    else:
        ylabel = 'crop production [t]'

    # read statistics dataframes for each crop type:
    stats_df = [pd.read_csv(input_dir/ f'stats_country_binned_self_{cr}.csv',
                             encoding="ISO-8859-1", header=0) for cr in crop]

    fig_efc = list()
    ax_efc = list()

    # loop over crops and ggcms:
    for df in stats_df: # crop
        fig_efc.append(plt.figure(facecolor='w'))
        ax_efc.append(fig_efc[-1].add_subplot(1,1,1))
        # init dicts and lists for loop over ggcms:
        values = dict()
        values[gmt_bin] = dict()
        values[ref_bin] = dict()
        ggcm_stack = dict()
        for iggcm, ggcm in enumerate(ggcms): # ggcm
            values[gmt_bin][ggcm] = list()
            values[ref_bin][ggcm] = list()
            if df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin), country].size==0:
                continue
            for f in frequencies: # fill list with value per frequency:
                values[gmt_bin][ggcm].append(
                    df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{f}'), country].values[0]
                    )
                values[ref_bin][ggcm].append(
                    df.loc[(df.ggcm==ggcm) & (df.bin==ref_bin) & (df.stat==f'{f}'), country].values[0]
                    )
                if deviation_from is not None and deviation_from=='mean':
                    values[gmt_bin][ggcm][-1] = values[gmt_bin][ggcm][-1] - \
                        df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='mean'), country].values[0]
                    values[ref_bin][ggcm][-1] = values[ref_bin][ggcm][-1] - \
                        df.loc[(df.ggcm==ggcm) & (df.bin==ref_bin) & (df.stat=='mean'), country].values[0]
                if deviation_from is not None and deviation_from=='mean_rel':
                    values[gmt_bin][ggcm][-1] = 100 * (values[gmt_bin][ggcm][-1] / \
                        df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='mean'), country].values[0]) - 100
                    values[ref_bin][ggcm][-1] = 100 * values[ref_bin][ggcm][-1] / \
                        df.loc[(df.ggcm==ggcm) & (df.bin==ref_bin) & (df.stat=='mean'), country].values[0] - 100
            # plot ggcms
            ax_efc[-1].plot(x_array, values[gmt_bin][ggcm],linestyle=ls_list['gmt_bin'],
                            color=co.colors_ggcm[co.input_version][ggcm], label=ggcm, alpha=alpha_ggcm)
            ax_efc[-1].plot(x_array, values[ref_bin][ggcm],linestyle=ls_list['ref_bin'],
                            color=co.colors_ggcm[co.input_version][ggcm], label=None, alpha=alpha_ggcm)
            if iggcm==0:
                ggcm_stack[gmt_bin] = np.array(values[gmt_bin][ggcm])
                ggcm_stack[ref_bin] = np.array(values[ref_bin][ggcm])
            else:
                ggcm_stack[gmt_bin] = np.vstack((ggcm_stack[gmt_bin], np.array(values[gmt_bin][ggcm])))
                ggcm_stack[ref_bin] = np.vstack((ggcm_stack[ref_bin], np.array(values[ref_bin][ggcm])))
         # plot ggcm mean:
        ax_efc[-1].plot(x_array, np.mean(ggcm_stack[gmt_bin], axis=0), linestyle=ls_list['gmt_bin'],
                        color='k', label=f'GGCM mean {gmt_bin:0.1f}-{gmt_bin+.5:0.1f}', alpha=1.)
        ax_efc[-1].plot(x_array, np.mean(ggcm_stack[ref_bin], axis=0), linestyle=ls_list['ref_bin'],
                        color='k', label=f'GGCM mean {ref_bin:0.1f}-{ref_bin+.5:0.1f}', alpha=1.)
        if log_x:
            ax_efc[-1].set_xscale('log')
        ax_efc[-1].set_xticks(x_array.astype(int))
        ax_efc[-1].set_xticklabels(x_array.astype(int))

        ax_efc[-1].legend()
    #    ax_efc[-1].set_xlim((-.5, idb+.5))
    #    ax_efc[-1].set_xticks(np.arange(idb+1))
    #    ax_efc[-1].set_xticklabels(np.sort(stats_df.bin.unique()))
        ax_efc[-1].grid('on', alpha=.5)
        if rp:
            ax_efc[-1].set_xlabel('Return period [years]')
        else:
            ax_efc[-1].set_xlabel('Frequency')
        ax_efc[-1].set_ylabel(ylabel)
        crop_str = crop[0]
        for cr in crop[1:]:
            crop_str + f' + {cr}'
        ax_efc[-1].set_title(f'EFC: {crop_str} - {country}')
        # save:
        if not (save_dir is None):
            fig_efc[-1].savefig(save_dir / f'EFC_{crop_str}_{country}_{deviation_from}_gmt_bin_{gmt_bin:0.2f}_ref_bin_{ref_bin:0.2f}.pdf', \
                            dpi=600, facecolor='w', edgecolor='w', \
                            orientation='portrait', papertype=None, format='pdf', \
                            transparent=True, bbox_inches=None, pad_inches=0.1, \
                            frameon=None, metadata=None)
    return fig_efc, ax_efc

def stats_country_bin_plots(stats_df, crop, combis=None, ylabel=None,
                            statistics_list=None, marker_list=None, ls_list=None, country='GLB'):
    """Plots based on alternative (exclusice) optional parameters"""
    # TODO : group small bins with less than 30 years?
    
    if not statistics_list: statistics_list = list()
    if not marker_list: marker_list = list()
    if not ls_list: ls_list = list()
    if combis:
        if 'mean+-' in combis:
            statistics_list.append('mean')
            statistics_list.append('mean+std')
            statistics_list.append('mean-std')
            marker_list.append('o')
            marker_list.append('_')
            marker_list.append('_')
            ls_list += ['-', '', '']
            ylabel = 'Crop production [t/y]'
        if 'median+-' in combis:
            statistics_list.append(.5)
            statistics_list.append(.25)
            statistics_list.append(.75)
            marker_list.append('o')
            marker_list.append('_')
            marker_list.append('_')
            ls_list += ['-', '', '']
            ylabel = 'Crop production [t/y]'
    # f_b = plt.figure(facecolor='w')
    # ax_b = f_b.add_subplot(1,1,1)
    values = dict()

    for ggcm in stats_df.ggcm.unique():
        label_ = ggcm
        values[ggcm] = dict()
        for ids, stat in enumerate(statistics_list):
            values[ggcm][stat] = list()
            for idb, gmt_bin in enumerate(np.sort(stats_df.bin.unique())):
                if idb or ids:
                    label_ = None
                if isinstance(stat, str) and '+' in stat:
                    val = stats_df.loc[(stats_df.ggcm==ggcm) & (stats_df.bin==gmt_bin) & (stats_df.stat==stat.split('+')[0]), country].values[0] +\
                          stats_df.loc[(stats_df.ggcm==ggcm) & (stats_df.bin==gmt_bin) & (stats_df.stat==stat.split('+')[1]), country].values[0]
                elif isinstance(stat, str) and '-' in stat:
                    val = stats_df.loc[(stats_df.ggcm==ggcm) & (stats_df.bin==gmt_bin) & (stats_df.stat==stat.split('-')[0]), country].values[0] -\
                          stats_df.loc[(stats_df.ggcm==ggcm) & (stats_df.bin==gmt_bin) & (stats_df.stat==stat.split('-')[1]), country].values[0]
                else:
                    val = stats_df.loc[(stats_df.ggcm==ggcm) & (stats_df.bin==gmt_bin) & (stats_df.stat==stat), country].values[0]
                # ax_b.plot(idb, val, marker=marker_list[ids], color=co.colors_ggcm[co.input_version][ggcm], label=label_, alpha=1.)
                values[ggcm][stat].append(val)

    f_b = plt.figure(facecolor='w')
    ax_b = f_b.add_subplot(1,1,1)
    for ggcm in stats_df.ggcm.unique():
        label_ = ggcm
        for ids, stat in enumerate(statistics_list):
            if ids: label_ = None
            ax_b.plot(np.arange(idb+1), values[ggcm][stat], marker=marker_list[ids],linestyle=ls_list[ids], color=co.colors_ggcm[co.input_version][ggcm], label=label_, alpha=1.)
    ax_b.legend()
    ax_b.set_xlim((-.5, idb+.5))
    ax_b.set_xticks(np.arange(idb+1))
    ax_b.set_xticklabels(np.sort(stats_df.bin.unique()))
    ax_b.set_xlabel('GMT bin (lower bound) [K]')
    ax_b.set_ylabel(ylabel)
    ax_b.set_title(f'{crop}: {str(statistics_list)}')

    return f_b, ax_b

def stats_country_binned_rel2bin(impact_comb, reference_bin=None, statistics_list=None):
    """Compute statistics per country and globally for a given combined country-impact dataframe
    Statistics are relative to the statistics of a reference_bin of the same GGCM.
    To answer question: what is the percentile of years below a certain percentile value from reference

    TODO: not finished, .... is it really necessary?

    Parameters:
        impact_comb (df): DataFrame with country impacts combined by GMT bin

    Optional Parameters:
        reference_bin (object): 
        statistics_list(list): c.f. crop_config.statistics_list (default)
    """
    if not reference_bin: reference_bin = co.reference_bin
    if not statistics_list: statistics_list = co.statistics_list_ref
    # 1. use stats_country_binned_self to get reference statistics
    stats_ref = stats_country_binned_self(impact_comb, statistics_list=statistics_list)
    stats_ref = stats_ref.loc[stats_ref.bin==reference_bin]
    # 2. get statistics relative to reference statistics of reference_bin and same ggcm
    return None

def single_stat_mat(mat, stat):
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
    if not statistics_list: statistics_list = co.statistics_list_mat
    if mask_zeros: mask = single_stat_mat(imp_mat, 'max==0')
    stat_array_out_list = list()
    for idx, stat in enumerate(statistics_list):
        # filter by bin, ggcm, and stat and calculate statistic
        stat_array_out_list.append(single_stat_mat(imp_mat, stat))
        if mask_zeros: stat_array_out_list[-1][mask] = np.nan
    return stat_array_out_list, statistics_list

def stats_mat_binned_rel2bin(imp_mat, imp_mat_ref, statistics_list=None, mask_zeros=True):
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
            print(f'gmt_bin {gmt_bin}-{gmt_bin+0.5}: N={imp_mats_binned_tmp.shape[0]}')
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
