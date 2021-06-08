#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:44:37 2020

@author: eberenzs
"""

import os
import numpy as np
import pandas as pd
from collections import OrderedDict 
# import copy
import matplotlib.pyplot as plt

import isimip3b_crop_config as co
import isimip3b_gmt_binning as gmt_binning

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
    cols = stats_df.columns.tolist()[1:4] + stats_df.columns.tolist()[-1:] + stats_df.columns.tolist()[0:1] + stats_df.columns.tolist()[7:-1]
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


def plot_lines_with_intervals(x, y_lines, y_fill=None, y_ens=None, colors=None, linestyles=None,
                              labels=None, title=None, xlabel=None, ylabel=None, fill=True,
                              alpha_fill=.15, alpha_ens=.22, linestyles_interval=None,
                              figsize=None, x_type=int, log_x=True, log_y=False,
                              xticks=None, ylims=None, yticks=None, grid_on=False,
                              in_labels=None,
                              in_label_keys=None, hlines=None, label_str=None):
    """Generic plot of lines with intervals
    
    Parameters:
        x (list or array)
        y_lines (list of lists or arrays): inner lists need to have same length as x
        y_fill(list of lists of lists or arrays): most inner lists need to have same length as x
        y_ens(list of lists of lists or arrays): most inner lists need to have same length as x
        colors (list): same length as number of lists or arrays in y_lines
    """
    if label_str is None:
        label_str = '%f'
    print(label_str)
    if linestyles is None:
        linestyles = ['-' for y in y_lines]
    if linestyles_interval is None: linestyles_interval = linestyles
    if colors is None:
        colors = co.colors8_heating[0:len(linestyles)]
    if figsize is None: # tuple: (wdth, hght)
        fig = plt.figure(facecolor='w')
    else:
        fig = plt.figure(facecolor='w', figsize=figsize)
    ax = fig.add_subplot(1,1,1)
    if hlines:
        for hline in hlines:
            ax.axhline(y=hline, c='k', ls=':', lw=.5, alpha=.77, zorder=1)
    for i, y in enumerate(y_lines):
        if 'perc' in label_str:
            label_ = label_str % (100*labels[i])
        else:
            label_ = label_str % (labels[i])
        ax.plot(x, y, linestyle=linestyles[i], color=colors[i], label=label_, alpha=1, zorder=5)
        if fill and y_fill:
            ax.fill_between(x, y_fill[i][0], y_fill[i][-1], facecolor=colors[i],
                            label=None, alpha=alpha_fill, zorder=3)
        elif y_fill:
            for y_ in y_fill[i]:
                ax.plot(x, y_fill, color=colors[i], linestyle=linestyles_interval[i],
                        label=None, alpha=alpha_ens, zorder=3)
        if y_ens:
            for i_, y_ in enumerate(y_ens[i]):
                ax.plot(x, y_, color=colors[i], linestyle=linestyles_interval[i],
                        label=None, alpha=alpha_ens, zorder=4)
                if in_labels:

                    ax.text(x[-5-3*i], y_[-5-3*i], f'({in_label_keys[i_]})', color=colors[i], zorder=7, label=None)
                    ax.text(1.1*x[-5-3*i], y_[-5-3*i], in_labels[in_label_keys[i_]], color=colors[i], zorder=7, label=None)
                    # ax.text(x[-2-2*i], y_[-2-2*i], f'({in_label_keys[i_]})', color=colors[i], zorder=7, label=label_)
                    # ax.text(1.2*x[-2-2*i], y_[-2-2*i], in_labels[in_label_keys[i_]], color=colors[i], zorder=7, label=None)
    if log_x:
        ax.set_xscale('log')
    if log_y:
        ax.set_yscale('log')
    if xticks is None:
        xticks = np.array(x)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks.astype(x_type))
    if ylims:
        ax.set_ylim(np.min(ylims), np.max(ylims))
    if yticks:
        ax.set_yticks(np.array(yticks))

    ax.grid(grid_on)
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return fig, ax


def plot_PR_change_per_bin_country(input_dir=co.out_dir / 'Stats_countries' / 'overlapping',
                                   crop=None, countries=None, percentiles=None,
                                   remove_small_samples=True, log_x=False, log_y=False,
                                   gmt_bins=None, ref_bin=0.5, max_bin=5.,
                                   ggcms=None, in_labels=True,
                                   x_ticks=None,
                                   ylims=(-1, 21), yticks=[0, 1, 2, 5, 10, 15, 20],
                                   save_dir=co.plot_dir, save_plot='pdf'):
    abc = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']

    if save_plot and (not save_dir):
        save_dir = co.plot_dir
    if percentiles is None:
        percentiles = [1/20, 1/40]
    if ggcms is None: ggcms = co.ggcms[co.input_version]
    ggcm_abc = dict()
    for iggcm, ggcm in enumerate(ggcms): # ggcm
        ggcm = ggcm.lower()
        ggcm_abc[abc[iggcm]] = ggcm
    if isinstance(crop, str):
        crop = [crop]
    elif crop is None:
        crop = co.crop_types
    if countries is None:
        countries = co.in_dir / 'main_producer.csv'
    if (isinstance(countries, os.PathLike) or isinstance(countries, str)) and os.path.isfile(countries):
        main_producers = pd.read_csv(countries, encoding="ISO-8859-1", header=0)
        countries = []
        for cr in crop:
            countries += list(main_producers.loc[main_producers.crop==cr, 'country'])#[:5]
        countries = ['GLB'] + list(OrderedDict.fromkeys(countries))
        print(countries)
    elif isinstance(countries, str):
        countries = [countries]
    print(crop)
    if log_y:
        hlines=[0.]
    else:
        hlines=[1.]
    colors = [co.colors8_heating[0],  co.colors8_heating[3], co.colors8_heating[2]]

    stats_df = [pd.read_csv(input_dir/ f'stats_country_binned_rel2bin_{10*ref_bin:02.0f}_{cr}.csv', # stats_country_binned_rel2bin_0_mai.csv
                             encoding="ISO-8859-1", header=0) for cr in crop]
    if gmt_bins is None:
        gmt_bins = np.sort(stats_df[0]['bin'].unique())
        gmt_bins = gmt_bins[gmt_bins >= ref_bin]
        gmt_bins = gmt_bins[gmt_bins <= max_bin]
    if x_ticks is None:
        x_ticks = np.arange(ref_bin, max_bin+.001, .5)
    x_array = gmt_bins
    fig_pr = list()
    for country in countries: # loop over countries
        print(country)
        # loop over crops and ggcms:
        for i_crop, df in enumerate(crop): # crop
            df = stats_df[i_crop]
            print(crop[i_crop])
            # init dicts and lists for loop over ggcms:
            values = list()
            means = list()
            medians = list()
            min_max = list()
            iqr = list()
            for pc in percentiles: # loop over percentiles / frequencies
                values.append(list())
                pc_stack = dict()
                i_pc = 0
                in_label_keys = list()
                for ggcm in ggcms: # ggcm
                    ggcm = ggcm.lower()
                    if df.loc[(df.ggcm==ggcm) & (df.stat==f'{pc}'), country].size==0:
                        continue
                    values[-1].append(list())
                    if in_labels:
                        in_label_keys.append(list(ggcm_abc.keys())[list(ggcm_abc.values()).index(ggcm)])
                    for gmt_bin in gmt_bins: # fill list with value per gmt_bin:
                        if df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{pc}'), country].size==0:
                            continue
                        sample_size = df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{pc}'), 'N_years'].values[0]
                        if remove_small_samples and sample_size < 1/pc: # sample size lower than return period? --> NaN
                            values[-1][-1].append(np.nan)
                        else: 
                            values[-1][-1].append( # probability ratio: percentile divided by ref. pc
                                df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{pc}'), country].values[0] / pc
                                )
                    if i_pc==0:
                        pc_stack[pc] = np.array(values[-1][-1])
                    else:
                        pc_stack[pc] = np.vstack((pc_stack[pc], np.array(values[-1][-1])))
                    i_pc += 1

                means.append(np.nanmean(pc_stack[pc], axis=0))
                medians.append(np.nanmedian(pc_stack[pc], axis=0))
                min_max.append((np.nanmin(pc_stack[pc], axis=0), np.nanmax(pc_stack[pc], axis=0)))
                iqr.append((np.nanquantile(pc_stack[pc], .25, axis=0), np.nanquantile(pc_stack[pc], .75, axis=0)))

            if not in_labels:
                ggcm_abc = False
            fig_pr.append(
                plot_lines_with_intervals(x_array, medians, y_fill=iqr, y_ens=values, labels=percentiles,
                                          fill=True, x_type=float, colors=colors,
                                          xlabel='Global mean warming (°C)',
                                          ylabel='PR',
                                          title=f'{country} - {crop[i_crop]}',
                                          log_x=log_x, log_y=log_y,
                                          xticks=x_ticks, ylims=ylims, yticks=yticks,
                                          alpha_ens=.2, in_labels=ggcm_abc, in_label_keys=in_label_keys,
                                          hlines=hlines, label_str='below %2.1fth perc.')[0]
                )
            if save_plot:
                fig_pr[-1].savefig(save_dir / f'ProbRatio_{crop[i_crop]}_{country}_ref_bin_{ref_bin:0.2f}_ggcmlabels_{in_labels}.{save_plot}', \
                            dpi=600, facecolor='w', edgecolor='w', \
                            orientation='portrait', papertype=None, format=save_plot, \
                            transparent=True, bbox_inches=None, pad_inches=0.1, \
                            frameon=None, metadata=None)
    return fig_pr, ggcm_abc



def plot_RP_per_bin_country(input_dir=co.out_dir / 'Stats_countries' / 'overlapping',
                            crop=None, countries=None,
                            frequencies=None, remove_small_samples=True, log_x=True,
                            log_y=False, rp=True, gmt_bins=None,
                            ggcms=None, deviation_from=None, save_dir=co.plot_dir, save_plot='pdf',
                            rel2ref=False, plot_ens=False, in_labels=False):
    """ relative to same bin?"""
    abc = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']
    if deviation_from is None: deviation_from = 'mean_rel'
    if gmt_bins is None:
        gmt_bins = [0.5, 2, 4.]
    if frequencies is None:
        frequencies = [1/100, 1/80, 1/60, 1/50, 1/40, 1/30, 1/20, 1/10, 1/5, 1/3, 1/2]
        # [1/100, 1/40] # [1/100, 1/80, 1/60, 1/40, 1/20, 1/10, 1/5, 1/2]
    if rp:
        frequencies.reverse()
        x_array = np.array(frequencies)**(-1)
    else: x_array = frequencies
    if ggcms is None: ggcms = co.ggcms[co.input_version]
    ggcm_abc = dict()
    for iggcm, ggcm in enumerate(ggcms): # ggcm
        ggcm = ggcm.lower()
        ggcm_abc[abc[iggcm]] = ggcm
    if isinstance(crop, str):
        crop = [crop]
    elif crop is None:
        crop = co.crop_types
    if countries is None:
        countries = co.in_dir / 'main_producer.csv'
    if (isinstance(countries, os.PathLike) or isinstance(countries, str)) and os.path.isfile(countries):
        main_producers = pd.read_csv(countries, encoding="ISO-8859-1", header=0)
        countries = []
        for cr in crop:
            countries += list(main_producers.loc[main_producers.crop==cr, 'country'])#[:5]
        countries = ['GLB'] + list(OrderedDict.fromkeys(countries))
        print(countries)
        # countries = countries[0:6]
    elif isinstance(countries, str):
        countries = [countries]
    print(crop)
    colors = [co.colors8_heating_bins[int(gmt_bin)] for gmt_bin in gmt_bins]

    stats_df = [pd.read_csv(input_dir/ f'stats_country_binned_self_{cr}.csv',
                             encoding="ISO-8859-1", header=0) for cr in crop]
    fig_rp = list()
    for country in countries: # loop over countries
        print(country)
        # loop over crops and ggcms:
        for i_crop, df in enumerate(crop): # crop
            df = stats_df[i_crop]
            print(crop[i_crop])
            # init dicts and lists for loop over ggcms:
            values = list()
            means = list()
            medians = list()
            min_max = list()
            iqr = list()
            twothirds = list()
            for gmt_bin in gmt_bins: # loop over GMT bins
                values.append(list())
                ggcm_stack = dict()
                iggcm = 0
                in_label_keys = list()
                for ggcm in ggcms: # ggcm
                    ggcm = ggcm.lower()
                    if df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin), country].size==0:
                        continue
                    values[-1].append(list())
                    if in_labels:
                        in_label_keys.append(list(ggcm_abc.keys())[list(ggcm_abc.values()).index(ggcm)])
                    for f in frequencies: # fill list with value per frequency:
                        sample_size = df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{f}'), 'N_years'].values[0]
                        if remove_small_samples and sample_size < 1/f: # sample size lower than return period? --> NaN
                            values[-1][-1].append(np.nan)
                        else: 
                            values[-1][-1].append(
                                df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{f}'), country].values[0]
                                )
                            if deviation_from is not None and 'mean' in deviation_from:
                                if rel2ref:
                                    ref_values = df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bins[0]) & (df.stat=='mean'), country].values[0]
                                else:
                                    ref_values = df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='mean'), country].values[0]
                            elif deviation_from is not None and 'median' in deviation_from:
                                if rel2ref:
                                    ref_values = df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bins[0]) & (df.stat=='0.5'), country].values[0]
                                else:
                                    ref_values = df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='0.5'), country].values[0]
    
                            if deviation_from is not None and '_rel' in deviation_from:
                                values[-1][-1][-1] = 100 * (values[-1][-1][-1] / ref_values) - 100
                            elif deviation_from is not None:
                                values[-1][-1][-1] = values[-1][-1][-1] - ref_values

                    # # plot ggcms
                    # ax_efc[-1].plot(x_array, values[gmt_bin][ggcm],linestyle=ls_list['gmt_bin'],
                    #                 color=co.colors_ggcm[co.input_version][ggcm], label=ggcm, alpha=alpha_ggcm)
                    # ax_efc[-1].plot(x_array, values[ref_bin][ggcm],linestyle=ls_list['ref_bin'],
                    #                 color=co.colors_ggcm[co.input_version][ggcm], label=None, alpha=alpha_ggcm)
                    if iggcm==0:
                        ggcm_stack[gmt_bin] = np.array(values[-1][-1])
                    else:
                        ggcm_stack[gmt_bin] = np.vstack((ggcm_stack[gmt_bin], np.array(values[-1][-1])))
                    iggcm += 1
                means.append(np.nanmean(ggcm_stack[gmt_bin], axis=0))
                medians.append(np.nanmedian(ggcm_stack[gmt_bin], axis=0))
                min_max.append((np.nanmin(ggcm_stack[gmt_bin], axis=0), np.nanmax(ggcm_stack[gmt_bin], axis=0)))
                iqr.append((np.nanquantile(ggcm_stack[gmt_bin], .25, axis=0), np.nanquantile(ggcm_stack[gmt_bin], .75, axis=0)))
                twothirds.append((np.nanquantile(ggcm_stack[gmt_bin], 1/6, axis=0), np.nanquantile(ggcm_stack[gmt_bin], 5/6, axis=0)))

            # fig_rp.append(
            #     plot_lines_with_intervals(x_array, means, min_max, labels=gmt_bins, fill=True, x_type=float, colors=colors, 
            #                               title=f'{country} - {crop[i_crop]}')
            #     )

            if not in_labels:
                ggcm_abc = False
            if not plot_ens:
                values=None
            fig_rp.append(
                plot_lines_with_intervals(x_array, medians, y_fill=iqr, # y_fill=iqr, # y_fill=twothirds
                                          y_ens=values, labels=gmt_bins, fill=True, x_type=int, colors=colors, 
                                          title=f'{country} - {crop[i_crop]}',
                                          xlabel='Return period (years)',
                                          ylabel='$\Delta$ Crop production (%)',
                                          log_x=log_x, log_y=log_y, alpha_ens=.2, in_labels=ggcm_abc,
                                          in_label_keys=in_label_keys,
                                          hlines=[0], label_str='%1.1f°C')[0]
                )
            if save_plot:
                if rel2ref:
                    fig_file_str = f'ReturnP_{crop[i_crop]}_{country}_ggcmlabels_{in_labels}_rel2ref.{save_plot}'
                else:
                    fig_file_str = f'ReturnP_{crop[i_crop]}_{country}_ggcmlabels_{in_labels}.{save_plot}'
                fig_rp[-1].savefig(save_dir / fig_file_str, \
                            dpi=600, facecolor='w', edgecolor='w', \
                            orientation='portrait', papertype=None, format=save_plot, \
                            transparent=True, bbox_inches=None, pad_inches=0.1, \
                            frameon=None, metadata=None)
    return fig_rp, ggcm_abc



def plot_stats_country_bin_EFC(input_dir=co.stats_cntry_dir, crop=None, country=None,
                               frequencies=None, rp=True, log_x=False,
                               gmt_bin=2.0, ref_bin=0.5, ggcms=None, deviation_from=None,
                               save_dir = None):

    if not frequencies:
        frequencies=[1/60, 1/40, 1/20, 1/10, 1/5, 1/2] # [1/100, 1/80, 1/60, 1/40, 1/20, 1/10, 1/5, 1/2]
    if ggcms is None: ggcms = co.ggcms[co.input_version]
    if rp:
        frequencies.reverse()
        x_array = np.array(frequencies)**(-1)
    else: x_array = frequencies
    if isinstance(crop, str):
        crop = [crop]
    elif crop is None:
        crop = co.crop_types
    if country is None:
        country = 'GLB'

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
            ggcm = ggcm.lower()
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
            print(stat)
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
