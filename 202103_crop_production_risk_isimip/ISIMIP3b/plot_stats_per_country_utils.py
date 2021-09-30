#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 10:11:04 2021

@author: eberenzs
"""
import os
import numpy as np
import pandas as pd
from collections import OrderedDict 

import matplotlib.pyplot as plt

import isimip3b_crop_config as co




def plot_ax_binned_impacts(df, fig=None, ax=None, subplot_i=(1,1,1),
                           gmt_bins_col=None, colors_per_bin=None,
                           hlines=None, title='', ylims=None, ylim_min=False,
                           marker='.', figsize=None, country=None,
                           ylabel='', xlabel='Year', legend=True):
    if country is None:
        country = 'GLB'
    if fig is None:
        if figsize is None: # tuple: (wdth, hght)
            fig = plt.figure(facecolor='w')
        else:
            fig = plt.figure(facecolor='w', figsize=figsize)
    if ax is None:
        ax = fig.add_subplot(subplot_i[0], subplot_i[1], subplot_i[2])
    if gmt_bins_col is None:
        gmt_bins_col = [0.5, 2, 4]
    if colors_per_bin is None:
        colors_per_bin = {0.5: co.colors8_heating[0],
                      2: co.colors8_heating[2],
                      4: co.colors8_heating[4],}
        
    for gmt_bin in df.bin.unique():
        if not (gmt_bin in gmt_bins_col):
            colors_per_bin[gmt_bin] = '#636363' # grey

    for gmt_bin in np.sort(df.bin.unique()):
        if gmt_bin in gmt_bins_col:

            if hlines:
                xx = [df.loc[df.bin==gmt_bin, 'Year'].min(), df.loc[df.bin==gmt_bin, 'Year'].max()]
                for hline in hlines:
                    if isinstance(hline, float):
                        yy = [df.loc[df.bin==gmt_bin, country].quantile(q=hline), df.loc[df.bin==gmt_bin, country].quantile(q=hline)]
                        ls = ':'
                        label_hline = f'{100*hline:1.1f}th pctl.'
                    elif hline == 'mean':
                        yy = [df.loc[df.bin==gmt_bin, country].mean(), df.loc[df.bin==gmt_bin, country].mean()]
                        ls = '-'
                        label_hline = hline
                    elif hline == 'mean-10%':
                        ls = '--'
                        yy = [.9*df.loc[df.bin==gmt_bin, country].mean(), .9*df.loc[df.bin==gmt_bin, country].mean()]
                        label_hline = hline
                    if gmt_bin == gmt_bins_col[0] and hline==0.025:
                        ax.axhline(y=yy[0], c=colors_per_bin[gmt_bin], ls='-', lw=3,
                            alpha=.5, zorder=21+int(gmt_bin), label=None)
                    else:
                        ax.plot(xx, yy, c=colors_per_bin[gmt_bin], ls='-', lw=3,
                                alpha=.5, zorder=20+int(gmt_bin), label=None)
                    if gmt_bin == gmt_bins_col[0]:
                        label = label_hline
                    else:
                        label = None
                    if gmt_bin == gmt_bins_col[0] and hline==0.025:
                        ax.axhline(y=yy[0], c='k', ls=ls, lw=1,
                            alpha=1, zorder=21+int(gmt_bin), label=label)
                    else:
                        ax.plot(xx, yy, c='k', ls=ls, lw=1,
                            alpha=1, zorder=21+int(gmt_bin), label=label)
            ax.scatter(df.loc[df.bin==gmt_bin, 'Year'], (-1)*df.loc[df.bin==gmt_bin, country], c=colors_per_bin[gmt_bin],
                       s=40, marker=marker, alpha=.99, label=r'$\Delta$GMT=' + f'{gmt_bin:1.1f}°C', zorder=-1)
            ax.scatter(df.loc[df.bin==gmt_bin, 'Year'], df.loc[df.bin==gmt_bin, country], c=colors_per_bin[gmt_bin],
                       s=10, marker=marker, alpha=.667, label=None, zorder=10+int(gmt_bin))


        else:
            ax.scatter(df.loc[df.bin==gmt_bin, 'Year'], df.loc[df.bin==gmt_bin, country], c=colors_per_bin[gmt_bin],
                       s=10, marker=marker, alpha=.667, label=None, zorder=3)


    ax.set_xlim(1850, 2100)
    if ylims is None:
        ylims = (df[country].min(), df[country].max())
    if ylim_min:
        ax.set_ylim(ylims[0], ylims[1])
    else:
        ax.set_ylim(0, ylims[1])
    if legend:
        #ax.legend()
        
        handles, labels = plt.gca().get_legend_handles_labels()
        labels_sort, handles_sort = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))

        order = [0, 1, 2, 4, 5, 3]
        #ax.legend(handles_sort, labels_sort, ncol=2)
        ax.legend([handles_sort[idx] for idx in order],[labels_sort[idx] for idx in order], ncol=2)
        #ax.legend(ncol=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    return fig, ax


def plot_binned_impacts(input_dir=co.out_dir, fn_str=None, crop=None, bin_colors=None,
                        hlines=None, hlabels=None, country=None, per_ggcm=True,
                        marker=None, figsize=None, detrended=False, factor=1e6,
                        ylabel='Crop production (mio. t/y)', xlabel='Year',
                        save_dir=co.plot_dir, save_plot='png'):
    if hlines is None:
        hlines = ['mean', 'mean-10%', 0.025]
    if country is None:
        country = 'GLB'
    if fn_str is None:
        if detrended:
            fn_str = 'country_impact_t_binned_%s_None_Detrended.csv' # 'country_impact_t_binned_{crop}_{irr}_{subdir}.csv'
        else:
            fn_str = 'country_impact_t_binned_%s_None_.csv'
    if isinstance(crop, str):
        crop = [crop]
    elif crop is None:
        crop = co.crop_types + ['combi-tons']
    if marker is None:
        marker = '.'
    for icr, cr in enumerate(crop):
        df = pd.read_csv(input_dir / (fn_str % (cr)), encoding="ISO-8859-1", header=0)
        df[country] = df[country] / factor
        if ('combi' in cr) and co.input_version=='ISIMIP3b':
            ggcms_sel = np.sort(co.ggcms['ISIMIP3b_allcrops'])
        else:
            ggcms_sel = np.sort(df.ggcm.unique())
        if per_ggcm:
            abc = 'abcdefghijklmnopr'
            if len(ggcms_sel) < 10:
                suplotshape = (3,3)
                figsize = (12,10)
            elif len(ggcms_sel) == 10:
                suplotshape = (4,3)
                figsize = (12,15)
            else:
                suplotshape = (4,3)
                figsize = (12,15)
            fig = plt.figure(facecolor='w', figsize=(12,10))
            for iggcm, ggcm in enumerate(ggcms_sel):
                subplot_i = (suplotshape[0],suplotshape[1],iggcm+1)
                if iggcm in [0, subplot_i[1], 2*subplot_i[1]]:
                    ylabel_ggcm = ylabel
                else:
                    ylabel_ggcm = None
                ylims = (df[country].min(), df[country].max())
                ax = fig.add_subplot(subplot_i[0], subplot_i[1], subplot_i[2])
                df_sel = df.loc[df.ggcm==ggcm]

                if iggcm == 0:
                    df_combi = df.loc[df.ggcm==ggcm]
                else:
                    df_combi = df_combi.append(df.loc[df.ggcm==ggcm])
                _, _ = plot_ax_binned_impacts(df_sel, ax=ax, fig=fig, subplot_i=subplot_i,
                                               hlines=hlines, title=f'({abc[iggcm]}) {ggcm}',
                                               marker=marker, figsize=figsize, country=country,
                                               ylabel=ylabel_ggcm, xlabel=None, legend=False,
                                               ylims=ylims, ylim_min=True)
            if save_plot:
                fig.savefig(save_dir / f'imp_scatter_{cr}_{country}_ref_bin_05_GGCMs_detrended_{detrended}.{save_plot}', \
                            dpi=300, facecolor='w', edgecolor='w', \
                            orientation='portrait', papertype=None, format=save_plot, \
                            transparent=True, bbox_inches=None, pad_inches=0.1, \
                            frameon=None, metadata=None)
        # combi all plots:
        subplot_i = (1,1,1)
        if (not 'combi' in cr) or (not co.input_version=='ISIMIP3b'):
            df_combi = df
        fig, ax = plot_ax_binned_impacts(df_combi, fig=None, subplot_i=subplot_i,
                                         hlines=hlines,
                                         marker=marker, figsize=figsize, country=country,
                                         ylabel=ylabel, xlabel=xlabel)
        if save_plot:
            fig.savefig(save_dir / f'imp_scatter_{cr}_{country}_ref_bin_05_0ALLGGCM_detrended_{detrended}.{save_plot}', \
                                dpi=300, facecolor='w', edgecolor='w', \
                                orientation='portrait', papertype=None, format=save_plot, \
                                transparent=True, bbox_inches=None, pad_inches=0.1, \
                                frameon=None, metadata=None)

    return fig, ax


def plot_lines_with_intervals(x, y_lines, y_fill=None, y_ens=None, colors=None, linestyles=None,
                              labels=None, title=None, xlabel=None, ylabel=None, fill=True,
                              alpha_fill=.15, alpha_ens=.22, linestyles_interval=None,
                              figsize=None, x_type=int, log_x=True, log_y=False,
                              xticks=None, ylims=None, yticks=None, grid_on=False,
                              in_labels=None,
                              in_label_keys=None, hlines=None, label_str=None,
                              ylabel_right=False):
    """Generic plot of lines with intervals
    
    Parameters:
        x (list or array)
        y_lines (list of lists or arrays): inner lists need to have same length as x
        y_fill(list of lists of lists or arrays): most inner lists need to have same length as x
        y_ens(list of lists of lists or arrays): most inner lists need to have same length as x
        colors (list): same length as number of lists or arrays in y_lines
    """
    if label_str is None:
        label_str = ['%f']
    elif isinstance(label_str, str):
        label_str = [label_str]
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
        if len(label_str) == len(y_lines):
            label_str_single = label_str[i]
        elif len(label_str) == 1:
            label_str_single = label_str[0]
        else:
            label_str_single = '%f'

        if 'perc' in label_str_single:
            label_ = label_str_single % (100*labels[i])
        elif 'mean' in label_str_single:
            label_ = label_str_single % (1+labels[i])
        else:
            label_ = label_str_single % (labels[i])
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
                        linewidth = 1,
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
    if isinstance(yticks, list) or isinstance(yticks, np.ndarray):
        ax.set_yticks(np.array(yticks))

    ax.grid(grid_on)
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(labeltop=False, labelright=ylabel_right, right=True, left=True, labelleft=True)

    ax.set_title(title)
    return fig, ax


def plot_PR_change_per_bin_country(input_dir=co.out_dir / 'Stats_countries' / 'overlapping',
                                   metric = None,
                                   crop=None, countries=None, percentiles=None,
                                   mean_deviations = None,
                                   remove_small_samples=True, log_x=False, log_y=False,
                                   gmt_bins=None, ref_bin=0.5, max_bin=4.25,
                                   ggcms=None, in_labels=True, colors=None,
                                   x_ticks=None,
                                   ylims=(-1, 21), yticks=[0, 1, 2, 5, 10, 15, 20],
                                   save_dir=co.plot_dir, save_plot='pdf',
                                   ylabel_right=True):
    abc = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']

    if metric is None: 
        metric = 'PR' # 'PR': probability ratio, 'PP': probability in percentage, 'PD': probability decimal
    if save_plot and (not save_dir):
        save_dir = co.plot_dir
    if percentiles is None:
        percentiles = [1/40]
    if mean_deviations is None:
        mean_deviations = [-0.1]
    labels = percentiles + mean_deviations
    if isinstance(crop, str):
        crop = [crop]
    elif crop is None:
        crop = co.crop_types + ['combi-kcal']
    if ggcms is None:
        if co.input_version=='ISIMIP3b':
            ggcms_sel = np.sort(co.ggcms['ISIMIP3b_allcrops'])
        else:
            ggcms_sel = co.ggcms[co.input_version]
        ggcms = co.ggcms[co.input_version]
    ggcm_abc = dict()
    for iggcm, ggcm in enumerate(ggcms): # ggcm
        ggcm = ggcm.lower()
        ggcm_abc[abc[iggcm]] = ggcm
    if countries is None:
        countries = co.in_dir / 'main_producer.csv'
    if (isinstance(countries, os.PathLike) or isinstance(countries, str)) and os.path.isfile(countries):
        main_producers = pd.read_csv(countries, encoding="ISO-8859-1", header=0)
        countries = []
        for cr in crop:
            if not 'kcal' in cr:
                countries += list(main_producers.loc[main_producers.crop==cr, 'country'])#[:5]
        countries = ['GLB'] + list(OrderedDict.fromkeys(countries))
        print(countries)
    elif isinstance(countries, str):
        countries = [countries]
    print(crop)
    if metric == 'PR' and log_y:
        hlines=[0.]
    elif metric == 'PR':
        hlines=[0., 1.]
    else:
        hlines = None
    if colors is None:
        colors = [co.colors8_heating[3], co.colors8_heating[0],  co.colors8_heating[2]]

    if percentiles:
        stats_rel2bin_df = [pd.read_csv(input_dir/ f'stats_country_binned_rel2bin_{10*ref_bin:02.0f}_{cr}.csv', # stats_country_binned_rel2bin_05_mai.csv
                             encoding="ISO-8859-1", header=0) for cr in crop]
    else:
        stats_rel2bin_df = None
    if mean_deviations:
        stats_self_df = [pd.read_csv(input_dir/ f'stats_country_binned_self_{cr}.csv', # stats_country_binned_self_whe.csv
                             encoding="ISO-8859-1", header=0) for cr in crop]
    else:
        stats_self_df = None
    if gmt_bins is None:
        gmt_bins = np.sort(stats_rel2bin_df[0]['bin'].unique())
        gmt_bins = gmt_bins[gmt_bins >= ref_bin]
        gmt_bins = gmt_bins[gmt_bins <= max_bin]
    if x_ticks is None:
        x_ticks = np.arange(ref_bin, max_bin+.001, .5)
    x_array = gmt_bins
    fig_pr = list()
    for country in countries: # loop over countries
        print(country)
        # loop over crops and ggcms:
        for i_crop, cr in enumerate(crop): # crop
            if 'combi' in cr:
                ggcms_loop = ggcms_sel
            else:
                ggcms_loop = ggcms
            print(crop[i_crop])
            # init dicts and lists for loop over ggcms:
            values = list()
            means = list()
            medians = list()
            min_max = list()
            iqr = list()
            label_str=list()
            linestyles = list()
            for pc in percentiles: # loop over percentiles / frequencies
                label_str.append('below %2.1fth perc.')
                linestyles.append('-')
                df = stats_rel2bin_df[i_crop]
                values.append(list())
                pc_stack = dict()
                i_pc = 0
                in_label_keys = list()
                for ggcm in ggcms_loop: # ggcm
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
                        elif metric == 'PR': 
                            values[-1][-1].append( # PR probability ratio: percentile divided by ref. pc
                                df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{pc}'), country].values[0] / pc
                                )
                        elif metric == 'PP': 
                            values[-1][-1].append( # PR probability ratio: percentage * 100
                                100 * df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{pc}'), country].values[0]
                                )
                        elif metric == 'PD': 
                            values[-1][-1].append( # PR probability ratio: percentage * 100
                                df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{pc}'), country].values[0]
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

            for md in mean_deviations: # loop over mean - deviations
                label_str.append('below %1.1f*mean')
                linestyles.append('--')
                if isinstance(md, float):
                    md = "('count_rel_mean', %1.1f)" % (md)
                df = stats_self_df[i_crop]
                values.append(list())
                pc_stack = dict()
                i_pc = 0
                in_label_keys = list()
                for ggcm in ggcms_loop: # ggcm
                    ggcm = ggcm.lower()
                    if df.loc[(df.ggcm==ggcm) & (df.stat==f'{md}'), country].size==0:
                        continue
                    values[-1].append(list())
                    if in_labels:
                        in_label_keys.append(list(ggcm_abc.keys())[list(ggcm_abc.values()).index(ggcm)])
                    for gmt_bin in gmt_bins: # fill list with value per gmt_bin:
                        if df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{md}'), country].size==0:
                            continue
                        sample_size = df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{md}'), 'N_years'].values[0]
                        if remove_small_samples and sample_size < 1/percentiles[0]: # sample size lower than first return period? --> NaN
                            values[-1][-1].append(np.nan)
                        elif metric == 'PR': 
                            values[-1][-1].append( # PR probability ratio: frequency of bin divided by reference bin frequency
                                df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{md}'), country].values[0] / \
                                    df.loc[(df.ggcm==ggcm) & (df.bin==ref_bin) & (df.stat==f'{md}'), country].values[0]
                                )
                        elif metric == 'PP': 
                            values[-1][-1].append( # PR probability in percentage: percentage * 100
                                100 * df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{md}'), country].values[0]
                                )
                        elif metric == 'PD': 
                            values[-1][-1].append( # PR probability (decimal): percentage
                                df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat==f'{md}'), country].values[0]
                                )
                    if i_pc==0:
                        pc_stack[md] = np.array(values[-1][-1])
                    else:
                        pc_stack[md] = np.vstack((pc_stack[md], np.array(values[-1][-1])))
                    i_pc += 1

                means.append(np.nanmean(pc_stack[md], axis=0))
                medians.append(np.nanmedian(pc_stack[md], axis=0))
                min_max.append((np.nanmin(pc_stack[md], axis=0), np.nanmax(pc_stack[md], axis=0)))
                iqr.append((np.nanquantile(pc_stack[md], .25, axis=0), np.nanquantile(pc_stack[md], .75, axis=0)))

            if not in_labels:
                ggcm_abc = False
            y_labels = {'PR': 'PR',
                        'PP': 'Probability [%]',
                        'PD': 'Probability'}
            fig_pr.append(
                plot_lines_with_intervals(x_array, medians, y_fill=iqr, y_ens=values, labels=labels,
                                          fill=True, x_type=float, colors=colors, linestyles=linestyles,
                                          xlabel='Global mean warming (°C)',
                                          ylabel=y_labels[metric],
                                          title=f'{country} - {cr}',
                                          log_x=log_x, log_y=log_y,
                                          xticks=x_ticks, ylims=ylims, yticks=yticks,
                                          alpha_ens=.2, in_labels=ggcm_abc, in_label_keys=in_label_keys,
                                          hlines=hlines, label_str=label_str,
                                          ylabel_right=ylabel_right)[0]
                )
            if save_plot:
                fig_pr[-1].savefig(save_dir / f'Prob_{metric}_{cr}_{country}_ref_bin_{ref_bin:0.2f}_ggcmlabels_{in_labels}.{save_plot}', \
                            dpi=300, facecolor='w', edgecolor='w', \
                            orientation='portrait', papertype=None, format=save_plot, \
                            transparent=True, bbox_inches=None, pad_inches=0.1, \
                            frameon=None, metadata=None)
    return fig_pr, ggcm_abc



def plot_CV_per_bin_country(input_dir=co.out_dir / 'Stats_countries' / 'overlapping',
                                   metric = None,
                                   crop=None, countries=None,
                                   remove_small_samples=30, log_x=False, log_y=False,
                                   gmt_bins=None, ref_bin=0.5, max_bin=4.25,
                                   ggcms=None, in_labels=True, colors=None,
                                   x_ticks=None,
                                   ylims=(0, .2), yticks=[0, .05, .1, .15, .2],
                                   save_dir=co.plot_dir, save_plot='pdf',
                                   ylabel_right=True):
    """Coefficient of Variation"""
    abc = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']

    if metric is None: 
        metric = 'CV' # 'Coefficient of Variation'
    if save_plot and (not save_dir):
        save_dir = co.plot_dir
    labels = [metric]
    if ggcms is None:
        if co.input_version=='ISIMIP3b':
            ggcms_sel = np.sort(co.ggcms['ISIMIP3b_allcrops'])
        else:
            ggcms_sel = co.ggcms[co.input_version]
        ggcms = co.ggcms[co.input_version]
    ggcm_abc = dict()
    for iggcm, ggcm in enumerate(ggcms): # ggcm
        ggcm = ggcm.lower()
        ggcm_abc[abc[iggcm]] = ggcm
    if isinstance(crop, str):
        crop = [crop]
    elif crop is None:
        crop = co.crop_types + ['combi-kcal']
    if countries is None:
        countries = co.in_dir / 'main_producer.csv'
    if (isinstance(countries, os.PathLike) or isinstance(countries, str)) and os.path.isfile(countries):
        main_producers = pd.read_csv(countries, encoding="ISO-8859-1", header=0)
        countries = []
        for cr in crop:
            if not 'kcal' in cr:
                countries += list(main_producers.loc[main_producers.crop==cr, 'country'])#[:5]
        countries = ['GLB'] + list(OrderedDict.fromkeys(countries))
        print(countries)
    elif isinstance(countries, str):
        countries = [countries]
    print(crop)
    hlines=[0.]
    if colors is None:
        colors = [co.colors8_heating[0], co.colors8_heating[2],  co.colors8_heating[3]]

    stats_self_df = [pd.read_csv(input_dir/ f'stats_country_binned_self_{cr}.csv', # stats_country_binned_self_whe.csv
                             encoding="ISO-8859-1", header=0) for cr in crop]

    if gmt_bins is None:
        gmt_bins = np.sort(stats_self_df[0]['bin'].unique())
        gmt_bins = gmt_bins[gmt_bins >= ref_bin]
        gmt_bins = gmt_bins[gmt_bins <= max_bin]
    if x_ticks is None:
        x_ticks = np.arange(ref_bin, max_bin+.001, .5)
    x_array = gmt_bins
    fig_cv = list()
    for country in countries: # loop over countries
        print(country)
        # loop over crops and ggcms:
        for i_crop, cr in enumerate(crop): # crop
            if 'combi' in cr:
                ggcms_loop = ggcms_sel
            else:
                ggcms_loop = ggcms
            print(crop[i_crop])
            # init dicts and lists for loop over ggcms:
            values = list()
            means = list()
            medians = list()
            min_max = list()
            iqr = list()
            label_str=list()
            linestyles = list()

            label_str.append('%s')
            linestyles.append('-')
            df = stats_self_df[i_crop]
            # values.append(list())
            pc_stack = dict()
            i_ggcm = 0
            in_label_keys = list()
            for ggcm in ggcms_loop: # ggcm
                ggcm = ggcm.lower()
                if df.loc[(df.ggcm==ggcm) & (df.stat=='mean'), country].size==0:
                    continue
                values.append(list())
                if in_labels:
                    in_label_keys.append(list(ggcm_abc.keys())[list(ggcm_abc.values()).index(ggcm)])
                for gmt_bin in gmt_bins: # fill list with value per gmt_bin:
                    if df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='mean'), country].size==0:
                        continue
                    sample_size = df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='mean'), 'N_years'].values[0]
                    if remove_small_samples and sample_size < remove_small_samples: # sample size lower than return period? --> NaN
                        values[-1].append(np.nan)
                    elif metric == 'CV' and df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='mean'), country].values[0] > 0: 
                        values[-1].append( # std divided by mean
                            df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='std'), country].values[0] / df.loc[(df.ggcm==ggcm) & (df.bin==gmt_bin) & (df.stat=='mean'), country].values[0]
                            )
                    else:
                        continue

                    
                if i_ggcm==0:
                    pc_stack[metric] = np.array(values[-1])
                else:
                    pc_stack[metric] = np.vstack((pc_stack[metric], np.array(values[-1])))
                i_ggcm += 1

            means.append(np.nanmean(pc_stack[metric], axis=0))
            medians.append(np.nanmedian(pc_stack[metric], axis=0))
            min_max.append((np.nanmin(pc_stack[metric], axis=0), np.nanmax(pc_stack[metric], axis=0)))
            iqr.append((np.nanquantile(pc_stack[metric], .25, axis=0), np.nanquantile(pc_stack[metric], .75, axis=0)))

            print(values)
            print(medians)
            print(iqr)

            if not in_labels:
                ggcm_abc = False
            fig_cv.append(
                plot_lines_with_intervals(x_array, medians, y_fill=iqr, y_ens=[values], labels=labels,
                                          fill=True, x_type=float, colors=colors, linestyles=linestyles,
                                          xlabel='Global mean warming (°C)',
                                          ylabel='CV',
                                          title=f'{country} - {cr}',
                                          log_x=log_x, log_y=log_y,
                                          xticks=x_ticks, ylims=ylims, yticks=yticks,
                                          alpha_ens=.2, in_labels=ggcm_abc, in_label_keys=in_label_keys,
                                          hlines=hlines, label_str=label_str,
                                          ylabel_right=ylabel_right)[0]
                )
            if save_plot:
                fig_cv[-1].savefig(save_dir / f'{metric}_{cr}_{country}_ref_bin_{ref_bin:0.2f}_ggcmlabels_{in_labels}.{save_plot}', \
                            dpi=300, facecolor='w', edgecolor='w', \
                            orientation='portrait', papertype=None, format=save_plot, \
                            transparent=True, bbox_inches=None, pad_inches=0.1, \
                            frameon=None, metadata=None)
    return fig_cv, ggcm_abc


def plot_RP_per_bin_country(input_dir=co.out_dir / 'Stats_countries' / 'overlapping',
                            crop=None, countries=None,
                            frequencies=None, remove_small_samples=True, log_x=True,
                            log_y=False, rp=True, gmt_bins=None, ylims=None,
                            ggcms=None, deviation_from=None, save_dir=co.plot_dir, save_plot='pdf',
                            rel2ref=False, plot_ens=False, in_labels=False,
                            title_str=None):
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
    if ggcms is None:
        if co.input_version=='ISIMIP3b':
            ggcms_sel = np.sort(co.ggcms['ISIMIP3b_allcrops'])
        else:
            ggcms_sel = co.ggcms[co.input_version]
        ggcms = co.ggcms[co.input_version]
    ggcm_abc = dict()
    for iggcm, ggcm in enumerate(ggcms): # ggcm
        ggcm = ggcm.lower()
        ggcm_abc[abc[iggcm]] = ggcm
    if isinstance(crop, str):
        crop = [crop]
    elif crop is None:
        crop = co.crop_types + ['combi-kcal']
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
        for i_crop, cr in enumerate(crop): # crop
            df = stats_df[i_crop]
            if 'combi' in cr:
                ggcms_loop = ggcms_sel
            else:
                ggcms_loop = ggcms
            print(crop[i_crop])
            # init dicts and lists for loop over ggcms:
            values = list() # values[i_bin][i_ggcm][i_frequency]
            means = list()
            medians = list()
            min_max = list()
            iqr = list()
            twothirds = list()
            ggcm_stack = dict()
            for gmt_bin in gmt_bins: # loop over GMT bins
                values.append(list())
                iggcm = 0
                in_label_keys = list()
                for ggcm in ggcms_loop: # ggcm
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

            for gmt_bin in gmt_bins: # loop over GMT bins
                means.append(np.mean(ggcm_stack[gmt_bin], axis=0))
                medians.append(np.median(ggcm_stack[gmt_bin], axis=0))
                min_max.append((np.min(ggcm_stack[gmt_bin], axis=0), np.max(ggcm_stack[gmt_bin], axis=0)))
                iqr.append((np.quantile(ggcm_stack[gmt_bin], .25, axis=0), np.quantile(ggcm_stack[gmt_bin], .75, axis=0)))
                twothirds.append((np.quantile(ggcm_stack[gmt_bin], 1/6, axis=0), np.quantile(ggcm_stack[gmt_bin], 5/6, axis=0)))

            # fig_rp.append(
            #     plot_lines_with_intervals(x_array, means, min_max, labels=gmt_bins, fill=True, x_type=float, colors=colors, 
            #                               title=f'{country} - {crop[i_crop]}')
            #     )

            if not in_labels:
                ggcm_abc = False
            if not plot_ens:
                values=None
            if title_str is None:
                title = f'{country} - {crop[i_crop]}'
            else:
                title = title_str
            fig_rp.append(
                plot_lines_with_intervals(x_array, medians, y_fill=iqr, # y_fill=iqr, # y_fill=twothirds
                                          y_ens=values, labels=gmt_bins, fill=True, x_type=int, colors=colors, 
                                          title=title,
                                          xlabel='Return period (years)',
                                          ylabel='$\Delta$ Crop production (%)',
                                          log_x=log_x, log_y=log_y, ylims=ylims,
                                          alpha_ens=.2, in_labels=ggcm_abc,
                                          in_label_keys=in_label_keys,
                                          hlines=[0], label_str='%1.1f°C')[0]
                )
            if save_plot:
                if rel2ref:
                    fig_file_str = f'ReturnP_{crop[i_crop]}_{country}_ggcmlabels_{in_labels}_rel2ref_ens_{plot_ens}.{save_plot}'
                else:
                    fig_file_str = f'ReturnP_{crop[i_crop]}_{country}_ggcmlabels_{in_labels}_ens_{plot_ens}.{save_plot}'
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
    if ggcms is None:
        if co.input_version=='ISIMIP3b':
            ggcms_sel = np.sort(co.ggcms['ISIMIP3b_allcrops'])
        else:
            ggcms_sel = co.ggcms[co.input_version]
        ggcms = co.ggcms[co.input_version]
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
    for icr, df in enumerate(stats_df): # crop
        if 'combi' in crop[icr]:
            ggcms_loop = ggcms_sel
        else:
            ggcms_loop = ggcms
        fig_efc.append(plt.figure(facecolor='w'))
        ax_efc.append(fig_efc[-1].add_subplot(1,1,1))
        # init dicts and lists for loop over ggcms:
        values = dict()
        values[gmt_bin] = dict()
        values[ref_bin] = dict()
        ggcm_stack = dict()
        for iggcm, ggcm in enumerate(ggcms_loop): # ggcm
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

    if co.input_version=='ISIMIP3b' and 'combi' in crop:
        ggcms_sel = np.sort(co.ggcms['ISIMIP3b_allcrops'])
    else:
        ggcms_sel = stats_df.ggcm.unique()

    for ggcm in np.sort(ggcms_sel):
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
