#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:44:37 2020

functions for plotting world maps with crop statistics etc.

@author: eberenzs
"""


import cartopy.crs as ccrs
import cartopy.util as cutil
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import matplotlib.pyplot as plt
import matplotlib.colors as mpl_colors
from matplotlib.colors import from_levels_and_colors
from matplotlib.ticker import LogLocator, LogFormatterSciNotation as LogFormatter

import numpy as np
import pandas as pd

import seaborn as sns
import xarray as xr

from pathlib import Path

import mplotutils as mpu
# https://github.com/mathause/mplotutils
# i.e. pip install https://github.com/mathause/mplotutils/archive/v0.2.0.zip
import isimip3b_crop_config as co


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
    if isinstance(ggcm, str):
        df_self_return = df_self.loc[df_self.ggcm==ggcm]
        df_rel2bin_return = df_rel2bin.loc[df_rel2bin.ggcm==ggcm]
    elif isinstance (ggcm, list) or isinstance(ggcm, np.ndarray):
        for i, crop_model in enumerate(ggcm):
            if i == 0:
                df_self_return = df_self.loc[df_self.ggcm==crop_model]
                df_rel2bin_return = df_rel2bin.loc[df_rel2bin.ggcm==crop_model]
            else:
                df_self_return = df_self_return.append(df_self.loc[df_self.ggcm==crop_model])
                df_rel2bin_return = df_rel2bin_return.append(df_rel2bin.loc[df_rel2bin.ggcm==crop_model])
    else:
        df_self_return = df_self
        df_rel2bin_return = df_rel2bin
    return df_self_return, df_rel2bin_return, fn_list

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

def load_mat_stat_combi(fn_str_end='.csv',
                        fn_str_start='stats',
                        ggcm=co.ggcms[co.input_version][0],
                        crop=co.crop_types[-1],
                        soc=co.soc,
                        co2=co.co2_fert,
                        stat=co.statistics_list_mat_ref[0],
                        gmt_bin=1.5,
                        ref_bin=0.5,
                        directory=None,
                        coord_exp=None,
                        calc_inverse=True
                       ):
    fn_list = list()
    fn_list.append(f'{fn_str_start}_self_{ggcm}_{soc}_{co2}_{crop}_{stat}_bin_{ref_bin:1.1f}{fn_str_end}') # ref_self
    fn_list.append(f'{fn_str_start}_self_{ggcm}_{soc}_{co2}_{crop}_{stat}_bin_{gmt_bin:1.1f}{fn_str_end}') # bin_self
    fn_list.append(f'{fn_str_start}_rel_{ggcm}_{soc}_{co2}_{crop}_{stat}_bin_{ref_bin:1.1f}_ref_{ref_bin:1.1f}{fn_str_end}') # ref_rel
    fn_list.append(f'{fn_str_start}_rel_{ggcm}_{soc}_{co2}_{crop}_{stat}_bin_{gmt_bin:1.1f}_ref_{ref_bin:1.1f}{fn_str_end}') # bin_rel

    ds_dict = dict()
    try:
        ds_dict['ref_self'] = load_mat_stat_single(fn_list[0], directory=directory,
                                                   coord_exp=coord_exp, calc_inverse=calc_inverse)
    except FileNotFoundError:
        ds_dict['ref_self'] = None
    try:
        ds_dict['bin_self'] = load_mat_stat_single(fn_list[1], directory=directory,
                                                   coord_exp=coord_exp, calc_inverse=calc_inverse)
    except FileNotFoundError:
        ds_dict['bin_self'] = None
    try:
        ds_dict['ref_rel']  = load_mat_stat_single(fn_list[2], directory=directory,
                                                   coord_exp=coord_exp, calc_inverse=calc_inverse)
    except FileNotFoundError:
        ds_dict['ref_rel'] = None
    try:
        ds_dict['bin_rel']  = load_mat_stat_single(fn_list[3], directory=directory,
                                                   coord_exp=coord_exp, calc_inverse=calc_inverse)
    except FileNotFoundError:
        ds_dict['bin_rel'] = None
    return ds_dict, fn_list


def get_colormap_pr(levels=None, pal=None, extend=None):
    if not levels:
        levels = np.array([0, 0.5, 0.8, 1., 1.2, 1.5, 2., 2.5, 3, 5, 7, 10, 15, 20, 40]) # 14 intervals
    if not extend:
        extend = 'neither'
    if not pal:
        pal = [
            # '#313695', '#4575b4',
            '#74add1', '#abd9e9', '#e0f3f8', # 3 colors, 0. bis 1.
            '#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026', # 9 colors bis 1. bis 15.
            '#630915', '#8e025c', # 2 colors, 15 bis 40.
        ]
    cmap, norm = from_levels_and_colors(levels, pal, extend=extend)
    return cmap, norm, levels, pal

def get_colormap_pr_thiery(levels=None, pal=None, extend=None):
    if not levels:
        levels = np.array([1/10, 1/8, 1/6, 1/4, 1/2, 1, 2, 4, 6, 8, 10]) # 10 intervals
    if not extend:
        extend = 'both'
    if not pal:
        pal = ['#003572', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0', # 5x blue, darkest: '#053061', 
         '#f7f7f7', '#f7f7f7', # 2x grey
         '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f'] # 5x red"""
        """pal = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', # blue
         '#f7f7f7', # grey
         '#fddbc7', '#f4a582', '#d6604d', '#b2182b'] # red"""
    cmap, norm = from_levels_and_colors(levels, pal, extend=extend)
    tick_labels = ["1/10", "1/8", "1/6", "1/4", "1/2", "1", "2", "4", "6", "8", "10"]
    return cmap, norm, levels, pal, tick_labels

def plot_mat_stat_single_old(ds, figsize=(10,10), title=None, inverse=False, replace_zeros=False, cmap=None):
    if cmap is None: cmap = 'YlOrRd'
    if title is None: title = ''
    f, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=figsize)
    ax.coastlines()
    if replace_zeros:
        # replace all values equal to 0 with np.nan
        ds = ds.where(ds['statistic'] != 0)
    if inverse:
        h = ax.pcolormesh(ds.lon, ds.lat, ds.inverse, cmap=cmap)  
    else:
        h = ax.pcolormesh(ds.lon, ds.lat, ds.statistic, cmap=cmap)
    mpu.colorbar(h, ax, orientation='horizontal', aspect=30)
    ax.set_extent([-180, 180, -62, 78], ccrs.PlateCarree())
    # ax.set_global()
    ax.set_title(title)
    return f, ax, h

#class MF(LogFormatter):
#    def set_locs(self, locs=None):
#        labels=[1,2,5,10]
#        self._sublabels = set(labels)

def plot_mat_stat_single(ds, figsize=None, title=None, levels=None, diverging=False, cmap=None, inverse=False, replace_zeros=False, log=False, stipple=None, hatches=None, bbox=None, cbar_label=None, variable_name='statistic'):
    if cbar_label is None:
        cbar_label = ''
    if bbox is None:
        (lonmin, latmin, lonmax, latmax) = co.bbox
    else:
        (lonmin, latmin, lonmax, latmax) = bbox
    ds = ds.sel(lon=slice(lonmin, lonmax), lat=slice(latmax, latmin))
    if stipple is not None:
        stipple = stipple.sel(lon=slice(lonmin, lonmax), lat=slice(latmax, latmin))
    if hatches is None:
        hatches=['...', '']
    if title is None: title=''
    if cmap is None or (not cmap=='PR'):
        if diverging:
            if cmap is None: cmap='RdBu_r'
            if (levels is None) or levels=='abs_80':
                if inverse:
                    max_ = np.max([np.abs(np.min(ds.inverse.values)), np.abs(np.max(ds.inverse.values))])
                else:
                    max_ = np.max([np.abs(np.min(ds[variable_name].values)), np.abs(np.max(ds[variable_name].values))])
                levels = np.arange(-max_, max_, np.int((max_+max_)/5)) # np.arange(-1., 1., 0.2)
        else:
            if cmap is None: cmap='YlOrRd'
            if levels is None or levels=='abs_80':
                if levels=='abs_80':
                    max_fac = .8
                else:
                    max_fac = 1 
                if inverse:
                    levels = np.arange(ds.inverse.min(), ds.inverse.max()*max_fac, np.int((ds.inverse.max()*max_fac-ds.inverse.min())/10)) # np.arange(0., 1., 0.1)
                else:
                    levels = np.arange(ds[variable_name].min(), ds[variable_name].max()*max_fac, np.int((ds[variable_name].max()*max_fac-ds[variable_name].min())/10)) # np.arange(0., 1., 0.1)  

    if cmap == 'PR':
        cmap, norm, levels, _ = get_colormap_pr(levels=levels)
        ticks = levels
    else:
        ticks=None
        cmap, norm = mpu.from_levels_and_cmap(levels, cmap, extend='both')
    if log: norm = mpl_colors.LogNorm(vmin=levels.min(), vmax=levels.max())
    LON, LAT = mpu.infer_interval_breaks(ds.lon, ds.lat, clip=True)

    f, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=figsize)       
    ax.set_facecolor('w') # '#d9d9d9''#f0f0f0') # paint background (no data) grey
    ax.add_feature(cfeature.LAND, facecolor='#e0e0e0', edgecolor='none', lw=0, zorder=.9)
    ax.coastlines()
    if replace_zeros:
        ds = ds.where(ds[variable_name] != 0) # replace all values equal to 0 with np.nan
    if inverse:
        h = ax.pcolormesh(LON, LAT, ds.inverse, cmap=cmap, norm=norm)
    else:
        h = ax.pcolormesh(LON, LAT, ds[variable_name], cmap=cmap, norm=norm)
    if log:
        cbar = mpu.colorbar(h, ax, orientation='horizontal', aspect=30, norm=norm)#, ticks=levels)
    else:
        cbar = mpu.colorbar(h, ax, orientation='horizontal', aspect=30, ticks=ticks)
    if not (stipple is None):
        levels_stipple = [0, .5, 1]
        h = ax.contourf(stipple.lon, stipple.lat, stipple.statistic, levels=levels_stipple, hatches=hatches, colors='none',)
    cbar.ax.set_xlabel(cbar_label)#, rotation=270)
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    
    # set the ticks
    lon = np.arange(-180, 181, 60)
    lat = np.arange(-50, 55, 25)
    ax.set_xticks(lon, crs=ccrs.PlateCarree());
    ax.set_yticks(lat, crs=ccrs.PlateCarree());

    # format the ticks as e.g 60°W
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    
    
    ax.set_title(title)
    return f, ax, h
    
def xr_concat_stat(objs, dim, agg_stat=None, divide=0, share_conistent=.7):
    dataset = None
    if agg_stat is None: agg_stat='median'
    # print(objs)
    stack_xr = xr.concat(objs, dim=dim)
    # print(stack_xr)
    if agg_stat=='median':
        dataset = stack_xr.median(dim=dim)
    elif agg_stat=='mean':
        dataset = stack_xr.mean(dim=dim)
    elif agg_stat=='std':
        dataset = stack_xr.std(dim=dim)
    elif agg_stat=='min':
        dataset = stack_xr.min(dim=dim)
    elif agg_stat=='max':
        dataset = stack_xr.max(dim=dim)
    elif agg_stat=='sign_consistent':
        positive = stack_xr.where(stack_xr >= divide).count(dim=dim) >= share_conistent * stack_xr.count(dim=dim)
        negative = stack_xr.where(stack_xr < divide).count(dim=dim) >= share_conistent * stack_xr.count(dim=dim)
        positive*1 + negative*1
        dataset = positive*1 + negative*1
    return dataset


def plot_map_of_colored_countries(country_array, data_array, levels=None, colors=None, extent=None, figsize=(10,10),
                                 title=None, cbar_label=None, extend=None, hatches=None, tick_labels=None):
    """
    data_array and country_array need to correspond.
    levels need to have length len(colors)+1
    """
    if colors is None:
        # 12 colors
        colors = ['#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', '#67001f', '#ae017e', '#dd1c77']
    if levels is None:
        # 12 levels
        levels = np.array([0, 5e4, 1e5, 2.5e5, 5e5, 7.5e5, 1e6, 2.5e6, 5e6, 7.5e6, 10e6, 12.5e6])
    if extent is None:
        extent = [-180, 180, -62, 78]
    if extend is None:
        extend = 'max'
 
    cmap, norm = from_levels_and_colors(levels, colors, extend=extend)


    shpfilename = shpreader.natural_earth(resolution='110m',
                                              category='cultural',
                                              name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    countries_nat_earth = reader.records()

    # init world map plot
    fig_map, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()),
                               figsize=figsize)#, constrained_layout=True)
    ax.set_facecolor('w') # '#d9d9d9''#f0f0f0')
    # paint background (no data) grey:
    ax.add_feature(cfeature.LAND, facecolor='#e0e0e0', edgecolor='none', lw=0, zorder=.9)
    # ax.coastlines()
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    countries_included = list()
    countries_excluded = list()
    for country in countries_nat_earth:
        if country.attributes['ADM0_A3'] in country_array:
            idx = country_array.tolist().index(country.attributes['ADM0_A3'])
        elif country.attributes['ISO_N3'] in country_array:
            idx = country_array.tolist().index(country.attributes['ISO_N3'])
            # set color based on position of value per country in levels
        else:
            countries_excluded.append(country.attributes['ADM0_A3'])
            idx = None
            # print('Skipping country: %s' %(country.attributes['ADM0_A3']))
            continue
        if np.isnan(data_array[idx]):
            idx = None
            # print('Skipping country: %s' %(country.attributes['ADM0_A3']))
            continue            
        if extend == 'max':
            if data_array[idx]>=max(levels):
                color = colors[-1]
            else:
                color = colors[(data_array[idx]>=levels).sum()-1]
        elif extend == 'both':
            color = colors[(data_array[idx]>=levels).sum()]
        elif extend == 'neither':
            color = colors[(data_array[idx]>=levels).sum()-1]
        if (isinstance(hatches, list) or isinstance(hatches, np.ndarray)) and hatches[idx]:
            ax.add_geometries([country.geometry], ccrs.PlateCarree(), linewidth=0.3, \
                              edgecolor='#737373', facecolor=color,
                              # label=country.attributes['ADM0_A3'],
                              zorder=3, hatch='...'
                             )
        else:
            ax.add_geometries([country.geometry], ccrs.PlateCarree(), linewidth=0.3, \
                              edgecolor='#737373', facecolor=color,
                              # label=country.attributes['ADM0_A3'],
                              zorder=3,
                             )
        """if (isinstance(hatches, list) or isinstance(hatches, np.ndarray)) and hatches[idx]:
                ax.add_geometries([country.geometry], ccrs.PlateCarree(), linewidth=0, \
                          edgecolor=None, facecolor=None, fill=False,
                          hatch='...', zorder=43,
                         )"""
        countries_included.append(country_array[idx])
    
    # sneak in a colorbar
    nrows = 2
    ncols = 2
    Z = np.arange(nrows * ncols).reshape(nrows, ncols)*1e5
    x = np.arange(ncols + 1)
    y = -np.arange(2 + 1) - extent[-1] -1  # -np.arange(nrows + 1) + extent[-1] 
    h = ax.pcolormesh(x, y, Z, cmap=cmap, norm=norm, zorder=0.5, alpha=1)
    cbar = mpu.colorbar(h, ax, orientation='horizontal', aspect=30, ticks=levels)
    if tick_labels:
        cbar.ax.set_xticklabels(tick_labels)
    cbar.ax.set_xlabel(cbar_label)#, rotation=270)

    # set the ticks
    lon = np.arange(-180, 181, 60)
    lat = np.arange(-50, 55, 25)
    ax.set_xticks(lon, crs=ccrs.PlateCarree());
    ax.set_yticks(lat, crs=ccrs.PlateCarree());

    # format the ticks as e.g 60°W
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    
    ax.set_title(title)
    plt.tight_layout()
    plt.show()
    return fig_map, ax, countries_included, countries_excluded