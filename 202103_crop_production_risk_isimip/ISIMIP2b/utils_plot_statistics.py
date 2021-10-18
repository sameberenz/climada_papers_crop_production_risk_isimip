#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:44:37 2020

@author: eberenzs
"""


import cartopy.crs as ccrs
import cartopy.util as cutil
import cartopy.feature as cfeature

import matplotlib.pyplot as plt
import matplotlib.colors as mpl_colors
from matplotlib.ticker import LogLocator, LogFormatterSciNotation as LogFormatter
import numpy as np
import pandas as pd

import seaborn as sns
import xarray as xr

from pathlib import Path

import mplotutils as mpu
# https://github.com/mathause/mplotutils
import crop_config as co


def load_mat_stat_single(fn, directory=co.stats_dir, coord_exp=None, replace_zeros=False):

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
                        coord_exp=None
                       ):
    fn_list = list()
    fn_list.append(f'{fn_str_start}_self_{ggcm}_{soc}_{co2}_{crop}_{stat}_bin_{ref_bin:1.1f}{fn_str_end}') # ref_self
    fn_list.append(f'{fn_str_start}_self_{ggcm}_{soc}_{co2}_{crop}_{stat}_bin_{gmt_bin:1.1f}{fn_str_end}') # bin_self
    fn_list.append(f'{fn_str_start}_rel_{ggcm}_{soc}_{co2}_{crop}_{stat}_bin_{ref_bin:1.1f}_ref_{ref_bin:1.1f}{fn_str_end}') # ref_rel
    fn_list.append(f'{fn_str_start}_rel_{ggcm}_{soc}_{co2}_{crop}_{stat}_bin_{gmt_bin:1.1f}_ref_{ref_bin:1.1f}{fn_str_end}') # bin_rel

    ds_dict = dict()
    try:
        ds_dict['ref_self'] = load_mat_stat_single(fn_list[0], directory=directory, coord_exp=coord_exp)
    except FileNotFoundError:
        ds_dict['ref_self'] = None
    try:
        ds_dict['bin_self'] = load_mat_stat_single(fn_list[1], directory=directory, coord_exp=coord_exp)
    except FileNotFoundError:
        ds_dict['bin_self'] = None
    try:
        ds_dict['ref_rel']  = load_mat_stat_single(fn_list[2], directory=directory, coord_exp=coord_exp)
    except FileNotFoundError:
        ds_dict['ref_rel'] = None
    try:
        ds_dict['bin_rel']  = load_mat_stat_single(fn_list[3], directory=directory, coord_exp=coord_exp)
    except FileNotFoundError:
        ds_dict['bin_rel'] = None
    return ds_dict, fn_list


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

def plot_mat_stat_single(ds, figsize=(10,10), title=None, levels=None, diverging=False, cmap=None, inverse=False, replace_zeros=False, log=False, stipple=None):
    if title is None: title=''
    if diverging:
        if cmap is None: cmap='RdBu_r'
        if (levels is None) or levels=='abs_80':
            if inverse:
                max_ = np.max([np.abs(np.min(ds.inverse.values)), np.abs(np.max(ds.inverse.values))])
            else:
                max_ = np.max([np.abs(np.min(ds.statistic.values)), np.abs(np.max(ds.statistic.values))])
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
                levels = np.arange(ds.statistic.min(), ds.statistic.max()*max_fac, np.int((ds.statistic.max()*max_fac-ds.statistic.min())/10)) # np.arange(0., 1., 0.1)  
    f, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=figsize)
    ax.coastlines()
    cmap, norm = mpu.from_levels_and_cmap(levels, cmap, extend='both')
    if log: norm = mpl_colors.LogNorm(vmin=levels.min(), vmax=levels.max())
    LON, LAT = mpu.infer_interval_breaks(ds.lon, ds.lat, clip=True)
    if replace_zeros:
        ds = ds.where(ds['statistic'] != 0) # replace all values equal to 0 with np.nan
    if inverse:
        h = ax.pcolormesh(LON, LAT, ds.inverse, cmap=cmap, norm=norm)
    else:
        h = ax.pcolormesh(LON, LAT, ds.statistic, cmap=cmap, norm=norm)
    if log:
        cbar = mpu.colorbar(h, ax, orientation='horizontal', aspect=30, norm=norm)#, ticks=levels)
    else:
        mpu.colorbar(h, ax, orientation='horizontal', aspect=30)
    if not (stipple is None):
        levels_stipple = [-1, -0.001, 1]
        h = ax.contourf(stipple.lon, stipple.lat, stipple.statistic, levels=levels_stipple, hatches=['...', ''], colors='none',)
    

    ax.set_extent([-180, 180, -62, 78], ccrs.PlateCarree())
    ax.set_title(title)
    return f, ax, h
    
def xr_concat_stat(objs, dim, agg_stat=None):
    if agg_stat is None: agg_stat='mean'
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
        dataset = np.sign(stack_xr.max(dim=dim)*stack_xr.min(dim=dim))
    return dataset