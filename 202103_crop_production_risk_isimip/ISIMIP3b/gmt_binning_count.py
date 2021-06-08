#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 10:02:10 2021

@author: eberenzs
"""
import isimip3b_gmt_binning as gmt_binning

df, bin_count = gmt_binning.make_delta_gmt_df()

yr_eoc = (2079, 2099)

means = dict()

for gcm in df.gcm.unique():
    means[gcm] = dict()
    print(gcm)
    for scenario in ('ssp126', 'ssp585'):
        
        print(scenario)
        df_tmp = df.loc[(df.gcm==gcm) & (df.experiment==scenario) & (df.year>=yr_eoc[0])]
        means[gcm][scenario] = df_tmp.dGMT.mean()