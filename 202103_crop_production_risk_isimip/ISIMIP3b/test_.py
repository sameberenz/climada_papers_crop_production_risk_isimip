#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 18:08:08 2021

@author: eberenzs
"""
import os
import numpy as np
import pandas as pd
from collections import OrderedDict
import copy

gmt_bins = [0.5, 2, 4]
ggcms_loop = ['m1', 'm2', 'm3']
frequencies = [.5, .2, .1, .01]

datain = [[2,4,6,10], [3,3,9,9], []]
values = list() # values[i_bin][i_ggcm][i_frequency]
means = list()
medians = list()
min_max = list()
iqr = list()
twothirds = list()

ggcm_stack = dict() ## place ouside loop

for gmt_bin in gmt_bins: # loop over GMT bins
    values.append(list())
    
    iggcm = 0
    in_label_keys = list()
    for ggcm in ggcms_loop: # ggcm
        ggcm = ggcm.lower()

        values[-1].append(list())

        for f in frequencies: # fill list with value per frequency:
            
            if ggcm=='m3' and f==0.01 and gmt_bin==4:
                values[-1][-1].append(np.nan)
            else: 
                values[-1][-1].append((gmt_bin/f)*(iggcm+1))

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

for iggcm, gmt_bin in enumerate(gmt_bins): # loop over GMT bins
    if iggcm==0:
        ggcm_stack_nan = ggcm_stack[gmt_bin]
    else:
        ggcm_stack_nan += ggcm_stack[gmt_bin]
for gmt_bin in gmt_bins: # loop over GMT bins
    ggcm_stack[gmt_bin] = ggcm_stack[gmt_bin][:, ~np.isnan(ggcm_stack_nan).any(axis=0)]
    means.append(np.nanmean(ggcm_stack[gmt_bin], axis=0))
    medians.append(np.nanmedian(ggcm_stack[gmt_bin], axis=0))
    min_max.append((np.nanmin(ggcm_stack[gmt_bin], axis=0), np.nanmax(ggcm_stack[gmt_bin], axis=0)))
    iqr.append((np.nanquantile(ggcm_stack[gmt_bin], .25, axis=0), np.nanquantile(ggcm_stack[gmt_bin], .75, axis=0)))
    twothirds.append((np.nanquantile(ggcm_stack[gmt_bin], 1/6, axis=0), np.nanquantile(ggcm_stack[gmt_bin], 5/6, axis=0)))