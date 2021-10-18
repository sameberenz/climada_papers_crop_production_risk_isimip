#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:39:21 2020

@author: eberenzs
"""
import impact_statistics
import crop_config as co


save_dir = co.res_dir / 'efc_plots'

crops = co.crop_types
countries = ['GLB', 'USA', 'FRA', 'MAR', 'KAZ', 'IRN', 'CHN']
deviations = [None, 'mean', 'mean_rel']

for crop in crops:
    for country in countries:
        for deviation in deviations:
            impact_statistics.plot_stats_country_bin_EFC(crop=crop, country=country,
                                                         deviation_from=deviation,
                                                         save_dir=save_dir)
