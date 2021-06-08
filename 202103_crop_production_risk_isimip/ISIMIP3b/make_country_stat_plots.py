#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 16:37:07 2020


Make plots per country for
production deviation per return periods (RP),
Probability ratio (PR),
probability (PP),
Coefficient of Variation (CV).


@author: Samuel Eberenz
"""
import numpy as np

import isimip3b_crop_config as co
#import isimip3b_stats_plot as stats_plot
import plot_stats_per_country_utils as stats_plot

# set input dir:
input_dir = co.out_dir / 'Stats_countries'

# set directory to save resulting plots to:
save_dir = co.plot_dir

# select which plot types to make
impact_data_plot = 0
pr_plots_crops = 0 # y : Probability Ratio, x : global warming level (GMT bin)
pp_plots_crops = 0 # y : Probability, x : GMT bin
rp_plots_crops = 1 # y : Delta Production, x : Return Period
cv_plots_crops = 0 # y : Coefficient of variation, x : GMT_bin

# crop type; can also be list. Default: None --> all 
crop = None # ['combi-kcal'] # None # ['ric'] # None # ['ric'] # ['whe'] # None

# countries (list); None --> main producers
countries = None # ['GLB', 'ARG', 'USA'] # None # ['GLB'] # ['CHE'] # None # ['GLB', 'USA', 'CHN']

max_bin = 4.5 # max. GMT bin on x axis

title_str = None # default: None

if impact_data_plot:
    stats_plot.plot_binned_impacts(detrended=False, country='GLB', per_ggcm=True, crop=crop)
    stats_plot.plot_binned_impacts(detrended=True, country='GLB', per_ggcm=True, crop=crop)

if pr_plots_crops:
    # no GGCM labels:
    stats_plot.plot_PR_change_per_bin_country(input_dir=input_dir,
                                       crop=crop, countries=countries, percentiles=None,
                                       remove_small_samples=True, log_x=False, log_y=0,
                                       gmt_bins=None, ref_bin=0.5, max_bin=max_bin,
                                       ggcms=None, save_dir=save_dir, in_labels=False,
                                       ylims=(-1, 16.5), yticks=[0, 1, 2, 4, 6, 8, 10, 12, 14, 16])
     # GGCM labels
    stats_plot.plot_PR_change_per_bin_country(input_dir=input_dir,
                                       crop=crop, countries=countries, percentiles=None,
                                       remove_small_samples=True, log_x=False, log_y=0,
                                       gmt_bins=None, ref_bin=0.5, max_bin=max_bin+.5,
                                       ggcms=None, save_dir=save_dir, in_labels=True,
                                       ylims=(-1, 21))


if pp_plots_crops:
    # no GGCM labels:
    stats_plot.plot_PR_change_per_bin_country(input_dir=input_dir, metric='PP',
                                       crop=crop, countries=countries, percentiles=None,
                                       remove_small_samples=True, log_x=False, log_y=0,
                                       gmt_bins=None, ref_bin=0.5, max_bin=max_bin,
                                       ggcms=None, save_dir=save_dir, in_labels=False,
                                       ylims=(0, 100), yticks=np.arange(0,110,10))
     # GGCM labels
    stats_plot.plot_PR_change_per_bin_country(input_dir=input_dir, metric='PP',
                                       crop=crop, countries=countries, percentiles=None,
                                       remove_small_samples=True, log_x=False, log_y=0,
                                       gmt_bins=None, ref_bin=0.5, max_bin=max_bin+.5,
                                       ggcms=None, save_dir=save_dir, in_labels=True,
                                       ylims=(0, 100), yticks=np.arange(0,110,10))

if rp_plots_crops:
    #frequencies = [1/60, 1/50, 1/40, 1/30, 1/20, 1/10, 1/5, 1/3, 1/2]
    frequencies = [1/100, 1/80, 1/60, 1/50, 1/40, 1/30, 1/20, 1/10, 1/5, 1/3, 1/2]
    stats_plot.plot_RP_per_bin_country(input_dir=input_dir, countries=countries, crop=crop,
                                frequencies=frequencies, remove_small_samples=True, log_x=True,
                                log_y=False, rp=True, gmt_bins=[0.5, 2.0, 4.0],
                                ggcms=None, deviation_from=None, save_dir=save_dir, save_plot='pdf',
                                rel2ref=True, plot_ens=True, in_labels=False,
                                ylims=(-30,30), title_str=title_str)

    stats_plot.plot_RP_per_bin_country(input_dir=input_dir, countries=countries, crop=crop,
                                frequencies=frequencies, remove_small_samples=True, log_x=True,
                                log_y=False, rp=True, gmt_bins=[0.5, 2.0, 4.0],
                                ggcms=None, deviation_from=None, save_dir=save_dir, save_plot='pdf',
                                rel2ref=True, plot_ens=False, in_labels=False,
                                ylims=(-30,30), title_str=title_str)

    stats_plot.plot_RP_per_bin_country(input_dir=input_dir, countries=countries, crop=crop,
                                frequencies=frequencies, remove_small_samples=True, log_x=True,
                                log_y=False, rp=True, gmt_bins=[0.5, 2.0, 4.0],
                                ggcms=None, deviation_from=None, save_dir=save_dir, save_plot='pdf',
                                rel2ref=False, plot_ens=True, in_labels=False,
                                ylims=(-25,2.5), title_str=title_str)

    stats_plot.plot_RP_per_bin_country(input_dir=input_dir, countries=countries, crop=crop,
                                frequencies=frequencies, remove_small_samples=True, log_x=True,
                                log_y=False, rp=True, gmt_bins=[0.5, 2.0, 4.0],
                                ggcms=None, deviation_from=None, save_dir=save_dir, save_plot='pdf',
                                rel2ref=False, plot_ens=False, in_labels=False,
                                ylims=(-25,2.5), title_str=title_str)


if cv_plots_crops:
    stats_plot.plot_CV_per_bin_country(input_dir=input_dir,
                                       crop=crop, countries=countries,
                                       remove_small_samples=30, log_x=False, log_y=0,
                                       gmt_bins=None, ref_bin=0.5, max_bin=max_bin,
                                       ggcms=None, save_dir=save_dir, in_labels=False)
