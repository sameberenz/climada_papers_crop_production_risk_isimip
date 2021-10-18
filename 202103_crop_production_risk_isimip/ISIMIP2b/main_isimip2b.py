#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 11:29:39 2020
@author: eberenzs

main_isimip2b.py

calls the following modules:
    1. crop_config (co):
        from config script, get paths, options, parameters, etc.
    2. climada_wrappers (cw):
        a) init hazard, hist_mean, and hazard
        b) calculate impact per country
        c) calculate complete impact sets <-- requires to assign hazard to exposure etc.
    3. impact_postprocess:
        a) sort impact sets by GMT bins
        b) aggregate per country / global
            --> GMT bin / time series
        c) combine impact sets by crop, GGCM, scenario, etc., adjust frequency
    4. impact_statistics:
        a) RP per grid cell
        b) RP per country + global
        --> reference hist and/or same
    5. compare_and_plot:
        a) comparison of statistics between different scnearios, GMT bins, GGCMs, ...
        b) plots / tables:
            i) world maps
            ii) EFC curves
            iii) compariosn plots(?))
            iv) compare imported and direct risk statistics for countriues
    6. risk_io:
        a) imported risk statistics per country

"""

import os
from scipy import sparse

import climada_wrappers as cw
import gmt_binning
import impact_statistics
import crop_config as co
"""1. import functional modules & config from config.py """

"""2a. init hazard, hist_mean, and exposure set from all files provided in input folders"""
haz_filenames = None
if co.init_exp or co.init_haz:
    print('\n ')
    print('Loading or initiating hazard and exposures...')
    for crop_irr_combi in co.haz_in_sub_dir:
        if not os.listdir(co.haz_in_sub_dir[crop_irr_combi]):
            continue
        print(f'{crop_irr_combi} ...')
        haz_filenames, _, exp_filenames, _ = cw.init_hazs_and_exps(input_dir=co.haz_in_sub_dir[crop_irr_combi])
# for plotting: exp_sets[0].plot_hexbin(pop_name=False, ignore_zero=True)

"""2c full impact sets"""
if co.run_calc_impact_sets and not co.load_impact_mats:
    print('\n Computing full impact sets...')
    imp_filenames, impact_sets, haz_df, combi_df, coord_exp, event_name_df = cw.calc_impact_sets(normalize=co.normalize_exposure, haz_filenames=None, ignore_irr=True)
elif co.run_calc_impact_sets and co.load_impact_mats:
    print('\n Computing full impact sets...')
    imp_filenames, _, _, _, _, _ = cw.calc_impact_sets(normalize=co.normalize_exposure, haz_filenames=None, ignore_irr=True, return_all_mats=False)

if co.load_impact_mats:
    print('Load impact matrices and compute statistics:')
    for idc, crop_type in enumerate(co.crop_types):
        print(f'\n {crop_type}... ')
        if not idc:
            imp_mat_filenames, impact_sets, haz_df, combi_df, coord_exp, \
                event_name_df = cw.import_imp_mat(input_dir=co.imp_mat_dir,
                                                  mat_dir = co.imp_mat_crop_dir[crop_type],
                                                  filenames=co.impact_mats_filename_list)
        else:
            imp_mat_filenames, impact_sets, _, combi_df, _, _ = cw.import_imp_mat(input_dir=co.imp_mat_dir,
                                                                           mat_dir = co.imp_mat_crop_dir[crop_type],
                                                                           filenames=co.impact_mats_filename_list)

        #  add cp exposure mean to each event
        impact_sets_total = cw.imp_mat_add_mean_cp(imp_mat_filenames, impact_sets) # or imp_mats if co.load_impact_mats
        # combine: sum irrigation types? via 'both_irr_id' and event_name
        impact_sets_total, combi_id, event_name_df = cw.imp_mat_combine_irr(imp_mat_filenames, impact_sets_total, combi_df, event_name_df)
        event_name_df = cw.event_name_df_add_columns(event_name_df)
        event_name_df = gmt_binning.bin_imp_events(event_name_df)
        imp_mat_statistics_binned = impact_statistics.stats_mat_binned_wrapper(event_name_df, impact_sets_total, combi_id, gmt_bins=co.gmt_bins_stat, i_ref_bin=co.gmt_bins_stat_i_ref_bin)
        if co.combine_crop_mats:
            event_name_df.loc[(event_name_df.crop==crop_type) & (event_name_df.irr=='noirr')].to_csv(co.imp_mat_crop_dir['combi'] / f'event_names_combi_{crop_type}.csv', index=False)
            for i_combi, both_irr_id_ in enumerate(combi_id):
                sparse.save_npz(co.imp_mat_crop_dir['combi'] / f'imp_full_combi_id_{int(both_irr_id_):02.0f}.npz', impact_sets_total[i_combi])

    if co.combine_crop_mats:
        combined_imp_mats_total, cc_id_list, combi_df, event_name_df_crop_combi = cw.import_imp_mat_combine_crops(combi_df, input_dir=co.imp_mat_crop_dir['combi'])
        imp_mat_statistics_binned_crops_combined = impact_statistics.stats_mat_binned_wrapper(event_name_df_crop_combi, combined_imp_mats_total, cc_id_list, gmt_bins=co.gmt_bins_stat, i_ref_bin=co.gmt_bins_stat_i_ref_bin, crops_combi=True)

elif co.run_calc_impact_sets:
    impact_sets_total = cw.imp_mat_add_mean_cp(imp_mat_filenames, impact_sets) # or imp_mats if co.load_impact_mats
    # combine: sum irrigation types? via 'both_irr_id' and event_name
    impact_sets_total, combi_id, event_name_df = cw.imp_mat_combine_irr(imp_mat_filenames, impact_sets_total, combi_df, event_name_df)
    event_name_df = cw.event_name_df_add_columns(event_name_df)
    event_name_df = gmt_binning.bin_imp_events(event_name_df)
    # at this point we got impact matrices of total crop production per event, combined
    # over gcm, scenario, irr; still differentiated effectively by ggcm and crop (and soc + co2)
    # from within each combined full imp set, group into gmt bins based on year, scenario, gcm:
    event_name_df = gmt_binning.bin_imp_events(event_name_df)
    imp_mat_statistics_binned = impact_statistics.stats_mat_binned_wrapper(event_name_df, impact_sets_total, combi_id, gmt_bins=co.gmt_bins_stat, i_ref_bin=co.gmt_bins_stat_i_ref_bin)
    # result: spare matrix per combined_full_imp and bin + list of event names per group
#    impact_mats_binned = gmt_binning.bin_imp_mat(impact_sets_irr_combi, combi_df)
    # statistics..... reduce dimension to grid for each statistic


"""2b. Calculate impact = crop production.
    Optional: exposure is normalized with FAO crop production statistics"""
if co.run_calc_country_impacts:
    print('\n ')
    print('Computing crop production impacts per country...')
    baseline_exposure_per_country = cw.calc_country_impacts(normalize=co.normalize_exposure,
                                                           haz_filenames=haz_filenames)

if co.calc_country_statistics:
    print('\n ')
    print('Computing crop production statistics per country...\n')
    print('Input dir:')
    print(co.impact_dir)
    print('output dir:')
    print(co.stats_cntry_dir)
    country_stats = dict()
    
    #country_impact_combined_dict = dict()
    #impact_in_filenames_dict = dict()
    #country_impact_combined_relp_dict = dict()
    #impact_in_filenames_relp_dict = dict()

    for crop in co.crop_types:
        print(crop)
        print('GMT binning...')
        # country_impact_combined_dict[crop], impact_in_filenames_dict[crop] = gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=True)
        # country_impact_combined_relp_dict[crop], impact_in_filenames_relp_dict[crop] = gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=False)

        country_impact_combined = gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=True)[0]




        if country_impact_combined.size > 0:
            print('Statistics...\n')
            country_stats[crop] = impact_statistics.stats_country_binned_self(country_impact_combined, crop=crop)
            if co.make_plots:
                fig_quantiles_GLB, ax_quantiles_GLB = impact_statistics.stats_country_bin_plots(country_stats[crop], crop, statistics_list=[.025, .25, .5, .75, .975], marker_list=['','','','',''], ls_list=[':', '--', '-', '--', ':'], country='GLB')
                fig_mean10_GLB, ax_mean10_GLB = impact_statistics.stats_country_bin_plots(country_stats[crop], crop, statistics_list=[('count_rel_mean', -.1)], marker_list=[''], ls_list=['-'], ylabel='quantile below', country='GLB')
                figstd_GLB, ax_mean10_GLB = impact_statistics.stats_country_bin_plots(country_stats[crop], crop, statistics_list=['std'], marker_list=[''], ls_list=['-'], ylabel='Standard Deviation [t/y]', country='GLB')
        else:
            print('[No data, no statistics]\n')


"""
# TODO
    # x 1. rename  cw.rename_events(haz_sets, haz_filenames)
    # x 2. append all with same CROP TYPE, IRR TYPE, SOC, SCEN, GGCM (?), but not forcing GCM
    # 0 3. adjust frequency(?)
    # x 4. calc impact
    # > 5. add baaseline exp to impact?
    # > 6. sort impact by bins?
    # 7. statistics
"""