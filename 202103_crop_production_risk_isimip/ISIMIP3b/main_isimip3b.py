#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main_isimip3b.py
Created on Mon Sep 21 11:29:39 2020
@author: Samuel Eberenz

Main wrapper python script to load analysis configuration and run data crunching and
analysis starting from raw GGCM simulations (netCDF format) and FAO statistics.

Imports and uses functions and variables from tlocal python scripts.
Runs the following steps (each can be switched on /off via boolean variables in isimip3b_crop_config):
    1. isimip3b_crop_config (co):
        setting paths, options, parameters, etc.;
        users can manipulate settings in isimip3b_crop_config.
    2. climada_wrappers (cw, requires CLIMADA repository):
        a) init hazard, hist_mean, and exposure from gridded crop yield simulations
            in input/Hazard folder
        b) calculate full gridded impact sets using CLIMADA
        c) calculate impact per country; detrend over time (optional))
    3. bin country-level impacts per GMT bin and compute crop production statistics
        a) country-level crop production statistics in tons per year
        b) country-level crop production statistics in kcal per year
"""

import os
from scipy import sparse

"""1. import functional modules & config from config.py """
import isimip3b_crop_config as co
if co.init_exp or co.init_haz or co.run_calc_impact_sets or co.load_impact_mats or \
    co.run_calc_country_impacts or co.calc_country_statistics or co.calc_country_statistics_rel2bin or \
    co.save_country_impact_binned or co.calc_country_statistics_kcal:
    print('importing climada wrappers (requires CLIMADA)')
    import isimip3b_climada_wrappers as cw
import isimip3b_gmt_binning as gmt_binning
import isimip3b_impact_statistics as impact_statistics

if co.save_results_table:
    from make_results_table_countries import save_results_table


"""2a. init hazard, hist_mean, and exposure set from all files provided in input folders"""
haz_filenames = None
if co.init_exp or co.init_haz:
    print('\n ')
    print('Loading or initiating hazard and exposures...')

    for crop_irr_combi in co.haz_in_sub_dir:
        if not os.listdir(co.haz_in_sub_dir[crop_irr_combi]):
            continue
        print(f'Init hazard: {crop_irr_combi} ...')
        haz_filenames, _, _, _ = cw.init_hazs_and_exps(input_dir=co.haz_in_sub_dir[crop_irr_combi], init_exp=False, init_haz=co.init_haz)
    if co.init_exp:
        print('Init Exposures ...')
        _, _, exp_filenames, _ = cw.init_hazs_and_exps(input_dir=co.haz_in_sub_dir[crop_irr_combi], init_exp=True, init_haz=False)
# for plotting: exp_sets[0].plot_hexbin(pop_name=False, ignore_zero=True)

"""2b full impact sets: calculate full gridded impact sets using CLIMADA"""
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

        #  add baseline exposure mean to each event:
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
        combined_imp_mats_total, cc_id_list, combi_df, event_name_df_crop_combi = cw.import_imp_mat_combine_crops(
            combi_df, input_dir=co.imp_mat_crop_dir['combi'])
        imp_mat_statistics_binned_crops_combined = impact_statistics.stats_mat_binned_wrapper(
            event_name_df_crop_combi, combined_imp_mats_total, cc_id_list, gmt_bins=co.gmt_bins_stat,
            i_ref_bin=co.gmt_bins_stat_i_ref_bin, crops_combi=True)

elif co.run_calc_impact_sets:
    impact_sets_total = cw.imp_mat_add_mean_cp(imp_mat_filenames, impact_sets) # or imp_mats if co.load_impact_mats
    # imp_mat_combine_irr: sum over irrigation types by 'both_irr_id' and event_name
    impact_sets_total, combi_id, event_name_df = cw.imp_mat_combine_irr(imp_mat_filenames, impact_sets_total, combi_df, event_name_df)
    event_name_df = cw.event_name_df_add_columns(event_name_df)
    event_name_df = gmt_binning.bin_imp_events(event_name_df)
    # at this point we got impact matrices of total crop production per event, combined
    # over GCM, scenario, irr; still differentiated effectively by GGCM and crop (and soc + co2)
    # from within each combined full imp set, group into gmt bins based on year, scenario, GCM:
    event_name_df = gmt_binning.bin_imp_events(event_name_df)
    imp_mat_statistics_binned = impact_statistics.stats_mat_binned_wrapper(event_name_df, impact_sets_total, combi_id, gmt_bins=co.gmt_bins_stat,
                                                                           i_ref_bin=co.gmt_bins_stat_i_ref_bin)
    # result: spare matrix per combined_full_imp and bin + list of event names per group
    # impact_mats_binned = gmt_binning.bin_imp_mat(impact_sets_irr_combi, combi_df)
    # statistics..... reduce dimension to grid for each statistic


"""2c. calculate impact per country ( = yearly crop production)
    Optional: exposure is normalized with FAO crop production statistics
    per country"""
if co.run_calc_country_impacts:
    print('\n ')
    print('Computing crop production impacts per country...')
    for irr in co.irr_type_list:
        if irr is None:
            print('\n combined irrigation:')
            combine_exposure_irr = False
        elif '_fullexp' in irr:
            print('\n ' + irr)
            combine_exposure_irr = True
            irr = irr.split('_')[0]
        else:
            print('\n ' + irr)
            combine_exposure_irr = False

        cw.calc_country_impacts(normalize=co.normalize_exposure, haz_filenames=None, irr=irr,
                                combine_exposure_irr=combine_exposure_irr)

    if co.detrend_country_impacts:
        # detrend country-level impacts over time (2nd order polynomial)
        print('\n ')
        print('Detrending crop production impacts per country...')
        for irr in co.irr_type_list: # [None]: # [None, 'noirr', 'firr', 'noirr_fullexp', 'firr_fullexp']:
            print(f'irr: {irr}')
            if irr is None:
                impact_dir = co.impact_dir
            elif irr == 'noirr':
                impact_dir = co.impact_dir_noirr
            elif irr == 'firr':
                impact_dir = co.impact_dir_firr
            elif irr == 'noirr_fullexp':
                impact_dir = co.impact_dir_noirr_fullexp
            elif irr == 'firr_fullexp':
                impact_dir = co.impact_dir_firr_fullexp
            output_dir = impact_dir / 'Detrended'
            output_dir.mkdir(parents=True, exist_ok=True)
            impact_statistics.polynomial_detrending_impact_multiple(input_dir=impact_dir,
                                                                    output_dir=output_dir,
                                                                    order=2, save_prefix='')
if co.detrend_country_impacts:
    result_subfolders = ['', 'Detrended']
else:
    result_subfolders = ['']


"""3. Bin country-level impacts per GMT bin and compute crop production statistics"""
if co.calc_country_statistics or co.calc_country_statistics_rel2bin or co.calc_country_statistics_kcal:
    print('\n ')
    print('Computing crop production statistics per country...\n')
    print('Input dir:')
    print(co.impact_dir)
    print('output dir:')
    print(co.stats_cntry_dir)

for subdir in result_subfolders: # loop over detrended / not detrended
    print(subdir)

    """3a. TONS: country-level crop production statistics in tons per year"""
    if co.calc_country_statistics or co.calc_country_statistics_rel2bin or co.save_country_impact_binned:
        if subdir == '': # not detrended
            country_impact_combined_nodetrend_dict = dict()
        for irr in co.irr_type_list: # loop over irrigation types.
            # example options: [None], [None, 'noirr', 'firr', 'noirr_fullexp', 'firr_fullexp']
            # Default None combines irrigated and non-irrigated impacts
            print(f'irr: {irr}')
            if subdir == '':
                country_impact_combined_nodetrend_dict[irr] = dict()
            if irr is None:
                impact_dir = co.impact_dir
            elif irr == 'noirr':
                impact_dir = co.impact_dir_noirr
            elif irr == 'firr':
                impact_dir = co.impact_dir_firr
            elif irr == 'noirr_fullexp':
                impact_dir = co.impact_dir_noirr_fullexp
            elif irr == 'firr_fullexp':
                impact_dir = co.impact_dir_firr_fullexp

            for crop in co.crop_types_tons: # loop over crop types
                print(crop)
                if 'combi' in crop:
                    baseline_exp_file_path = cw.country_impact_combine_crops(imp_dir=impact_dir / subdir, crops=None,
                                                                         unit=crop.split('-')[1], baseline_exp_file_path=impact_dir)
                else:
                    baseline_exp_file_path = impact_dir

                print('GMT binning...')
                if co.overlapping: # overlapping GMT bins
                    stats_cntry_dir_tmp = co.stats_cntry_dir / 'overlapping'
                else: # distinct GMT bins
                    stats_cntry_dir_tmp = co.stats_cntry_dir

                if ('etrended' in subdir) and co.detrend_country_impacts_mean_corr:
                    country_impact_combined_nodetrend = country_impact_combined_nodetrend_dict[irr][crop]
                else:
                    country_impact_combined_nodetrend = None
                if co.overlapping: # overlapping GMT bins
                    for i_bin, gmt_bins_ol in enumerate(co.gmt_bins_overlapping):
                        # loop over GMT bin groups (overlapping)
                        print(gmt_bins_ol)
                        if i_bin == 0:
                            country_impact_combined = gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=True,
                                                                                 baseline_exp_file_path=baseline_exp_file_path,
                                                                                 gmt_bins=gmt_bins_ol, drop_nobin=True,
                                                                                 impact_dir=impact_dir / subdir,
                                                                                 country_impact_combi_mean_ref=country_impact_combined_nodetrend)[0]
                        else:
                            country_impact_combined = country_impact_combined.append(gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=True,
                                                                                 baseline_exp_file_path=baseline_exp_file_path,
                                                                                 gmt_bins=gmt_bins_ol, drop_nobin=True,
                                                                                 impact_dir=impact_dir / subdir,
                                                                                 country_impact_combi_mean_ref=country_impact_combined_nodetrend)[0],
                                                                                 ignore_index=True)
                else: # no overlapping GMT bins
                    country_impact_combined = gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=True,
                                                                         baseline_exp_file_path=baseline_exp_file_path, #impact_dir / subdir,
                                                                         impact_dir=impact_dir / subdir,
                                                                         country_impact_combi_mean_ref=country_impact_combined_nodetrend)[0]
                if co.save_country_impact_binned:
                    print('Saving country_impact_combined to')
                    print(co.out_dir / f'country_impact_t_binned_{crop}_{irr}_{subdir}.csv')
                    country_impact_combined.to_csv(co.out_dir / f'country_impact_t_binned_{crop}_{irr}_{subdir}.csv', index=False)
                if subdir == '':
                    country_impact_combined_nodetrend_dict[irr][crop] = country_impact_combined
                if irr:
                    stats_cntry_dir_tmp = stats_cntry_dir_tmp / irr / subdir
                else:
                    stats_cntry_dir_tmp = stats_cntry_dir_tmp / subdir
                stats_cntry_dir_tmp.mkdir(parents=True, exist_ok=True)
                if country_impact_combined.size > 0:
                    if co.calc_country_statistics:
                        print('\n Statistics...\n')
                        # country_stats[crop] =
                        impact_statistics.stats_country_binned_self(country_impact_combined,
                                                                    crop=crop,
                                                                    out_dir=stats_cntry_dir_tmp)
                    if co.calc_country_statistics_rel2bin:
                        print('\n Statistics country rel2bin...\n')
                        # country_stats_rel2bin[crop] =
                        impact_statistics.stats_country_binned_rel2bin(country_impact_combined,
                                                                       crop=crop,
                                                                       out_dir=stats_cntry_dir_tmp)
                else:
                    print('[No data, no statistics]\n')

    """3b. KCAL: country-level crop production statistics in kcal per year"""
    if co.calc_country_statistics_kcal:
        if subdir == '':
            country_impact_combined_kcal_nodetrend_dict = dict()
        for irr in co.irr_type_list: # [None]: #[None, 'noirr_fullexp', 'firr_fullexp', 'noirr', 'firr']:
            print(f'irr: {irr}')
            if subdir == '':
                country_impact_combined_kcal_nodetrend_dict[irr] = dict()
            if irr is None:
                impact_dir = co.impact_dir
            elif irr == 'noirr_fullexp':
                impact_dir = co.impact_dir_noirr_fullexp
            elif irr == 'firr_fullexp':
                impact_dir = co.impact_dir_firr_fullexp
            elif irr == 'noirr':
                impact_dir = co.impact_dir_noirr
            elif irr == 'firr':
                impact_dir = co.impact_dir_firr

            for crop in co.crop_types_kcal: # ['mai-kcal']
                print(crop)
                if crop.split('-')[0] == 'combi':
                    crops_ = None
                elif crop.split('-')[0] in co.crop_types:
                    crops_ = [crop.split('-')[0]]
                else:
                    raise ValueError(f'invalid crop {crop}')
                # impacts are combined and  baseline exposure value are added to impact values.
                baseline_exp_file_path = cw.country_impact_combine_crops(imp_dir=impact_dir / subdir, crops=crops_,
                                                                         unit=crop.split('-')[1],
                                                                         baseline_exp_file_path=impact_dir)
                print(baseline_exp_file_path)
                # impact_dir = impact_dir / subdir # sub directory for detrended impact data
                print('GMT binning...')

                if co.overlapping: # overlapping GMT bins
                    stats_cntry_dir_tmp = co.stats_cntry_dir / 'overlapping'
                else:
                    stats_cntry_dir_tmp = co.stats_cntry_dir

                if ('etrended' in subdir) and co.detrend_country_impacts_mean_corr:
                    country_impact_combined_nodetrend = country_impact_combined_kcal_nodetrend_dict[irr][crop]
                else:
                    country_impact_combined_nodetrend = None
                if co.overlapping: # overlapping GMT bins
                    for i_bin, gmt_bins_ol in enumerate(co.gmt_bins_overlapping):
                        print(gmt_bins_ol)
                        if i_bin == 0:
                            country_impact_combined = gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=True,
                                                                                 baseline_exp_file_path=baseline_exp_file_path,
                                                                                 gmt_bins=gmt_bins_ol, drop_nobin=True,
                                                                                 impact_dir=impact_dir / subdir,
                                                                                 country_impact_combi_mean_ref=country_impact_combined_nodetrend)[0]
                        else:
                            country_impact_combined = country_impact_combined.append(gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=True,
                                                                                 baseline_exp_file_path=baseline_exp_file_path,
                                                                                 gmt_bins=gmt_bins_ol, drop_nobin=True,
                                                                                 impact_dir=impact_dir / subdir,
                                                                                 country_impact_combi_mean_ref=country_impact_combined_nodetrend)[0],
                                                                                     ignore_index=True)
                else: # no overlapping GMT bins
                    country_impact_combined = gmt_binning.bin_country_impact(crop, co2=co.co2_fert, absolute_production=True,
                                                                         baseline_exp_file_path=baseline_exp_file_path,
                                                                         impact_dir=impact_dir / subdir,
                                                                         country_impact_combi_mean_ref=country_impact_combined_nodetrend)[0]
                print(country_impact_combined.shape)
                print(country_impact_combined.bin.unique())
                if subdir == '':
                    country_impact_combined_kcal_nodetrend_dict[irr][crop] = country_impact_combined
                if irr:
                    stats_cntry_dir_tmp = stats_cntry_dir_tmp / irr / subdir
                else:
                    stats_cntry_dir_tmp = stats_cntry_dir_tmp / subdir
                stats_cntry_dir_tmp.mkdir(parents=True, exist_ok=True)
                if country_impact_combined.size > 0:
                    if True:
                        print('\n Statistics...\n')
                        # country_stats[crop] = i
                        impact_statistics.stats_country_binned_self(country_impact_combined,
                                                                    crop=crop,
                                                                    out_dir=stats_cntry_dir_tmp)
                    if True:
                        print('\n Statistics country rel2bin...\n')
                        # country_stats_rel2bin[crop] =
                        impact_statistics.stats_country_binned_rel2bin(country_impact_combined,
                                                                       crop=crop,
                                                                       out_dir=stats_cntry_dir_tmp)
                else:
                    print('[No data, no statistics]\n')

if co.save_results_table:
    if co.overlapping and co.detrend_country_impacts:
        stats_dir = co.stats_cntry_dir / 'overlapping' / 'Detrended'
    elif co.overlapping:
        stats_dir = co.stats_cntry_dir / 'overlapping'
    else:
        stats_dir = co.stats_cntry_dir
    print(stats_dir)

    results_df, baseline_exp_df = save_results_table(ref_bin=co.reference_bin,
                                                     stats_dir=stats_dir)
