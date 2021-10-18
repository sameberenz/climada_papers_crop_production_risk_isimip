#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 18:08:01 2020
@author: carmensteinmann

Write a .csv table containing FAO data per country (ISO-3 alpha country code).
"""

import pandas as pd
import os
import numpy as np
from iso3166 import countries as iso_cntry
import climada.util.coordinates as co

data_path = '/Volumes/Steinplatte/ETH/Arbeit/01_WCR/ISIMIP/ISIMIP_2b/Input/Exposure'

FAO_FILE = "FAOSTAT_data_producer_prices.csv"
FAO_FILE2 = "FAOSTAT_data_production_quantity.csv"

filename = os.path.join(data_path, FAO_FILE2)

#convert_FAO_data(filename)

def convert_FAO_data(filename, var='Value'):
    """Write a .csv table containing FAO data per country (ISO-3 alpha country code)
        Parameters:
            filename (str): file to be modified
            var (str): Varible to be extracted from the FAO file (the year, crop type and country
                                                                  are included by default)
    """
    fao = pd.read_csv(filename)
    fao_countrycode = getattr(fao, 'Area Code').values
    fao_crops = fao.Item.values.tolist()
    fao_year = fao.Year.values
    fao_value = getattr(fao, var).values
    
    #convert FAO country code into ISO numeric-3
    fao_country = co.country_faocode2iso(fao_countrycode)
    
    #arrays of all countries, crops and years contained in the input data 
    single_countries = np.unique(fao_country)
    single_crops = np.unique(fao_crops).tolist()
    single_years = np.unique(fao_year)
    
    #country names for column headers
    country_names = list()
    for country_idx, country_name in enumerate(single_countries):
        if (country_name !=0) and (country_name != 891):
            country_names.append(iso_cntry.get(int(country_name)).alpha3)
        else:
            country_names.append(int(country_name))
    
    #all crop and year combinations as the first two data columns
    new_years = np.zeros([len(single_years)*len(single_crops)], dtype=int)
    new_crops = list()
    idx = 0
    for _, crop in enumerate(single_crops):
        for _, year in enumerate(single_years):
            new_crops.append(crop)
            new_years[idx] = year
            idx = idx+1
    
    #FAO data in the new structure
    values = np.empty([len(new_years),len(single_countries)])
    values[:] = np.NaN
    for country_idx, country_name in enumerate(single_countries[1:],1):
        for crop_idx, crop_name in enumerate(single_crops):
            for year_idx, year in enumerate(single_years):
                idx = np.where((np.asarray(fao_country) == country_name) & \
                               (np.asarray(fao_crops) == crop_name) & \
                               (fao_year == year))[0]
                idx_new = np.where((np.asarray(new_crops) == crop_name) & \
                                   (new_years == year))[0]
                if len(idx) >= 1:
                    values[idx_new, country_idx] = fao_value[idx]
                   
    #save data to .csv file
    impact_filename = filename.split('.')[0]+'_converted.csv'
    data = {'Crop': new_crops, 'Year': new_years}
    dataframe1 = pd.DataFrame(data)
    dataframe2 = pd.DataFrame(values, columns=country_names)
    result = pd.concat([dataframe1, dataframe2], axis=1, sort=False)
    result.to_csv(impact_filename)