#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 14:19:23 2021

@author: eberenzs
"""



from pathlib import Path

import numpy as np
import xarray as xr
import pandas as pd
from geopandas import GeoDataFrame

from matplotlib import pyplot as plt

import climada.util.coordinates as u_coord
from climada.entity.exposures.crop_production import CropProduction

BBOX = [-180, -85, 180, 85]  # [Lon min, lat min, lon max, lat max]

data_dir = Path('/Users/eberenzs/Documents/Projects/climada_papers_crop_production_risk_isimip/202103_crop_production_risk_isimip/data/ISIMIP3b')
filename_mirca = 'cultivated_area_MIRCA_GGCMI.nc4'
filepath_mirca = data_dir / 'input' / 'Exposure' / filename_mirca
filename_spam_ray = 'spam_ray_yields.nc4'
filepath_spam_ray = data_dir / 'input' / 'Exposure' / filename_spam_ray

crop_idx = dict()
crop_idx['mirca'] = {'mai': 0,
                     'whe': 1,
                     'soy': 2,
                     'ric': 3,
                     }
crop_idx['spam_ray'] = {'mai': 0,
                        'whe': 1,
                        'soy': 3,
                        'ric': 2,
                         }

exp_mai = CropProduction()
exp_mai.set_from_area_and_yield_nc4('mai', 1, 1,
                                    input_dir=data_dir / 'input' / 'Exposure',
                                    bbox = (-140, 0, -40, 80))
exp_mai.plot_scatter(pop_name=False)
val = exp_mai.gdf.value
exp_mai.gdf['value'] = exp_mai.gdf.region_id
exp_mai.plot_scatter(pop_name=False)

"""
### import MIRCA data:
# MIRCA: Cultivated area [ha] for crop types: 0: maize, 1: wheat, 2: soybean, 3: rice
# Dataset is opened and data within the bbox extends is extracted
data_set = xr.open_dataset(filepath_mirca, decode_times=False)

[lonmin, latmin, lonmax, latmax] = BBOX
area = data_set.sel(lon=slice(lonmin, lonmax), lat=slice(latmax, latmin))
# plot:
area.isel(crop=0)['cultivated area all'].plot.contour()

# The latitude and longitude are set; the region_id is determined
lon, lat = np.meshgrid(area.lon.values, area.lat.values)
gdf = GeoDataFrame()
gdf['latitude'] = lat.flatten()
gdf['longitude'] = lon.flatten()
gdf['region_id'] = u_coord.get_country_code(gdf.latitude, gdf.longitude)

### import SPAM Ray data:
# SPAM/Ray: Yields [t/ha/y] for crop types: 0: maize, 1: wheat, 2: rice, 3: soybean
# Dataset is opened and data within the bbox extends is extracted
data_set = xr.open_dataset(filepath_spam_ray, decode_times=False)

yields = data_set.sel(lon=slice(lonmin, lonmax), lat=slice(latmax, latmin))
# plot:
#yields.isel(crop=4)['yield.tot'].plot.contour()




for key in crop_idx['spam_ray'].keys():
    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    area.isel(crop=crop_idx['mirca'][key])['cultivated area all'].plot.contour(ax=ax)

    ax = fig.add_subplot(2, 1, 2)
    yields.isel(crop=crop_idx['spam_ray'][key])['yield.tot'].plot.contour(ax=ax)
    fig.suptitle(key)
"""


