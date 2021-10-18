#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 14:38:17 2021

@author: eberenzs
"""
from pathlib import Path
import shapely
from climada.hazard.relative_cropyield import RelativeCropyield



haz_dict = dict()

haz = RelativeCropyield()
var = 'ric-firr'
filename = f'gepic_gfdl-esm2m_ewembi_historical_2005soc_co2_yield-{var}_global_annual_1861_2005.nc'
haz.set_raster([str(Path(Path("/Users/eberenzs/Documents/Projects/crop_analysis/data/ISIMIP2b/input/Hazard") / var, filename))],
               band=[1,2], geometry=list([shapely.geometry.box(0, 0, 90, 70)]))
haz.plot_intensity(event=0)

haz2 = RelativeCropyield()
var = 'ric-noirr'
filename = f'gepic_gfdl-esm2m_ewembi_historical_2005soc_co2_yield-{var}_global_annual_1861_2005.nc'
haz2.set_raster([str(Path(Path("/Users/eberenzs/Documents/Projects/crop_analysis/data/ISIMIP2b/input/Hazard") / var, filename))],
               band=[1,2], geometry=list([shapely.geometry.box(0, 0, 90, 70)]))
haz2.plot_intensity(event=0)

haz3 = haz2
haz3.intensity = haz.intensity / haz2.intensity
haz2.plot_intensity(event=0)
"""
for var in ['ric-firr', 'ric-noirr', 'whe-firr', 'whe-noirr']:
    
    haz_dict[var] = RelativeCropyield()
    path = Path("/Users/eberenzs/Documents/Projects/crop_analysis/data/ISIMIP2b/input/Hazard") / var
    fn = f'gepic_gfdl-esm2m_ewembi_historical_2005soc_co2_yield-{var}_global_annual_1861_2005.nc'
    haz_dict[var].set_from_isimip_netcdf(input_dir=path, filename=fn, bbox=None,
                               yearrange=(2001, 2003))
"""