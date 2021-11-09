#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 10:31:50 2021

@author: eberenzs
"""

from cartopy.io import shapereader
import geopandas
from shapely.geometry import MultiPolygon

import climada.util.coordinates as u_coord

def get_admin1_shape(country_iso3a, admin1_name):
    """
    returns shape obejct with shape of admin1 region for given country and region

    Parameters
    ----------
    country_iso3a : str
        Country code as ISO3alpha, e.g., 'FRA' for France, 'DEU' for Germany, 'USA', ...
        Find all codes online at https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3
    admin1_name : str
        Name of admin 1 region you want to extract. If it does not match exactly,
        None value is returned and a list of available region names is printed.

    Returns
    -------
    Polygon
        (Multi-)Polygon shape instance for admin 1 region.

    """

    # extract admin1 infos and shapes for country using climada.util.coordinates:
    admin1_info, admin1_shapes = u_coord.get_admin1_info(country_iso3a)
    admin1_info = admin1_info[country_iso3a]
    admin1_shapes = admin1_shapes[country_iso3a]
    # create list of admin 1 region names in given country:
    admin1_names = [record['name'] for record in admin1_info]

    # try to match given admin1_name and return shape:
    for idx, name in enumerate(admin1_names):
        if admin1_names[idx]==admin1_name:
            geoseries = geopandas.GeoSeries(admin1_shapes[idx])
            if geoseries.size == 1:
                return geoseries[0]
            else:
                return MultiPolygon(list(geoseries))
            
    # if no match was found, return None and print list of regions available:
    print(f'Could not match admin1_name {admin1_name} to region names available. Use one of the following names instead:\n')
    print(admin1_names)
    return None

def get_admin1_names(country_iso3a):
    """
    returns list of admin1 regions for given country

    Parameters
    ----------
    country_iso3a : str
        Country code as ISO3alpha, e.g., 'FRA' for France, 'DEU' for Germany, 'USA', ...
        Find all codes online at https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3

    Returns
    -------
    list
        list of strings with region names

    """

    # extract admin1 infos and shapes for country using climada.util.coordinates:
    admin1_info, admin1_shapes = u_coord.get_admin1_info(country_iso3a)
    admin1_info = admin1_info[country_iso3a]
    admin1_shapes = admin1_shapes[country_iso3a]
    # create list of admin 1 region names in given country:
    return [record['name'] for record in admin1_info]

