#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 15:34:00 2020

@author: eberenzs
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path
import datetime as dt
import h5py


local_path = Path(os.path.dirname(__file__))
DATA_DIR = local_path.parent / 'data' / 'ISIMIP2b_country_impacts' / 'Impacts_normalized'
RES_DIR = local_path.parent / 'results'

COL4 = ['#e41a1c','#377eb8','#4daf4a','#984ea3']

impact_csv_name = "impact_%s_%s_%s_%s.csv" # % (ggcm, gcm, scenario, crop)

