#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 10:26:17 2020

@author: eberenzs
"""
import numpy as np

x = np.array([[1, 2, 3, 4, 5, 6], [11, 12, 13, 14, 15, 16], [21, 22, 23, 24, 25, 26], [31, 32, 33, 34, 35, 36]])
print(x.shape)
print(x)
names1 = np.array(['A', 'B', 'C', 'D'])
names2 = np.array(['D', 'B', 'A', 'C'])
names_map = [np.where(names2 == name)[0][0] for i, name in enumerate(names1)]
print(names2[names_map])
x2 = x[names_map, :] # Sort the rows and then the columns
print(x2.shape)
print(x2)