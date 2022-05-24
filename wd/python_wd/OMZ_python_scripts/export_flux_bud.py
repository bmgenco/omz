#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:52:42 2022

first attempt at hdf file imported from: http://spg-satdata.ucsd.edu/CC4km/

from that website:
 
'the datasets are on a grid of 540 (width) x 417 (height) 
 with approximately 4000 m step in HDF4 format.  
the upper-left corner (lat, lon) is 45N; -140E; 
the lower-right corner is 30.03597N; -115.5454E
Chla values in each pixel of unsigned byte are log10-scaled and can be 
calculated from the pixel value (PV) as: 
Chl (mg m-3) = 10^(0.015 * PV - 2.0) i.e. 10 to the power of 0.015 * PV - 2.0. 
Pixel values 0 and 1 (black in Fig. 1) and 255 (white in Fig. 1) 
are considered invalid and must be excluded from any statistics.  
PV = 1 is used for coastline and has to be excluded too. 
The annotation (color bar) is written into the dataset. 
When reading with python or Matlab the unsigned byte (uint8) variable is read as
signed byte (int8, values from -127 to 128)
and not as unsigned byte (values from 0 to 255) and 
values over 128 become negative.
A simple fix is to add 256 if the signed pixel value is negative. 
Net Primary Production (NPP) was calculated according to the 
modified VGPM algorithm (Kahru et al. 2009). 
Pixel values are signed 2-byte integers in mg C m-2 day-1 with values
-32767 or 32767 meaning no data.''

to do:
1) convert fake dims to lat  lon


@author: brandon
"""
import os
import mpl_toolkits

#from mpl_toolkits.basemap import Basemap # no longer compatible
#https://github.com/matplotlib/basemap/issues/494
import cartopy

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np  

#import pyhdf

from pyhdf.SD import SD, SDC
from natsort import natsorted

import netCDF4 as nc

#import re #regular expressions 



# file import and slection

data_d='/home/brandon/vestawd/omz/data/POC_Flux_data_sets/4km_derived_daily_export_flux/ef_2018/2018'
os.chdir(data_d)

files=os.listdir(data_d)
#since number of files = 365 can pick date by indexing this list. need nat_sore
files=natsorted(files)

t0=166
t0=t0-1

window=files[(t0-31):(t0+90)] #change window size as needed
del files

#in future use regeular expressions

#lopp here:

i=40
file_name=window[i]

hdf = SD(file_name, SDC.READ)

#print(hdf.datasets())
#print(hdf.attributes())

data=hdf.select('EF_2018_175')

lat = hdf.select('fakeDim0')
latitude = lat[:,:]
lon = hdf.select('fakeDim1')
longitude = lon[:,:]

hdf.end()

m.pcolormesh(x, y, data)
