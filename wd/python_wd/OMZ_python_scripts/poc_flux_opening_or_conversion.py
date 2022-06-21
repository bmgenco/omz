#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:51:21 2022

@author: brandon
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:52:42 2022

first attempt at hdf file imported from: http://spg-satdata.ucsd.edu/CC4km/


@author: brandon
"""
import os
import mpl_toolkits

#from mpl_toolkits.basemap import Basemap # no longer compatible
#https://github.com/matplotlib/basemap/issues/494
#import cartopy

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np  

#import pyhdf

#from pyhdf.SD import SD, SDC
from natsort import natsorted

#import netCDF4 as nc
#import re #regular expressions 

import xarray as xr

from datetime import datetime

# file import and slection

# data_d="/home/brandon/vestawd/omz/data/POC_Flux_data_sets/annual_vgpm_NPP_hdf4"
# out_d="/home/brandon/vestawd/omz/data/POC_Flux_data_sets/annual_vgpm_NPP_netcdf"

# data_d="/home/brandon/vestawd/omz/data/POC_Flux_data_sets/annual_vgpm_modis_chl_hdf4"
# out_d="/home/brandon/vestawd/omz/data/POC_Flux_data_sets/annual_vgpm_modis_chl_netcdf"

#data_d="/home/brandon/vestawd/omz/data/POC_Flux_data_sets/annual_vgpm_modis_sst_hdf4"
data_d="/media/brandon/8160add2-78ae-4cdd-a1c5-792e889fb0b6/home/brandon/Desktop/temp"


out_d="/home/brandon/vestawd/omz/data/POC_Flux_data_sets/annual_vgpm_modis_sst_netcdf"



os.chdir(data_d)

files=os.listdir(data_d)
#since number of files = 365 can pick date by indexing this list. need nat_sore
files=natsorted(files)

# adding corrds

# max/min values derived from xyz file
min=(5/60)

lon=np.arange(-179.95833, 179.95836, min)
lat=np.arange(-89.958336, 89.958336, min)

t1= datetime.now()

print(t1)

# saves hdf as nc for use in r
os.chdir(data_d)
for i in range(len(files)):
    file_name=files[i]
    file = xr.open_dataset(file_name, engine='netcdf4')
    file = file.rename_dims({"fakeDim0":"lat", "fakeDim1":"lon"})
    file = file.assign_coords({"lat":lat, "lon":lon})
    size=len(file_name)
    new_file_name=(file_name[:size - 3] +"nc")
    file.to_netcdf(path=out_d+"/"+new_file_name)
   
t2= datetime.now()
print(t2)

print("elapsed time", t2-t1)




#print(hdf.datasets())
#print(hdf.attributes())

#data=hdf.select('EF_2018_175')

#lat = hdf.select('fakeDim0')
#latitude = lat[:,:]
#lon = hdf.select('fakeDim1')
#longitude = lon[:,:]

#hdf.end()

#m.pcolormesh(x, y, data)


