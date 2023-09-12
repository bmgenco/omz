#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 07:03:33 2023

@author: adamstoer
"""

import os
import pandas as pd
import xarray as xr
import urllib.request as request
#import urllib3.request as request
from contextlib import closing
import shutil
#import matplotlib.pyplot as plt

#root_dir = '/Users/adamstoer/Synced Documents/Data/'
root_dir ='/home/brandon/europawd/bgcargo_workshop/data/'

data_lst_df = pd.read_csv('https://data-argo.ifremer.fr/etc/argo_sprof_index.txt', 
                          skiprows =  8)
data_lst_df['wmo'] = [int(flle_loc.split('/')[1]) for flle_loc in data_lst_df.file.tolist()]

float_in_use = 6901472
sprof_wmo_file = data_lst_df.loc[data_lst_df['wmo']==6901472]['file'].values[0]
sprof_ftp_file = 'ftp://ftp.ifremer.fr/ifremer/argo/dac/' + sprof_wmo_file
sprof_loc_file = root_dir + str(float_in_use) + '_Sprof.nc'


os.chdir(root_dir) # deposit metadata download here
if (os.path.isfile(sprof_loc_file) == False) or (os.path.getsize(sprof_loc_file)==0):
    print('Downloading Sprof File')
    with closing(request.urlopen(sprof_ftp_file)) as r:
        with open(sprof_ftp_file.split('/')[-1], 'wb') as f:
            shutil.copyfileobj(r, f)
            
ds = xr.open_dataset(sprof_loc_file)
df = ds.to_dataframe().reset_index()


"""
# use sergey's code'

# MAKE A MAP
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import proplot as pplt

land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face')
lakes_50m = cfeature.NaturalEarthFeature('physical', 'lakes', '50m',
                                        edgecolor='face')


fig, axes = pplt.subplots([1], proj={1:'pcarree'}, width=6, height = 4,
                                   sharex = False, sharey = False, grid = False)
axes[0].set_adjustable('datalim')
#bounds = larry_df.dissolve(by='SID').geometry.bounds
#axes[0].set_ylim(bounds.miny.values[0]-5,bounds.maxy.values[0]+5)
#axes[0].set_xlim(bounds.minx.values[0]-5,bounds.maxx.values[0]+5)
axes[0].set_facecolor('cyan0')
axes[0].add_feature(land_50m, edgecolor='black',lw = 0.5, alpha = 0.8, facecolor = 'darkolivegreen')
axes[0].add_feature(cfeature.LAKES, edgecolor='black', facecolor = 'cyan0',zorder=14,lw = 0.5)

'''
axes[0].plot(track_df.LON,track_df.LAT, 
             color = 'black', transform = ccrs.PlateCarree(),
             label = 'Hurricane Maria')
'''

axes[0].plot(df['LONGITUDE'], 
             df['LATITUDE'], 
             markersize = 1, marker = 'o', label = 'My Float')

axes[0].legend(facecolor = 'white', facealpha = 1, ncol = 2)
plt.show()
"""