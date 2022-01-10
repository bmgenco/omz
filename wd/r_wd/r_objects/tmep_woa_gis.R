
# https://pjbartlein.github.io/REarthSysSci/netCDF.html#get-coordinate-including-time-variables
# micromoles_per_kilogram ->(µmol/kg)

# mole per litre (mol* L−1) sometimes denoted by M.
# https://www.nodc.noaa.gov/OC5/WOD/wod18-notes.html

# 1 ml/l of O2 is approximately 43.570 µmol/kg 
# (assumes a molar volume of O2 of 22.392 l/mole and a constant seawater potential density of 1025 kg/m3).
# The conversion of units on a per-volume 
#(e.g., per liter) to a per-mass (e.g., per kilogram) basis assumes a constant seawater potential density of 1025 kg/m3
# so..

value = 20/1.025



library("ncdf4")
library('raster')
library("sf")


wd<-"/home/brandon/vestawd/omz/wd/r_wd"



robj<-"./r_objects"
fig<-"../../figures"
gis_data<-"../../data/gis_data"
data_d<-"../../data/2018_data"
ocean_color<-"../../data/ocean_color_bud"
bathy_d<-"../../data/bathy" # edit this


ctdata<-"../../data/2018_data/ctd/data"
fdata<-"../../data/2018_data/fluor/get.rvdata.us/cruise/OC1806A/fileset/129998/OC1806A/129998/data" # Cruise flow fluorometer
tdata<-"../../data/2018_data/temp/OC1806A/130001/data" # Cruise flow Temp
wdata<-"../../data/2018_data/anem/OC1806A/130002/data" # Cruise wind/anemometer
gps_data<-"../../data/2018_data/gnss/data" # Cruise gps raw
# 
# setwd(wd)
# setwd(data_d)
# setwd("./ODV_exported")
# woa<-read_tsv("isosurface_data.txt", skip=1) %>% select(., 'Depth [m]  at  Oxygen [ml/l]=0.45 @ Depth [m]  at  Oxygen [ml/l]=0.45=first', 'Longitude', "Latitude")
# names(woa)<-c("depth", "lon", "lat")
# woa$lon<-(360-woa$lon)*-1
# woa.dataframe<-woa
setwd(wd)

# create a loop, function here

x<-nc_open("../../data/WOA/woa18_all_o00_01.nc")
o<-ncvar_get(x, "o_an")
lon<-ncvar_get(x, "lon")
lat<-ncvar_get(x, "lat")
depth<-ncvar_get(x, "depth")
time<-ncvar_get(x, "time")
