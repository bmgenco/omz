# to mak kml from carious thingngs



rm(list=ls())
# relative directories
wd<-"/home/brandon/callistowd/omz/wd/r_wd"
robj<-"r_objects"
fig<-"../../figures"
gis_data<-"../../data/gis_data"
data_d<-"../../data/2018_data"
# data_d<-"../../data"
ocean_color<-"../../data/ocean_color_bud"
bathy_d<-"../../data/bathy" # edit this

f.ipak <- function(pkg){
  # loads packages, quietly, given by a vector of package names e.g., pkg<-c("ggplot", "tidyverse")
  # will install  packages listed , and their dependencies, if needed.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, quiet=T, verbose = F)
  sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts=F)
}

# packages<-c("stars", "gridExtra", "cowplot", "sf", "ggspatial","stringr","sp", "rgdal",  "rgeos", "raster","readr" ,"tidyverse", "ggplot2", "lubridate",  "ggthemes", "data.table", "reshape2", "RColorBrewer", "marmap", "extrafont", "oce", "MODIS", "measurements")  # set packages here

# tidyverse: is all of the following packages: ggplot2, dplyr, tidyr, readr, purr, tibble, sringr, & forcats.

packages<-c("sp", "rgdal",  "rgeos", "raster", "readr", "tidyverse", "lubridate",  "ggthemes",  "sf", "cmocean", "ncdf4", "RNetCDF",  "plot3D", "tidync", "devtools", "stars", "ncmeta", "maps", "oce","data.table","raster", "fasterize", "RStoolbox", "scales", "purrr", "marmap", "cowplot")

# "HURDAT",

f.ipak(packages)


setwd(wd)
setwd(robj)

olaf_profiles<-readRDS("olaf_profiles_list.R")

h<-readRDS("hurdat.R")
h$DateTime<-with_tz(h$DateTime, tz="UTC")
storm_name<-"OLAF"
# print(filter(h, Name == storm_name) %>% select(., Key) %>% unique(.))
temp<-filter(h, Name == storm_name & year(DateTime) >= 2021) # Knew year from Mike's email
temp<-filter(h, Key == "EP152021" & year(DateTime) >= 2021)
h.pts<- sf::st_as_sf(temp, coords = c("Lon","Lat"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
conv<-1.94384 # convert to m/s
h.pts$max_sustained_wind_m_s<-h.pts$Wind/conv
olaf.pts<-h.pts
rm(temp, h.pts)

stations_mapped<-readRDS("OC1806A_mapped_ctd_stations.R")
bud_fig1<-readRDS("bud_fig1.R")