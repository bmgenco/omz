rm(list=ls())

# wd<-"/home/brandon/vestawd/omz/wd/r_wd/r_scripts"
# wd<-"/home/brandon/vestawd/omz/wd/r_wd/"



# wd<-"/home/brandon/europawd/omz/wd/r_wd/"
wd<-"/home/brandon/callistowd/omz/wd/r_wd"

robj<-"r_objects"
setwd(wd)

f.ipak <- function(pkg){
  # loads packages, quietly, given by a vector of package names e.g., pkg<-c("ggplot", "tidyverse")
  # will install  packages listed , and their dependencies, if needed.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, quiet=T, verbose = F)
  sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts=F)
}

packages<-c("sp", "rgdal",  "rgeos", "raster", "readr", "tidyverse", "lubridate",  "ggthemes",  "sf", "cmocean", "ncdf4", "RNetCDF",  "plot3D", "tidync", "devtools", "stars", "ncmeta", "maps", "oce",
            "data.table","raster", "fasterize", "RStoolbox", "scales", "purrr", "marmap", "gsw", "R.utils","Matrix")

f.ipak(packages)

setwd(wd)
# path_code = "r_scripts/submodules/BGC-ARGO_R_WORKSHOP/"

source(paste0(path_code, "initialize_argo.R"))
source(paste0(path_code, "try_download.R"))
source(paste0(path_code, "do_download.R"))
source(paste0(path_code, "download_float.R"))
source(paste0(path_code, "download_multi_floats.R"))
source(paste0(path_code, "check_dir.R"))
source(paste0(path_code, "get_var_name_units.R"))
source(paste0(path_code, "select_profiles.R"))
source(paste0(path_code, "load_float_data.R"))
source(paste0(path_code, "plot_trajectories.R"))
source(paste0(path_code, "get_lon_lat_lims.R"))
source(paste0(path_code, "show_trajectories.R"))
source(paste0(path_code, "do_pause.R"))
source(paste0(path_code, "depth_interp.R"))
source(paste0(path_code, "calc_auxil.R"))
source(paste0(path_code, "get_multi_profile_mean.R"))
source(paste0(path_code, "show_profiles.R"))
source(paste0(path_code, "plot_profiles.R"))
source(paste0(path_code, "show_sections.R"))
source(paste0(path_code, "plot_sections.R"))
source(paste0(path_code, "plot_profiles.R"))
rm(f.ipak, packages, path_code)

setwd(robj)
stations_mapped<-readRDS("OC1806A_mapped_ctd_stations.R")
e<-extent(stations_mapped)

#does alrger extent matter
e[1]=-125.00
e[2]=-100.00
e[3]=10.00
e[4]=30.00

load("init_argo.RData")

h<-readRDS("hurdat.R")
h$DateTime<-with_tz(h$DateTime, tz="UTC")
storm_name<-"OLAF"
print(filter(h, Name == storm_name) %>% select(., Key) %>% unique(.))
temp<-filter(h, Name == storm_name & year(DateTime) >= 2021) # Knew year form Mike's email
h.pts<- sf::st_as_sf(temp, coords = c("Lon","Lat"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
conv<-1.94384 # convert to m/s
h.pts$max_sustained_wind_m_s<-h.pts$Wind/conv
rm(h, temp )

st_date<-min(h.pts$DateTime)
ed_date<-max(h.pts$DateTime)

storm_profiles<-select_profiles(lon_lim=c(e[1],e[2]), lat_lim=c(e[3], e[4]), start_date=date(st_date), end_date=date(ed_date), outside="none", sensor=NULL)

# download_multi_floats(storm_profiles$floats)
# data<-load_float_data(storm_profiles$floats)
# data<-data$Data



## profiles 


### functions

f.profile.select<-function(storm_profiles, days_around){ 
profile.windows<-vector(mode="list", length=length(storm_profiles$floats))

for(i in 1:length(profile.windows)){  
float<-storm_profiles$floats[i] 
float_idx <-which(Float$wmoid==float)
prof_ids = c(Float$prof_idx1[float_idx]:Float$prof_idx2[float_idx])
dates = Sprof$date[prof_ids] 

overpass<-dates[which.min(abs(dates-st_date)):which.min(abs(dates-ed_date))]

window<-c((which.min(abs(dates-st_date))-days_around):(which.min(abs(dates-ed_date))+days_around))
profiles<-prof_ids[window]
TIME<-dates[window]
t<-data.frame(TIME, profiles)
profile.windows[[i]]$data<-t
profile.windows[[i]]$overpass<-overpass
profile.windows[[i]]$window<-window
profile.windows[[i]]$bufer_number<-days_around
names(profile.windows)[i]<-storm_profiles$floats[i]  
}
return(profile.windows)
rm(t)
}
f.spacetime_agro<-function(storm_profiles){
  floats.sp<-vector(mode="list", length=length(storm_profiles$floats))
  for(i in 1:length(floats.sp)){
    d<-load_float_data(storm_profiles$floats[i])
    d<-d$Data[[1]] %>% lapply(., as.vector) %>% as.data.frame() %>% distinct() %>% na.omit()
    d<-d %>% filter(., DIRECTION=="A") %>% st_as_sf(., coords = c("LONGITUDE", "LATITUDE" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    # d<-d %>% st_as_sf(., coords = c("LONGITUDE", "LATITUDE" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    floats.sp[[i]]<-d
    names(floats.sp)[i]<-storm_profiles$floats[i]
  }
  return(floats.sp)
  rm(d)
}
f.subset<-f.subset<-function(profile.windows, floats.sp){
  window.sp<-vector(mode="list", length=length(floats.sp))
  for(i in 1:length(floats.sp)){
    x<-profile.windows[[i]]$data$TIME %>% as_datetime(.)
    y<-floats.sp[[i]]
    y$TIME<-y$TIME %>% as_datetime(.)
    window<-c(which.min(abs(y$TIME-x[1]))):(which.min(abs(y$TIME-x[length(x)])))
    window.sp[[i]]<-y[window,]
    names(window.sp)[i]<-storm_profiles$floats[i]
  }
  return(window.sp)
  
}



profile.windows<-f.profile.select((storm_profiles), days_around =5)
floats.sp<-f.spacetime_agro(storm_profiles)
window.sp<-f.subset(profile.windows, floats.sp)

olaf_profiles<-list(profile.windows, floats.sp, window.sp)
names(olaf_profiles)<-c("window", "floats.sp", "window.sp")

setwd(wd)
# setwd(robj)
# saveRDS(profile.windows, "olaf_argo_profiles_temporal.R")
# saveRDS(floats.sp, "olaf_argo_floats_spatial.R")
# saveRDS(floats.sp, "olaf_argo_floats_spatial.R")

# saveRDS(olaf_profiles, "olaf_profiles_list.R")

## bah
## updates 2023-07-31
setwd(wd)
setwd(robj)
olaf_profiles<-readRDS("olaf_profiles_list.R")


x<-olaf_profiles$window.sp
z<-x$`6903093` %>% filter(., CYCLE_NUMBER %in% c(48, 49, 50, 51, 52))
e<-extent(z)
q<-select_profiles(lon_lim=c(e[1],e[2]), lat_lim=c(e[3], e[4]), start_date=(z$TIME[1]-1), end_date=(z$TIME[length(z$TIME)]+1), outside="none", sensor=NULL)
d<-load_float_data(q$floats)
x<-d$Data$F6903093

library(purrr)
x<-x %>%  keep(~min(.$TIME) >= z$TIME[1] & max(.$TIME) <= z$TIME[length(z$TIME)])

keep(~min(.$range_value) >= 0 & max(.$range_value) <= 100)

d$Mdata



setwd(wd)
setwd(robj)


plot_profiles(profile_ids=q$profiles, variables='DOXY')

load_float_data(float_ids=as.numeric(q$floats), float_profs=q$profiles)

show_sections()

