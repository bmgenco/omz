#set_up

rm(list=ls())
# r_scriptwd<-getwd()
# wd<-substring(r_scriptwd,1,nchar(r_scriptwd) -10) # hardcoding replace here
wd<-"/home/brandon/vestawd/omz/wd/r_wd"
setwd(wd)

# relative directories
robj<-"r_objects"
fig<-"../../figures"
gis_data<-"../../data/gis_data"
data_d<-"../../data/2018_data"
# data_d<-"../../data"
ocean_color<-"../../data/ocean_color_bud"
bathy_d<-"../../data/bathy" # edit this

getwd()
f.ipak <- function(pkg){
  # loads packages, quietly, given by a vector of package names e.g., pkg<-c("ggplot", "tidyverse")
  # will install  packages listed , and their dependencies, if needed.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, quiet=T, verbose = F)
  sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts=F)
}
packages<-c("plyr","sp", "rgdal",  "rgeos", "raster", "readr", "tidyverse", "lubridate",  "ggthemes",  "sf", "cmocean", "ncdf4", "RNetCDF",  "plot3D", "tidync", "devtools", "stars", "ncmeta", "maps", "oce", "data.table", "fasterize", "RStoolbox", "scales", "purrr", "HURDAT" , "gsw", "R.utils","Matrix")
# "HURDAT",
f.ipak(packages)


# functions
setwd(wd)

path<-"r_scripts/submodules/OneArgo-R/"

source(paste0(path, "get_lon_lat_time.R"))
source(paste0(path, "initialize_argo.R"))
source(paste0(path, "load_float_data.R"))
source(paste0(path, "select_profiles.R"))
source(paste0(path, "show_profiles.R"))
source(paste0(path, "show_sections.R"))
source(paste0(path, "show_time_series.R"))
source(paste0(path, "show_trajectories.R"))

path = "r_scripts/submodules/OneArgo-R/auxil/"
file_path<-paste0(wd,"/",path)
files<-list.files(file_path, pattern ="\\.R$")

for(i in 1:length(files)){
  source(paste0(path, files[i]))
}

rm(path, file_path, files)

# my functions
f.select_tc<-function(x, key_id){
  x<-x$tc
  x<-x[which(names(x) == key_id)]
  return(x)
}

# data

data<-readRDS("r_objects/ch2_seq_1_data.R")
#objects
h<-data$selected_zone_hurdat$h
function_variables<-data$function_variables
# hurdat 
h.pts<-data$selected_zone_hurdat$h.pts
h.tracks<-data$selected_zone_hurdat$h.tracks


list.100<-data$tc_search_radi_profiles$km_100
list.200<-data$tc_search_radi_profiles$km_200
list.500<-data$tc_search_radi_profiles$km_500

# get meta data
initialize_argo()

# olaf
key_id<-"EP152021"
tc.pts<-filter(h.pts, Key ==key_id )
conv<-1.94384 # convert to m/s
tc.pts$max_sustained_wind_m_s<-tc.pts$Wind/conv

f.select_tc_specific_profiles<-function(list.500, key_id){

tc.p<-f.select_tc(x=list.500, key_id) %>%.[[1]] %>%st_as_sf()
e<-extent(tc.p)%>%as.vector(.)
q<-select_profiles(lon_lim=c(e[1],e[2]), lat_lim=c(e[3], e[4]), start_date = as.character(min(tc.p$time)-1), end_date =as.character(max(tc.p$time)+1),  outside="none", sensor=NULL)
}

#convert tc.p into frame




# specifc float
x<-tc.p$EP152021 %>% filter(. , float == "6903093")
float.x<-as.numeric(x$float[1])

floats_2<-unique(tc.p$EP152021$float) %>%as.numeric(.) %>%.[-1]

# show_time_series(float_ids = float.x, plot_depth =50)
show_time_series(float_ids = floats_2, plot_depth =90)
show_trajectories(float_ids = floats_2, )



min(tc.p$EP152021$)





