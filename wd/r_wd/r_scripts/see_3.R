
### set up ####
rm(list=ls())

r_scriptwd<-getwd()
# wd<-substring(r_scriptwd,1,nchar(r_scriptwd) -10) # hardcoding replace here
wd<-"/home/brandon/callistowd/omz/wd/r_wd/"


setwd(wd)
# relative directories
robj<-"r_objects"
fig<-"../../figures"
gis_data<-"../../data/gis_data"
output<-"../../output"
data_d<-"../../data"

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

path<-"r_scripts/submodules/OneArgo-R/"

source(paste0(path, "get_lon_lat_time.R"))
source(paste0(path, "initialize_argo.R"))
source(paste0(path, "load_float_data.R"))
source(paste0(path, "select_profiles.R"))
source(paste0(path, "show_profiles.R"))
source(paste0(path, "show_sections.R"))
source(paste0(path, "show_time_series.R"))
source(paste0(path, "show_trajectories.R"))

# auxilary functions
path = "r_scripts/submodules/OneArgo-R/auxil/"
file_path<-paste0(wd,"/",path)
files<-list.files(file_path, pattern ="\\.R$")

for(i in 1:length(files)){
  source(paste0(path, files[i]))
}

# packages<-c("gsw", "R.utils","Matrix")
# f.ipak(packages)
# 
rm(path, file_path, files)

### Functions ####
f.select_name<-function(x, key_id){
  x<-x$tc
  x<-x[which(names(x) == key_id)]
  return(x)}

### Data ####

data<-readRDS("r_objects/ch2_seq_1_data.R")  
ep_doxy_floats<-readRDS("r_objects/ep_doxy_floats.rds")

setwd(gis_data)
coastline<-st_read("./ne_50m_land/ne_50m_land.shp") 
setwd(wd)

function_variables<-data$function_variables
# hurdat objects
h<-data$selected_zone_hurdat$h
h.pts<-data$selected_zone_hurdat$h.pts
h.tracks<-data$selected_zone_hurdat$h.tracks

# match_up

list.100<-data$tc_search_radi_profiles$km_100
list.200<-data$tc_search_radi_profiles$km_200
list.500<-data$tc_search_radi_profiles$km_500

rm(data)


# TC 1 FLECIA EP082009

key_id<-"EP082009"
p<-f.select_name(x=list.500, key_id)
tc<-filter(h, Key == key_id )



#### profiles ###

p1<-show_profiles(float_ids = p$EP082009$float[1], float_profs = p$EP082009$profile[1])
p2<-show_profiles(float_ids = p$EP082009$float[1], float_profs = p$EP082009$profile[2])
p3<-show_profiles(float_ids = p$EP082009$float[1], float_profs = p$EP082009$profile[3])
p4<-show_profiles(float_ids = p$EP082009$float[1], float_profs = p$EP082009$profile[4])
p5<-show_profiles(float_ids = p$EP082009$float[1], float_profs = p$EP082009$profile[5])
all<-show_profiles(float_ids = p$EP082009$float[1], float_profs = p$EP082009$profile)


show_sections(float_ids = p$EP082009$float[1], float_profs =as.array(p$EP082009$profile))

