# script header

## Load One argo scripts
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


# set packages
packages<-c("plyr","sp", "rgdal",  "rgeos", "raster", "readr", "tidyverse", "lubridate",  "ggthemes",  "sf", "cmocean", "ncdf4", "RNetCDF",  "plot3D", "tidync",
"devtools", "stars", "ncmeta", "maps", "oce", "data.table", "fasterize", "RStoolbox", "scales", "purrr", "HURDAT" , "gsw", "R.utils","Matrix", "cowplot", "purrr")
# "HURDAT",