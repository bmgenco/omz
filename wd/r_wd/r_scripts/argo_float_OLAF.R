rm(list=ls())

# wd<-"/home/brandon/vestawd/omz/wd/r_wd/r_scripts"
wd<-"/home/brandon/europawd/omz/wd/r_wd/"
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
path_code = "r_scripts/submodules/BGC-ARGO_R_WORKSHOP/"

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


rm(f.ipak, packages, path_code)

# here

# run code in arfo float OLAF..


setwd(robj)
load("init_argo.RData")
# initialize_argo()
float<-"6903093"
 

float_idx <-which(Float$wmoid==float)
prof_ids = c(Float$prof_idx1[float_idx]:Float$prof_idx2[float_idx])

dates = Sprof$date[prof_ids] 
# edit here for dates I want, then slect profiles for plottinf functions and laoding nc data

h.pts$DateTime

st_date<-min(h.pts$DateTime)
ed_date<-max(h.pts$DateTime)

window<-dates[dates <= ed_date] %>% .[dates >= st_date]
which(dates = st_date)


sensors = unique(Sprof$sens[prof_ids])
load_float_data <-(float_ids=float)

# plot_trajectories(float)

WMO<-as.numeric(float)
# plot_trajectories()
setwd("..")

float_file = nc_open(paste0("Profiles/", WMO,"_Sprof.nc"))
data = load_float_data(float_ids= WMO, # specify WMO number
                       variables=c('PSAL','TEMP', 'DOXY'))

show_trajectories(float_ids=WMO)


show_profiles( profile_ids=WMO, 
               variables=c('TEMP', 'DOXY'),
               type="floats", # given IDs refer to the floats
               obs='on', # 'on' shows points on the profile at which each measurement was made
               raw="yes" # show the unadjusted data 
)

