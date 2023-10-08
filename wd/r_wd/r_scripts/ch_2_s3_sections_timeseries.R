#### set_up ####

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
# ocean_color<-"../../data/ocean_color_bud"bathy_d<-"../../data/bathy" # edit this
my_functions<-"r_scripts/submodules/BG_functions/"

# load functions and default packges

source(paste0(my_functions, "ch2_functions.R"))
source(paste0(my_functions, "ch2_header.R"))

f.ipak(packages)
## get metadata ##
initialize_argo()


#### data ####
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



# olaf
key_id<-"EP152021"
tc.pts<-f.single_storm_h.pts(key_id = key_id, h.pts)
tc.p<-list.500$tc$EP152021 # can make a function
#co .ld be used to find profiels in time and sapce not asscoaited with hurricnae
# tc.p<-f.select_tc_specific_profiles(list.x = list.500, key_id = key_id) 

# specific float
x<-tc.p%>% filter(. , float == "6903093")
s<-f.reformat_profiles_for_ONE_Argo(x)


all.floats<-f.reformat_profiles_for_ONE_Argo(tc.p)

# picking tjus t to floats aof intertest
# floats_2<-unisue(tc.p$EP152021$float) %>%as.numeric(.) %>%.[-1]



# trying porfiles
data.argo<-load_float_data(float_ids =s$float_ids, float_profs = s$float_profs)

x<-data.argo$Data$F6903093
# y< as.argo(time=x$TIME, longitude = x$LONGITUDE, latitude =  x$LATITUDE, salinity = x$)

#### Plotting ####
show_sections(float_ids = all.floats$float_ids[3])

show_sections(float_ids = all.floats$float_ids[2], max_depth =150, plot_mld = 2)

# float 94
show_sections(float_ids = all.floats$float_ids[1], float_profs = all.floats$float_profs[1], max_depth =150, plot_mld = 2)


# show_time_series(float_ids = float.x, plot_depth =50)
show_time_series(float_ids = s$float_ids, float_profs = s$float_profs, plot_depth =90)
plot.traj<-show_trajectories(float_ids = all.floats$float_ids, float_profs = all.floats$float_profs, return_ggplot = T )
## add pgplot
plot.tc_and_trag<-p

# messing with sections
show_sections(float_ids = s$float_ids, float_profs = s$float_profs[[1]])

show_sections(float_ids = s$float_ids[2], float_profs = s$float_profs[2] )

show_trajectories(float_ids = s$float_ids[2], float_profs = s$float_profs[2])
t.50<-show_time_series(float_ids = s$float_ids[2], float_profs = s$float_profs[2], plot_depth =50)
t.60<-show_time_series(float_ids = s$float_ids[2], float_profs = s$float_profs[2], plot_depth =60)
t.70<-show_time_series(float_ids = s$float_ids[2], float_profs = s$float_profs[2], plot_depth =70)
t.80<-show_time_series(float_ids = s$float_ids[2], float_profs = s$float_profs[2], plot_depth =80)
t.90<-show_time_series(float_ids = s$float_ids[2], float_profs = s$float_profs[2], plot_depth =90)
show_profiles(float_ids = s$float_ids[2], float_profs = s$float_profs[2])

### oce object




