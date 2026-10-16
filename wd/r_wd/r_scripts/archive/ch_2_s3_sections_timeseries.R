#### set_up ####

rm(list=ls())
# r_scriptwd<-getwd()
# wd<-substring(r_scriptwd,1,nchar(r_scriptwd) -10) # hardcoding replace here
# wd<-"/home/brandon/vestawd/omz/wd/r_wd"
wd<-"/home/brandon/callistowd/omz/wd/r_wd"
# wd<-"/home/brandon/europawd/omz/wd/r_wd"

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
# load("r_objects/init_argo.RData")


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

rm(data)

#load and creatw dta objects




#### olaf ####
key_id<-"EP152021"
tc.pts<-f.single_storm_h.pts(key_id = key_id, h.pts)
tc.p<-list.500$tc$EP152021 # can make a function
#co .ld be used to find profiels in time and sapce not asscoaited with hurricnae
# tc.p<-f.select_tc_specific_profiles(list.x = list.500, key_id = key_id) 


#load and creatw dta objects




# specific float
# x<-tc.p%>% filter(. , float == "6903094")
x<-tc.p%>% filter(. , float == "5905068")

# x<-tc.p%>% filter(. , float == "5905068")

s<-f.reformat_profiles_for_ONE_Argo(x)




all.floats<-f.reformat_profiles_for_ONE_Argo(tc.p)

# olaf.profiles<-load_float_data(float_ids = all.floats$float_ids, float_profs = all.floats$float_profs)



# picking tjus t to floats aof intertest
# floats_2<-unisue(tc.p$EP152021$float) %>%as.numeric(.) %>%.[-1]



# trying porfiles
# data.argo<-load_float_data(float_ids =s$float_ids, float_profs = s$float_profs)

# x<-data.argo$Data$F6903093
# y< as.argo(time=x$TIME, longitude = x$LONGITUDE, latitude =  x$LATITUDE, salinity = x$)

#### Plotting olaf ####
show_profiles(float_ids = s$float_ids[1], float_profs=s$float_profs[[1]][2], raw='yes')

show_sections(float_ids = s$float_ids[1], float_profs = s$float_profs[1], max_depth=200, raw='no', plot_mld=2)

show_sections(float_ids = all.floats$float_ids[3])

show_sections(float_ids = all.floats$float_ids[1], max_depth =200, plot_mld = 2, raw='no')

# float 94
show_sections(float_ids = all.floats$float_ids[1], float_profs = all.floats$float_profs[1], max_depth =150, plot_mld = 2, raw='no')


# show_time_series(float_ids = float.x, plot_depth =50)
show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[[1]], plot_depth =90)
plot.traj<-show_trajectories(float_ids = all.floats$float_ids, float_profs = all.floats$float_profs, return_ggplot = T )
# farting around with ggplot build



tc.p## add pgplot
olaf.track<-filter(h.tracks, Key== key_id) %>% st_as_sf(
ggplot()+geom_sf(data=olaf.track)+plot.traj

# messing with sections
show_sections(float_ids = s$float_ids[1], float_profs = s$float_profs[1])

show_sections(float_ids = s$float_ids[2], float_profs = s$float_profs[2] )

show_trajectories(float_ids = s$float_ids[1], float_profs = s$float_profs[1])
t.50<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =50)
dev.copy2pdf(file="t1_50.pdf")
t.60<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =60)
t.70<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =70)
t.80<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =80)
t.90<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =90)
t.100<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =100)
show_profiles(float_ids = s$float_ids[1], float_profs = s$float_profs[1])


t.50<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =50, variables=c("DOXY", "TEMP"))
t.60<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =60, variables=c("DOXY", "TEMP", "CHLA"))
t.70<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =70, variables=c("DOXY", "TEMP", "CHLA"))
t.80<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =80, variables=c("DOXY", "TEMP")
t.90<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =90, variables=c("DOXY", "TEMP", "CHLA"))
t.100<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =100, variables=c("DOXY", "TEMP", "CHLA"))                                             

### oce object













#### Geneive ####
key_id<-"EP122020"
tc.pts<-f.single_storm_h.pts(key_id = key_id, h.pts)
tc.p<-list.500$tc$EP122020
s<-f.reformat_profiles_for_ONE_Argo(tc.p)
Geneive.profiles<-load_float_data(float_ids = s$float_ids, float_profs = s$float_profs)

### plotting genieve
t.50<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =50)
t.60<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =60)
t.70<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =70)
t.80<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =80)
t.90<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =90)
t.100<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =100)
show_profiles(float_ids = s$float_ids[1], float_profs = s$float_profs[1])

show_sections(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_mld = 2, max_depth = 200)




### kay ####

# did I plot the worng porfiles?

key_id<-"EP122022"
tc.pts<-f.single_storm_h.pts(key_id = key_id, h.pts)
tc.p<-list.500$tc$EP122021 # can make a function
#co .ld be used to find profiels in time and sapce not asscoaited with hurricnae
# tc.p<-f.select_tc_specific_profiles(list.x = list.500, key_id = key_id) 

# specific float
x<-tc.p%>% filter(. , float == "6903093")
s<-f.reformat_profiles_for_ONE_Argo(x)
all.floats<-f.reformat_profiles_for_ONE_Argo(tc.p)

kay.profiles<-load_float_data(float_ids = all.floats$float_ids, float_profs = all.floats$float_profs)

sb<-s$float_profs[[1]][1:4]
show_profiles(float_ids = s$float_ids[1], float_profs = sb, variables=c("TEMP", "CHLA", "DOXY"))
sd<-s$float_profs[[1]][5:6]
show_profiles(float_ids = s$float_ids[1], float_profs = 54, variables=c("TEMP", "CHLA", "DOXY"))
sa<-s$float_profs[[1]][7:10]
show_profiles(float_ids = s$float_ids[1], float_profs = sa, variables=c("TEMP", "CHLA", "DOXY"))
show_sections(float_ids = s$float_ids[1], float_profs = s$float_profs[1])

x<-tc.p%>% filter(. , float == "5905068")
s<-f.reformat_profiles_for_ONE_Argo(x)

### plotting genieve
t.50<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =50)
t.60<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =60)
t.70<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =70)
t.80<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =80)
t.90<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =90)
t.100<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =100)
show_profiles(float_ids = s$float_ids[1], float_profs = s$float_profs[1])

show_profiles(float_ids = s$float_ids[1], float_profs = s$float_profs[1][1:3], variables=c("TEMP", "CHLA", "DOXY"))
t.50<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =50, variables=c("TEMP", "CHLA"))
t.60<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =60, variables=c("TEMP", "CHLA"))
t.70<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =70, variables=c( "TEMP", "CHLA"))
t.80<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =80, variables=c("TEMP", "CHLA"))
t.90<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =90, variables=c( "TEMP", "CHLA"))
t.100<-show_time_series(float_ids = s$float_ids[1], float_profs = s$float_profs[1], plot_depth =100, variables=c("TEMP", "CHLA"))         

### saveing porifles

profiles.data<-list(olaf.profiles, kay.profiles, geneive.profiles)
names(profiles.data)<-c("olaf.profiles", "kay.profiles", "geneive.profiles")
setwd(robj)
saveRDS(profiles.data, "profiles_data.rds")

setwd(robj)
test<-readRDS("profiles_data.rds")
