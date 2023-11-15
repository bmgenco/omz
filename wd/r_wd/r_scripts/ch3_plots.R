#### setup ####

rm(list=ls())
wd<-"/home/brandon/callistowd/omz/wd/r_wd"
# wd<-"/home/brandon/vestawd/omz/wd/r_wd"


# relative directories

robj<-"r_objects"
fig<-"../../figures"
gis_data<-"../../data/gis_data"
output<-"../../output"
data_d<-"../../data"
cor_d<-"../../data/argo/coriolis_dac"
my_functions<-"r_scripts/submodules/BG_functions/"

setwd(wd)

source(paste0(my_functions, "ch2_functions.R"))
source(paste0(my_functions, "ch2_header.R"))

### data 1 ####
data<-readRDS("r_objects/ch2_seq_1_data.R")

h<-data$selected_zone_hurdat$h
function_variables<-data$function_variables
# hurdat 
h.pts<-data$selected_zone_hurdat$h.pts
h.tracks<-data$selected_zone_hurdat$h.tracks

list.100<-data$tc_search_radi_profiles$km_100
list.200<-data$tc_search_radi_profiles$km_200
list.500<-data$tc_search_radi_profiles$km_500

profile.pts<-data$selected_zone_float_metadata$df.profile_points

rm(data)

### Olaf Profiles  ####


setwd(wd)
setwd(cor_d)
files<-list.files()

# plot
t<-read.argo(files[4])
u<-read.argo(files[5])
v<-read.argo(files[6])
w<-read.argo(files[7])
x<-read.argo(files[8])
y<-read.argo(files[9])
z<-read.argo(files[10])
a<-read.argo(files[11])

olaf.sec<-as.section(list(t,u,v,w,x,y,z,a))

## ploting sections

plot(olaf.sec, which="oxygen2", ztype="image", xtype="time", ylim=c(80, 0), xlab=NULL)


### ploting: ### old plots

oce::plotProfile(t, xtype = "oxygen", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = 0, col="blue", mar=c(0,3.1,3.1,0.2)) 
oce::plotProfile(t, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "dashed", col="blue", add=T)
