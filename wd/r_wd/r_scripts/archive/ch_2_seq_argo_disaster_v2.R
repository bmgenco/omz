#to do 

# - change default server to https://data-argo.ifremer.fr/dac/coriolis/6903093/profiles/ 


#### attmep 1: ####

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
# lapply(packages, require, character.only = TRUE)
# rm(f.ipak, packages)

setwd(wd)
print(paste0("Current Working Directory is ", getwd()))

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

packages<-c("gsw", "R.utils","Matrix")
f.ipak(packages)
rm(f.ipak, packages, path_code)


## bah

setwd(wd)
setwd("../../data/argo_temp")
directory<-list.files()
setwd(directory)
file<-list.files()
x<-read.csv(file)


# old
setwd(wd)
setwd(robj)
load("init_argo.RData")

olaf_profiles<-readRDS("olaf_profiles_list.R")


# x<-olaf_profiles$window.sp
# z<-x$`6903093` %>% filter(., CYCLE_NUMBER %in% c(48, 49, 50, 51, 52))
# e<-extent(z)

x<-olaf_profiles$window
z<-x$`6903093`$data

# z<-x$`6903093` %>% filter(., CYCLE_NUMBER %in% c(48, 49, 50, 51, 52))
# e<-extent(z)
# q<-select_profiles(lon_lim=c(e[1],e[2]), lat_lim=c(e[3], e[4]), start_date=(z$TIME[1]-1), end_date=(z$TIME[length(z$TIME)]+1), outside="none", sensor=NULL)
# float_n<-as.numeric(q$floats[1])

float_n<-6903093
profiles<-as.character(z$profiles)
data = load_float_data(float_ids = float_n float_profs = profiles )

float_idx <-which(Float$wmoid==q$floats[1])

prof_ids = c(Float$prof_idx1[float_idx]:Float$prof_idx2[float_idx])
v<-which(prof_ids==q$profiles)


trajectory = show_trajectories(float_ids=q$floats, return_ggplot=TRUE, prof_ids = q$profiles)

plot_profiles(profile_ids=q$profiles[1], variables=c('TEMP'), type="profiles", qc_flags=c(0,9))


q$floats

setwd(wd)
setwd(robj)
setwd("Profiles")

profiles<-as.character(q$profiles)


data=load_float_data(float_ids =float_n)

show_profiles( profile_ids=q$profiles, \
               variables = "DOXY"
               type="profiles", # given IDs refer to the profiles
               per_float=0,  # show all profiles in one plot:
               obs='on' # plot a marker at each observation
               
)

show_profiles(profile_ids=q$floats[1], 
              variables=c('PSAL','DOXY'),
              type="floats", # given IDs refer to the floats
              obs='on',# 'on' shows points on the profile at which
              #  each measurement was made
)





#### attemp 2: ####

rm(list=ls())
# relative directories
wd<-"/home/brandon/callistowd/omz/wd/r_wd"
robj<-"r_objects"
fig<-"../../figures"
argoFloat<-"../../data/argo/argoFloat"
#coriolis
argocor<-"../../data/argo/coriolis_dac"
setwd(wd)
setwd(argoFloat)
argofloat_d<-getwd()
setwd(wd)
setwd(argocor)
cor_d<-getwd()

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

packages<-c("sp", "rgdal",  "rgeos", "raster", "readr", "tidyverse", "lubridate",  "ggthemes",  "sf", "cmocean", "ncdf4", "RNetCDF",  "plot3D", "tidync", "devtools", "stars", "ncmeta", "maps", "oce","data.table","raster", "fasterize", "RStoolbox", "scales", "purrr", "marmap", "cowplot", "argoFloats", "ggplotify")

f.ipak(packages)
setwd(wd)

# options(argoFloats.server="usgodae")
# indexAll<-getIndex(destdir = argofloat_d)
setwd(robj)
# olaf_profiles<-readRDS("olaf_profiles_list.R")

# x<-olaf_profiles$window$`6903093`
# z<-olaf_profiles$window.sp$`6903093`

## download vi argoFloat:
# index<-subset(indexAll, ID="6903093")
# index<-subset(index, cycle=z$CYCLE_NUMBER)
# 
# st<-as.character(date(x$data$TIME[1]))
# en<-as.character(date(x$data$TIME[12]))
# olaf<-index%>%subset(., time=list(from=st, to=en))
# setwd(wd)
# 
# 
# # profiles<-getProfiles(index, destdir = argo_d, age=0)
# # argos <-readProfiles(profiles, destdir = argo_d)
# test<-argos[["oxygen"]]

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






# s<-as.section(c(u,v,w,x,y,z))

# z<-tidync(files[1])

# plot(x@data$pressure, x@data$oxygen)
# 
# data<-x@data
# i<-which(names(x@data)=="oxygenAdjusted")
# # a<-"[µmol/kg]"
# names(x@data)[i]<-"Oxygen"

# names(x@data)[i]<-expression("Oxygen ["*µmol/kg*"]")
# oce::plotProfile(x, xtype = "oxygenAdjusted", ytype="depth", ylim=c(200,0),xlim=c(5,250), xlab= expression("Oxygen ["*µmol/kg*"]"))


"oxygenAdjustedError"


#### errors wit plot size record plot methpd ####

mar=c(0,3.1,3.1,0.15), cex.main=.25, cex.axis=.25, cex.lab=.25, cex=.25

dev.off()
par(mfrow=c(1,1), cex=.75)

# dev.new(width=2 ,height=6, noRStudioGD = F, unit="in")
oce::plotProfile(t, xtype = "oxygen", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = 0, col="blue", mar=c(0,3.1,3.1,0.2)) # hack for x lab
oce::plotProfile(t, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "dashed", col="blue", add=T)
oce::plotProfile(u, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "dotted", col="blue", add=T)
oce::plotProfile(v, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "solid", col="blue", add=T)
bo<-recordPlot()


dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(w, xtype = "oxygen", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = 0, col="blue",  mar=c(0,3.1,3.1,0.2))
oce::plotProfile(w, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "dotted", col="blue", add=T)
oce::plotProfile(x, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "solid", col="blue", add=T)
do<-recordPlot()

dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(y, xtype = "oxygen", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = 0, col="blue",  mar=c(0,3.1,3.1,0.2)) # hack for x lab
oce::plotProfile(y, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "dashed", col="blue", add=T)
oce::plotProfile(z, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "dotted", col="blue", add=T)
oce::plotProfile(a, xtype = "oxygenAdjusted", ytype="depth", ylim=c(175,0),xlim=c(5,250), lty = "solid", col="blue", add=T)
ao<-recordPlot()

## Temp

dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(t, xtype = "temperature", ytype="depth", ylim=c(175,0), xlim=c(5,35), lty = "dashed", col="red",  mar=c(0,3.1,3.1,0.2))
oce::plotProfile(u, xtype = "temperature", ytype="depth", ylim=c(175,0), xlim=c(5,35), lty = "dotted", col="red", add=T)
oce::plotProfile(v, xtype = "temperature", ytype="depth", ylim=c(175,0),  xlim=c(5,35), lty = "solid", col="red", add=T)
bt<-recordPlot()

dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(w, xtype = "temperature", ytype="depth", ylim=c(175,0), xlim=c(5,35), lty = "dotted", col="red",  mar=c(0,3.1,3.1,0.2))
oce::plotProfile(x, xtype = "temperature", ytype="depth", ylim=c(175,0), xlim=c(5,35), lty = "solid", col="red", add=T)
dt<-recordPlot()

dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(y, xtype = "temperature", ytype="depth", ylim=c(175,0), xlim=c(5,35), lty = "dashed", col="red",  mar=c(0,3.1,3.1,0.2))
oce::plotProfile(z, xtype = "temperature", ytype="depth", ylim=c(175,0), xlim=c(5,35), lty = "dotted", col="red", add=T)
oce::plotProfile(a, xtype = "temperature", ytype="depth", ylim=c(175,0), xlim=c(5,35), lty = "solid", col="red", add=T)
at<-recordPlot()

## CHl-a
dev.off()
par(mfrow=c(1,1), cex=.75)

# dev.new(width=2 ,height=6, noRStudioGD = F, unit="in")
oce::plotProfile(t, xtype = "chlorophyllA", ytype="depth", ylim=c(175,0), xlim=c(0,1), lty = 0, col="green", mar=c(0,3.1,3.1,0.2)) # hack for x lab
oce::plotProfile(t, xtype = "chlorophyllAAdjusted", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dashed", col="green", add=T)
oce::plotProfile(u, xtype = "chlorophyllAAdjusted", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dotted", col="green", add=T)
oce::plotProfile(v, xtype = "chlorophyllAAdjusted", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "solid", col="green", add=T)
bc<-recordPlot()


dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(w, xtype = "chlorophyllA", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = 0, col="green",  mar=c(0,3.1,3.1,0.2))
oce::plotProfile(w, xtype = "chlorophyllAAdjusted", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dotted", col="green", add=T)
oce::plotProfile(x, xtype = "chlorophyllAAdjusted", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "solid", col="green", add=T)
dc<-recordPlot()

dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(y, xtype = "chlorophyllA", ytype="depth", ylim=c(175,0),xlim=c(0,2), lty = 0, col="green",  mar=c(0,3.1,3.1,0.2)) # hack for x lab
oce::plotProfile(y, xtype = "chlorophyllAAdjusted", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dashed", col="green", add=T)
oce::plotProfile(z, xtype = "chlorophyllAAdjusted", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dotted", col="green", add=T)
oce::plotProfile(a, xtype = "chlorophyllAAdjusted", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "solid", col="green", add=T)
ac<-recordPlot()


## CDOM
par(mfrow=c(1,1), cex=.75)

# dev.new(width=2 ,height=6, noRStudioGD = F, unit="in")
oce::plotProfile(t, xtype = "CDOM", ytype="depth", ylim=c(175,0), xlim=c(0,1), lty = 0, col="brown", mar=c(0,3.1,3.1,0.2)) # hack for x lab
oce::plotProfile(t, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dashed", col="brown", add=T)
oce::plotProfile(u, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dotted", col="brown", add=T)
oce::plotProfile(v, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "solid", col="brown", add=T)
bcd<-recordPlot()


dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(w, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = 0, col="brown",  mar=c(0,3.1,3.1,0.2))
oce::plotProfile(w, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dotted", col="brown", add=T)
oce::plotProfile(x, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "solid", col="brown", add=T)
dcd<-recordPlot()

dev.off()
par(mfrow=c(1,1), cex=.75)
# dev.new(width=.5 ,height=4, noRStudioGD = F, unit="in")
oce::plotProfile(y, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,2), lty = 0, col="brown",  mar=c(0,3.1,3.1,0.2)) # hack for x lab
oce::plotProfile(y, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dashed", col="brown", add=T)
oce::plotProfile(z, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "dotted", col="brown", add=T)
oce::plotProfile(a, xtype = "CDOM", ytype="depth", ylim=c(175,0),xlim=c(0,1), lty = "solid", col="brown", add=T)
acd<-recordPlot()

setwd(wd)

plot_list<-list(bo, do, ao, bt, dt, at, bc, dc, ac, bcd, dcd, acd)
setwd(wd)
setwd(robj)

saveRDS(plot_list, "oflaf_profiles.R")

setwd(wd)
setwd(fig)
dev.off()
pdf("olaf_profiles.pdf", height =9, width=8)
plot_grid(bo, do, ao, bt, dt, at, nrow=2, ncol=3,
         labels = c('a', 'b', 'c', 'd', 'e', 'f'), align = "h", scale=.8, greed=T, axis="r")

dev.off()



setwd(wd)
setwd(fig)
dev.off()
pdf("chla_cdom.pdf", height =9, width=8)
plot_grid(bc, dc, ac, bcd, dcd, acd, nrow=2, ncol=3,
          labels = c('a', 'b', 'c', 'd', 'e', 'f'), align = "h", scale=.8, greed=T, axis="r")
dev.off()

### adding letter after the fact


mtext("My Multiplot Title", side = 3, line = - 2,outer = TRUE)


# oce::plotProfile(x, xtype = "Oxygen", ytype="depth", ylim=c(200,0),xlim=c(5,250))
# Here
# colors before druiinf after
# lines date

#see here:



plot(x, which="oxygenAdjusted", ztype="image", xtype="time")



, col=ramp[1], ylim=c(500,0),xlim=c(5,260))

plot(argos, which="profile", col="red")

summary(argos)




setwd(wd)
setwd(argo_d)

files<-list.files()
file<-files[6]

t<-tidync(file)



#### attemp 3 ####
setwd(wd)
setwd("../../data/argo_temp")
directory<-list.files()
setwd(directory)
file<-list.files()
x<-read.csv(file)


float_idx <-which(Float$wmoid==float_n

                  