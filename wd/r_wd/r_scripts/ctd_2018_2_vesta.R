# dissertation messing about CTD data 2018 v2 
# machine: vesta
# as far a
# new code for bud rewrite papper 2022-10-01

# plotting depends on package oce
# see https://dankelley.github.io/oce/
# and: https://bscheng.com/2017/10/08/oceanographic-data/

#################################################################
#                         TO DO                                 #
# 1: explore use of "useSmoothScatter"
# 2: write functions to simplify
# 3: see function options:
# https://dankelley.github.io/oce/reference/plot-ctd-method.html
# 4: plot sections
# 5: om section maps rmeove line ad numbers
# 5: fix casts with incorrect staion ids

#################################################################
rm(list=ls())


### install packages ###
f.ipak <- function(pkg){
  
  # loads packages, quietly, given by a vector of package names e.g., pkg<-c("ggplot", "tidyverse")
  # will install  packages listed , and their dependencies, if needed.
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, quiet=T, verbose = F)
  sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts=F)
}
# packages<-c("oce", "marmap", "sf", "raster", "lubridate", "ggplot2")  # set packages here
packages<-c("oce", 'sp', "sf", "rgeos", "lubridate", "tidyverse", "data.table")
f.ipak(packages)

## set directory
# wd<-"/home/brandon/vestawd/omz/r_wd"
# robj<-"/home/brandon/vestawd/omz/r_wd/r_objects"
# fig<-"/home/brandon/vestawd/omz/r_wd/figures"
# gis_data<-"/home/brandon/vestawd/omz/gis_data"
# data_d<-"/home/brandon/vestawd/omz/2018_data/"
# ctdata<-"/home/brandon/vestawd/omz/2018_data/ctd/data"


## new directoires:
wd<-"/home/brandon/vestawd/omz/wd/r_wd"
robj<-"/home/brandon/vestawd/omz/wd/r_wd/r_objects"
fig<-"/home/brandon/vestawd/omz/figures"
gis_data<-"/home/brandon/vestawd/omz/data/gis_data"
data_d<-"/home/brandon/vestawd/omz/data/2018_data/"
ctdata<-"/home/brandon/vestawd/omz/data/2018_data/ctd/data"

setwd(wd)
#list.files(ctdata)

#### functions ####
f.read<-function(ctdata){
  files<-list.files(ctdata)
  data<-vector(mode="list", length=length(files))
  
  for(i in 1:length(files)){
    # add if else statement
    data[[i]]<-read.ctd(file.path(ctdata, files[]),)  
   # trim downcast
    x<-ctdTrim(x, method = "downcast")
    }
return(data)}

#temp
# x<-ctdTrim(x, method = "downcast")
# plotProfile(x)
# plot(station4, which=c(1,2,3,5), type="l", span=150)


#### new code ####
setwd(robj)
d.sec<-readRDS("OC1806A_ctd_sections.R")

d35<-d.sec$st3.5
d2<-d.sec$st2
d4<-d.sec$st4




####  create stations section - make a function - upcasts ####

files<-list.files(ctdata)
## station 1 - casts ##


## station 2 - casts ##
c9<-read.ctd(file.path(ctdata, files[89]),) %>% ctdTrim(., method = "upcast")
c10<-read.ctd(file.path(ctdata, files[3]),) %>% ctdTrim(., method = "upcast")
c11<-read.ctd(file.path(ctdata, files[5]),) %>% ctdTrim(., method = "upcast")
c12<-read.ctd(file.path(ctdata, files[7]),) %>% ctdTrim(., method = "upcast")
c13<-read.ctd(file.path(ctdata, files[9]),) %>% ctdTrim(., method = "upcast")
c14<-read.ctd(file.path(ctdata, files[11]),) %>% ctdTrim(., method = "upcast")
c15<-read.ctd(file.path(ctdata, files[13]),) %>% ctdTrim(., method = "upcast")
c16<-read.ctd(file.path(ctdata, files[15]),) %>% ctdTrim(., method = "upcast")
c17<-read.ctd(file.path(ctdata, files[17]),) %>% ctdTrim(., method = "upcast")

#section
u2<-as.section(list(c9,c10,c11,c12,c13,c14,c15,c16,c17))

#c30 is mis labelledas station 3

## station 3.5 - casts ##
c30<-read.ctd(file.path(ctdata, files[47]),) %>% ctdTrim(., method = "upcast")
c31<-read.ctd(file.path(ctdata, files[49]),) %>% ctdTrim(., method = "upcast")
c32<-read.ctd(file.path(ctdata, files[51]),) %>% ctdTrim(., method = "upcast")
c33<-read.ctd(file.path(ctdata, files[53]),) %>% ctdTrim(., method = "upcast")
c34<-read.ctd(file.path(ctdata, files[55]),) %>% ctdTrim(., method = "upcast")
#section
u35<-as.section(list(c30,c31,c32,c33,c34))

## station 4 - casts ##
c35<-read.ctd(file.path(ctdata, files[57]),) %>% ctdTrim(., method = "upcast")
c36<-read.ctd(file.path(ctdata, files[59]),) %>% ctdTrim(., method = "upcast")
c37<-read.ctd(file.path(ctdata, files[61]),) %>% ctdTrim(., method = "upcast")
c38<-read.ctd(file.path(ctdata, files[63]),) %>% ctdTrim(., method = "upcast")
c39<-read.ctd(file.path(ctdata, files[65]),) %>% ctdTrim(., method = "upcast")
c40<-read.ctd(file.path(ctdata, files[69]),) %>% ctdTrim(., method = "upcast")
#section
u4<-as.section(list(c35,c36,c37,c38,c39,c40))


setwd(robj)

x<-vector(mode = "list", length =6)
x[[1]]<-d2
x[[2]]<-d35
x[[3]]<-d4
x[[4]]<-u2
x[[5]]<-u35
x[[6]]<-u4

saveRDS(x, "2018_ctd_comparison_up_v_downcasts.R")

####
plot(u35, which="oxygen2", ztype="image", xtype="time", ylim=c(0,200))



###  OLD cod ####

####  create stations section - make a function -downcasts ####

files<-list.files(ctdata)
## station 1 - casts ##
c1<-read.ctd(file.path(ctdata, files[1]),) %>% ctdTrim(., method = "downcast")
c2<-read.ctd(file.path(ctdata, files[23]),) %>% ctdTrim(., method = "downcast")
c3<-read.ctd(file.path(ctdata, files[45]),) %>% ctdTrim(., method = "downcast")
c4<-read.ctd(file.path(ctdata, files[67]),) %>% ctdTrim(., method = "downcast")
c5<-read.ctd(file.path(ctdata, files[81]),) %>% ctdTrim(., method = "downcast")
c6<-read.ctd(file.path(ctdata, files[83]),) %>% ctdTrim(., method = "downcast")
c7<-read.ctd(file.path(ctdata, files[85]),) %>% ctdTrim(., method = "downcast")
c8<-read.ctd(file.path(ctdata, files[87]),) %>% ctdTrim(., method = "downcast")
#section
st1<-as.section(list(c1,c2,c3,c4,c5,c6,c7,c6,c8))

## station 2 - casts ##
c9<-read.ctd(file.path(ctdata, files[89]),) %>% ctdTrim(., method = "downcast")
c10<-read.ctd(file.path(ctdata, files[3]),) %>% ctdTrim(., method = "downcast")
c11<-read.ctd(file.path(ctdata, files[5]),) %>% ctdTrim(., method = "downcast")
c12<-read.ctd(file.path(ctdata, files[7]),) %>% ctdTrim(., method = "downcast")
c13<-read.ctd(file.path(ctdata, files[9]),) %>% ctdTrim(., method = "downcast")
c14<-read.ctd(file.path(ctdata, files[11]),) %>% ctdTrim(., method = "downcast")
c15<-read.ctd(file.path(ctdata, files[13]),) %>% ctdTrim(., method = "downcast")
c16<-read.ctd(file.path(ctdata, files[15]),) %>% ctdTrim(., method = "downcast")
c17<-read.ctd(file.path(ctdata, files[17]),) %>% ctdTrim(., method = "downcast")

#section
st2<-as.section(list(c9,c10,c11,c12,c13,c14,c15,c16,c17))

## station 3 - casts ##
c18<-read.ctd(file.path(ctdata, files[19]),) %>% ctdTrim(., method = "downcast")
c19<-read.ctd(file.path(ctdata, files[21]),) %>% ctdTrim(., method = "downcast")
c20<-read.ctd(file.path(ctdata, files[25]),) %>% ctdTrim(., method = "downcast")
c21<-read.ctd(file.path(ctdata, files[27]),) %>% ctdTrim(., method = "downcast")
c22<-read.ctd(file.path(ctdata, files[29]),) %>% ctdTrim(., method = "downcast")
c23<-read.ctd(file.path(ctdata, files[31]),) %>% ctdTrim(., method = "downcast")
c24<-read.ctd(file.path(ctdata, files[33]),) %>% ctdTrim(., method = "downcast")
c25<-read.ctd(file.path(ctdata, files[35]),) %>% ctdTrim(., method = "downcast") 
c26<-read.ctd(file.path(ctdata, files[37]),) %>% ctdTrim(., method = "downcast") 
c27<-read.ctd(file.path(ctdata, files[39]),) %>% ctdTrim(., method = "downcast") 
c28<-read.ctd(file.path(ctdata, files[41]),) %>% ctdTrim(., method = "downcast") 
c29<-read.ctd(file.path(ctdata, files[43]),) %>% ctdTrim(., method = "downcast")

#section
st3<-as.section(list(c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29))
#c30 is mis labelledas station 3

## station 3.5 - casts ##
c30<-read.ctd(file.path(ctdata, files[47]),) %>% ctdTrim(., method = "downcast")
c31<-read.ctd(file.path(ctdata, files[49]),) %>% ctdTrim(., method = "downcast")
c32<-read.ctd(file.path(ctdata, files[51]),) %>% ctdTrim(., method = "downcast")
c33<-read.ctd(file.path(ctdata, files[53]),) %>% ctdTrim(., method = "downcast")
c34<-read.ctd(file.path(ctdata, files[55]),) %>% ctdTrim(., method = "downcast")
#section
st3.5<-as.section(list(c30,c31,c32,c33,c34))

## station 4 - casts ##
c35<-read.ctd(file.path(ctdata, files[57]),) %>% ctdTrim(., method = "downcast")
c36<-read.ctd(file.path(ctdata, files[59]),) %>% ctdTrim(., method = "downcast")
c37<-read.ctd(file.path(ctdata, files[61]),) %>% ctdTrim(., method = "downcast")
c38<-read.ctd(file.path(ctdata, files[63]),) %>% ctdTrim(., method = "downcast")
c39<-read.ctd(file.path(ctdata, files[65]),) %>% ctdTrim(., method = "downcast")
c40<-read.ctd(file.path(ctdata, files[69]),) %>% ctdTrim(., method = "downcast")
#section
st4<-as.section(list(c35,c36,c37,c38,c39,c40))

## station 5 - casts ##
c41<-read.ctd(file.path(ctdata, files[71]),) %>% ctdTrim(., method = "downcast")
c42<-read.ctd(file.path(ctdata, files[73]),) %>% ctdTrim(., method = "downcast")
c43<-read.ctd(file.path(ctdata, files[75]),) %>% ctdTrim(., method = "downcast")
c44<-read.ctd(file.path(ctdata, files[77]),) %>% ctdTrim(., method = "downcast")
c45<-read.ctd(file.path(ctdata, files[79]),) %>% ctdTrim(., method = "downcast")
#section
st5<-as.section(list(c41,c42,c43,c44,c45))




#### larger sections ####

# all
all<-as.section(list(c1,c2,c3,c4,c5,c6,c7,c6,c9,c10,c11,c12,c13,c14,c15,c16,c17,
c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,
c38,c39,c40,c41,c42,c43,c44,c45))

# close in time and space 2-4
# sec.time<-as.section(list(c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,
# c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39,c40))

# triangle 2-3.5
sec.tri<-as.section(list(c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,
c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34))

# SE/NW omz "Y axis" stations 2, 3.5, 5 # come back here
sec.y<-as.section(list(c9,c10,c11,c12,c13,c14,c15,c16,c17,c30,c31,c32,c33,c34,c41,c42,c43,c44,c45))

# EN/WS omz "x axis" stations 1 3.5 3
sec.x<-as.section(list(c1,c2,c3,c4,c5,c6,c7,c8,c30,c31,c32,c33,c34,c18,c19,c20,c21,c22,c23,c24,
c25,c26,c27,c28,c29))





#### creating spatail dtaframe of centorids ####

namez<-c("time", "longitude", "latitude")

one<-cbind.data.frame(st1@metadata$time, st1@metadata$longitude, st1@metadata$latitude)
names(one)<-namez

two<-cbind.data.frame(st2@metadata$time, st2@metadata$longitude, st2@metadata$latitude)
names(two)<-namez

three<-cbind.data.frame(st3@metadata$time, st3@metadata$longitude, st3@metadata$latitude)
names(three)<-namez

three_five<-cbind.data.frame(st3.5@metadata$time, st3.5@metadata$longitude, st3.5@metadata$latitude)
names(three_five)<-namez

four<-cbind.data.frame(st4@metadata$time, st4@metadata$longitude, st4@metadata$latitude)
names(four)<-namez

five<-cbind.data.frame(st5@metadata$time, st5@metadata$longitude, st5@metadata$latitude)
names(five)<-namez

ctd_stations_list<-vector(mode = "list", length =6)
ctd_stations_list[[1]]<-one
ctd_stations_list[[2]]<-two
ctd_stations_list[[3]]<-three
ctd_stations_list[[4]]<-three_five
ctd_stations_list[[5]]<-four
ctd_stations_list[[6]]<-five

f.time<-function(x){
x$time<-with_tz(x$time, tz="US/Mountain")
return(x)}
ctd_stations_list<-lapply(ctd_stations_list, f.time)
names(ctd_stations_list)<-c("st1", "st2", "st3", "st3.5", "st4", "st5")




f.cent<-function(x){
y<-SpatialPoints(coords=cbind(x$longitude, x$latitude), proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
x<-gCentroid(y)
x<-cbind(x@coords[1], x@coords[2])
return(x)}

cent<-lapply(ctd_stations_list,f.cent)
cent<-data.frame(matrix(unlist(cent), nrow=length(cent), byrow=TRUE))
cent<-cbind(as.vector(c("1", "2", "3", "3.5", "4", "5" )), cent)
names(cent)<-c("station", "longitude", "latitude")

f.mean_time<-function(x){
x<-with_tz(x$time, tzone="UTC")%>%as.numeric%>%mean
x<-with_tz(with_tz(as_datetime(x), tzone="UTC"), tzone="US/Mountain")
return(x)}

x<-unlist(lapply(ctd_stations_list, f.mean_time))
x<-with_tz(with_tz(as_datetime(x), tzone="UTC"), tzone="US/Mountain")
cent$time<-x

cent_ctd<-st_as_sf(cent, coords = c("longitude", "latitude"),
          crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# load flow through start stop times
setwd(data_d)
stations<-read_csv("2018_time_and_station.csv")
names(stations)<-c("date", "hour", "station")
stations$station<-as.factor(stations$station)
stations[10,2]<-"0400"
stations$time<-mdy_hm(paste(stations$date, stations$hour, tz="UTC"))
stations$time<-with_tz(stations$time, tz="US/Mountain")
stations<-dplyr::select(stations, time, station)

stations<-dplyr::filter(stations, station != "Steam" & station != "Unknown")
stations<-droplevels(stations)
tz(stations$time)<-"US/Mountain"

start<-stations %>% group_by(station) %>% filter(time== min(time)) %>%
  slice(1) %>% ungroup()
start<-as.data.frame(start)
names(start)<-c("start_time", "station")

end<-stations %>% group_by(station) %>% filter(time== max(time)) %>%
  slice(1) %>% ungroup()
end<-as.data.frame(end)
names(end)<-c("end_time", "station")

x<-merge(start, end, by="station")
x$middle_time<-((x$end_time-x$start_time)/2 +x$start_time)

### clear and next step
rm(list=ls()[! ls() %in% c("x", "cent",  "robj", "fig", "gis_data", "wd", "stations", "data_d")])
stations_centroid_ctd<-cent
station_time<-x
rm(x, cent)

station_time<-with_tz(station_time, tzone="US/Mountain")
setwd(robj)
full<-readRDS('data.R')
pos<-readRDS("2018_gps.R")

full$hour<-round_date(full$time, unit="hour" )
match<-full %>% distinct(hour,  .keep_all = TRUE)
rm(full)

f.time_join<-function(f, t){
  #for joining data by "nearest" time/date in seconds
  #requires "data.table" package
  setDT(f)
  setDT(t)
  setkey(f, time)
  setkey(t, time)
  f<-f[t, roll='nearest']
}

pos<-f.time_join(pos, match)
pos<-select(pos, time, lat, lon)

f.time_join2<-function(pos, data){
  setDT(data)
  setDT(pos)
  setkey(data, time)
  setkey(pos, time)
  data<-data[pos, roll='nearest']
}

#o.names<-names(station_time)

# hard coding, needs to be a function....

names(station_time)[4]<-"time"
station_time<-f.time_join2(station_time,pos)
names(station_time)<-c("middle_time", "middle_lat", "middle_lon", names(station_time)[4:5], "time")
names(station_time)[5]<-"time"
station_time<-f.time_join2(station_time,pos)
names(station_time)<-c("start_time", "start_lat", "start_lon", names(station_time)[4:7], "time")
station_time<-f.time_join2(station_time,pos)
names(station_time)<-c("end_time", "end_lat", "end_lon", names(station_time)[4:10])

flowthrough_stations<-select(station_time, station, start_time, start_lat, start_lon, middle_time, middle_lat, middle_lon, end_time, end_lat, end_lon)
rm(station_time)


### fixing hardcoding issue
backup<-station_time

start<-dplyr::select(station_time, station, start_time)
names(start)[2]<-"time"
mid<-dplyr::select(station_time, station, middle_time)
names(mid)[2]<-"time"
end<-dplyr::select(station_time, station, end_time)
names(end)[2]<-"time"

start<-f.time_join2_flow(start,pos)
names(start)<-c("station", "start_time", "start_lat", "start_lon")

mid<-f.time_join2_flow(mid,pos)
names(mid)<-c("station", "middle_time", "middle_lat", "middle_lon")

end<-f.time_join2_flow(end,pos)
names(station_time)<-c("station", "end_time", "end_lat", "end_lon")



#### saving objects ###
#sp_list<-vector(mode="list", length=6)
stations_positions_list<-vector(mode="list", length=3)

read_me<-('list of stations of dataframes of positions and times from 2018 cruise.list object[[2]] "flowthrough_stations" is start, time and  middle position for flowthrough data. list object[[3]]  "stations_centroid_ctd" is centroid of ctd casts, and time at that position')

stations_positions_list[[1]]<-read_me
stations_positions_list[[2]]<-flowthrough_stations
stations_positions_list[[3]]<-stations_centroid_ctd
names(stations_positions_list)<-c("read_me", "flowthrough", "ctd_centroids")


rm(list=ls()[! ls() %in% c("stations_positions_list",  "robj", "fig", "gis_data", "wd", "data_d")])
setwd(robj)
saveRDS(stations_positions_list, "stations_positions_list.R")

#### Fix script here to redo ctd data analysis ####

####loadinng other data ####

setwd(robj)

temp<-readRDS("temp.R")
cl<-temp[[1]]
test<-temp[[2]]
mapped<-temp[[3]]
sites<-temp[[4]]
rm(temp)


#### marmap bathymetry ----
# stationbathy<-getNOAA.bathy(lon1=-120, lon2 = -100, lat1=15, lat2=26)
# odd white lines in final image?

# convert in gmt tools
# gmt grdconvert SRTM15+V2.1.nc -R-120/-100/14/26 omz.nc

setwd("/home/brandon/vestawd/wd_ideal_error/data")
bathy.stations<-raster("omz_stations.nc", varname="z")
bathy.stations<-as.bathy(bathy.stations)
setwd(robj)
saveRDS(bathy.stations, file="bathy.stations.rda")

setwd(robj)
bathy<-readRDS("bathy.rda")

ramp<-c("#2D3184", "#254289", "#1C518F", "#115F96", "#046C9C", "#0078A2", "#0084A7",
        "#068FAB", "#189AAF", "#28A4B3", "#45B6B8", "#54BEBA", "#62C5BC", "#70CCBD", "#7DD2BF", 
        "#8BD8C0", "#97DDC1", "#A3E2C3", "#AFE5C4", "#BAE9C6", "#C4EBC8", "#CDEECA", "#D6F0CD",
        "#DEF1D0", "#E4F2D3", "#EAF2D6", "#EEF2DA", "#F2F2DE", "#F3F1E4")

sites<-sf:::as_Spatial(sites)
labels<-layer(sp.text(coordinates(sites), txt = sites$station, pos = 1))

setwd(fig)
bmp("omz_stations_and_bathymetery.bmp", res=300, height=4800, width=4800)
plot(bathy.stations, bpal=ramp,  image=T, deep=c(-6500,0), shallow=c(-100,0),
     step=c(200,0), lwd=c(0.4,0.8), lty=c(1,1))
scaleBathy(bathy.stations, deg=1, x="topright")
points(sites, pch=21, col="orange",bg=col2alpha("red",.7),cex=7)
points(sites, pch=21, col="white", bg="white",cex=5)
text(x=coordinates(sites)[1:6], y=coordinates(sites)[7:12],labels=sites$station, cex=1.5)
dev.off()

#### cruise plots ----
par(mfrow=c(2,2))
par(oma=c(0,0,2,0))

# density plots
plot(all, which=3, ztype="image", xtype="time")
plot(all, which=3, ztype="image", xtype="track")
plot(all, which=3, ztype="image", xtype="latitude")
plot(all, which=3, ztype="image", xtype="longitude")
title("Density: Sigma-theta", outer=T)

# salinity plots
plot(all, which="salinity", ztype="image", xtype="time")
plot(all, which="salinity", ztype="image", xtype="track")
plot(all, which="salinity", ztype="image", xtype="latitude")
plot(all, which="salinity", ztype="image", xtype="longitude")
title("Salinity: PSU", outer=T)

# temperature plots
plot(all, which="temperature", ztype="image", xtype="time")
plot(all, which="temperature", ztype="image", xtype="track")
plot(all, which="temperature", ztype="image", xtype="latitude")
plot(all, which="temperature", ztype="image", xtype="longitude")
title("Temperature: C", outer=T)

# oxygen2 plots
plot(all, which="oxygen2", ztype="image", xtype="time")
plot(all, which="oxygen2", ztype="image", xtype="track")
plot(all, which="oxygen2", ztype="image", xtype="latitude")
plot(all, which="oxygen2", ztype="image", xtype="longitude")
title("Oxygen: µmol/kg", outer=T)

# fluorescence plots
plot(all, which="fluorescence", ztype="image", xtype="time")
plot(all, which="fluorescence", ztype="image", xtype="track")
plot(all, which="fluorescence", ztype="image", xtype="latitude")
plot(all, which="fluorescence", ztype="image", xtype="longitude")
title("Fluorescence: mg/m^3", outer=T)

# PAR plots
plot(all, which="par", ztype="image", xtype="time")
plot(all, which="par", ztype="image", xtype="track")
plot(all, which="par", ztype="image", xtype="latitude")
plot(all, which="par", ztype="image", xtype="longitude")
title("PAR:  µmol/m2/s ?", outer=T)
#double check units here

## section plots
par(mfrow=c(2,4))
par(oma=c(0,0,2,0))
# Clump
plot(sec.time, drawIsobaths=T, which="map")
plot(sec.time, which=3, ztype="image", xtype="time")
plot(sec.time, which=3, ztype="image", xtype="longitude")
plot(sec.time, which=3, ztype="image", xtype="latitude")
plot(sec.time, which="salinity", ztype="image", xtype="time")
plot(sec.time, which="temperature", ztype="image", xtype="time")
plot(sec.time, which="oxygen2", ztype="image", xtype="time")
plot(sec.time, which="fluorescence", ztype="image", xtype="time")
title("Central: Stations 2-4 - Time plots", outer=T)

##### map of stations #####
# slides per station

m.coast <- ggplot()+ geom_sf(data =cl, size=2)+ theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme_bw()+ theme(text = element_text(size = 10)) 

m.c<-m.coast+ theme(legend.title = element_blank())+
  ggtitle("Stations") + theme(plot.title = element_text(hjust = 0.5))+
  geom_sf(data =sitesa, shape =2, size =3) + geom_sf_text(aes(label = station), hjust=3, data= sitesa, size = 3) +
  geom_sf(data =sitesb, shape =2, size =3) + geom_sf_text(aes(label = station), vjust=2, data= sitesb, size = 3) +
  geom_sf(data =sitesc, shape =2, size =3) + geom_sf_text(aes(label = station), hjust=1.75, data= sitesc, size = 3)+
  geom_sf(data =sitesd, shape =2, size =3) + geom_sf_text(aes(label = station), hjust=2.5, data= sitesd, size = 3) +
  geom_sf(data =sitese, shape =2, size =3) + geom_sf_text(aes(label = station), hjust=3, data= sitese,size = 3)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 

par(mfrow=c(2,1))

#cast times plotted
times<-all@metadata$time
bud$lat<-as.numeric(str_extract_all(bud$lat, '[:digit:][:digit:][:punct:][:digit:]', simplify = T))
times
times<-mdy_hms(times)
times<-with_tz(times, tz="US/Mountain")


par(mfrow=c(2,3))
#station 1
plot(st1, drawIsobaths=T, which="map")
plot(st1, which=3, ztype="image", xtype="time")
plot(st1, which="salinity", ztype="image", xtype="time")
plot(st1, which="temperature", ztype="image", xtype="time")
plot(st1, which="oxygen2", ztype="image", xtype="time")
plot(st1, which="fluorescence", ztype="image", xtype="time")

#station 2
plot(st2, drawIsobaths=T, which="map")
plot(st2, which=3, ztype="image", xtype="time")
plot(st2, which="salinity", ztype="image", xtype="time")
plot(st2, which="temperature", ztype="image", xtype="time")
plot(st2, which="oxygen2", ztype="image", xtype="time")
plot(st2, which="fluorescence", ztype="image", xtype="time")

#station 3
plot(st3, drawIsobaths=T, which="map")
plot(st3, which=3, ztype="image", xtype="time")
plot(st3, which="salinity", ztype="image", xtype="time")
plot(st3, which="temperature", ztype="image", xtype="time")
plot(st3, which="oxygen2", ztype="image", xtype="time")
plot(st3, which="fluorescence", ztype="image", xtype="time")

# station 3.5
plot(st3.5, drawIsobaths=T, which="map")
plot(st3.5, which=3, ztype="image", xtype="time")
plot(st3.5, which="salinity", ztype="image", xtype="time")
plot(st3.5, which="temperature", ztype="image", xtype="time")
plot(st3.5, which="oxygen2", ztype="image", xtype="time")
plot(st3.5, which="fluorescence", ztype="image", xtype="time")

# station 4
plot(st4, drawIsobaths=T, which="map")
plot(st4, which=3, ztype="image", xtype="time")
plot(st4, which="salinity", ztype="image", xtype="time")
plot(st4, which="temperature", ztype="image", xtype="time")
plot(st4, which="oxygen2", ztype="image", xtype="time")
plot(st4, which="fluorescence", ztype="image", xtype="time")

# station 5
plot(st5, drawIsobaths=T, which="map")
plot(st5, which=3, ztype="image", xtype="time")
plot(st5, which="salinity", ztype="image", xtype="time")
plot(st5, which="temperature", ztype="image", xtype="time")
plot(st5, which="oxygen2", ztype="image", xtype="time")
plot(st5, which="fluorescence", ztype="image", xtype="time")


### Mega time plots ####

# ToDO add surface tracks map at station

par(mfrow=c(2,6))
#station 1
plot(c5, drawIsobaths=T, which="map")
plot(st1, which=3, ztype="image", xtype="time")
plot(st1, which="salinity", ztype="image", xtype="time")
plot(st1, which="temperature", ztype="image", xtype="time")
plot(st1, which="oxygen2", ztype="image", xtype="time")
plot(st1, which="fluorescence", ztype="image", xtype="time")

#station 2
plot(c13, drawIsobaths=T, which="map")
plot(st2, which=3, ztype="image", xtype="time")
plot(st2, which="salinity", ztype="image", xtype="time")
plot(st2, which="temperature", ztype="image", xtype="time")
plot(st2, which="oxygen2", ztype="image", xtype="time")
plot(st2, which="fluorescence", ztype="image", xtype="time")

#station 3
plot(c24, drawIsobaths=T, which="map")
plot(st3, which=3, ztype="image", xtype="time")
plot(st3, which="salinity", ztype="image", xtype="time")
plot(st3, which="temperature", ztype="image", xtype="time")
plot(st3, which="oxygen2", ztype="image", xtype="time")
plot(st3, which="fluorescence", ztype="image", xtype="time")

# station 3.5
plot(c32, drawIsobaths=T, which="map")
plot(st3.5, which=3, ztype="image", xtype="time")
plot(st3.5, which="salinity", ztype="image", xtype="time")
plot(st3.5, which="temperature", ztype="image", xtype="time")
plot(st3.5, which="oxygen2", ztype="image", xtype="time")
plot(st3.5, which="fluorescence", ztype="image", xtype="time")

# station 4
plot(c38, drawIsobaths=T, which="map")
plot(st4, which=3, ztype="image", xtype="time")
plot(st4, which="salinity", ztype="image", xtype="time")
plot(st4, which="temperature", ztype="image", xtype="time")
plot(st4, which="oxygen2", ztype="image", xtype="time")
plot(st4, which="fluorescence", ztype="image", xtype="time")

# station 5
plot(c43, drawIsobaths=T, which="map")
plot(st5, which=3, ztype="image", xtype="time")
plot(st5, which="salinity", ztype="image", xtype="time")
plot(st5, which="temperature", ztype="image", xtype="time")
plot(st5, which="oxygen2", ztype="image", xtype="time")
plot(st5, which="fluorescence", ztype="image", xtype="time")


#### maps of spatial distributions ####
par(mfrow=c(2,3))
plot(st1, drawIsobaths=T, which="map")
title("Station 1")
plot(st2, drawIsobaths=T, which="map")
title("Station 2")
plot(st3, drawIsobaths=T, which="map")
title("Station 3")
plot(st3.5, drawIsobaths=T, which="map")
title("Station 3.5")
plot(st4, drawIsobaths=T, which="map")
title("Station 4")
plot(st5, drawIsobaths=T, which="map")
title("Station 5")

#### profiles per station ####

par(mfrow=c(2,4))
par(oma=c(0,0,2,0))
plotProfile(c1, ytype="depth")
plotProfile(c2, ytype="depth")
plotProfile(c3, ytype="depth")
plotProfile(c4, ytype="depth")
plotProfile(c5, ytype="depth")
plotProfile(c6, ytype="depth")
plotProfile(c7, ytype="depth")
plotProfile(c8, ytype="depth")
title("Station 1", outer=T)

par(mfrow=c(3,3))
par(oma=c(0,0,2,0))
plotProfile(c9, ytype="depth")
plotProfile(c10, ytype="depth")
plotProfile(c11, ytype="depth")
plotProfile(c12, ytype="depth")
plotProfile(c13, ytype="depth")
plotProfile(c14, ytype="depth")
plotProfile(c15, ytype="depth")
plotProfile(c16, ytype="depth")
plotProfile(c17, ytype="depth")
title("Station 2", outer=T)

par(mfrow=c(3,4))
par(oma=c(0,0,2,0))
plotProfile(c18, ytype="depth")
plotProfile(c19, ytype="depth")
plotProfile(c20, ytype="depth")
plotProfile(c21, ytype="depth")
plotProfile(c22, ytype="depth")
plotProfile(c23, ytype="depth")
plotProfile(c24, ytype="depth")
plotProfile(c25, ytype="depth")
plotProfile(c26, ytype="depth")
plotProfile(c27, ytype="depth")
plotProfile(c28, ytype="depth")
plotProfile(c29, ytype="depth")
title("Station 3", outer =T)

par(mfrow=c(2,3))
par(oma=c(0,0,2,0))
plotProfile(c30, ytype="depth")
plotProfile(c31, ytype="depth")
plotProfile(c32, ytype="depth")
plotProfile(c33, ytype="depth")
plotProfile(c34, ytype="depth")
title("Station 3.5", outer=T)

par(mfrow=c(2,3))
par(oma=c(0,0,2,0))
plotProfile(c35, ytype="depth")
plotProfile(c36, ytype="depth")
plotProfile(c37, ytype="depth")
plotProfile(c38, ytype="depth")
plotProfile(c39, ytype="depth")
plotProfile(c40, ytype="depth")
title("Station 4", outer=T)

par(mfrow=c(2,3))
par(oma=c(0,0,2,0))
plotProfile(c41, ytype="depth")
plotProfile(c42, ytype="depth")
plotProfile(c43, ytype="depth")
plotProfile(c44, ytype="depth")
plotProfile(c45, ytype="depth")
title("Station 5", outer =T)

#### cruise  surface plots 0-200 M####

#time
par(mfrow=c(4,2))
plot(all, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(all, which="temperature", ztype="image", xtype="time", ylim=c(0,200))
plot(all, which="oxygen2", ztype="image", xtype="time", ylim=c(0,200))
plot(all, which="fluorescence", ztype="image", xtype="time", ylim=c(0,200))

#track

plot(all, which=3, ztype="image", xtype="track", ylim=c(0,200))
plot(all, which="temperature", ztype="image", xtype="track", ylim=c(0,200))
plot(all, which="oxygen2", ztype="image", xtype="track", ylim=c(0,200))
plot(all, which="fluorescence", ztype="image", xtype="track", ylim=c(0,200))

#latitude

plot(all, which=3, ztype="image", xtype="latitude", ylim=c(0,200))
plot(all, which="temperature", ztype="image", xtype="latitude", ylim=c(0,200))
plot(all, which="oxygen2", ztype="image", xtype="latitude", ylim=c(0,200))
plot(all, which="fluorescence", ztype="image", xtype="latitude", ylim=c(0,200))

#longitude

plot(all, which=3, ztype="image", xtype="longitude", ylim=c(0,200))
plot(all, which="temperature", ztype="image", xtype="longitude", ylim=c(0,200))
plot(all, which="oxygen2", ztype="image", xtype="longitude", ylim=c(0,200))
plot(all, which="fluorescence", ztype="image", xtype="longitude", ylim=c(0,200))

#### station surface plots 0-200 M####


par(mfrow=c(2,5))
#station 1

plot(st1, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(st1, which="salinity", ztype="image", xtype="time", ylim=c(0,200))
plot(st1, which="temperature", ztype="image", xtype="time", ylim=c(0,200))
plot(st1, which="oxygen2", ztype="image", xtype="time", ylim=c(0,200))
plot(st1, which="fluorescence", ztype="image", xtype="time",ylim=c(0,200))

#station 2
plot(st2, which=3, ztype="image", xtype="time",ylim=c(0,200))
plot(st2, which="salinity", ztype="image", xtype="time",ylim=c(0,200))
plot(st2, which="temperature", ztype="image", xtype="time",ylim=c(0,200))
plot(st2, which="oxygen2", ztype="image", xtype="time",ylim=c(0,200))
plot(st2, which="fluorescence", ztype="image", xtype="time", ylim=c(0,200))

#station 3
plot(st3, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(st3, which="salinity", ztype="image", xtype="time", ylim=c(0,200))
plot(st3, which="temperature", ztype="image", xtype="time", ylim=c(0,200))
plot(st3, which="oxygen2", ztype="image", xtype="time",ylim=c(0,200))
plot(st3, which="fluorescence", ztype="image", xtype="time",ylim=c(0,200))

# station 3.5

plot(st3.5, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(st3.5, which="salinity", ztype="image", xtype="time", ylim=c(0,200))
plot(st3.5, which="temperature", ztype="image", xtype="time", ylim=c(0,200))
plot(st3.5, which="oxygen2", ztype="image", xtype="time", ylim=c(0,200))
plot(st3.5, which="fluorescence", ztype="image", xtype="time", ylim=c(0,200))

# station 4
plot(st4, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(st4, which="salinity", ztype="image", xtype="time", ylim=c(0,200))
plot(st4, which="temperature", ztype="image", xtype="time", ylim=c(0,200))
plot(st4, which="oxygen2", ztype="image", xtype="time", ylim=c(0,200))
plot(st4, which="fluorescence", ztype="image", xtype="time", ylim=c(0,200))

# station 5
plot(st5, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(st5, which="salinity", ztype="image", xtype="time", ylim=c(0,200))
plot(st5, which="temperature", ztype="image", xtype="time", ylim=c(0,200))
plot(st5, which="oxygen2", ztype="image", xtype="time", ylim=c(0,200))
plot(st5, which="fluorescence", ztype="image", xtype="time", ylim=c(0,200))

#### other section surface plots 0-220m ####

# triangle
par(mfrow=c(2,4))
plot(sec.tri, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(sec.tri, which=3, ztype="image", xtype="track", ylim=c(0,200))
plot(sec.tri, which=3, ztype="image", xtype="longitude", ylim=c(0,200))
plot(sec.tri, which=3, ztype="image", xtype="latitude", ylim=c(0,200))

plot(sec.tri, which="salinity", ztype="image", xtype="time", ylim=c(0,200))
plot(sec.time, which="temperature", ztype="image", xtype="time", ylim=c(0,200))
plot(sec.time, which="oxygen2", ztype="image", xtype="time", ylim=c(0,200))
plot(sec.time, which="fluorescence", ztype="image", xtype="time", ylim=c(0,200))
#title("Central: Stations 2-4 - Time plots", outer=T)

# y axis

par(mfrow=c(2,4))
plot(sec.y, which="map")
plot(sec.y, which=3, ztype="image", xtype="distance", ylim=c(0,200))
plot(sec.y, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(sec.y, which=3, ztype="image", xtype="longitude", ylim=c(0,200))

plot(sec.y, which="salinity", ztype="image", xtype="time", ylim=c(0,200))
plot(sec.y, which="temperature", ztype="image", xtype="time", ylim=c(0,200))
plot(sec.y, which="oxygen2", ztype="image", xtype="time", ylim=c(0,200))
plot(sec.y, which="fluorescence", ztype="image", xtype="time", ylim=c(0,200))

# x axis
par(mfrow=c(2,4))
plot(sec.x, which="map")
plot(sec.x, which=3, ztype="image", xtype="distance", ylim=c(0,200))
plot(sec.x, which=3, ztype="image", xtype="time", ylim=c(0,200))
plot(sec.x, which=3, ztype="image", xtype="latitude", ylim=c(0,200))

plot(sec.x, which="salinity", ztype="image", xtype="distance", ylim=c(0,200))
plot(sec.x, which="temperature", ztype="image", xtype="distance", ylim=c(0,200))
plot(sec.x, which="oxygen2", ztype="image", xtype="distance", ylim=c(0,200))
plot(sec.x, which="fluorescence", ztype="image", xtype="distance", ylim=c(0,200))

