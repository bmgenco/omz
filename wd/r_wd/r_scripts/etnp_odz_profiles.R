#  crusie profiles from entmp_od

# to do:
# ADD SKQ2016 to the Station 3.5 plot, 
# remove SKQ2016 and TT66 from the Station 4 plot.   

rm(list=ls())

#### directories ####

wd<-"/home/brandon/vestawd/omz/wd/r_wd"
# wd<-"/home/brandon/callistowd/omz/wd/r_wd"
setwd(wd)  
robj<-"r_objects"
fig<-"../../figures"
other_crusies<-"../../data/additonal_cruise_profiles/ENTP_ODZ_time_series/"




#### packages ####
f.ipak <- function(pkg){
  # loads packages, quietly, given by a vector of package names e.g., pkg<-c("ggplot", "tidyverse")
  # will install  packages listed , and their dependencies, if needed.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, quiet=T, verbose = F)
  sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts=F)
}

packages<-c("sp", "rgdal",  "rgeos", "raster", "readr", "tidyverse", "lubridate",  
            "ggthemes",  "sf", "cmocean",   "plot3D", "tidync", "devtools", 
            "stars", "ncmeta", "maps", "oce", "data.table", "fasterize", "RStoolbox", "scales", "purrr", "nngeo", "readxl", "stringr")

# install.packages("remotes")
# remotes::install_github("decisionpatterns/stringr.tools")

f.ipak(packages)
library("stringr.tools")

### load staions centorids:
setwd(wd)
setwd(robj)
stations_positions_list<-readRDS("OC1806A_stations_positions_list.R")
stations<-sf::st_as_sf(stations_positions_list$ctd_centroids, coords =c("longitude", "latitude"),  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
rm(stations_positions_list)

### our profiles ####


#### entnp ODZ ####
setwd(wd)
setwd(other_crusies)

# data_old<-read_delim("eOMP/Inputs/All 8 Cruises_adjusted_omp.csv", col_names = T, delim =",")
# data_older<-read_excel("Integration calculations/allCruises_updated_plus_PotDens_20220405.xlsx", col_names = T)
# data2<-read.ctd("eOMP/Inputs/All 8 Cruises_adjusted_omp.csv")

data<-read_delim("All 8 Cruises_adjusted including KM1919 .csv", skip=3, col_names = T, delim =",")
data<-Filter(function(x)!all(is.na(x)), data) # reove emopt data columns
metadata<-read_xlsx("ETNP ODZ file metadata.xlsx") %>% select(., Year, 'Start month', Shorthand)

# hack
metadata$month<-c("02","01","12","03","11","12","03","09")
metadata$Year<-as.character(metadata$Year)
metadata$month<-str_prefix(metadata$month, "-")
metadata$date<-str_postfix(metadata$Year, metadata$month )
metadata<-select(metadata, date, Shorthand)
names(metadata)<-c("date", "cruise")



## for use in remote working stations:
setwd(wd)
setwd(robj)
data<-readRDS("etnp_odz_vertical_profiles.R")

############################  notes on data: ############################ 
# "this dataset set O2 = 0 when O2 < 10 μmol kg-1 and NO2 > 0.1 μmol kg-1."
#
#

cid<-unique(data$Cruise)
cid<-as.vector(na.omit(cid))

#fix discrepanices  between metadatshorthand and cid
metadata[which(metadata$cruise =="P_18"),2]<-"P-18"
metadata[which(metadata$cruise =="RR1804"),2]<-"RR_1804"
metadata[which(metadata$cruise =="Clivar"),2]<-"CLIVAR"

metadata$date<-str_prefix(metadata$date, " ")


cruises<-vector(mode="list", length = length(cid))
names(cruises)<-cid
for(i in 1:length(cid)){
cruises[[i]]<-as.data.frame(filter(data, Cruise ==cid[[i]]  ))}


# f.closest<-function(x,y){
# #This vers does not return correct number profoeils from cruise  RR_1804.. becasue of sligt varitions in lat coordiantes for this cruise/station/pt
# q<-st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# z<-st_nn(y, q, k = 1, returnDist = F, sparse =T)
# nearest<-st_coordinates(q[unlist(z),])
# x<-filter(x, Latitude == nearest[2], Longitude == nearest[1])
# }

# update to find station 3
f.closest<-function(x,y){
  pts<-unique(x[c("Latitude","Longitude")])
  q<-st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  z<-st_nn(y, q, k = 1, returnDist = F, sparse =T)
  nearest<-st_coordinates(q[unlist(z),])
  z<-filter(x, Latitude == nearest[2], Longitude == nearest[1])
  q<-unique(z$Station)
  x<-filter(x, Station==q)
}

# add dates to name, a bit of a hack..
for(i in 1:length(cruises)){
  x<-names(cruises)[i]
  y<-which(metadata$cruise==x)
  z<-str_prefix(metadata$date[y], metadata$cruise[y])
  names(cruises)[i]<-z}

names(cruises)[which(names(cruises)=="RR_1804 2018-03")]<-"RR1804 2018-03"
names(cruises)[which(names(cruises)=="P-18 2016-11")]<-"P1816 2016-11"


y<-stations[5,]$geometry
odz_st4<-lapply(cruises, f.closest, y=y)


y<-stations[4,]$geometry
odz_st3.5<-lapply(cruises, f.closest, y=y)


#### testinf for stations etc
# i<-8
# test<-as.data.frame(filter(data, Cruise ==cid[[i]]))
# view(test)
# view(odz_st4[[i]])

#### as. ctd

# saveRDS(data, "etnp_odz_vertical_profiles.R")

setwd(wd)
setwd(robj)
d.sec<-readRDS("OC1806A_ctd_sections.R")
d4<-d.sec$st4


t1<-d4@data$station[1][[1]]
t2<-d4@data$station[2][[1]]
t3<-d4@data$station[3][[1]]
t4<-d4@data$station[4][[1]]
t5<-d4@data$station[5][[1]]
t6<-d4@data$station[6][[1]]

f.ctd<-function(x){
names(x)<-tolower(names(x))
x<-as.oce(x) %>% as.ctd(.,)}

y<-lapply(odz_st4, f.ctd)

names(y)

y<-within(y, rm('TT66 1972-02', 'SKQ2016 2016-12'))



ramp<-c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffcc33", "#a65628", "#f781bf", "#575757")
ramp<-ramp[-7]
ramp<-ramp[-7]
dm<-400
# removing objects from list


setwd(wd)
setwd(fig)
pdf("Station_4.pdf", height = 10, width = 7)

# plotProfile(y$TT66, xtype = "oxygen", ytype="depth", col=ramp[1], ylim=c(500,0),xlim=c(5,260))
plotProfile(y$'WOCE 1994-01', xtype = "oxygen", ytype="depth", col=ramp[1],ylim=c(dm,0), xlim=c(5,260))
plotProfile(y$'CLIVAR 2007-12', xtype = "oxygen", ytype="depth", col=ramp[2], ylim=c(dm,0),add=T)
plotProfile(y$'TN278 2012-03', xtype = "oxygen", ytype="depth", col=ramp[3], ylim=c(dm,0), add=T)
# plotProfile(y$SKQ2016, xtype = "oxygen", ytype="depth", col=ramp[5], ylim=c(500,0), add=T)
plotProfile(y$`P1816 2016-11`, xtype = "oxygen", ytype="depth", col=ramp[4], ylim=c(dm,0), add=T)
plotProfile(y$'RR1804 2018-03', xtype = "oxygen", ytype="depth", col=ramp[5], ylim=c(dm,0), add=T)
plotProfile(y$'KM1919 2019-09', xtype = "oxygen", ytype="depth", col=ramp[6], ylim=c(dm,0), add=T)
plotProfile(t1, xtype = "oxygen2", ytype="depth", col=ramp[7], ylim=c(dm,0), add=T)
plotProfile(t2, xtype = "oxygen2", ytype="depth", col=ramp[7], ylim=c(dm,0), add=T)
plotProfile(t3, xtype = "oxygen2", ytype="depth", col=ramp[7], ylim=c(dm,0), add=T)
plotProfile(t4, xtype = "oxygen2", ytype="depth", col=ramp[7], ylim=c(dm,0), add=T)
plotProfile(t5, xtype = "oxygen2", ytype="depth", col=ramp[7], ylim=c(dm,0), add=T)
plotProfile(t6, xtype = "oxygen2", ytype="depth", col=ramp[7], ylim=c(dm,0), add=T)
legend(175, 300, legend=c(names(y), "Station 4"),col=ramp, lty=1, cex=0.95)
dev.off()

# station 3.5

y<-lapply(odz_st3.5, f.ctd)
ramp<-c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#575757")

# resorted names and color patterm

t<-names(y)
n<-c(t[2], t[3], t[4], t[6], t[7], t[8], t[5], t[1])

plotProfile(y$'WOCE 1994-01', xtype = "oxygen", ytype="depth", col=ramp[1],ylim=c(dm,0), xlim=c(5,260))
plotProfile(y$'CLIVAR 2007-12', xtype = "oxygen", ytype="depth", col=ramp[2], ylim=c(dm,0),add=T)
plotProfile(y$'TN278 2012-03', xtype = "oxygen", ytype="depth", col=ramp[3], ylim=c(dm,0), add=T)
# plotProfile(y$SKQ2016, xtype = "oxygen", ytype="depth", col=ramp[5], ylim=c(500,0), add=T)
plotProfile(y$`P1816 2016-11`, xtype = "oxygen", ytype="depth", col=ramp[4], ylim=c(dm,0), add=T)
plotProfile(y$'RR1804 2018-03', xtype = "oxygen", ytype="depth", col=ramp[5], ylim=c(dm,0), add=T)
plotProfile(y$'KM1919 2019-09', xtype = "oxygen", ytype="depth", col=ramp[7], ylim=c(dm,0), add=T)

plotProfile(y$'SKQ2016 2016-12', xtype = "oxygen", ytype="depth", col=ramp[8], ylim=c(dm,0),add=T)
plotProfile(y$'TT66 1972-02', xtype = "oxygen", ytype="depth", col=ramp[6], ylim=c(dm,0),add=T)

plotProfile(t1, xtype = "oxygen2", ytype="depth", col=ramp[9], ylim=c(dm,0), add=T)
plotProfile(t2, xtype = "oxygen2", ytype="depth", col=ramp[9], ylim=c(dm,0), add=T)
plotProfile(t3, xtype = "oxygen2", ytype="depth", col=ramp[9], ylim=c(dm,0), add=T)
plotProfile(t4, xtype = "oxygen2", ytype="depth", col=ramp[9], ylim=c(dm,0), add=T)
plotProfile(t5, xtype = "oxygen2", ytype="depth", col=ramp[9], ylim=c(dm,0), add=T)
plotProfile(t6, xtype = "oxygen2", ytype="depth", col=ramp[9], ylim=c(dm,0), add=T)
legend(175, 300, legend=c(n, "Station 3.5"),col=ramp, lty=1, cex=0.95)


## other stufff... ( not sure)

lat<-lapply(y, function(x){x<-unlist(unique(x@data$latitude)[[1]])})
lon<-lapply(y, function(x){x<-unlist(unique(x@data$longitude)[[1]])})
st<-lapply(y, function(x){x<-unlist(unique(x@data$station))})

text("Station 4")

plotProfile(d4, xtype = "oxygen2", ytype="depth", col=ramp[7], add=T)

plotProfile(d4, xtype = "oxygen", ytype="depth", col=ramp[9])
plot(d4, which="oxygen2", ztype="image", xtype="time")


# Making average


ml<-length(t1@data$pressure)

l2=ml-length(t2@data$pressure)
a2<-(rep(NA, l2))
p2<-append(t2@data$pressure, a2)
rm(l2,a2)

l3=ml-length(t3@data$pressure)
a3<-(rep(NA, l3))
p3<-append(t3@data$pressure, a3)
rm(l3,a3)

l4=ml-length(t4@data$pressure)
a4<-(rep(NA, l4))
p4<-append(t4@data$pressure, a4)
rm(l4,a4)

l5=ml-length(t5@data$pressure)
a5<-(rep(NA, l5))
p5<-append(t5@data$pressure, a5)
rm(l5,a5)

l6=ml-length(t6@data$pressure)
a6<-(rep(NA, l6))
p6<-append(t6@data$pressure, a6)
rm(l6,a6)



p<-apply(cbind(t1@data$pressure, p2, p3, p4, p5, p6), MARGIN=1, FUN=median, na.rm = FALSE)

plotProfile(d4, xtype = "oxygen", ytype="depth", col=ramp[7], ylim=c(200,0))

