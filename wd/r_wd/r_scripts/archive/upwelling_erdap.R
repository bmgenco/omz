# see here: 
# https://oceanview.pfeg.noaa.gov/products/upwelling/bakun

library(plotdap)


#### set directories #####

wd<-"/home/brandon/vestawd/omz/wd/r_wd"
robj<-"/home/brandon/vestawd/omz/wd/r_wd/r_objects"
fig<-"/home/brandon/vestawd/omz/wd/r_wd/figures"
gis_data<-"/home/brandon/vestawd/omz/data/gis_data"
data_d<-"/home/brandon/vestawd/omz/data/2018_data"


f.upwell <- function(ektrx, ektry, coast_angle) {
  pi <- 3.1415927
  degtorad <- pi/180.
  alpha <- (360 - coast_angle) * degtorad
  s1 <- cos(alpha)
  t1 <- sin(alpha)
  s2 <- -1 * t1
  t2 <- s1
  perp <- (s1 * ektrx) + (t1 * ektry)
  para <- (s2 * ektrx) + (t2 * ektry)
  return(perp/10)
}

# min(data$lat)
# max(data$lat)
# min(data$lon)
# max(data$lon)

#data from:
#https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdlasFnTran6.html

setwd(data_d)
files<-list.files()
x<-read.csv2(files[6], header = F, skip = 2, sep=",")
metadata<-read.table(files[6], header = F, nrows=2, sep=",")
names(x)<-as.vector(metadata[1,], mode="character")

