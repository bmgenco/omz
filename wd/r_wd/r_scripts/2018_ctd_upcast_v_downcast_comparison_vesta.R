# duplicate of ctd code from"ctd_2018_2_vesta.R"
#ctd comparison

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
# rm(list=ls())


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


## vesta directoires:
wd<-"/home/brandon/vestawd/omz/wd/r_wd"
robj<-"/home/brandon/vestawd/omz/wd/r_wd/r_objects"
fig<-"/home/brandon/vestawd/omz/figures"
gis_data<-"/home/brandon/vestawd/omz/data/gis_data"
data_d<-"/home/brandon/vestawd/omz/data/2018_data/"
ctdata<-"/home/brandon/vestawd/omz/data/2018_data/ctd/data"


#### load r object
setwd(robj)
x<-readRDS("2018_ctd_comparison_up_v_downcasts.R")
x[[1]]->d2
x[[2]]->d35
x[[3]]->d4
x[[4]]->u2
x[[5]]->u35
x[[6]]->u4


#### combined 

files<-list.files(ctdata)

## station 2 - casts ##
c9<-read.ctd(file.path(ctdata, files[89]),) 
c10<-read.ctd(file.path(ctdata, files[3]),) 
c11<-read.ctd(file.path(ctdata, files[5]),) 
c12<-read.ctd(file.path(ctdata, files[7]),) 
c13<-read.ctd(file.path(ctdata, files[9]),) 
c14<-read.ctd(file.path(ctdata, files[11]),) 
c15<-read.ctd(file.path(ctdata, files[13]),) 
c16<-read.ctd(file.path(ctdata, files[15]),) 
c17<-read.ctd(file.path(ctdata, files[17]),) 

#section
st2<-as.section(list(c9,c10,c11,c12,c13,c14,c15,c16,c17))

#c30 is mis labelledas station 3

## station 3.5 - casts ##
c30<-read.ctd(file.path(ctdata, files[47]),) 
c31<-read.ctd(file.path(ctdata, files[49]),) 
c32<-read.ctd(file.path(ctdata, files[51]),) 
c33<-read.ctd(file.path(ctdata, files[53]),) 
c34<-read.ctd(file.path(ctdata, files[55]),) 
#section
st35<-as.section(list(c30,c31,c32,c33,c34))

## station 4 - casts ##
c35<-read.ctd(file.path(ctdata, files[57]),) 
c36<-read.ctd(file.path(ctdata, files[59]),) 
c37<-read.ctd(file.path(ctdata, files[61]),) 
c38<-read.ctd(file.path(ctdata, files[63]),) 
c39<-read.ctd(file.path(ctdata, files[65]),) 
c40<-read.ctd(file.path(ctdata, files[69]),) 
#section
st4<-as.section(list(c35,c36,c37,c38,c39,c40))

#### plotting ###
plot(c13, which="oxygen2")
p
plot(st35, which="oxygen2", depth)
plot(u35, which="oxygen2")


#data combination:
up_35 <- unlist(lapply(u35, function(x) x[['oxygen2']]))
y<-ctdDecimate(u35)

u35_1<-cbind(u35@data$station[[1]]@data$oxygen2,u35@data$station[[1]]@data$depth)  
u35_1<-as.data.frame(u35_1)
names(u35_1)<-c("o2", "depth")

u35_2<-cbind(u35@data$station[[2]]@data$oxygen2,u35@data$station[[2]]@data$depth)  
u35_2<-as.data.frame(u35_2)
names(u35_2)<-c("o2", "depth")

u35_3<-cbind(u35@data$station[[3]]@data$oxygen2,u35@data$station[[3]]@data$depth)  
u35_3<-as.data.frame(u35_3)
names(u35_3)<-c("o2", "depth")

u35_4<-cbind(u35@data$station[[4]]@data$oxygen2,u35@data$station[[4]]@data$depth)  
u35_4<-as.data.frame(u35_4)
names(u35_4)<-c("o2", "depth")

u35_5<-cbind(u35@data$station[[5]]@data$oxygen2,u35@data$station[[5]]@data$depth)  
u35_5<-as.data.frame(u35_5)
names(u35_5)<-c("o2", "depth")


d35_1<-cbind(d35@data$station[[1]]@data$oxygen2,d35@data$station[[1]]@data$depth)  
d35_1<-as.data.frame(d35_1)
names(d35_1)<-c("o2", "depth")

d35_2<-cbind(d35@data$station[[2]]@data$oxygen2,d35@data$station[[2]]@data$depth)  
d35_2<-as.data.frame(d35_2)
names(d35_2)<-c("o2", "depth")

d35_3<-cbind(d35@data$station[[3]]@data$oxygen2,d35@data$station[[3]]@data$depth)  
d35_3<-as.data.frame(d35_3)
names(d35_3)<-c("o2", "depth")

d35_4<-cbind(d35@data$station[[4]]@data$oxygen2,d35@data$station[[4]]@data$depth)  
d35_4<-as.data.frame(d35_4)
names(d35_4)<-c("o2", "depth")

d35_5<-cbind(d35@data$station[[5]]@data$oxygen2,d35@data$station[[5]]@data$depth)  
d35_5<-as.data.frame(d35_5)
names(d35_5)<-c("o2", "depth")

par(mfrow=c(2,3))
plot(u35_1, xlim=c(0, 250), ylim=c(75,0), col="red");par(new=TRUE)
plot(d35_1,xlim=c(0, 250), ylim=c(75, 0), col ="blue", main="station 3.5")

plot(u35_2, xlim=c(0, 250), ylim=c(75,0), col="red");par(new=TRUE)
plot(d35_2,xlim=c(0, 250), ylim=c(75, 0), col ="blue", main="red = upcast")

plot(u35_3,xlim=c(0, 250), ylim=c(75,0), col="red");par(new=TRUE)
plot(d35_3,xlim=c(0, 250), ylim=c(75, 0), col ="blue", main="blue = downcast")

plot(u35_4,xlim=c(0, 250), ylim=c(75,0), col="red");par(new=TRUE)
plot(d35_4,xlim=c(0, 250), ylim=c(75, 0), col ="blue")

plot(u35_5,xlim=c(0, 250), ylim=c(75,0), col="red");par(new=TRUE)
plot(d35_5,xlim=c(0, 250), ylim=c(75, 0), col ="blue")
