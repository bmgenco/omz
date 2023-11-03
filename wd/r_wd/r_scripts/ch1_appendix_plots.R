# ploting is easier in r code, run r workbook tio get dat aobjects

#### setup ####

rm(list=ls())
# wd<-"/home/brandon/callistowd/omz/wd/r_wd"
wd<-"/home/brandon/vestawd/omz/wd/r_wd"


# relative directories
robj<-"r_objects"
fig<-"../../figures"
gis_data<-"../../data/gis_data"
output<-"../../output"

setwd(wd)

# packages 
f.ipak <- function(pkg){
  # loads packages, quietly, given by a vector of package names e.g., pkg<-c("ggplot", "tidyverse")
  # will install  packages listed , and their dependencies, if needed.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, quiet=T, verbose = F)
  sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts=F)
}
packages<-c("sp", "rgdal",  "rgeos", "raster", "readr", "tidyverse", "lubridate",  "ggthemes",  "sf", "cmocean", "ncdf4", "RNetCDF",  "plot3D", "tidync", "devtools", "stars", "ncmeta", "maps", "oce", "data.table", "fasterize", "RStoolbox", "scales", "purrr", "HURDAT", "ggpattern", "zoo", "xts", "dsa",  "tsibble", "feasts", "tsbox", "rsoi", "cowplot")

f.ipak(packages)

setwd(wd)
d.sec<-readRDS("r_objects/OC1806A_ctd_sections.R")
d35<-d.sec$st3.5
d4<-d.sec$st4

d3<-d.sec$st3


### multiple casts station 3 ####
c1<-d3@data$station[[1]]

c2<-d3@data$station[[2]]
c3<-d3@data$station[[3]]
c4<-d3@data$station[[4]]
c5<-d3@data$station[[5]]
c6<-d3@data$station[[6]]

c7<-d3@data$station[[7]]

c8<-d3@data$station[[8]]
c9<-d3@data$station[[9]]
c10<-d3@data$station[[10]]
c11<-d3@data$station[[11]]
c12<-d3@data$station[[12]]


plotProfile(c1, xtype="oxygen2",ylim=c(160,10), col="blue", ytype= "depth")
mtext("Cast 1", side=3, line=3)
pa<-recordPlot()
plotProfile(c6, xtype="oxygen2",ylim=c(160,10), col="blue", ytype= "depth")
mtext("Cast 6", side=3, line=3)
pb<-recordPlot()
plotProfile(c12, xtype="oxygen2",ylim=c(160,10), col="blue", ytype= "depth")
mtext("Cast 12", side=3, line=3)
pc<-recordPlot()


plotProfile(c1, xtype="beamAttenuation",ylim=c(160,10), col ="black", ytype= "depth")
pd<-recordPlot()
plotProfile(c6, xtype="beamAttenuation",ylim=c(160,10), col ="black", ytype= "depth")
pe<-recordPlot()
plotProfile(c12, xtype="beamAttenuation",ylim=c(160,10), col ="black", ytype= "depth")
pf<-recordPlot()


plotProfile(c1, xtype="fluorescence",ylim=c(160,10), col ="green", ytype="depth")
pg<-recordPlot()
plotProfile(c6, xtype="fluorescence",ylim=c(160,10), col ="green", ytype="depth")
ph<-recordPlot()
plotProfile(c12, xtype="fluorescence",ylim=c(160,10), col ="green", ytype="depth")
pi<-recordPlot()



plotProfile(c1, xtype="density",  ytype= "depth", col="brown", ylim=c(160,10))
pj<-recordPlot()
plotProfile(c6, xtype="density", ytype= "depth", col="brown", ylim=c(160,10))
pk<-recordPlot()
plotProfile(c12, xtype="density", ytype= "depth", col="brown", ylim=c(160,10))
pl<-recordPlot()

dev.off()
setwd(wd)
setwd(fig)

pdf("station_3_depth_TS.pdf", height=16, width=12)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, pj, pk, pl, nrow=4, ncol=3,
          labels ="auto", align = "h", scale=.82, greed=T, axis="r")
dev.off()

#### staion 3 density ts ####



plotProfile(c1, xtype="oxygen2",ylim=c(27,24), col="blue", ytype= "sigmaTheta")
mtext("Cast 1", side=3, line=3)
pa<-recordPlot()
plotProfile(c6, xtype="oxygen2",ylim=c(27,24), col="blue", ytype= "sigmaTheta")
mtext("Cast 6", side=3, line=3)
pb<-recordPlot()
plotProfile(c12, xtype="oxygen2",ylim=c(27,24), col="blue", ytype= "sigmaTheta")
mtext("Cast 12", side=3, line=3)
pc<-recordPlot()


plotProfile(c1, xtype="beamAttenuation",ylim=c(27,24), col ="black", ytype= "sigmaTheta")
pd<-recordPlot()
plotProfile(c6, xtype="beamAttenuation",ylim=c(27,24), col ="black", ytype= "sigmaTheta")
pe<-recordPlot()
plotProfile(c12, xtype="beamAttenuation",ylim=c(27,24), col ="black", ytype= "sigmaTheta")
pf<-recordPlot()


plotProfile(c1, xtype="fluorescence",ylim=c(27,24), col ="green", ytype="sigmaTheta")
pg<-recordPlot()
plotProfile(c6, xtype="fluorescence",ylim=c(27,24), col ="green", ytype="sigmaTheta")
ph<-recordPlot()
plotProfile(c12, xtype="fluorescence",ylim=c(27,24), col ="green", ytype="sigmaTheta")
pi<-recordPlot()

dev.off()
setwd(wd)
setwd(fig)



pdf("station_3_density_TS.pdf", height=12, width=12)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, nrow=3, ncol=3,
          labels ="auto", align = "h", scale=.82, greed=T, axis="r")
dev.off()



### station 3 ocenagprahic plots

plotProfile(c6, xtype="density+N2", ylim=c(180,10),  ytype="depth")
pa<-recordPlot()
plotProfile(c6, xtype="salinity+temperature", ylim=c(180,10))
pb<-recordPlot()
plotProfile(c6, xtype="spice", ylim=c(180,10),  ytype="depth")
pc<-recordPlot()
plotProfile(c6, xtype="Rrho", ylim=c(180,10),  ytype="depth")
pd<-recordPlot()

dev.off()
setwd(wd)
setwd(fig)

pdf("station_3_ocenao.pdf", height=12, width=12)
plot_grid(pa, pb, pc, pd, nrow=2, ncol=2,
          labels ="auto", align = "h", scale=.82, greed=T, axis="r")
dev.off()





### multiple casts station 3.5 ####
c1<-d35@data$station[[1]]

c2<-d35@data$station[[2]]
c3<-d35@data$station[[3]]
c4<-d35@data$station[[4]]
c5<-d35@data$station[[5]]

plotProfile(c1, xtype="oxygen2",ylim=c(80,10), col="blue", ytype= "depth")
mtext("Cast 1", side=3, line=3)
pa<-recordPlot()
plotProfile(c2, xtype="oxygen2",ylim=c(80,10), col="blue", ytype= "depth")
mtext("Cast 2", side=3, line=3)
pb<-recordPlot()
plotProfile(c3, xtype="oxygen2",ylim=c(80,10), col="blue", ytype= "depth")
mtext("Cast 3", side=3, line=3)
pc<-recordPlot()
plotProfile(c4, xtype="oxygen2",ylim=c(80,10), col="blue", ytype= "depth")
mtext("Cast 4", side=3, line=3)
pd<-recordPlot()
plotProfile(c5, xtype="oxygen2",ylim=c(80,10), col="blue", ytype= "depth")
mtext("Cast 5", side=3, line=3)
pe<-recordPlot()
plotProfile(c1, xtype="beamAttenuation",ylim=c(80,10), col ="black", ytype= "depth")
pf<-recordPlot()
plotProfile(c2, xtype="beamAttenuation",ylim=c(80,10), col ="black", ytype= "depth")
pg<-recordPlot()
plotProfile(c3, xtype="beamAttenuation",ylim=c(80,10), col ="black", ytype= "depth")
ph<-recordPlot()
plotProfile(c4, xtype="beamAttenuation",ylim=c(80,10), col ="black", ytype= "depth")
pi<-recordPlot()
plotProfile(c5, xtype="beamAttenuation",ylim=c(80,10), col ="black", ytype= "depth")
pj<-recordPlot()
plotProfile(c1, xtype="fluorescence",ylim=c(80,10), col ="green", ytype="depth")
pk<-recordPlot()
plotProfile(c2, xtype="fluorescence",ylim=c(80,10), col ="green", ytype="depth")
pl<-recordPlot()
plotProfile(c3, xtype="fluorescence",ylim=c(80,10), col ="green", ytype="depth")
pm<-recordPlot()
plotProfile(c4, xtype="fluorescence",ylim=c(80,10), col ="green", ytype="depth")
pn<-recordPlot()
plotProfile(c5, xtype="fluorescence",ylim=c(80,10), col ="green", ytype="depth")
po<-recordPlot()

dev.off()

pdf("station_35_depth_TS.pdf", height=12, width=16)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, pj, pk, pl, pm, pn, po, nrow=3, ncol=5,
          labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o'), align = "h", scale=.82, greed=T, axis="r")
dev.off()


#### ploting DCM  TS station 3.5 ####


dev.off()



plotProfile(c1, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26.5,25))
mtext("Cast 1", side=3, line=3)
pa<-recordPlot()
plotProfile(c2, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26.5,25))
mtext("Cast 2", side=3, line=3)
pb<-recordPlot()
plotProfile(c3, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26.5,25))
mtext("Cast 3", side=3, line=3)
pc<-recordPlot()
plotProfile(c4, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26.5,25))
mtext("Cast 4", side=3, line=3)
pd<-recordPlot()
plotProfile(c5, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26.5,25))
mtext("Cast 5", side=3, line=3)
pe<-recordPlot()
plotProfile(c1, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26.5,25))
pf<-recordPlot()
plotProfile(c2, xtyp="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26.5,25))
pg<-recordPlot()
plotProfile(c3, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26.5,25))
ph<-recordPlot()
plotProfile(c4, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26.5,25))
pi<-recordPlot()
plotProfile(c5, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26.5,25))
pj<-recordPlot()
plotProfile(c1, xtype="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26.5,25))
pk<-recordPlot()
plotProfile(c2, xtyp="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26.5,25))
pl<-recordPlot()
plotProfile(c3, xtype="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26.5,25))
pm<-recordPlot()
plotProfile(c4, xtype="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26.5,25))
pn<-recordPlot()
plotProfile(c5, xtype="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26.5,25))
po<-recordPlot()

dev.off()

setwd(wd)
setwd(fig)

pdf("station_35_density_TS.pdf", height=12, width=16)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, pj, pk, pl, pm, pn, po, nrow=3, ncol=5,
          labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o'), align = "h", scale=.82, greed=T, axis="r")
dev.off()



#### Plotting depth all casta station 4 ####

plotProfile(c1, xtype="oxygen2",  ytype= "depth", col="blue", ylim=c(80,10))
mtext("Cast 1", side=3, line=3)
pa<-recordPlot()
plotProfile(c2, xtype="oxygen2",  ytype= "depth", col="blue", ylim=c(80,10))
mtext("Cast 2", side=3, line=3)
pb<-recordPlot()
plotProfile(c3, xtype="oxygen2",  ytype= "depth", col="blue", ylim=c(80,10))
mtext("Cast 3", side=3, line=3)
pc<-recordPlot()
plotProfile(c4, xtype="oxygen2",  ytype= "depth", col="blue", ylim=c(80,10))
mtext("Cast 4", side=3, line=3)
pd<-recordPlot()
plotProfile(c5, xtype="oxygen2",  ytype= "depth", col="blue", ylim=c(80,10))
mtext("Cast 5", side=3, line=3)
pe<-recordPlot()
plotProfile(c6, xtype="oxygen2",  ytype= "depth", col="blue", ylim=c(80,10))
mtext("Cast 6", side=3, line=3)
pf<-recordPlot()

plotProfile(c1, xtype="beamAttenuation",  ytype= "depth", col="black", ylim=c(80,10))
pg<-recordPlot()
plotProfile(c2, xtyp="beamAttenuation",  ytype= "depth", col="black", ylim=c(80,10))
ph<-recordPlot()
plotProfile(c3, xtype="beamAttenuation",  ytype= "depth", col="black", ylim=c(80,10))
pi<-recordPlot()
plotProfile(c4, xtype="beamAttenuation",  ytype= "depth", col="black", ylim=c(80,10))
pj<-recordPlot()
plotProfile(c5, xtype="beamAttenuation",  ytype= "depth", col="black", ylim=c(80,10))
pk<-recordPlot()
plotProfile(c6, xtype="beamAttenuation",  ytype= "depth", col="black", ylim=c(80,10))
pl<-recordPlot()


plotProfile(c1, xtype="fluorescence",  ytype= "depth", col="green", ylim=c(80,10))
pm<-recordPlot()
plotProfile(c2, xtype="fluorescence", ytype= "depth", col="green", ylim=c(80,10))
pn<-recordPlot()
plotProfile(c3, xtype="fluorescence",ytype= "depth", col="green", ylim=c(80,10))
po<-recordPlot()
plotProfile(c4, xtype="fluorescence",  ytype= "depth", col="green", ylim=c(80,10))
pp<-recordPlot()
plotProfile(c5, xtype="fluorescence",  ytype= "depth", col="green", ylim=c(80,10))
pq<-recordPlot()
plotProfile(c6, xtype="fluorescence",  ytype= "depth", col="green", ylim=c(80,10))
pr<-recordPlot()


dev.off()

setwd(wd)
setwd(fig)

pdf("station_4_all_depth.pdf", height=12, width=19.2)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, pj, pk, pl, pm, pn, po, pp, pq, pr, nrow=3, ncol=6,
          labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r'), align = "h", scale=.82, greed=T, axis="r")
dev.off()



plotProfile(c1, xtype="density",  ytype= "depth", col="brown", ylim=c(80,10))
ps<-recordPlot()
plotProfile(c2, xtype="density", ytype= "depth", col="brown", ylim=c(80,10))
pt<-recordPlot()
plotProfile(c3, xtype="density", ytype= "depth", col="brown", ylim=c(80,10))
pu<-recordPlot()
plotProfile(c4, xtype="density",  ytype= "depth", col="brown", ylim=c(80,10))
pv<-recordPlot()
plotProfile(c5, xtype="density",  ytype= "depth", col="brown", ylim=c(80,10))
pw<-recordPlot()
plotProfile(c6, xtype="density",  ytype= "depth", col="brown", ylim=c(80,10))
px<-recordPlot()




dev.off()

setwd(wd)
setwd(fig)

pdf("station_4_all_depth_with_density.pdf", height=16, width=19.2)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, pj, pk, pl, pm, pn, po, pp, pq, pr, ps, pt, pu, pv, pw, px, nrow=4, ncol=6,
          labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x'), align = "h", scale=.82, greed=T, axis="r")
dev.off()



#### ploting DCM  TS station 4 ####




## need to adjust here

c1<-d4@data$station[[1]]
c2<-d4@data$station[[2]]
c3<-d4@data$station[[3]]
c4<-d4@data$station[[4]]
c5<-d4@data$station[[5]]
c6<-d4@data$station[[6]]

# plotting

plotProfile(c1, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26,24))
mtext("Cast 1", side=3, line=3)
pa<-recordPlot()
plotProfile(c2, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26,24))
mtext("Cast 2", side=3, line=3)
pb<-recordPlot()
plotProfile(c3, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26,24))
mtext("Cast 3", side=3, line=3)
pc<-recordPlot()
plotProfile(c4, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26,24))
mtext("Cast 4", side=3, line=3)
pd<-recordPlot()
plotProfile(c5, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26,24))
mtext("Cast 5", side=3, line=3)
pe<-recordPlot()
plotProfile(c6, xtype="oxygen2",  ytype="sigmaTheta", col="blue", ylim=c(26,24))
mtext("Cast 6", side=3, line=3)
pf<-recordPlot()

plotProfile(c1, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26,24))
pg<-recordPlot()
plotProfile(c2, xtyp="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26,24))
ph<-recordPlot()
plotProfile(c3, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26,24))
pi<-recordPlot()
plotProfile(c4, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26,24))
pj<-recordPlot()
plotProfile(c5, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26,24))
pk<-recordPlot()
plotProfile(c6, xtype="beamAttenuation",  ytype="sigmaTheta", col="black", ylim=c(26,24))
pl<-recordPlot()


plotProfile(c1, xtype="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26,24))
pm<-recordPlot()
plotProfile(c2, xtype="fluorescence", ytype="sigmaTheta", col="green", ylim=c(26,24))
pn<-recordPlot()
plotProfile(c3, xtype="fluorescence",ytype="sigmaTheta", col="green", ylim=c(26,24))
po<-recordPlot()
plotProfile(c4, xtype="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26,24))
pp<-recordPlot()
plotProfile(c5, xtype="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26,24))
pq<-recordPlot()
plotProfile(c6, xtype="fluorescence",  ytype="sigmaTheta", col="green", ylim=c(26,24))
pr<-recordPlot()




dev.off()

setwd(wd)
setwd(fig)

pdf("station_4_density_TS_.pdf", height=12, width=19.2)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, pj, pk, pl, pm, pn, po, pp, pq, pr, nrow=3, ncol=6,
          labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r'), align = "h", scale=.82, greed=T, axis="r")
dev.off()



### sig

### multiple casts station 4

#### plotting dpeth station 4 all ####


c1<-d4@data$station[[1]]
c2<-d4@data$station[[2]]
c3<-d4@data$station[[3]]
c4<-d4@data$station[[4]]
c5<-d4@data$station[[5]]
c6<-d4@data$station[[6]]


plotProfile(c1, xtype="oxygen2",  ytype="depth", col="blue", ylim=c(80,10))
mtext("Cast 1", side=3, line=3)
pa<-recordPlot()
plotProfile(c2, xtype="oxygen2",  ytype="depth", col="blue", ylim=c(80,10))
mtext("Cast 2", side=3, line=3)
pb<-recordPlot()
plotProfile(c3, xtype="oxygen2",  ytype="depth", col="blue", ylim=c(80,10))
mtext("Cast 3", side=3, line=3)
pc<-recordPlot()
plotProfile(c4, xtype="oxygen2",  ytype="depth", col="blue", ylim=c(80,10))
mtext("Cast 4", side=3, line=3)
pd<-recordPlot()
plotProfile(c5, xtype="oxygen2",  ytype="depth", col="blue", ylim=c(80,10))
mtext("Cast 5", side=3, line=3)
pe<-recordPlot()
plotProfile(c6, xtype="oxygen2",  ytype="depth", col="blue", ylim=c(80,10))
mtext("Cast 6", side=3, line=3)
pf<-recordPlot()

plotProfile(c1, xtype="beamAttenuation",  ytype="depth", col="black", ylim=c(80,10))
pg<-recordPlot()
plotProfile(c2, xtyp="beamAttenuation",  ytype="depth", col="black", ylim=c(80,10))
ph<-recordPlot()
plotProfile(c3, xtype="beamAttenuation",  ytype="depth", col="black", ylim=c(80,10))
pi<-recordPlot()
plotProfile(c4, xtype="beamAttenuation",  ytype="depth", col="black", ylim=c(80,10))
pj<-recordPlot()
plotProfile(c5, xtype="beamAttenuation",  ytype="depth", col="black", ylim=c(80,10))
pk<-recordPlot()
plotProfile(c6, xtype="beamAttenuation",  ytype="depth", col="black", ylim=c(80,10))
pl<-recordPlot()

plotProfile(c1, xtype="fluorescence",  ytype="depth", col="green", ylim=c(80,10))
pm<-recordPlot()
plotProfile(c2, xtype="fluorescence", ytype="depth", col="green", ylim=c(80,10))
pn<-recordPlot()
plotProfile(c3, xtype="fluorescence",ytype="depth", col="green", ylim=c(80,10))
po<-recordPlot()
plotProfile(c4, xtype="fluorescence",  ytype="depth", col="green", ylim=c(80,10))
pp<-recordPlot()
plotProfile(c5, xtype="fluorescence",  ytype="depth", col="green", ylim=c(80,10))
pq<-recordPlot()
plotProfile(c6, xtype="fluorescence",  ytype="depth", col="green", ylim=c(80,10))
pr<-recordPlot()


plotProfile(c1, xtype="density",  ytype="depth", col="brown", ylim=c(80,10))
ps<-recordPlot()
plotProfile(c2, xtype="density", ytype="depth", col="brown", ylim=c(80,10))
pt<-recordPlot()
plotProfile(c3, xtype="density",ytype="depth", col="brown", ylim=c(80,10))
pu<-recordPlot()
plotProfile(c4, xtype="density",  ytype="depth", col="brown", ylim=c(80,10))
pv<-recordPlot()
plotProfile(c5, xtype="density",  ytype="depth", col="brown", ylim=c(80,10))
pw<-recordPlot()
plotProfile(c6, xtype="density",  ytype="depth", col="brown", ylim=c(80,10))
px<-recordPlot()

dev.off()
setwd(wd)
setwd(fig)


test<-list(ps, pt, pu, pv, pw, px)

pltz<-list()
pltz[[1]]<-ps
pltz[[2]]<-pt
pltz[[3]]<-pu

pdf("station_4_depth_TS_v2.pdf", height=16, width=16)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, pj, pk, pl, pm, pn, po, pp, pq, pr, ps, pt, pu, pv, pw, px, nrow=4, ncol=6,
labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x'), align = "h", scale=.82, greed=T, axis="r")
dev.off()


plot_grid(pltz)

library(grid)
grid.arrange(grobs =pltz)


dev.off()
par(mfrow = c(2, 2))
plotProfile(c3, xtype="density+N2", ylim=c(80,10),  ytype="depth")
plotProfile(c3, xtype="salinity+temperature", ylim=c(80,10))
plotProfile(c3, xtype="spice", ylim=c(80,10),  ytype="depth")
plotProfile(c3, xtype="RrhoSF", ylim=c(80,10),  ytype="depth")

plot_1<-recordPlot()
par(mfrow = c(2, 2))
plotProfile(c4, xtype="density+N2", ylim=c(80,10),  ytype="depth")
plotProfile(c4, xtype="salinity+temperature", ylim=c(80,10))
plotProfile(c4, xtype="spice", ylim=c(80,10),  ytype="depth")
plotProfile(c4, xtype="RrhoSF", ylim=c(80,10),  ytype="depth")

plot_2<-recordPlot()


dev.off()
par(mfrow = c(1, 6))
plotProfile(c1, xtype="density",  ytype= "depth", col="brown", ylim=c(80,10))
mtext("Cast 1", side=4, line=.5)
plotProfile(c2, xtype="density", ytype= "depth", col="brown", ylim=c(80,10))
mtext("Cast 2", side=4, line=.5)
plotProfile(c3, xtype="density", ytype= "depth", col="brown", ylim=c(80,10))
mtext("Cast 3", side=4, line=.5)
plotProfile(c4, xtype="density",  ytype= "depth", col="brown", ylim=c(80,10))
mtext("Cast 4", side=4, line=.5)
plotProfile(c5, xtype="density",  ytype= "depth", col="brown", ylim=c(80,10))
mtext("Cast 5", side=4, line=.5)
plotProfile(c6, xtype="density",  ytype= "depth", col="brown", ylim=c(80,10))
mtext("Cast 6", side=4, line=.5)

plot_3<-recordPlot()


setwd(wd)
setwd(fig)

dev.off()
pdf("station_4_density_subplots_c_v2.pdf", height=7, width =14)
plot_3
dev.off()

