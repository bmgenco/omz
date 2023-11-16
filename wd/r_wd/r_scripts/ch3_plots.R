#### setup ####

rm(list=ls())
# wd<-"/home/brandon/callistowd/omz/wd/r_wd"
wd<-"/home/brandon/vestawd/omz/wd/r_wd"


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


t<-read.argo(files[1])
u<-read.argo(files[2])
v<-read.argo(files[3])
w<-read.argo(files[4])
x<-read.argo(files[5])
y<-read.argo(files[6])
z<-read.argo(files[7])
a<-read.argo(files[8])

olaf.sec<-as.section(list(t,u,v,w,x,y,z,a))

## ploting sections

plot(olaf.sec, which="BBP700", ztype="image", xtype="time", ylim=c(80, 0))


### ploting   old plots: ####
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

saveRDS(plot_list, "oflaf_profiles_v2.R")

setwd(wd)
setwd(fig)
dev.off()
pdf("olaf_profiles_v2.pdf", height =9, width=8)
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

# not cdom vut attenuatio


### testing density ####


plotProfile(w, xtype="oxygenAdjusted",ylim=c(27,24), col="blue", ytype= "sigmaTheta")
mtext("Cast 1", side=3, line=3)
pa<-recordPlot()
plotProfile(c6, xtype="oxygen2",ylim=c(27,24), col="blue", ytype= "sigmaTheta")
mtext("Cast 6", side=3, line=3)
pb<-recordPlot()
plotProfile(c12, xtype="oxygen2",ylim=c(27,24), col="blue", ytype= "sigmaTheta")
mtext("Cast 12", side=3, line=3)
pc<-recordPlot()

