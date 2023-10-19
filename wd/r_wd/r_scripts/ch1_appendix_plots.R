# ploting is easier in r code, run r workbook tio get dat aobjects



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

pdf("station_35_density_TS_v2.pdf", height=12, width=16)
plot_grid(pa, pb, pc, pd, pe, pf, pg, ph, pi, pj, pk, pl, pm, pn, po, nrow=3, ncol=5,
          labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o'), align = "h", scale=.82, greed=T, axis="r")
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
