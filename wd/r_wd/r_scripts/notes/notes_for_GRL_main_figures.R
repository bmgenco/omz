# notes and scraps from GRL Main figures
 plotting works
data(coastlineWorld)
## make a colormap that looks like the website
col <- colorRampPalette(c("purple", "#00007F", "blue",
                          "#007FFF", "cyan","#7FFF7F",
                        "yellow", "#FF7F00", "red", "#7F0000"))

#so many color ramps....
# col<- colorRampPalette(c("#ffffff", "#f7fcf5","#ffffe5","#f7fcb9","#d9f0a3",
# "#addd8e","#78c679","#41ab5d","#238443","#006837","#004529", "#41b6c4", "#1d91c0", "#225ea8","#253494","#081d58"))



## green
col<- colorRampPalette(c(
  "#DBF7D8", "#D4F1CF", "#CDEDC8", "#C7E8C0", "#C0E4B9", "#B8DFB1", "#B2DBAA", "#ABD7A3", "#A4D39C", "#9DCF95", "#95CA8E", 
  "#8EC787", "#87C381", "#7FBF7A", "#78BC74", "#6FB76C", "#67B466", "#5EB060", "#54AD5A", "#4AAA55", "#3EA64F", "#2FA34B", 
  "#1D9F49", "#0E9B49", "#059649", "#039149", "#048C49", "#068849", "#098348", "#0D7E47", "#117946", "#147445", "#167044", 
  "#186B42", "#1A6741", "#1B613F", "#1C5D3D", "#1D583B", "#1E5438", "#1E4F36", "#1E4A33", "#1E4631", "#1E412E", "#1D3D2B", 
  "#1C3928", "#1B3424", "#1A2F21", "#192B1E", "#17271B", "#152317"))


### green and purple

col<- colorRampPalette(c("#ffffff",
                         "#DBF7D8", "#D4F1CF", "#CDEDC8", "#C7E8C0", "#C0E4B9", "#B8DFB1", "#B2DBAA", "#ABD7A3", "#A4D39C", "#9DCF95", "#95CA8E", 
                         "#8EC787", "#87C381", "#7FBF7A", "#78BC74", "#6FB76C", "#67B466", "#5EB060", "#54AD5A", "#4AAA55", "#3EA64F", "#2FA34B", 
                         "#1D9F49", "#0E9B49", "#059649", "#039149", "#048C49", "#068849", "#098348", "#0D7E47", "#117946", "#147445", "#167044", 
                         "#186B42", "#1A6741", "#1B613F", "#1C5D3D", "#1D583B", "#1E5438", "#1E4F36", "#1E4A33", "#1E4631", "#1E412E", "#1D3D2B", 
                         "#1C3928", "#1B3424", "#1A2F21", "#192B1E", "#17271B", "#152317", "#3E4A89", "#404588", 
                         "#423E85","#443882","#46327E", "#472C7A", "#482576","#481E6F", "#481769", "#471063", "#46085C","#440154"))  

# deltte every 4th green and prupkle
col<- colorRampPalette(c("#ffffff","#DBF7D8", "#D4F1CF", "#C7E8C0", "#C0E4B9", "#B8DFB1",  "#ABD7A3", "#A4D39C", "#9DCF95", "#95CA8E", 
                         "#8EC787", "#87C381", "#7FBF7A", "#78BC74", "#6FB76C", "#67B466", "#5EB060", "#54AD5A", "#4AAA55", "#3EA64F", "#2FA34B", 
                         "#1D9F49", "#0E9B49", "#059649", "#039149", "#048C49", "#068849", "#098348", "#0D7E47", "#117946", "#147445", "#167044", 
                         "#186B42", "#1A6741", "#1B613F", "#1C5D3D", "#1D583B", "#1E5438", "#1E4F36", "#1E4A33", "#1E4631", "#1E412E", "#1D3D2B", 
                         "#1C3928", "#1B3424", "#1A2F21", "#192B1E", "#17271B", "#152317", "#3E4A89", "#404588", 
                         "#423E85","#443882","#46327E", "#472C7A", "#482576","#481E6F", "#481769", "#471063", "#46085C","#440154"))  

col<-colorRampPalette(c("#ffffff", "#DBF7D8", "#C9EAC3" ,"#B8DFB0" ,"#A6D49E", "#93CA8C" ,"#80BF7A", "#6BB66A" ,"#54AD5A", "#35A44C", "#0B9A49", 
                        "#038D49", "#0B8148", "#147445","#196841", "#1C5D3D", "#1E5037", "#1E4530", "#1C3928", "#1A2E20", "#152317","#3E4A89", "#404588", 
                        "#423E85","#443882","#46327E", "#472C7A", "#482576","#481E6F", "#481769", "#471063", "#46085C","#440154"))


col <- colorRampPalette(c("#ffffff", "#f7fcf5","#ffffe5","#f7fcb9","#d9f0a3",
                          "#addd8e","#78c679","#41ab5d","#238443","#006837","#004529",
                          "#88419d", "#810f7c","#4d004b"))

# 
# col <- colorRampPalette(c("#f7fcf5","#ffffe5","#f7fcb9","#d9f0a3",
# "#addd8e","#78c679","#41ab5d","#238443","#006837","#004529"))


# col <- colorRampPalette(c("#f7fcf5", "#e5f5e0", "#c7e9c0", 
# "#a1d99b", "#74c476", "#41ab5d",
# "#238b45", "#006d2c", "#00441b"))


# col <- colorRampPalette(c( "#c7e9c0", 
# "#a1d99b", "#74c476", "#41ab5d",
# "#238b45", "#006d2c", "#00441b"))



imagep(lon, lat, log10(test), col=col, zlim=log10(c(0.01, 20)), missingColor=1, zlab='log10(Chl)')


polygon(coastlineWorld[['longitude']], coastlineWorld[['latitude']], col='grey')


##### testing 



f.tmp1<-function(x){
  y<-nc_open(x)
  chl<-ncvar_get(y, 'chlor_a')
  lon<-ncvar_get(y, 'lon')
  lat<-ncvar_get(y, 'lat')
  rownames(chl)<-lon
  colnames(chl)<-lat
  nc_close(y)
  rm(y)
  x<-chl
  return(x)
}
f.tmp2<-function(x){
  y<-nc_open(x)
  # start<-ncatt_get(y, attname = 'time_coverage_start', varid = 0, verbose=F)
  # start<-start$value
  # start<-ymd_hms(start)
  # start<-with_tz(start, tz="US/Mountain")
  
  end<-ncatt_get(y, attname = 'time_coverage_end', varid = 0, verbose=F)
  end<-end$value
  end<-ymd_hms(end)
  end<-with_tz(end, tz="US/Mountain")
  # time<-c(start, end, (end-start))
  day<-day(end)
  nc_close(y)
  # rm(y, start, end)
  rm(y, end)
  
  return(as.character(day))
  
  
  lon<-ncvar_get(y, 'lon')
  lat<-ncvar_get(y, 'lat')
  rownames(chl)<-lon
  colnames(chl)<-lat
  nc_close(y)
  rm(y)
  x<-chl
  return(x)
}
# f.tmp3<-function(a,t){
#   z<-(a+t)/2
#   return(z)}

#aqua
directory<-"/home/brandon/vestawd/omz/ocean_color_bud/aqua"
setwd(directory)
files<-list.files(directory)

a<-lapply(files, f.tmp1)
dates<-lapply(files, f.tmp2)
names(a)<-as.character(dates)
rm(files, directory, dates)

#terra
directory<-"/home/brandon/vestawd/omz/ocean_color_bud/terra"
setwd(directory)
files<-list.files(directory)

t<-lapply(files, f.tmp1)
dates<-lapply(files, f.tmp2)
names(t)<-as.character(dates)
rm(files, directory, dates)

lat<-rownames(a[[1]])
lon<-colnames(a[[1]])















#2021-06-07 rmeots ensing image notes

combined<-vector(mode="list", length=length(a))
for(i in 1:length(combined)){
combined[[i]]<-(a[[i]]+t[[i]])/2
}

for(i in 1:length(combined)){
  combined[[i]]<-ifelse(is.na(combined[[i]]),a[[i]],combined[[i]])
}

for(i in 1:length(combined)){
  combined[[i]]<-ifelse(is.na(combined[[i]]),t[[i]],combined[[i]])
}




t1<- matrix(NA, nrow=3, ncol = 3)
t1[1,1]<-1
t1[3,3]<-1


t2<- matrix(NA, nrow=3, ncol = 3)
t2[1,1]<-2
t2[3,1]<-2

test<-(t1+t2)/2

test[,]<
  
apply(test, fun=if(is.na(test)) test<-t1)

test[is.na(test)]==t1

## this is what I want

combined<-ifelse(is.na(combined),a,combined)
combined<-ifelse(is.na(combined),t,combined)



directory<-"/home/brandon/vestawd/omz/ocean_color_bud/aqua"
setwd(directory)

files<-list.files(directory)
y<-nc_open(files[1])


f.tmp1<-function(x){
  y<-nc_open(x)
  chl<-ncvar_get(y, 'chlor_a')
  lon<-ncvar_get(y, 'lon')
  lat<-ncvar_get(y, 'lat')
  rownames(chl)<-lon
  colnames(chl)<-lat
  nc_close(y)
  rm(y)
  x<-chl
  return(x)
}
f.tmp2<-function(x){
  y<-nc_open(x)
  # start<-ncatt_get(y, attname = 'time_coverage_start', varid = 0, verbose=F)
  # start<-start$value
  # start<-ymd_hms(start)
  # start<-with_tz(start, tz="US/Mountain")
  
  end<-ncatt_get(y, attname = 'time_coverage_end', varid = 0, verbose=F)
  end<-end$value
  end<-ymd_hms(end)
  end<-with_tz(end, tz="US/Mountain")
  # time<-c(start, end, (end-start))
  day<-day(end)
  nc_close(y)
  # rm(y, start, end)
  rm(y, end)
  
  return(as.character(day))
  
  
  lon<-ncvar_get(y, 'lon')
  lat<-ncvar_get(y, 'lat')
  rownames(chl)<-lon
  colnames(chl)<-lat
  nc_close(y)
  rm(y)
  x<-chl
  return(x)
}
a<-lapply(files, f.tmp1)
dates<-lapply(files, f.tmp2)
names(a)<-as.character(dates)




y<-nc_open(files[1])
t<-ncvar_get(y, 'time_coverage_start:')



directory<-"/home/brandon/vestawd/omz/ocean_color_bud/aqua"

read_rs<-function(directory){
  files<-list.files(directory)
  x<-vector(mode="list", length=length(files))
  
  for(i in 1:length(files)){
    x[[i]]<-nc_open(files[i])
    chl<-ncvar_get(x[[i]], 'chlor_a')
    lon<-ncvar_get(x[[i]], 'lon')
    lat<-ncvar_get(x[[i]], 'lat')
    rownames(chl)<-lon
    colnames(chl)<-lat
    x[[i]]<-(chl) 
  }
  return(x)}  


test<-read_rs(directory)  

for(i in 1:length(files)){
  a[[i]]<-read.delim(file.path(directory, files[i]), sep =",", header = F)  
}
return(data)}

#1
chl<-comp
if comp=NA comp<-a e

med<-
  comp<-matrix(NA, nrow=612, ncol=600)

new<-a if a=!NA and t=NA else 
  
  
  
  # creat 50 km bufferr
  
  st_buffer(station3.5, dist = 5000)
#https://geocompr.robinlovelace.net/reproj-geo-data.html

















setwd(robj)
f<-readRDS("2018_flourometer_raw.R")
t<-readRDS("2018_temperature_raw.R")

x<-f.time_join(f,t)
rm(f,t)

x<-subset(x, select=c(date, time, fluorometer, celsius))

# testing
y<-sample_n(x, 300)


# works:

f1<-ggplot(x, aes(x=time)) +
  geom_line(aes(y=wind), color=wcolor) +
  geom_line(aes(y=chlorophyll * scalef), color=fcolor) +
  geom_line(aes(y=celsius), color=tcolor) +
  xlab("Date") +
  ggtitle("2018 Flowthrough: Surface Temperature and chlorophyll") +
  scale_x_datetime(date_breaks=('2 days'))+
  #scale_y_continuous(name ="Temperature (Celsius)", sec.axis = sec_axis(~./scalef,name="Chlorophyll-a [μg/L]")) + 
  scale_y_continuous(name ="Chlorophyll-a [μg/L]", sec.axis = sec_axis(~./scalef,name="°C")) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.ticks.y.left = element_line(color = fcolor),
        axis.text.y.left = element_text(color = fcolor), 
        axis.title.y.left = element_text(color = fcolor)
  ) +
  
  theme(axis.ticks.y.right = element_line(color = tcolor),
        axis.text.y.right = element_text(color = tcolor), 
        axis.title.y.right = element_text(color = tcolor)
  )


# Attemp 1


f1<-ggplot(x, aes(x=time)) +
  geom_line(aes(y=wind), color=wcolor) +
  geom_line(aes(y=chlorophyll * scalef), color=fcolor) +
  geom_line(aes(y=celsius), color=tcolor) +
  xlab("Date") +
  ggtitle("2018 Flowthrough: Surface Temperature and chlorophyll") +
  scale_x_datetime(date_breaks=('2 days'))+
  #scale_y_continuous(name ="Temperature (Celsius)", sec.axis = sec_axis(~./scalef,name="Chlorophyll-a [μg/L]")) + 
  scale_y_continuous(name ="Chlorophyll-a [μg/L]", sec.axis = sec_axis(~./scale.chla,name="°C")) + 
  scale_y_continuous(sec.axis =sec_axis(~./x$wind, name="Wind speed (m/s)"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.ticks.y.left = element_line(color = fcolor),
        axis.text.y.left = element_text(color = fcolor), 
        axis.title.y.left = element_text(color = fcolor)
  ) +
  
  theme(axis.ticks.y.right = element_line(color = tcolor),
        axis.text.y.right = element_text(color = tcolor), 
        axis.title.y.right = element_text(color = tcolor)
  ) +

  theme(axis.ticks.y.right = element_line(color = wcolor),
        axis.text.y.right = element_text(color = wcolor), 
        axis.title.y.right = element_text(color = wcolor))

  



stations_positions_list
dim(stations_positions_list$flowthrough_stations)
unique(stations_positions_list$flowthrough_stations$start_time)
