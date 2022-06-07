rm(list=ls())
options(warn=-1)
t.start<-Sys.time()

#### directories ####

wd<-"/home/brandon/vestawd/omz/wd/r_wd"
setwd(wd)  
robj<-"r_objects"
fig<-"../../figures"
data_poc_all<-"../../data/POC_Flux_data_sets/annual_vgpm_NPP_netcdf"


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
"stars", "ncmeta", "maps", "oce", "data.table", "fasterize", "RStoolbox", "scales", "purrr")

# "ncdf4", "RNetCDF",

f.ipak(packages)


### stations ####
setwd(wd)
setwd(robj)
# 
stations_positions_list<-readRDS("OC1806A_stations_positions_list.R")

st<-stations_positions_list[[3]][1,cbind(2,3)]
st<-sf::st_as_sf(st, coords =c("longitude", "latitude"),  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
flat<-sf::st_transform(st, "+proj=aeqd +lat_0=20.50006 +lon_0=-106.5006" )
circle_1<-sf::st_buffer(flat, dist=100000)%>%sf::st_transform(., "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
rm(flat, st)

st<-stations_positions_list[[3]][2,cbind(2,3)]
st<-sf::st_as_sf(st, coords =c("longitude", "latitude"),  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
flat<-sf::st_transform(st, "+proj=aeqd +lat_0=16.49850 +lon_0=-107.2007" )
circle_2<-sf::st_buffer(flat, dist=100000)%>%sf::st_transform(., "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
rm(flat, st)

st<-stations_positions_list[[3]][3,cbind(2,3)]
st<-sf::st_as_sf(st, coords =c("longitude", "latitude"),  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
flat<-sf::st_transform(st, "+proj=aeqd +lat_0=15.99901 +lon_0=-108.5020" )
circle_3<-sf::st_buffer(flat, dist=100000)%>%sf::st_transform(., "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

st<-stations_positions_list[[3]][4,cbind(2,3)]
st<-sf::st_as_sf(st, coords =c("longitude", "latitude"),  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
flat<-sf::st_transform(st, "+proj=aeqd +lat_0=18.5005 +lon_0=-108.502" )
circle_3.5<-sf::st_buffer(flat, dist=100000)%>%sf::st_transform(., "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

st<-stations_positions_list[[3]][5,cbind(2,3)]
st<-sf::st_as_sf(st, coords =c("longitude", "latitude"),  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
flat<-sf::st_transform(st, "+proj=aeqd +lat_0=21.50061 +lon_0=-109.4999" )
circle_4<-sf::st_buffer(flat, dist=100000)%>%sf::st_transform(., "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

st<-stations_positions_list[[3]][6,cbind(2,3)]
st<-sf::st_as_sf(st, coords =c("longitude", "latitude"),  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
flat<-sf::st_transform(st, "+proj=aeqd +lat_0=24.69887 +lon_024.69887" )
circle_5<-sf::st_buffer(flat, dist=100000)%>%sf::st_transform(., "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

rm(st, flat, stations_positions_list)


#make a list for next steps

station_circle_list<-list(circle_1, circle_2, circle_3, circle_3.5, circle_4, circle_5 )
rm(circle_1, circle_2, circle_3, circle_3.5, circle_4, circle_5)

f.bbox<-function(x){
  z<-st_bbox(x)
  z<-data.frame(lon= c(z[[1]], z[[3]]), lat = c(z[[2]], z[[4]])) #%>% st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  return(z)
}
station_box_list<-lapply(station_circle_list, f.bbox)
names(station_box_list)<-c("st1", "st2", "st3", "st3.5", "st4", "st5")
rm(station_circle_list)


### functions ####

f.ef<-function(NPP){
  # units =mg C m−2 d−1
  EF =(0.284 * NPP + 9.75)
  EF= round(EF, 3)
  return(EF)
}

f.vgpm_cal<-function(vgpm){
  # VGPM-CAL = 10^(LOG10(VGPM)) − 0.1924. 
  NPP= 10^(log10(vgpm)) - 0.1924
  return(NPP)
}

f.grid_square_area2<-function(x) {
  x.area<-vector(mode="numeric", length = nrow(x)) 
  
  for(i in 1:length(x.area)){
    bbox<-st_bbox(x[[2]][[i]])
    box<-vector(mode="numeric", length=4)
    names(box)<-c("xmin","ymin","ymax", "xmax")
    box$ymin<-round(bbox$ymin[[1]]-(5/120),3)
    box$ymax<-round(bbox$ymax[[1]]+(5/120),3)
    box$xmin<-round(bbox$xmin[[1]] +-(5/120),3)
    box$xmax<-round(bbox$xmax[[1]] +(5/120),3)
    x.area[i]<- rgeos::bbox2SP(n = box$ymax[[1]], s = box$ymin[[1]], w = box$xmin[[1]], e = box$xmax[[1]],
                               proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% st_as_sf(.,) %>% st_area(.,) %>% round(., -1)
  }
  sum_grid_area<-sum(x.area)
  x$value<-x.area
  y<-vector(mode = "list", length=2)
  names(y)<-c("grid_area", "sum_grid_area")
  y$grid_area<-x
  y$sum_grid_area<-sum_grid_area
  
  return(y)}

f.net_ef<-function(x){
  # grid_area<-f.grid_square_area2(x.sf)
  # x$value<-f.vgpm_cal(x$value)
  x$value<-f.vgpm_cal(x$npp)
  x$value<-f.ef(x$value)
  x<-select(x, value, geometry)
  y<-f.grid_square_area2(x)
  EF<-(x*y$grid_area)
  sum_grid_area<<-y$sum_grid_area
  return(sum(EF$value))
}

f.ef_data_nc<-function(z){
  
  vgpm_integrated<-as.data.frame(matrix(ncol = 13, nrow =length(files)))
  names(vgpm_integrated)<-c("date",  "vgpm_median (mgC/m^2/day)", "vgpm_mean (mgC/m^2/day)", "cal_vgpm_median (mgC/m^2/day)", "cal_vgpm_mean (mgC/m^2/day)", "export_flux_median (mgC/m^2/day)",  "export_flux_mean (mgC/m^2/day)", "net_vgpm_npp (mgC/day)", "window_area (m^2)", "sum_grid_area (m^2)", "vgpm_variance (mgC/m^2/day)", "total_export_flux (mgC/day)", "mean_export_flux (mgC/m^2/day) for_window_area (total/sum_grid_area)" )
 
  for(i in 1:length(files)){
    file<-files[i]
    t<-tidync(file)%>% hyper_filter(lat = lat >= min(z$lat) & lat <= max(z$lat), lon = lon>= min(z$lon) & lon<= max(z$lon)) %>% hyper_tibble()
    t$npp[t$npp == -9999] <- NA
    x<-t$npp
    x.sf<-t%>%st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    
    bbox<-st_bbox(x.sf)
    box<-vector(mode="numeric", length=4)
    names(box)<-c("xmin","ymin","ymax", "xmax")
    box$ymin<-round(bbox$ymin[[1]]-(5/120),3)
    box$ymax<-round(bbox$ymax[[1]]+(5/120),3)
    box$xmin<-round(bbox$xmin[[1]] +-(5/120),3)
    box$xmax<-round(bbox$xmax[[1]] +(5/120),3)
    
    area <- rgeos::bbox2SP(n = box$ymax[[1]], s = box$ymin[[1]], w = box$xmin[[1]], e = box$xmax[[1]],
                           proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% st_as_sf(.,) %>% st_area(.,)
    rm(box, bbox)
    
    t<-str_split(file, '\\.')
    t<-t[[1]][2]
    t<-str_split(t, "")
    y<-c(t[[1]][1:4])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)   
    yd<-c(t[[1]][5:7])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)
    origin<-toString(c(y,"-01-01")) %>% gsub(",","",.) %>% gsub(" ","",.)
    
    # claulated values:
    net_ef<-f.net_ef(x.sf)
    names(x.sf)<-c("value", "geometry")
    y<-f.grid_square_area2(x.sf)
    q<-x.sf*y$grid_area
    net_npp<-sum(q$value)
    rm(q,z)
    
    vgpm_integrated[i,1]<-as.character(as.Date(as.numeric(yd), origin=origin)) # date  
    vgpm_integrated[i,2]<-median(x) # vgpm_median
    vgpm_integrated[i,3]<-mean(x) # vgpm_mean
    vgpm_integrated[i,4]<-f.vgpm_cal(vgpm_integrated[i,2])# cal_vgpm_median
    vgpm_integrated[i,5]<-f.vgpm_cal(vgpm_integrated[i,3])# cal_vgpm_mean
    vgpm_integrated[i,6]<-f.ef(vgpm_integrated[i,4]) # export_flux_median
    vgpm_integrated[i,7]<-f.ef(vgpm_integrated[i,5]) # export_flux_mean
    vgpm_integrated[i,8]<-net_npp # net_vgpm_npp (mgC/day)
    vgpm_integrated[i,9]<-as.numeric(area) # window_area (m^2) 
    vgpm_integrated[i,10]<-sum_grid_area # sum_grid_area (m^2)
    vgpm_integrated[i,11]<-var(x) # vgpm_variance
    vgpm_integrated[i,12]<-net_ef # total_export_flux (mgC/day)
    vgpm_integrated[i,13]<-vgpm_integrated[i,12]/vgpm_integrated[i,10] # mean export flux (mgC/m^2/day) for sample window
    
  }
  
  
  n<-names(vgpm_integrated)
  vgpm_integrated<-dplyr::select(vgpm_integrated, n[1], n[12], n[13], n[6], n[7], n[9], n[10], n[8],n[2], n[3],  n[11], n[4], n[5] )
  
  return(vgpm_integrated)}

f.grided_ef<-function(x){
  x$value<-f.vgpm_cal(x$npp)
  x$value<-f.ef(x$value)
  return(x)
}

f.ef_plots<-function(z) {
  
  ef_per_date<-vector(mode = "list", length = length(files))
  nl<-vector(mode = "character", length = length(files))
  for(i in 1:length(files)){
    file<-files[i]
    t<-str_split(file, '\\.')
    t<-t[[1]][2]
    t<-str_split(t, "")
    y<-c(t[[1]][1:4])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)   
    yd<-c(t[[1]][5:7])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)
    origin<-toString(c(y,"-01-01")) %>% gsub(",","",.) %>% gsub(" ","",.)
    nl[i]<-as.character(as.Date(as.numeric(yd), origin=origin))
  }
  
  names(ef_per_date)<-nl  
  
  for(i in 1:length(files)){
    
    file<-files[i]
    t<-read_delim(file, col_select = c(1:3))
    t$value[t$value == -9999] <- NA
    
    y<- subset(t, lat >= min(z$lat) & lat <= max(z$lat)  & lon >= min(z$lon) & lon <= max(z$lon),   select=c(value, lon, lat))
    rm(t)
    x.sf<-y %>%st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    rm(y)
    
    ef_per_date[[i]]<-f.grided_ef(x.sf)
    
  }
  
  return(ef_per_date)}

### old code  ####

setwd(wd)
setwd(data_poc_all)

files<-list.files()

z<-station_box_list[[1]]
st1<-f.ef_data_nc(z)

z<-station_box_list[[2]]
st2<-f.ef_data_nc(z)

z<-station_box_list[[3]]
st3<-f.ef_data_nc(z)

z<-station_box_list[[4]]
st3.5<-f.ef_data_nc(z)

z<-station_box_list[[5]]
st4<-f.ef_data_nc(z)

z<-station_box_list[[6]]
st5<-f.ef_data_nc(z)

data<-vector(mode="list", length=length(station_box_list))

names(data)<-names(station_box_list)
data$st1<-st1
data$st2<-st2
data$st3<-st3
data$st3.5<-st3.5
data$st4<-st4
data$st5<-st5
rm(st1, st2, st3, st3.5, st4, st5)
rm(station_box_list) 

setwd(wd)
setwd(robj)
saveRDS(data, "vgpm_integrated_all_years_202200607.R")

t.end<-Sys.time()
r1.time<-t.end-t.start

f<-(paste0("data run time is ", r1.time," and finshed at ", t.end ))
setwd(wd)

fileConn<-file("run_time.txt")
writeLines(vars, fileConn)
close(fileConn)


# 
# ## trying to read hdf4s ###
# setwd(wd)
# setwd(data_poc_all)
# files<-list.files()
# i<-1
# file<-files[i]
# 
# x<-tidync(file)
# x<- RNetCDF::open.nc(file)
# http://127.0.0.1:43852/help/library/ncdf4/html/print.ncdf4.html
# 
# 
# setwd("/home/brandon/vestawd/omz/data")
# file<-"vgpm_2002185.hdf"
# 
# ### getting xyz
# setwd(wd)
# data_poc<-"../../data/POC_Flux_data_sets/"
# setwd(data_poc)
# files<-list.files()
# file<-files[3]
# t<-read_delim(file, col_select = c(1:3))
# 
# z<-read.delim2(file, header=T, sep= " ")
# a=-179.958328 
# b=-179.875000



