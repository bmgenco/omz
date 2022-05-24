#POC notes


# Chl-a methods:


# ODZ method:sread in ODZ actual
* using nc files from [online repository](https://www.bco-dmo.org/dataset/865316)
* use [tidy nc package](https://ropensci.org/blog/2019/11/05/tidync/)
+[additional reference](https://uomresearchit.github.io/r-tidyverse-intro/04-dplyr/)

create  raster of number of observations from odz atlas:
  
  ```{r}
file<-(file.path(data_d, "odz_atlas/nc_depth.nc"))
x<-tidync(file) %>% hyper_tibble(select_var = "numObs")
rm(file)
y<-x%>% group_by(Longitude, Latitude) %>% summarise(sum_obs = sum(numObs)) %>% filter(., sum_obs >1) %>% as.data.frame(.) #return 2d dataframe of depth integareted observations

print(paste0("For cells with a minimum of one observation the ", summary(y)[3,3]))

# data->  sf:
numobs_pts.odz<- sf::st_as_sf(y, coords = c("Longitude","Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 

x<-rasterFromXYZ(y) # raster for mask:
rm(y)
crs(x) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
numobs_raster.odz<-x
rm(x)
```

# Chl-a methods:


# ODZ method:sread in ODZ actual
* using nc files from [online repository](https://www.bco-dmo.org/dataset/865316)
* use [tidy nc package](https://ropensci.org/blog/2019/11/05/tidync/)
+[additional reference](https://uomresearchit.github.io/r-tidyverse-intro/04-dplyr/)

create  raster of number of observations from odz atlas:
  
  ```{r}
file<-(file.path(data_d, "odz_atlas/nc_depth.nc"))
x<-tidync(file) %>% hyper_tibble(select_var = "numObs")
rm(file)
y<-x%>% group_by(Longitude, Latitude) %>% summarise(sum_obs = sum(numObs)) %>% filter(., sum_obs >1) %>% as.data.frame(.) #return 2d dataframe of depth integareted observations

print(paste0("For cells with a minimum of one observation the ", summary(y)[3,3]))

# data->  sf:
numobs_pts.odz<- sf::st_as_sf(y, coords = c("Longitude","Latitude"), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 

x<-rasterFromXYZ(y) # raster for mask:
rm(y)
crs(x) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
numobs_raster.odz<-x
rm(x)
```





## subset
copied from graphing , should make a  function....
```{r}
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


#make alsit for next steps

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



```


1) list of files
2) list of 




#bah


```{r}
z2<-st_bbox(circle_2)
z2<-data.frame(lon= c(z2[[1]], z2[[3]]), lat = c(z2[[2]], z2[[4]])) %>% st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# # tresting bounding box
# plot(st_geometry(z3), axes = TRUE)
# plot(st_geometry(circle_2), col = 'red', add = TRUE)


```



# old attempt:
r The upper-left corner (lat, lon) is 45N; -140E; the lower-right corner is 30.03597N; -115.5454E
```{r}
setwd(data_poc)
files<-list.files()
file<-files[175]
GDALinfo(file)



hdfImage <- readGDAL(paste("HDF4_EOS:EOS_GRID:", files[1], ":MODIS_Grid_16DAY_250m_500m_VI:250m 16 days EVI", sep = "")) 
> image(hdfImage) 
> writeGDAL(hdfImage, "/Users/name/Documents/R_HDF4/output/hdfImage.tif", drivername = "GTiff", type = "Float32", mvFlag = NA, options=NULL, copy_drivername = "GTiff", setStatistics=FALSE) 



x<-tidync(file) 
# %>% hyper_tibble(select_var = "numObs")


```




```{r}
setwd(data_poc)
files<-list.files()
i<-2
file<-files[i]
t<-read_delim(file, col_select = c(1:3))

z<-station_box_list[[4]]

```

# f.grid_square_area = function(x) {
#   x.area<-vector(mode="numeric", length = nrow(x)) 
#   
#   corner<-st_bbox(x)[c("xmin", "ymin")]  
#   d<-ceiling(sqrt(nrow(x)))
#   up<-x[(nrow(x)-(2*d-1)),1]
#   right<-x[(nrow(x)-(d-2)),1]
#   y.d<-(up[1,1]$geometry[[1]][2]-st_bbox(x)$ymin[[1]])/2
#   x.d<-(right[1,1]$geometry[[1]][1]-st_bbox(x)$xmin[[1]])/2
#  
# 
#   box<-st_bbox(x)[c("xmin", "ymin")]
#   box[1][[1]]<-box[1][[1]]-x.d
#   box[2][[1]]<-box[2][[1]]-y.d
#   
#   
#   y<-st_make_grid(x,  
#   n = c(ceiling(sqrt(nrow(x)))-1, ceiling(sqrt(nrow(x)))),
#   # n=c(22,23),
#   crs = if (missing(x)) NA_crs_ else st_crs(x),
#   offset = box,
#   what = "polygons",
#   square = TRUE,
#   flat_topped = FALSE)
#   
#   
#   for(i in 1:length(y)){
#   box<-st_bbox(y[i])
#   x.area[i]<- rgeos::bbox2SP(n = box$ymax[[1]], s = box$ymin[[1]], w = box$xmin[[1]], e = box$xmax[[1]],
#                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% st_as_sf(.,) %>% st_area(.,)
#   }
#  return(x.area)}



https://stackoverflow.com/questions/13442461/populating-a-data-frame-in-r-in-a-loop

#extra spatial crap
## functions:

```{r}
# useful links:
f.grid_square_area = function(x) {
  x.area<-vector(mode="numeric", length = nrow(x)) 
  
  
  
  for(i in 1:nrow(x)){
    i<-1
    box<-st_bbox(x[i,1])
    x.area[i]<- rgeos::bbox2SP(n = box$ymax[[1]], s = box$ymin[[1]], w = box$xmin[[1]], e = box$xmax[[1]],
                               proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% st_as_sf(.,) %>% st_area(.,)
  }
  return(x.area)}

corner<-st_bbox(x)[c("xmin", "ymin")]  
d<-ceiling(sqrt(nrow(x.sf)))
up<-x.sf[(nrow(x.sf)-(2*d-1)),1]
right<-x.sf[(nrow(x.sf)-(d-2)),1]
y.d<-up[1,1]$geometry[[1]][2]-st_bbox(x.sf)$ymin[[1]]
x.d<-right[1,1]$geometry[[1]][1]-st_bbox(x.sf)$xmin[[1]]

box<-st_bbox(x)[c("xmin", "ymin")]
box[1][[1]]<-box[1][[1]]+x.d
box[2][[1]]<-box[2][[1]]-y.d
n.offset = 
  
  # trial
  x<-as_tibble(st_make_grid(x.sf))
x<-as_tibble(x.sf, .rows=ceiling(sqrt(nrow(x.sf))), .columns=ceiling(sqrt(nrow(x.sf))) )


(st_bbox(x.sf)$ymax[[1]]-st_bbox(x.sf)$ymin[[1]])/ceiling(sqrt(nrow(x.sf)))

st_bbox(x.sf)$ymax[[1]]-x.sf[1,1]$geometry[[1]][2]


grid<-st_make_grid(x.sf,  
                   n = c(ceiling(sqrt(nrow(x.sf))), ceiling(sqrt(nrow(x.sf)))),
                   crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                   offset = st_bbox(x.sf)[c("xmin", "ymin")],
                   what = "polygons",
                   square = TRUE,
                   flat_topped = FALSE
)

test<-as.tibble(grid)

x<-st_make_grid(x.sf,
                n = nrow(x.sf),
                crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                what = "corners",
                square = TRUE,
                flat_topped = FALSE)


x<-st_make_grid(x.sf, offset = st_bbox(x.sf)[c("xmin", "ymin")],
                crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                what = "polygons",
                square = TRUE,
                flat_topped = FALSE)

x<-st
x<-st_make_grid(x.sf, what="corners", square = T)

plot(st_geometry(x))
plot(st_geometry(x.sf), add=TRUE)
```








```{r}
setwd(robj)
saveRDS(data, "vgpm_integrated_20220510.R")
```


```{r}
setwd(data_poc)
files<-list.files()

f.input<-function(files, z) {
  vgpm_integrated<-matrix(ncol = 7, nrow =length(files))
  names(vgpm_integrated)<-c("date", "vgpm_median", "vgpm_mean", "cal_vgpm_median", "cal_vgpm_mean", "export_flux_median", "export_flux_mean")
  
  for(i in 1:length(files)){
    
    file<-files[i]
    t<-read_delim(file, col_select = c(1:3))
    t$value[t$value == -9999] <- NA
    
    x<- subset(t, lat >= min(z$lat) & lat <= max(z$lat)  & lon >= min(z$lon) & lon <= max(z$lon),   select=c(value, lon, lat)) 
    # x.sf<-st1 %>%st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    t<-str_split(file, '\\.')
    
    t<-t[[1]][2]
    t<-str_split(t, "")
    y<-c(t[[1]][1:4])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)   
    yd<-c(t[[1]][5:7])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)
    origin<-toString(c(y,"-01-01")) %>% gsub(",","",.) %>% gsub(" ","",.)
    
    vgpm_integrated[i,1]<-as.Date(as.numeric(yd), origin=origin)
    vgpm_integrated[i,2]<-median(x)
    vgpm_integrated[i,3]<-mean(x)
    vgpm_integrated[i,4]<-f.f.vgpm_cal(vgpm_integrated[i,2])
    vgpm_integrated[i,5]<-f.f.vgpm_cal(vgpm_integrated[i,5])
    vgpm_integrated[i,6]<-f.ef(vgpm_integrated[i,4])
    vgpm_integrated[i,7]f.ef(vgpm_integrated[i,5])
  }
  
  return(vgpm_integrated)}

data<-lapply(station_box_list, f.input)



#function, FOR LOOP


t<-read_delim(file, col_select = c(1:3))

t$value[t$value == -9999] <- NA

# should be a loop or function: 

x<- subset(t, lat >= min(z$lat) & lat <= max(z$lat)  & lon >= min(z$lon) & lon <= max(z$lon),   select=c(value, lon, lat)) 
x.sf<-st1 %>%st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


# t_mapped<-st_as_sf(t, coords = c("lon", "lat"),
#                  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


rm(z, t)
rm(circle_1, circle_2, circle_3, circle_3.5, circle_4, circle_5)

#### grepping options: hardcoding b/c I am lazy.... ####


# grepl("^acs_265.", files[i])
# dig<-y$value[1] %>% sub("^\\d\\.", "", .) %>% as.character %>% str_count # counts digits for rounding **improve** 
# sample.month <- sampling.day %>% round_date(., unit="month") %>% as.character %>% gsub('.{3}$', '', .) # "sampling.day" created by f.time_fix
# file<-str_subset(ts.files, pattern=sample.month)
# t<-file %>%sub("^\\d\\.", "", .)

# return file date from file

t<-str_split(file, '\\.')
t<-t[[1]][2]
t<-str_split(t, "")

y<-c(t[[1]][1:4])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)   
yd<-c(t[[1]][5:7])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)
origin<-toString(c(y,"-01-01")) %>% gsub(",","",.) %>% gsub(" ","",.)
t<-as.Date(as.numeric(yd), origin=origin)




```
## run fitst attempt:

#extra spatial crap
## functions:

```{r}
# useful links:
f.grid_square_area = function(x) {
  x.area<-vector(mode="numeric", length = nrow(x)) 
  
  
  
  for(i in 1:nrow(x)){
    i<-1
    box<-st_bbox(x[i,1])
    x.area[i]<- rgeos::bbox2SP(n = box$ymax[[1]], s = box$ymin[[1]], w = box$xmin[[1]], e = box$xmax[[1]],
                               proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% st_as_sf(.,) %>% st_area(.,)
  }
  return(x.area)}

corner<-st_bbox(x)[c("xmin", "ymin")]  
d<-ceiling(sqrt(nrow(x.sf)))
up<-x.sf[(nrow(x.sf)-(2*d-1)),1]
right<-x.sf[(nrow(x.sf)-(d-2)),1]
y.d<-up[1,1]$geometry[[1]][2]-st_bbox(x.sf)$ymin[[1]]
x.d<-right[1,1]$geometry[[1]][1]-st_bbox(x.sf)$xmin[[1]]

box<-st_bbox(x)[c("xmin", "ymin")]
box[1][[1]]<-box[1][[1]]+x.d
box[2][[1]]<-box[2][[1]]-y.d
n.offset = 
  
  # trial
  x<-as_tibble(st_make_grid(x.sf))
x<-as_tibble(x.sf, .rows=ceiling(sqrt(nrow(x.sf))), .columns=ceiling(sqrt(nrow(x.sf))) )


(st_bbox(x.sf)$ymax[[1]]-st_bbox(x.sf)$ymin[[1]])/ceiling(sqrt(nrow(x.sf)))

st_bbox(x.sf)$ymax[[1]]-x.sf[1,1]$geometry[[1]][2]


grid<-st_make_grid(x.sf,  
                   n = c(ceiling(sqrt(nrow(x.sf))), ceiling(sqrt(nrow(x.sf)))),
                   crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                   offset = st_bbox(x.sf)[c("xmin", "ymin")],
                   what = "polygons",
                   square = TRUE,
                   flat_topped = FALSE
)

test<-as.tibble(grid)

x<-st_make_grid(x.sf,
                n = nrow(x.sf),
                crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                what = "corners",
                square = TRUE,
                flat_topped = FALSE)


x<-st_make_grid(x.sf, offset = st_bbox(x.sf)[c("xmin", "ymin")],
                crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                what = "polygons",
                square = TRUE,
                flat_topped = FALSE)

x<-st
x<-st_make_grid(x.sf, what="corners", square = T)

plot(st_geometry(x))
plot(st_geometry(x.sf), add=TRUE)
```








```{r}
setwd(robj)
saveRDS(data, "vgpm_integrated_20220510.R")
```


```{r}
setwd(data_poc)
files<-list.files()

f.input<-function(files, z) {
  vgpm_integrated<-matrix(ncol = 7, nrow =length(files))
  names(vgpm_integrated)<-c("date", "vgpm_median", "vgpm_mean", "cal_vgpm_median", "cal_vgpm_mean", "export_flux_median", "export_flux_mean")
  
  for(i in 1:length(files)){
    
    file<-files[i]
    t<-read_delim(file, col_select = c(1:3))
    t$value[t$value == -9999] <- NA
    
    x<- subset(t, lat >= min(z$lat) & lat <= max(z$lat)  & lon >= min(z$lon) & lon <= max(z$lon),   select=c(value, lon, lat)) 
    # x.sf<-st1 %>%st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    t<-str_split(file, '\\.')
    
    t<-t[[1]][2]
    t<-str_split(t, "")
    y<-c(t[[1]][1:4])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)   
    yd<-c(t[[1]][5:7])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)
    origin<-toString(c(y,"-01-01")) %>% gsub(",","",.) %>% gsub(" ","",.)
    
    vgpm_integrated[i,1]<-as.Date(as.numeric(yd), origin=origin)
    vgpm_integrated[i,2]<-median(x)
    vgpm_integrated[i,3]<-mean(x)
    vgpm_integrated[i,4]<-f.f.vgpm_cal(vgpm_integrated[i,2])
    vgpm_integrated[i,5]<-f.f.vgpm_cal(vgpm_integrated[i,5])
    vgpm_integrated[i,6]<-f.ef(vgpm_integrated[i,4])
    vgpm_integrated[i,7]f.ef(vgpm_integrated[i,5])
  }
  
  return(vgpm_integrated)}

data<-lapply(station_box_list, f.input)



#function, FOR LOOP


t<-read_delim(file, col_select = c(1:3))

t$value[t$value == -9999] <- NA

# should be a loop or function: 

x<- subset(t, lat >= min(z$lat) & lat <= max(z$lat)  & lon >= min(z$lon) & lon <= max(z$lon),   select=c(value, lon, lat)) 
x.sf<-st1 %>%st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


# t_mapped<-st_as_sf(t, coords = c("lon", "lat"),
#                  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


rm(z, t)
rm(circle_1, circle_2, circle_3, circle_3.5, circle_4, circle_5)

#### grepping options: hardcoding b/c I am lazy.... ####


# grepl("^acs_265.", files[i])
# dig<-y$value[1] %>% sub("^\\d\\.", "", .) %>% as.character %>% str_count # counts digits for rounding **improve** 
# sample.month <- sampling.day %>% round_date(., unit="month") %>% as.character %>% gsub('.{3}$', '', .) # "sampling.day" created by f.time_fix
# file<-str_subset(ts.files, pattern=sample.month)
# t<-file %>%sub("^\\d\\.", "", .)

# return file date from file

t<-str_split(file, '\\.')
t<-t[[1]][2]
t<-str_split(t, "")

y<-c(t[[1]][1:4])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)   
yd<-c(t[[1]][5:7])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)
origin<-toString(c(y,"-01-01")) %>% gsub(",","",.) %>% gsub(" ","",.)
t<-as.Date(as.numeric(yd), origin=origin)




```
## run fitst attempt:



#extra spatial crap
## functions:

```{r}
# useful links:
f.grid_square_area = function(x) {
  x.area<-vector(mode="numeric", length = nrow(x)) 
  
  
  
  for(i in 1:nrow(x)){
    i<-1
    box<-st_bbox(x[i,1])
    x.area[i]<- rgeos::bbox2SP(n = box$ymax[[1]], s = box$ymin[[1]], w = box$xmin[[1]], e = box$xmax[[1]],
                               proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% st_as_sf(.,) %>% st_area(.,)
  }
  return(x.area)}

corner<-st_bbox(x)[c("xmin", "ymin")]  
d<-ceiling(sqrt(nrow(x.sf)))
up<-x.sf[(nrow(x.sf)-(2*d-1)),1]
right<-x.sf[(nrow(x.sf)-(d-2)),1]
y.d<-up[1,1]$geometry[[1]][2]-st_bbox(x.sf)$ymin[[1]]
x.d<-right[1,1]$geometry[[1]][1]-st_bbox(x.sf)$xmin[[1]]

box<-st_bbox(x)[c("xmin", "ymin")]
box[1][[1]]<-box[1][[1]]+x.d
box[2][[1]]<-box[2][[1]]-y.d
n.offset = 
  
  # trial
  x<-as_tibble(st_make_grid(x.sf))
x<-as_tibble(x.sf, .rows=ceiling(sqrt(nrow(x.sf))), .columns=ceiling(sqrt(nrow(x.sf))) )


(st_bbox(x.sf)$ymax[[1]]-st_bbox(x.sf)$ymin[[1]])/ceiling(sqrt(nrow(x.sf)))

st_bbox(x.sf)$ymax[[1]]-x.sf[1,1]$geometry[[1]][2]


grid<-st_make_grid(x.sf,  
                   n = c(ceiling(sqrt(nrow(x.sf))), ceiling(sqrt(nrow(x.sf)))),
                   crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                   offset = st_bbox(x.sf)[c("xmin", "ymin")],
                   what = "polygons",
                   square = TRUE,
                   flat_topped = FALSE
)

test<-as.tibble(grid)

x<-st_make_grid(x.sf,
                n = nrow(x.sf),
                crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                what = "corners",
                square = TRUE,
                flat_topped = FALSE)


x<-st_make_grid(x.sf, offset = st_bbox(x.sf)[c("xmin", "ymin")],
                crs = if (missing(x.sf)) NA_crs_ else st_crs(x.sf),
                what = "polygons",
                square = TRUE,
                flat_topped = FALSE)

x<-st
x<-st_make_grid(x.sf, what="corners", square = T)

plot(st_geometry(x))
plot(st_geometry(x.sf), add=TRUE)
```








```{r}
setwd(robj)
saveRDS(data, "vgpm_integrated_20220510.R")
```


```{r}
setwd(data_poc)
files<-list.files()

f.input<-function(files, z) {
  vgpm_integrated<-matrix(ncol = 7, nrow =length(files))
  names(vgpm_integrated)<-c("date", "vgpm_median", "vgpm_mean", "cal_vgpm_median", "cal_vgpm_mean", "export_flux_median", "export_flux_mean")
  
  for(i in 1:length(files)){
    
    file<-files[i]
    t<-read_delim(file, col_select = c(1:3))
    t$value[t$value == -9999] <- NA
    
    x<- subset(t, lat >= min(z$lat) & lat <= max(z$lat)  & lon >= min(z$lon) & lon <= max(z$lon),   select=c(value, lon, lat)) 
    # x.sf<-st1 %>%st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    t<-str_split(file, '\\.')
    
    t<-t[[1]][2]
    t<-str_split(t, "")
    y<-c(t[[1]][1:4])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)   
    yd<-c(t[[1]][5:7])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)
    origin<-toString(c(y,"-01-01")) %>% gsub(",","",.) %>% gsub(" ","",.)
    
    vgpm_integrated[i,1]<-as.Date(as.numeric(yd), origin=origin)
    vgpm_integrated[i,2]<-median(x)
    vgpm_integrated[i,3]<-mean(x)
    vgpm_integrated[i,4]<-f.f.vgpm_cal(vgpm_integrated[i,2])
    vgpm_integrated[i,5]<-f.f.vgpm_cal(vgpm_integrated[i,5])
    vgpm_integrated[i,6]<-f.ef(vgpm_integrated[i,4])
    vgpm_integrated[i,7]f.ef(vgpm_integrated[i,5])
  }
  
  return(vgpm_integrated)}

data<-lapply(station_box_list, f.input)



#function, FOR LOOP


t<-read_delim(file, col_select = c(1:3))

t$value[t$value == -9999] <- NA

# should be a loop or function: 

x<- subset(t, lat >= min(z$lat) & lat <= max(z$lat)  & lon >= min(z$lon) & lon <= max(z$lon),   select=c(value, lon, lat)) 
x.sf<-st1 %>%st_as_sf(., coords = c("lon", "lat" ), crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


# t_mapped<-st_as_sf(t, coords = c("lon", "lat"),
#                  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


rm(z, t)
rm(circle_1, circle_2, circle_3, circle_3.5, circle_4, circle_5)

#### grepping options: hardcoding b/c I am lazy.... ####


# grepl("^acs_265.", files[i])
# dig<-y$value[1] %>% sub("^\\d\\.", "", .) %>% as.character %>% str_count # counts digits for rounding **improve** 
# sample.month <- sampling.day %>% round_date(., unit="month") %>% as.character %>% gsub('.{3}$', '', .) # "sampling.day" created by f.time_fix
# file<-str_subset(ts.files, pattern=sample.month)
# t<-file %>%sub("^\\d\\.", "", .)

# return file date from file

t<-str_split(file, '\\.')
t<-t[[1]][2]
t<-str_split(t, "")

y<-c(t[[1]][1:4])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)   
yd<-c(t[[1]][5:7])%>% toString(.) %>% gsub(",","",.) %>% gsub(" ","",.)
origin<-toString(c(y,"-01-01")) %>% gsub(",","",.) %>% gsub(" ","",.)
t<-as.Date(as.numeric(yd), origin=origin)




```
## run fitst attempt:





