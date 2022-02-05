cell_list<-match
x<-cell_list %>% lapply(.,'[[',1 ) %>% lapply(., dim) %>%
  lapply(., function(x){return(x[1])}) %>% .[order(unlist(.), decreasing=TRUE)]
x<-names(x)
cell_list<-cell_list[x]
rm(x)

variables_used<-year_begin + spatial_box + (select.list(function_variables, sensor, fraction_fodz))
cell_list<-append(variables_used, cell_list)

if(is.na(function_variables$argo_sensor)==TRUE){
  cell_floats= select_profiles(lon_lim,lat_lim,start_date, end_date,
                               sensor=c("DOXY"), outside="both")} else{cell_floats= select_profiles(lon_lim,lat_lim,
                                                                                                    start_date, end_date,sensor=sensor_select, outside="both")}

spatial_box<-c(unname(st_bbox(y)$xmin), unname(st_bbox(y)$ymin), unname(st_bbox(y)$xmax), unname(st_bbox(y)$ymax))

if(is.na(function_variables$argo_sensor)==TRUE){
  cell_floats= select_profiles(lon_lim,lat_lim,start_date, end_date,
                               sensor=c("DOXY"), outside="both")} else{cell_floats= select_profiles(lon_lim,lat_lim,
start_date, end_date,sensor=sensor_select, outside="both")}

### older functions ####

test<-lengths(match)
x<-cell_list %>% lapply(.,'[[',1 ) %>% lapply(., dim) %>% 
lapply(., function(x){return(x[1])}) %>% .[order(unlist(.), decreasing=TRUE)]
x<-names(x)
cell_list<-cell_list[x]

List2<-List2[c("M2","M1")]
x<-match[[30]]
z<-x$profile_dates

z<-with_tz(z, tz="UTC")

y<-(time_length(diff(z, lag=1, difference=1), unit="days"))
summary(y) # outlier explains this

spatial_box<-c(-115,0, -100,23) #(xmin,ymin, xmax, ymax)

f.argo_by_database_w_selection_options<-function(database, year_begin, function_variables, spatial_box){
 
   ##############################################################################
  #                USER defined variables and defaults:                        #
  #                                                                            #
  # spatial_box = c(xmin,ymin, xmax, ymax). vector of 4 integers. nothing else #
  # function
  
  ##############################################################################
  
  #### default values if not defined: #### In FUTURE make spatail_box and year begin function variables
  if(missing(spatial_box)) {
    y<-database
  } else {
  y<-st_crop(database,  xmin=spatial_box[1], ymin=spatial_box[2], xmax=spatial_box[3], ymax=spatial_box[4])
  }
  if(missing(year_begin))  {
    year_begin<-year(min(database$DateTime)) %>% as.integer(.)
  }
  if(function_variables$fraction_fodz=NULL){function_variables$fraction_fodz<-0}
  if(function_variables$h.prior=NULL){function_variables$h.prior<-90}
  if(function_variables$h.post=NULL){function_variables$h.post<-90}
  
  #### Function : ####
  y<-filter(y, year(DateTime) >= year_begin, maxFODZ >= function_variables$fraction_fodz) 
  y<-y[order(y$grid_id, decreasing = F),]
  id<-unique(y$grid_id)
  function_list<-as.list(id)
  names(function_list)<-id
  
  f.get<-function(x, y){
    y.id<-x
    data<-filter(y, grid_id == y.id)
    
    z<-st_bbox(unique(data))
    lat_lim=c(z$ymin[[1]], z$ymax[[1]])
    lon_lim=c(z$xmin[[1]], z$xmax[[1]])
    start_date=date(min(y$DateTime)) - function_variables$h.prior
    end_date=date(max(y$DateTime)) + function_variables$h.post
    
    #cell_floats from bioargo submodule:
    cell_floats= select_profiles(lon_lim,lat_lim,start_date, end_date,
            sensor=function_variables$argo_sensor, outside="both" #  All floats that cross into time/space limits
    )
    #
    # 'outside', 'none' 'time' 'space' 'both': By default, only float profiles
    #           that are within both the temporal and spatial constraints are
    #           returned ('none'); specify to also maintain profiles outside
    #           the temporal constraints ('time'), spatial constraints
    #           ('space'), or both constraints ('both')
    # 
    # 'sensor', 'SENSOR_TYPE': By default, all floats within the lon/lat/time
    #           limits are considered. This option allows the selection by 
    #           sensor type. Available are: DOXY, CHLA, BBP700, 
    #           PH_IN_SITU_TOTAL, NITRATE, DOWN_IRRADIANCE380,
    #           DOWN_IRRADIANCE412, DOWN_IRRADIANCE490, DOWNWELLING_PAR
    #           (Currently, only one sensor type can be selected.)
    
    p<-vector(mode="list")
    p$data<-data
    p$profiles<-cell_floats$profiles
    p$floats<-cell_floats$floats
    p$profile_dates<-Sprof$date[cell_floats$profiles]
    
    # x<-match$'36762'
    # overpass_time<-x$data$DateTime
    # dates<-
    # 
    # 
    # # window function:
    # f.window<-function(p){
    #   
    # }
    # p$match<-f.wind(p)
    # 
    
    return(p)
    
  }
  
  cell_list<-lapply(function_list, y=y, f.get)
  rm(function_list)
  names(cell_list)<-id
  return(cell_list)

  }
f.argo_by_database<-function(database, year_begin){
  y<-filter(database, year(DateTime) >= year_begin) %>% .[order( .$grid_id, decreasing = F),]
  id<-unique(y$grid_id)
  # cell_list<-vector(mode="list", length=length(id))
  function_list<-as.list(id)
  names(function_list)<-id
  
  f.get<-function(x, y){
    y.id<-x
    data<-filter(y, grid_id == y.id)
    
    z<-st_bbox(unique(data))
    lat_lim=c(z$ymin[[1]], z$ymax[[1]])
    lon_lim=c(z$xmin[[1]], z$xmax[[1]])
    start_date=date(min(y$DateTime)) - 90
    end_date=date(max(y$DateTime)) + 90
    
    cell_floats= select_profiles(lon_lim,lat_lim,start_date, end_date,
                                 sensor=c(), 
                                 outside="both" #  All floats that cross into time/space limits
    )
    
    # 'DOXY'  
    p<-vector(mode="list")
    p$data<-data
    p$profiles<-cell_floats$profiles
    p$floats<-cell_floats$floats
    return(p)
    
  }
  
  
  cell_list<-lapply(function_list, y=y, f.get)
  rm(function_list)
  names(cell_list)<-id
  return(cell_list)
}
f.argo_by_database_w_clipFODZ<-function(database, year_begin, function_variables){
  y<-filter(database, year(DateTime) >= year_begin, maxFODZ >= function_variables$fraction_fodz) 
  y<-y[order(y$grid_id, decreasing = F),]
  id<-unique(y$grid_id)
  function_list<-as.list(id)
  names(function_list)<-id
  
  f.get<-function(x, y){
    y.id<-x
    data<-filter(y, grid_id == y.id)
    
    z<-st_bbox(unique(data))
    lat_lim=c(z$ymin[[1]], z$ymax[[1]])
    lon_lim=c(z$xmin[[1]], z$xmax[[1]])
    start_date=date(min(y$DateTime)) - 90
    end_date=date(max(y$DateTime)) + 90
    
    cell_floats= select_profiles(lon_lim,lat_lim,start_date, end_date,
                                 sensor=c("DOXY"), 
                                 outside="both" #  All floats that cross into time/space limits
    )
    
    # 'DOXY'  
    p<-vector(mode="list")
    p$data<-data
    p$profiles<-cell_floats$profiles
    p$floats<-cell_floats$floats
    return(p)
    
  }
  
  cell_list<-lapply(function_list, y=y, f.get)
  rm(function_list)
  names(cell_list)<-id
  return(cell_list)
}


f.argo_by_database<-function(database, year_begin){
  y<-filter(database, year(DateTime) >= year_begin) %>% .[order( .$grid_id, decreasing = F),]
  id<-unique(y$grid_id)
  # cell_list<-vector(mode="list", length=length(id))
  function_list<-as.list(id)
  names(function_list)<-id
  
  f.get<-function(x, y){
    y.id<-x
    data<-filter(y, grid_id == y.id)
    
    z<-st_bbox(unique(data))
    lat_lim=c(z$ymin[[1]], z$ymax[[1]])
    lon_lim=c(z$xmin[[1]], z$xmax[[1]])
    start_date=date(min(y$DateTime)) - 90
    end_date=date(max(y$DateTime)) + 90
    
    cell_floats= select_profiles(lon_lim,lat_lim,start_date, end_date,
                                 sensor=c(), 
                                 outside="both" #  All floats that cross into time/space limits
    )
    
    # 'DOXY'  
    p<-vector(mode="list")
    p$data<-data
    p$profiles<-cell_floats$profiles
    p$floats<-cell_floats$floats
    return(p)
    
  }
  
  cell_list<-lapply(function_list, y=y, f.get)
  rm(function_list)
  names(cell_list)<-id
  return(cell_list)
}
f.argo_by_database_w_clipFODZ<-function(database, summary_2d.odz, year_begin, function_variables){
  y<-filter(database, year(DateTime) >= year_begin) 
  z<-maxf %>% filter(., summary_2d.odz >= function_variables$fraction_fodz)
  y<-st_difference(y,z)
  rm(z)
  y<-y[order(y$grid_id, decreasing = F),]
  id<-unique(y$grid_id)
  # cell_list<-vector(mode="list", length=length(id))
  function_list<-as.list(id)
  names(function_list)<-id
  
  f.get<-function(x, y){
    y.id<-x
    data<-filter(y, grid_id == y.id)
    
    z<-st_bbox(unique(data))
    lat_lim=c(z$ymin[[1]], z$ymax[[1]])
    lon_lim=c(z$xmin[[1]], z$xmax[[1]])
    start_date=date(min(y$DateTime)) - 90
    end_date=date(max(y$DateTime)) + 90
    
    cell_floats= select_profiles(lon_lim,lat_lim,start_date, end_date,
                                 sensor=c("DOXY"), 
                                 outside="both" #  All floats that cross into time/space limits
    )
    
    # 'DOXY'  
    p<-vector(mode="list")
    p$data<-data
    p$profiles<-cell_floats$profiles
    p$floats<-cell_floats$floats
    return(p)
    
  }
  
  cell_list<-lapply(function_list, y=y, f.get)
  rm(function_list)
  names(cell_list)<-id
  return(cell_list)
}


#### combining datbase, needs to be improved create selection function ####
# https://community.rstudio.com/t/performing-a-full-join-on-sf-objects/43902/10
# https://r-spatial.github.io/sf/reference/st_make_grid.html

# older...
# 
# st_geometry(h.buf)<-"points"
# y<-st_join( x, h.buf, join = st_nearest_feature, left = T) # this keep geometry of x
# # function here
# # get boudning box for each entry.. # https://rdrr.io/cran/sf/man/st_bbox.html
# z<-y[order(y$DateTime, decreasing = T),]
# z1<-st_bbox(z[1,]$geometry)
# # Get coresponding date... maby form a difnfet set of joins... all dates and strom keys
# z2<-z[1,]$DateTime
# # convert to lubridate
###   select ####

# from Main _workshop script:
# lat_lim=c(45, 60)
# lon_lim=c(-150, -135)
# start_date="2008-01-01"
# end_date="2018-12-31"
# 
# OSP_data= select_profiles ( lon_lim,
#                             lat_lim,
#                             start_date,
#                             end_date,
#                             sensor=c(), # this selects only floats with nitrate sensors
#                             outside="both" #  All floats that cross into the time/space limits
# )
# 
# sensor=c('NITRATE'),

# x<-filter(database, year(DateTime) >= 2000) %>% .[order( .$weight, decreasing = T),]
# x.cells<-unique(x$cell)
# # x.grid<-unique()
# x.id<-x$grid_id[1]
# y<-filter(x, grid_id == x.id)
# z<-st_bbox(unique(y))




# lat_lim=c(z$ymin[[1]], z$ymin[[1]])
# lon_lim=c(z$xmin[[1]], z$xmax[[1]])
# start_date=date(min(y$DateTime)) - function_variables$h.post
# end_date=date(max(y$DateTime)) + function_variables$h.prior
# 
# 
# lat_lim=c(45, 60)
# lon_lim=c(-150, -135)

# Select profiles based on those limits with specified sensor (NITRATE)

# year_begin=1990

# match<-f.argo_by_database(database, year_begin<-2000)




f.buffer_select<-function(h.pts, function_variables ,cl){
  
  # y<-st_combine(st_geometry(cl))
  # y<-st_centroid(st_geometry(cl))
  # 
  # x<-sf:st_buffer(cl, dist=)
  # x<-st_difference(h.pts ,x)
  
  x<-h.pts
  radi.m<-function_variables$h.radi
  y<-st_geometry(x)
  
  f.circle<-function(y, radi.m){
    lat0<-y[[1]][2]
    lon0<-y[[1]][1]
    center.reproj<-paste0("+proj=aeqd +lat_0=", lat0, " ", "+lon_0=", lon0)
    flat<-sf::st_transform(y, center.reproj )
    circle<-sf::st_buffer(flat, dist=radi.m)%>%sf::st_transform(., "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 
    return(circle)
  }
  
  circle<-sf::st_as_sf(f.circle(y, radi.m))
  names(circle)<-"buffer"
  st_geometry(circle)<-"buffer"
  z<-cbind(x, circle)
  names(z)[names(z) == "geometry"] <- "points"
  st_geometry(z)<-"buffer"
  return(z)
  
}



plot(x,  main="Overlap: ODZ-Atlas (red/black) and Hurricane buffer (yellow)",legend=F )
plot(st_geometry(overlap.odz), col="black", add=T, legend=F)
plot(filtered_numobs.raster, col="red", add=T, legend=F)

# file not used needs improvement
# h<-read_csv2(file.path(data_d, "/all_omz_hurricanes/hurdat2-nepac-1949-2020-043021a.txt"),trim_ws=T, col_names = F, show_col_types=F)

# nc files: https://rpubs.com/boyerag/297592


#### Option 1 netcdf: ####
# x<-nc_open(file.path(data_d, "odz_atlas/nc_depth.nc"))
# lat<-ncvar_get(x, "Latitude")
# lon<-ncvar_get(x, "Longitude")
# z<-ncvar_get(x, "Depth")
# obs<-ncvar_get(x, "numObs")
# nc_close(x)
# 

# come back here:
# derivdepth<-ncvar_get(x, "maxDerivDepth")
# fodz<-ncvar_get(x, "fODZ")
# o2<-ncvar_get(x, "O2")


# test<-obs[,,1]
# 
# x<-dim()
# 
# scatter3D(obs[(1:256),,], obs[,(1:230),], obs[,,(1:50)])
# 
# x<-obs[(1:dim(obs)[1]),1,1]
# y<-obs[1,(1:dim(obs)[2]),1]
# z<-obs[1,1,(1:dim(obs)[3])]
# 
# scatter3D(x,y,z)
# 