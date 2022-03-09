#version1  
# based on " f.argo_by_database_w_selection_options" but edited accrodingly 
f.argo_by_hurdat_w_selection_options<-function(h.buf, year_begin, function_variables, spatial_box){
  
  #                USER defined variables and defaults:                                   #
  # year begin = year begin form hurricanes, profies will be function_variables$hprior b4 #                                                                               #
  # spatial_box = c(xmin,ymin, xmax, ymax). vector of 4 integers. nothing else            #
  # if function_variables$argo_sensor==NA THEN no sensor limits on profile return         #
  
  #### default values if not defined: #### In FUTURE make spatial_box and year begin function variables
  
  if(missing(year_begin)){
    year_begin<-year(min(h.buf$DateTime)) %>% as.integer(.)
    
  }
  if(missing(spatial_box)) {
    y<-h.buf
    spatial_box<-c(unname(st_bbox(h.buf)$xmin), unname(st_bbox(h.buf)$ymin), unname(st_bbox(h.buf)$xmax), unname(st_bbox(h.buf)$ymax))
  } else {y<-st_crop(h.buf,  xmin=spatial_box[1], ymin=spatial_box[2], xmax=spatial_box[3], ymax=spatial_box[4])}
  

  if(is.na(function_variables$h.prior)==TRUE){function_variables$h.prior<-90}
  if(is.na(function_variables$h.post)==TRUE){function_variables$h.post<-90}
  
  names(spatial_box)<-c("xmin", "ymin", "xmax", "ymax")
  
  #### Function : ####
  y<-filter(y, year(DateTime) >= year_begin)
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
    if(is.na(function_variables$argo_sensor)==TRUE) {cell_floats<-select_profiles(lon_lim,lat_lim,start_date, end_date, outside="both",
                                                                                  sensor=NULL)} else {cell_floats<-select_profiles(lon_lim,lat_lim,
                                                                                                                                   start_date, end_date, sensor=function_variables$argo_sensor, outside="both")}  
    
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
    
    data.select<-p$data %>% filter(., year(DateTime)>=year_begin) %>%
      select(., Key, Name, DateTime, lon, lat, Wind) %>% st_drop_geometry(.)%>%
      remove_rownames(.)
    
    p$events<-vector(mode="list", length=nrow(data.select))
    for(i in 1:length(p$events)){
      p$events[[i]]$data<-data.select[i,]}
    names(p$events)<-data.select$Key
    rm(data.select)
    
    p$profiles<-cell_floats$profiles
    p$floats<-cell_floats$floats
    p$profile_dates<-Sprof$date[cell_floats$profiles]
    
    f.closest<-function(x){
      f.distance<-function(y){
        pts<-SpatialPoints(coords=cbind(y$lon, y$lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), bbox = NULL)
        hpt<- SpatialPoints(coords=cbind(x$data$lon, x$data$lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), bbox = NULL)
        y$distance_km<-round(spDistsN1(pts, hpt, longlat = TRUE), 0)
        return(y)
      }
      
      before<-p$profile_dates[p$profile_dates<=x$data$DateTime]
      after<-p$profile_dates[p$profile_dates>=x$data$DateTime]
      
      if(length(before)==0){x$prior<-NA} else {
        
        before_window<-before[which(date(before) >= date(x$data$DateTime)-function_variables$h.prior)]
        
        if(length(before_window)==0){x$prior<-NA} else {
          
          x$prior<-data.frame(matrix(ncol=6, nrow=length(before_window))) # profiles withins user define function_
          # names(x$prior)<-c("wmo", "profile_index",  "DateTime", "lat", "lon", "distance_km", "days_prior")
          names(x$prior)<-c("profile_index",  "DateTime", "lat", "lon", "distance_km", "days_prior")
          x$prior$DateTime<-before_window
          before_p<-p$profiles[p$profile_dates<=x$data$DateTime]
          x$prior$profile_index<-before_p[which(date(before) >= date(x$data$DateTime)-function_variables$h.prior)]
          # need to figure out profile wmos
          # x$prior$wmo<-Sprof$wmo[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
          # change here for date/lubridate object if need be
          x$prior$days_prior<-as.integer((date(x$prior$DateTime) - date(x$data$DateTime)) *-1)
          x$prior$lat<-Sprof$lat[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
          x$prior$lon<-Sprof$lon[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
          rm(before,  before_p, before_window)
          x$prior<-f.distance(y=x$prior)
        }
      }
      
      if(length(after)==0){x$post<-NA} else {
        
        after_window<-after[which(date(after) <= date(x$data$DateTime)+function_variables$h.post)]
        
        if(length(after_window)==0){x$post<-NA} else{
          
          x$post<-data.frame(matrix(ncol=6, nrow=length(after_window)))
          # names(x$post)<-c("wmo", "profile_index","DateTime", "lat", "lon", "distance_km", "days_after")
          names(x$post)<-c("profile_index","DateTime", "lat", "lon", "distance_km", "days_after")
          x$post$DateTime<-after_window
          after_p<-p$profiles[p$profile_dates>=x$data$DateTime]
          x$post$profile_index<-after_p[which(date(after) <= date(x$data$DateTime)+function_variables$h.post)] 
          # need to figure out profile wmos
          # x$post$wmo<-Sprof$wmo[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
          # change here for date/lubridate object if need be
          x$post$days_after<-as.integer(date(x$post$DateTime)-date(x$data$DateTime))
          x$post$lat<-Sprof$lat[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
          x$post$lon<-Sprof$lon[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
          rm(after, after_p, after_window)
          x$post<-f.distance(y=x$post)
        }
      }
      
      return(x)}
    
    p$events<-lapply(p$events, f.closest)
    p$events<-p$events[!sapply(p$events,is.null)]
    return(p)
    
  }
  
  cell_list<-lapply(function_list, y=y, f.get)
  rm(function_list)
  names(cell_list)<-id
  f.rm<-function(x){if(length(x$profiles)==0){x<-NULL}
    return(x)}
  
  cell_list<-lapply(cell_list, f.rm)
  cell_list<-cell_list[!sapply(cell_list,is.null)]
  rm(f.get, f.rm)
  
  f.order<-function(cell_list){
    x<-cell_list %>% lapply(.,'[[',1 ) %>% lapply(., dim) %>% lapply(., function(x){return(x[1])})
    y<-x[order(as.vector(unlist(x)), decreasing=TRUE)]
    cell_list<-cell_list[names(y)]
    return(cell_list)}
  
  cell_list<-f.order(cell_list)
  
  # names(year_begin)<-"year_begin"
  input_variables<-function_variables %>% purrr::keep(names(.) %in% c("argo_sensor","fraction_fodz","h.prior","h.post"))
  input_variables$spatial_box<-spatial_box
  input_variables<-append(input_variables, year_begin, after=0)
  names(input_variables)[1]<-"year_begin"
  cell_list<-append(cell_list, list(input_variables), after=0)
  names(cell_list)[1]<-"input_variables"
  
  # function_variables<<-function_variables # would update function variables
  return(cell_list)
  
}


#version 3
f.argo_by_database_w_selection_options<-function(database, year_begin, function_variables, spatial_box){
  
  #                USER defined variables and defaults:                                   #
  # year begin = year begin form hurricanes, profies will be function_variables$hprior b4 #                                                                               #
  # spatial_box = c(xmin,ymin, xmax, ymax). vector of 4 integers. nothing else            #
  # if function_variables$argo_sensor==NA THEN no sensor limits on profile return         #
  
  #### default values if not defined: #### In FUTURE make spatial_box and year begin function variables
  
  if(missing(year_begin)){
    year_begin<-year(min(database$DateTime)) %>% as.integer(.)
    
  }
  if(missing(spatial_box)) {
    y<-database
    spatial_box<-c(unname(st_bbox(database)$xmin), unname(st_bbox(database)$ymin), unname(st_bbox(database)$xmax), unname(st_bbox(database)$ymax))
  } else {y<-st_crop(database,  xmin=spatial_box[1], ymin=spatial_box[2], xmax=spatial_box[3], ymax=spatial_box[4])}
  
  if(is.na(function_variables$fraction_fodz)==TRUE){function_variables$fraction_fodz<-0}
  if(is.na(function_variables$h.prior)==TRUE){function_variables$h.prior<-90}
  if(is.na(function_variables$h.post)==TRUE){function_variables$h.post<-90}
  
  names(spatial_box)<-c("xmin", "ymin", "xmax", "ymax")
  
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
    if(is.na(function_variables$argo_sensor)==TRUE) {cell_floats<-select_profiles(lon_lim,lat_lim,start_date, end_date, outside="both",
                                                                                  sensor=NULL)} else {cell_floats<-select_profiles(lon_lim,lat_lim,
                                                                                                                                   start_date, end_date, sensor=function_variables$argo_sensor, outside="both")}  
    
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
    
    data.select<-p$data %>% filter(., year(DateTime)>=year_begin) %>%
      select(., Key, Name, DateTime, lon, lat, Wind) %>% st_drop_geometry(.)%>%
      remove_rownames(.)
    
    p$events<-vector(mode="list", length=nrow(data.select))
      for(i in 1:length(p$events)){
      p$events[[i]]$data<-data.select[i,]}
    names(p$events)<-data.select$Key
    rm(data.select)
    
    p$profiles<-cell_floats$profiles
    p$floats<-cell_floats$floats
    p$profile_dates<-Sprof$date[cell_floats$profiles]
    
    f.closest<-function(x){
      f.distance<-function(y){
        pts<-SpatialPoints(coords=cbind(y$lon, y$lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), bbox = NULL)
        hpt<- SpatialPoints(coords=cbind(x$data$lon, x$data$lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), bbox = NULL)
        y$distance_km<-round(spDistsN1(pts, hpt, longlat = TRUE), 0)
        return(y)
      }
      
      before<-p$profile_dates[p$profile_dates<=x$data$DateTime]
      after<-p$profile_dates[p$profile_dates>=x$data$DateTime]
      
      if(length(before)==0){x$prior<-NA} else {
        
        before_window<-before[which(date(before) >= date(x$data$DateTime)-function_variables$h.prior)]
        
        if(length(before_window)==0){x$prior<-NA} else {
          
          x$prior<-data.frame(matrix(ncol=6, nrow=length(before_window))) # profiles withins user define function_
          # names(x$prior)<-c("wmo", "profile_index",  "DateTime", "lat", "lon", "distance_km", "days_prior")
          names(x$prior)<-c("profile_index",  "DateTime", "lat", "lon", "distance_km", "days_prior")
          x$prior$DateTime<-before_window
          before_p<-p$profiles[p$profile_dates<=x$data$DateTime]
          x$prior$profile_index<-before_p[which(date(before) >= date(x$data$DateTime)-function_variables$h.prior)]
          # need to figure out profile wmos
          # x$prior$wmo<-Sprof$wmo[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
          # change here for date/lubridate object if need be
          x$prior$days_prior<-as.integer((date(x$prior$DateTime) - date(x$data$DateTime)) *-1)
          x$prior$lat<-Sprof$lat[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
          x$prior$lon<-Sprof$lon[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
          rm(before,  before_p, before_window)
          x$prior<-f.distance(y=x$prior)
        }
      }
      
      if(length(after)==0){x$post<-NA} else {
        
        after_window<-after[which(date(after) <= date(x$data$DateTime)+function_variables$h.post)]
        
        if(length(after_window)==0){x$post<-NA} else{
          
          x$post<-data.frame(matrix(ncol=6, nrow=length(after_window)))
          # names(x$post)<-c("wmo", "profile_index","DateTime", "lat", "lon", "distance_km", "days_after")
          names(x$post)<-c("profile_index","DateTime", "lat", "lon", "distance_km", "days_after")
          x$post$DateTime<-after_window
          after_p<-p$profiles[p$profile_dates>=x$data$DateTime]
          x$post$profile_index<-after_p[which(date(after) <= date(x$data$DateTime)+function_variables$h.post)] 
          # need to figure out profile wmos
          # x$post$wmo<-Sprof$wmo[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
          # change here for date/lubridate object if need be
          x$post$days_after<-as.integer(date(x$post$DateTime)-date(x$data$DateTime))
          x$post$lat<-Sprof$lat[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
          x$post$lon<-Sprof$lon[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
          rm(after, after_p, after_window)
          x$post<-f.distance(y=x$post)
        }
      }
      
      return(x)}
    
    p$events<-lapply(p$events, f.closest)
    p$events<-p$events[!sapply(p$events,is.null)]
    return(p)
  
  }
  
  cell_list<-lapply(function_list, y=y, f.get)
  rm(function_list)
  names(cell_list)<-id
  f.rm<-function(x){if(length(x$profiles)==0){x<-NULL}
    return(x)}
  
  cell_list<-lapply(cell_list, f.rm)
  cell_list<-cell_list[!sapply(cell_list,is.null)]
  rm(f.get, f.rm)
  
  f.order<-function(cell_list){
    x<-cell_list %>% lapply(.,'[[',1 ) %>% lapply(., dim) %>% lapply(., function(x){return(x[1])})
    y<-x[order(as.vector(unlist(x)), decreasing=TRUE)]
    cell_list<-cell_list[names(y)]
    return(cell_list)}
  
  cell_list<-f.order(cell_list)
  
  # names(year_begin)<-"year_begin"
  input_variables<-function_variables %>% purrr::keep(names(.) %in% c("argo_sensor","fraction_fodz","h.prior","h.post"))
  input_variables$spatial_box<-spatial_box
  input_variables<-append(input_variables, year_begin, after=0)
  names(input_variables)[1]<-"year_begin"
  cell_list<-append(cell_list, list(input_variables), after=0)
  names(cell_list)[1]<-"input_variables"
  
  # function_variables<<-function_variables # would update function variables
  return(cell_list)
  
}







