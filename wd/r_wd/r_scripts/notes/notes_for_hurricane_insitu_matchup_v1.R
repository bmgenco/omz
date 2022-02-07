testing





#broken?
f.argo_by_database_w_selection_options<-function(database, year_begin, function_variables, spatial_box){
  
  #                USER defined variables and defaults:                                   #
  # year begin = year begin form hurricanes, profies will be function_variables$hprior b4 #                                                                               #
  # spatial_box = c(xmin,ymin, xmax, ymax). vector of 4 integers. nothing else            #
  # if function_variables$argo_sensor==NA THEN no sensor limits on profile return         #
  # profiles may only be index value will need to update mathc fo each use of Sprof
  
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
    
    p$profiles_postions ### need to do this for adjusting events data
    p$floats_slection ## tmeproary thing for matching floats to profiles??
    
    # f.window<-function(p){p}
    # p$individual_tracks<-vector()
    # p$match<-f.wind(p)
    
    
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
          names(x$prior)<-c("profile_index",  "DateTime", "lat", "lon", "distance_km", "days_prior")
          x$prior$DateTime<-before_window
          before_p<-p$profiles[p$profile_dates<=x$data$DateTime]
          x$prior$profile_index<-before_p[which(date(before) >= date(x$data$DateTime)-function_variables$h.prior)]
          #issue with wmo.. not sure about what wmo is
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
          #issue with wmo.. not sure about what wmo is
          x$post<-data.frame(matrix(ncol=6, nrow=length(after_window)))
          names(x$post)<-c("profile_index","DateTime", "lat", "lon", "distance_km", "days_after")
          x$post$DateTime<-after_window
          after_p<-p$profiles[p$profile_dates>=x$data$DateTime]
          x$post$profile_index<-after_p[which(date(after) <= date(x$data$DateTime)+function_variables$h.post)] 
          
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
  f.rm<-function(x){if(length(p$profiles)==0){x<-NULL} 
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




### older fucntion

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
    return(p)
    # f.window<-function(p){p}
    # p$individual_tracks<-vector()
    # p$match<-f.wind(p)
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





# 
# ## testing:
# p<-bud_match[[2]]
# year_begin<-bud_match$input_variables$year_begin
# 
# data.select<-p$data %>% filter(., year(DateTime)>=(year_begin-1)) %>% 
#   select(., Key, Name, DateTime, lon, lat) %>% st_drop_geometry(.)
# p$events<-vector(mode="list", length=nrow(data.select))
# names(p$events)<-data.select$Key
# 
# for(i in 1:length(p$events)){
#   p$events[[i]]$data<-data.select[i,]}
# rm(data.select)
# 
# 
# 
# x<-p$events[[18]]
# function_variables$h.prior<-90
# function_variables$h.post<-90
# 
# f.closest<-function(x){
# before<-p$profile_dates[p$profile_dates<=x$data$DateTime]
# after<-p$profile_dates[p$profile_dates>=x$data$DateTime]
# 
# before_window<-before[which(date(before) >= date(x$data$DateTime)-function_variables$h.prior)]
# after_window<-after[which(date(after) <= date(x$data$DateTime)+function_variables$h.post)] 
# 
# x$prior<-data.frame(matrix(ncol=8, nrow=length(before_window))) # profiles withins user define function_
# x$post<-data.frame(matrix(ncol=8, nrow=length(after_window)))
# 
# names(x$prior)<-c("float", "profile_index", "wmo", "DateTime", "lat", "lon", "distance_km", "days_prior")
# names(x$post)<-c("float", "profile_index", "wmo", "DateTime", "lat", "lon", "distance_km", "days_after")
# 
# x$prior$DateTime<-before_window
# x$post$DateTime<-after_window
# 
# before_p<-p$profiles[p$profile_dates<=x$data$DateTime]
# after_p<-p$profiles[p$profile_dates>=x$data$DateTime]
# 
# x$prior$profile_index<-before_p[which(date(before) >= date(x$data$DateTime)-function_variables$h.prior)]
# x$post$profile_index<-after_p[which(date(after) <= date(x$data$DateTime)+function_variables$h.post)] 
# 
# x$prior$wmo<-Sprof$wmo[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
# x$post$wmo<-Sprof$wmo[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
# 
# # change here for date/lubridate object if need be
# x$prior$days_prior<-as.integer((date(x$prior$DateTime) - date(x$data$DateTime)) *-1)
# x$post$days_after<-as.integer(date(x$post$DateTime)-date(x$data$DateTime))
# 
# x$prior$lat<-Sprof$lat[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
# x$prior$lon<-Sprof$lon[c(x$prior$profile_index[1]:x$prior$profile_index[length(x$prior$profile_index)])]
# 
# x$post$lat<-Sprof$lat[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
# x$post$lon<-Sprof$lon[c(x$post$profile_index[1]:x$post$profile_index[length(x$post$profile_index)])]
# 
# rm(before, after, before_p, before_window, after_p, after_window)
# 
# f.distance<-function(y){
#   pts<-SpatialPoints(coords=cbind(y$lon, y$lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), bbox = NULL)
#   hpt<- SpatialPoints(coords=cbind(x$data$lon, x$data$lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), bbox = NULL)
#   y$distance_km<-round(spDistsN1(pts, hpt, longlat = TRUE), 0)
#   return(y)
#   }
# 
# x$prior<-f.distance(y=x$prior)
# x$post<-f.distance(y=x$post)
# 
# # add if statements here
# 
# return(x)
# }
# 
# 
# 
# 
# 
# 
# #temp 
# 
# 
# 
# prior<-which(abs(x$data$DateTime - before) == min(abs(x$data$DateTime - before)))
# before_time<-p$profile_dates[prior]
# # overpass_time
# post<-which(abs(x$data$DateTime - after) == min(abs(x$data$DateTime - after)))+length(before)
# after_time<-p$profile_dates[post]
# 
# 
# 
# x$profiles_before<-before
# x$profiles_after<-before
# p<-bud_match[[2]]
# # year_begin<-bud_match[[1]]$year_begin
# 
# data.select<-p$data %>% filter(., year(DateTime)>=year_begin) %>%
#   select(., Key, Name, DateTime, lon, lat, Wind) %>% st_drop_geometry(.)%>%
#   remove_rownames(.)
# # data.select<-remove_rownames(data.select)
# # p<-vector(mode="list")
# names(p$events)<-data.select$Key
# p$events<-vector(mode="list", length=nrow(data.select))
# for(i in 1:length(p$events)){
# p$events[[i]]$data<-data.select[i,]}
# rm(data.select)
# 
# x<-
# f.closest(p$events){
#   x<-p$events
# }
# 
# 
# temp_list<-lapply(p f.closest)
# 
# 
# 
# # p$events$data<-split(data.select, seq(nrow(data.select)))
# 
# 
# 
# 
# # xbase::findInterval(overpass_time, x$profile_dates )
# 
# before<-x$profile_dates[x$profile_dates<=overpass_time]
# after<-x$profile_dates[x$profile_dates>=overpass_time]
# 
# 
# prior<-which(abs(overpass_time - before) == min(abs(overpass_time - before)))
# before_time<-x$profile_dates[prior]
# # overpass_time
# post<-which(abs(overpass_time - after) == min(abs(overpass_time - after)))+length(before)
# after_time<-x$profile_dates[post]
# 
# 
# ### adjust 'prior' and 'post' selection index
# 
# 
# which(abs(x$profile_dates - overpass_time) == (min(abs(x$profile_dates - overpass_time)) & overpass_time))
# before<-x$profiles_dates[
#   
# 
# before
# overpass_time
# after
# 
# profile<-match$'41272'$profiles[closest]
# 
# window<-match$'41272'$profiles[c((closest-5):(closest+5))]
# print(paste0("overpass: ", overpass_time, " profile: ", profile, " at ", date_closet))
# 
# 








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