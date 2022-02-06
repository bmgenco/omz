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
        
        data.select<-p$data %>% filter(., year(DateTime)>=year_begin) %>%
          select(., Key, Name, DateTime, lon, lat, Wind) %>% st_drop_geometry(.)%>%
          remove_rownames(.)
      
        names(p$events)<-data.select$Key
        p$events<-vector(mode="list", length=nrow(data.select))
        for(i in 1:length(p$events)){
          p$events[[i]]$data<-data.select[i,]}
        rm(data.selects)
        
        p$profiles<-cell_floats$profiles
        p$floats<-cell_floats$floats
        p$profile_dates<-Sprof$date[cell_floats$profiles]
        
        p$profiles_postions ### need to do this for adjusting events data
        p$floats_slection ## tmeproary thing for matching floats to profiles??
        
        # f.window<-function(p){p}
        # p$individual_tracks<-vector()
        # p$match<-f.wind(p)
        p$events<-lapply(p$events, f.closest)
        
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



data.select<-data %>% filter(., year(DateTime)>=(year_begin-1)) %>% 
  select(., Key, Name, DateTime, lon, lat) %>% st_drop_geometry(.)
p$events<-vector(mode="list", length=nrow(data.select))
names(p$events)<-data.select$Key

## testing:


q

x<-p$events[[18]]
function_variables$h.prior<-90
function_variables$h.post<-90
f.closest<-function(x){
  
  
  return(p$events)
 

# here


# how do i match float and fporfiles??

before<-p$profile_dates[p$profile_dates<=x$data$DateTime]
after<-p$profile_dates[p$profile_dates>=x$data$DateTime]

before_window<-before[which(date(before) >= date(x$data$DateTime)-function_variables$h.prior)]
after_window<-after[which(date(after) <= date(x$data$DateTime)+function_variables$h.post)] 

x$prior<-data.frame(matrix(ncol=5, nrow=length(before_window))) # profiles withins user define function_
x$post<-data.frame(matrix(ncol=5, nrow=length(after_window)))

names(x$prior)<-c("float", "profiles", "DateTime", "distance", "time_prior")
names(x$post)<-c("float", "profiles", "DateTime", "distance", "time_after")

x$prior$DateTime<-before_window
x$post$DateTime<-after_window
print(x)


# add if staemnets here

prior<-which(abs(x$data$DateTime - before) == min(abs(x$data$DateTime - before)))
before_time<-p$profile_dates[prior]
# overpass_time
post<-which(abs(x$data$DateTime - after) == min(abs(x$data$DateTime - after)))+length(before)
after_time<-p$profile_dates[post]



x$profiles_before<-before
x$profiles_after<-before
return(x)
}






#temp 

p<-bud_match[[2]]
# year_begin<-bud_match[[1]]$year_begin

data.select<-p$data %>% filter(., year(DateTime)>=year_begin) %>%
  select(., Key, Name, DateTime, lon, lat, Wind) %>% st_drop_geometry(.)%>%
  remove_rownames(.)
# data.select<-remove_rownames(data.select)
# p<-vector(mode="list")
names(p$events)<-data.select$Key
p$events<-vector(mode="list", length=nrow(data.select))
for(i in 1:length(p$events)){
p$events[[i]]$data<-data.select[i,]}
rm(data.select)

x<-
f.closest(p$events){
  x<-p$events
}


temp_list<-lapply(p f.closest)



# p$events$data<-split(data.select, seq(nrow(data.select)))




# xbase::findInterval(overpass_time, x$profile_dates )

before<-x$profile_dates[x$profile_dates<=overpass_time]
after<-x$profile_dates[x$profile_dates>=overpass_time]


prior<-which(abs(overpass_time - before) == min(abs(overpass_time - before)))
before_time<-x$profile_dates[prior]
# overpass_time
post<-which(abs(overpass_time - after) == min(abs(overpass_time - after)))+length(before)
after_time<-x$profile_dates[post]


### adjust 'prior' and 'post' selection index


which(abs(x$profile_dates - overpass_time) == (min(abs(x$profile_dates - overpass_time)) & overpass_time))
before<-x$profiles_dates[
  

before
overpass_time
after

profile<-match$'41272'$profiles[closest]

window<-match$'41272'$profiles[c((closest-5):(closest+5))]
print(paste0("overpass: ", overpass_time, " profile: ", profile, " at ", date_closet))


