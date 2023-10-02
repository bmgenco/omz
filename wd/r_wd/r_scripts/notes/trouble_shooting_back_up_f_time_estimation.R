####  testing ####

z<-all_storms_500km
z<-Map(cbind, z, Key= names(z))
z<-lapply(z, st_as_sf)
z<-z[[1]]
z.key<-unique(z$Key) 
y<-filter(h.pts, Key ==z.key)
w<-vector(mode = 'list', length=dim(z)[1])

test<-f.per_storm_t(z=z, h.pts=h.pts)






)



i<-5

rm(pt, pld, fld, pp, fp, np, end, first, w ,z, y2, i, pl, fl, direction, y)

##### function v1  ####

f.time_relative_overpass<-function(z, h.pts){
  ### internal functions ####
  f.line<-function(x, np) {
    l<-st_as_sf(rbind(np,x))
    l$multipoints<-st_combine(l$geometry)
    st_geometry(l)<-NULL
    l<-dplyr::select(l, DateTime,multipoints)
    l$geometry<-l$multipoints
    l<-select(l, -multipoints)
    l<-l[1,]
    st_geometry(l)<-l$geometry
    l<-st_cast(l, "LINESTRING")
    return(l)
  }
  f.overpass_time<-function(np, fp, pp, direction) {
    if(direction == "p") { 
      v<-(st_distance(np$geometry, pp$geometry))/as.numeric((difftime(np$DateTime, pp$DateTime, units="secs")))
      d<-st_distance(pt$overpass_point, np$geometry) 
      tt<-np$DateTime-(as.numeric(d)/as.numeric(v))
    } else if (direction == "f") {
      v<-(st_distance(np$geometry, fp$geometry))/as.numeric((difftime(np$DateTime, fp$DateTime, units="secs")))
      d<-st_distance(pt$overpass_point, np$geometry) 
      tt<-np$DateTime+(as.numeric(d)/as.numeric(v))
    }
    return(tt)
  }
  f.time_estimation<-function(y,z){
    
    z<-na.omit(z)
    w<-vector(mode = 'list', length=dim(z)[1])
    for(i in 1:length(w)) {
      pt<-z[i,]
      y$distance_m<-st_distance(y$geometry,pt$geometry) %>% as.vector(.)
      np<-filter(y, distance_m == min(y$distance_m))
      first<-y %>% .[1,]
      end<-y %>% .[length(y),]
      
      # np<-filter(y, distance_m == min(y$distance_m)) %>%dplyr::select(.,  DateTime, geometry, distance_m)
      # first<- y%>%dplyr::select(.,  DateTime, geometry, distance_m) %>% .[1,]
      # end<- y%>%dplyr::select(.,  DateTime, geometry, distance_m) %>% .[length(y),]
      
      if(is.null(np$Date)==T) {
        pt$overpass_time<-NA} else if (isTRUE(np$DateTime) && np$DateTime  == end$DateTime) {pt$overpass_time<-np$DateTime
        } else if (isTRUE(np$DateTime) && np$DateTime == first$DateTime) {
          pt$overpass_time<-np$DateTime} else {
            
            # pp<-y[(which(y$DateTime==np$DateTime) -1),] %>% dplyr::select(., DateTime, geometry, distance_m)
            # fp<-y[(which(y$DateTime==np$DateTime) +1),] %>% dplyr::select(., DateTime, geometry, distance_m)
            
            pp<-y[(which(y$DateTime==np$DateTime) -1),]
            fp<-y[(which(y$DateTime==np$DateTime) +1),]
            
            pl<-f.line(x=pp, np)
            fl<-f.line(x=fp, np)
            
            fld<-st_distance(pt$overpass_point, fl)
            pld<-st_distance(pt$overpass_point, pl)  
            
            
            if(fld < pld){
              direction <- "f"
            } else if (pld < fld){
              direction <-"p"}
            
            pt$overpass_time<-f.overpass_time(np=np, fp=fp, pp=pp, direction=direction) 
            if(is.na(pt$overpass_time) =="TRUE") {
              pt$day_difference<-NA} else {pt$day_difference<-as.numeric(date(pt$time)-date(pt$overpass_time))}
            pt<-dplyr::select(pt, -Key)   
            
            # if(isTRUE(pt)){
            #  w[[i]]<-pt 
            # } else {w[[i]]<-NULL}
            
            w[[i]]<-pt
          }
      
      
    }
    
    w<-as.data.frame(do.call(rbind, w))
    return(w)} 
  f.per_storm_t<-function(z, h.pts){
    z.key<-unique(z$Key) 
    y<-filter(h.pts, Key ==z.key) 
    z<-f.time_estimation(z=z,y=y)
    return(z)
  }  
  
  ### master function ####
  z<-Map(cbind, z, Key= names(z))
  z<-lapply(z, f.per_storm_t, h.pts = h.pts)
  z<-lapply(z, na.omit)
  return(z)
}

##### function v2 ####

f.time_relative_overpass<-function(z, h.pts){
  ### internal functions ####
  
  # f.custom_min <- function(x) {if (length(x)>0) min(x) else Inf}  
  f.line<-function(x, np) {
    x<-dplyr::select(x, DateTime, geometry)
    np<-dplyr::select(np, DateTime, geometry)
    l<-st_as_sf(rbind(np,x))
    l$multipoints<-st_combine(l$geometry)
    st_geometry(l)<-NULL
    l<-dplyr::select(l, DateTime, multipoints)
    l$geometry<-l$multipoints
    l<-select(l, -multipoints)
    l<-l[1,]
    st_geometry(l)<-l$geometry
    l<-st_cast(l, "LINESTRING")
    return(l)
  }
  f.overpass_time<-function(np, fp, pp, direction) {
    if(direction == "p") { 
      v<-(st_distance(np$geometry, pp$geometry))/(as.numeric((difftime(np$DateTime, pp$DateTime, units="secs"))))
      d<-st_distance(pt$overpass_point, np$geometry) 
      tt<-np$DateTime-(as.numeric(d)/(as.numeric(v)))
    } else if (direction == "f") {
      v<-(st_distance(np$geometry, fp$geometry))/as.numeric((difftime(np$DateTime, fp$DateTime, units="secs")))
      d<-st_distance(pt$overpass_point, np$geometry) 
      tt<-np$DateTime+(as.numeric(d)/as.numeric(v))
    }
    return(ymd_hms(tt))
  }
  f.time_estimation<-function(y,z){
    w<-vector(mode = 'list', length=dim(z)[1])
    for(i in 1:length(w)){
      pt<-z[i,]
      pt$overpass_time<-NA
      y2<-y
      y2$distance_m<-st_distance(y2$geometry,pt$overpass_point) %>% as.vector(.)
      
      np<-filter(y2, distance_m == min(y2$distance_m))
      first<-y2 %>% .[1,]
      end<-y2 %>% .[dim(y2)[1],]
      
      # np<-filter(y, distance_m == min(y$distance_m)) %>%dplyr::select(.,  DateTime, geometry, distance_m)
      # first<- y%>%dplyr::select(.,  DateTime, geometry, distance_m) %>% .[1,]
      # end<- y%>%dplyr::select(.,  DateTime, geometry, distance_m) %>% .[length(y),]
      
      # if(is.null(np$Date)==T) {
      #   pt$overpass_time<-NA} else if (isTRUE(np$DateTime) && np$DateTime  == end$DateTime) {pt$overpass_time<-np$DateTime
      #   } else if (isTRUE(np$DateTime) && np$DateTime == first$DateTime) {
      #     pt$overpass_time<-np$DateTime} else {
      
      
      
      # if (np$DateTime == end$DateTime) {pt$overpass_time<-np$DateTime
      #     } else if (np$DateTime == first$DateTime) {
      #       pt$overpass_time<-np$DateTime} else {
      
      pp<-y2[(which(y2$DateTime==np$DateTime) -1),]
      fp<-y2[(which(y2$DateTime==np$DateTime) +1),]
      
      pl<-f.line(x = pp, np)
      fl<-f.line(x = fp, np)
      
      fld<-st_distance(pt$overpass_point, fl)
      pld<-st_distance(pt$overpass_point, pl)  
      
      
      if(fld < pld){
        direction <- "f"
      } else if (pld < fld){
        direction <-"p"}
      
      pt$overpass_time<-f.overpass_time(np=np, fp=fp, pp=pp, direction=direction)
      # }
      # rm(pp, fp, pp, fp, pl, fl, fld, pld, np)} 
      
      pt$day_difference<-as.numeric(date(pt$time)-date(pt$overpass_time))    
      
      
      w[[i]]<-pt
      #rm(pt, y2)
    }
    
    w<-as.data.frame(do.call(rbind, w))
    return(w)} 
  f.per_storm_t<-function(z, h.pts){
    z.key<-unique(z$Key) 
    y<-filter(h.pts, Key ==z.key) 
    z<-f.time_estimation(z=z,y=y)
    return(z)
  }  
  
  ### master function ####
  z<-Map(cbind, z, Key= names(z))
  z<-lapply(z, st_as_sf)
  z<-lapply(z, f.per_storm_t, h.pts = h.pts)
  z<-lapply(z, na.omit)
  return(z)
}


#### bah ####