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