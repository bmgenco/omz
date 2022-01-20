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