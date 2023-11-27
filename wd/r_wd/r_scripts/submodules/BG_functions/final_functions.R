#                                                    to do:                                                   #
# (1) edit f.pressure to get local time of day for profile see:                                               #            
#  https://stackoverflow.com/questions/63257782/converting-utc-time-to-local-time-in-r-using-lat-and-long     #
#                                                                                                             #
###############################################################################################################

#1 per float
f.remove<-function(x){
  
x$s<-droplevels(x$s)
x$s<-as.numeric(x$s)
flts<-unique(x$float)
z<-vector(mode="list", length=length(flts))
for(i in 1:length(flts)){
  z[[i]]<-filter(x, float ==flts[i])
}
for(i in 1:length(z)){
  if(length(unique(z[[i]]$s))<=1) {z[[i]]<-NA} else {z[[i]]<-z[[i]]}
}

z<-z[!sapply(z, function(x) all(is.na(x)))]
return(z)


  
}

# match functions:

f.find_omz<-function(f){
  # if(all(is.na(f$DOXY_ADJUSTED))==T) {pick<-na.omit(f$DOXY)} else {pick<-na.omit(f$DOXY_ADJUSTED)}
  
  if(all(is.na(f$DOXY_ADJUSTED))==T) {pick<-f$DOXY} else {pick<-f$DOXY_ADJUSTED}
  
  # if(is.na(pick)==T) {o<-NA} else if(any(min(na.omit(pick))>20)==T) {o<-min(na.omit(pick))} else {o<-min(which(pick <=20))} 
  
  # if(all(is.na(pick)==T)) {o<-NA} else if(any(min(na.omit(pick))>20)==T) {o<-min(na.omit(pick))} else {o<-which.min(abs(pick-20))}

  # if(all(is.na(pick)==T)) {o<-NA} else {o<-which.min(abs(pick-20))}
  
  # only_pick<-na.omit(pick)
  if(all(is.na(pick)==T)) {o<-NA} else if(min(na.omit(pick)>=20)==T) {o<-which.min(abs(pick-20))} else {o<-which.min(pick[pick<=20])}

  # if(is.na(o)==T) {o<-NA} else {o<-o}
  if(is.na(o)==T) {p<-NA} else if(all(is.na(f$PRES_ADJUSTED))==T) {p<-f$PRES[o]} else {p<-f$PRES_ADJUSTED[o]}
  
  omz_stats<-list()
  
  if(is.na(o)==T) {omz_stats$o<-NA} else {omz_stats$o<-pick[o]}
  if(is.na(p)==T) {omz_stats$p<-NA} else {omz_stats$p<-p}
  
 
  return(omz_stats)
  }


f.get_profiles<-function(x){
  y<-load_float_data(float_ids= x$float, format = 'dataframe')
  y<-filter(y, CYCLE_NUMBER %in% x$profile)
  z<-list(x,y)
  names(z)<-c("tc_info", "profiles")
  return(z)
}



f.pressure_of_omz<-function(z){
  
  # add removal of profile here
  # add min value of moz chossen
  
  
  y<-z$tc_info
  x<-z$profiles
# origin <- lubridate::ymd_hms('1950-01-01 00:00:00')
# if(is.null(x$JULD)==T) {x$date_jd_utc<-NA} else {x$date_jd_utc<-origin + x$JULD * 3600*24}
  

  y$omz_p<-NA
# y$decimal_day_dif<-NA
  y$omz_value<-NA
  
  
  cyc<-unique(x$CYCLE_NUMBER)
 
  for(i in 1:length(cyc)){
   
    f<-filter(x, CYCLE_NUMBER==cyc[i])
# ti<-filter(y, profile==profile[i])
    omz_stats<-f.find_omz(f)
    
    # if(any(is.na(f$date_jd_utc)==T) {d<-ti$day_difference} else {d<-as.numeric(unique(f$date_jd_utc)-ti$overpass_time)}
# if(exists("date_jd_utc", f)==F) {d<-ti$day_difference} else {d<-as.numeric(unique(f$date_jd_utc)-ti$overpass_time)}
    # if(!exists(d)) {y$decimal_day_dif[i]<-NA} else {y$decimal_day_dif[i]<-d}
   
# if(exists('d')==F) {d<-ti$day_difference} else {d<-d}
# y$decimal_day_dif[i]<-d
    y$omz_value[i]<-omz_stats$o
    y$omz_p[i]<-omz_stats$p
    }
  
  
  return(y)
}


f.omz_change_stats<-function(y){

  y$mean<-mean(y$omz_p)
  y$median<-median(y$omz_p)
  y.p<-filter(y, day_difference < 0)
  y$prior_n<-length(y.p$omz_p)
  y$post_n<-length(y$omz_p)-y$prior_n
  y$prior_mean<-mean(y.p$omz_p)
  y$prior_min<-min(y.p$omz_p)
  y$change_prior<-NA
  for(i in 1:dim(y)[1]){
    y$change_prior[i]<-y$prior_mean[i]-y$omz_p[i]}
  y$var<-var(y$omz_p)
  y$max<-max(y$omz_p)
  y$min<-min(y$omz_p)
  y$range<-y$max-y$min
  return(y)  
}

# 
# origin <- lubridate::ymd_hms('1950-01-01 00:00:00')
# date_jd_utc<-origin + y$JULD * 3600*24



# f.pressure_of_omz<-function(z){
#   
#   y<-z$tc_info
#   x<-z$profiles
#   c<-unique(x$CYCLE_NUMBER)
#   z<-vector(mode="list", length=length(c))  
#   
#   for(i in 1:length(c)){
#     f<-filter(x, CYCLE_NUMBER==c[i])
#     p<-f.find_omz(f)
#     
#     
#     z[[i]]<-p}
#   
#   y$omz_p<-unlist(z)
#   return(y)
# }


f.final_match<-function(x){
# z<-f.get_profiles(x)
  
y<-f.pressure_of_omz(x)
y<-f.omz_change_stats(y)

# z<-list(x$profiles, y)
# names(z)<-c("profiles", "omz_change")
# return(z)    
return(y)
rm(x,y,z,o,p,i)
}  








#4 change

