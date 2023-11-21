



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
  if(all(is.na(f$DOXY_ADJUSTED))==T) {pick<-f$DOXY} else {pick<-f$DOXY_ADJUSTED}
  
  
  if(any(min(na.omit(pick))>20)==T) {o<-min(na.omit(pick))} else {o<-min(which(pick <=20))} 
  
  if(is.na(o)==T) {o<-1} else {o<-o}
  
  
  if(all(is.na(f$PRES_ADJUSTED))==T) {p<-f$PRES[o]} else {p<-f$PRES_ADJUSTED[o]}} 

f.get_profiles<-function(x){
  y<-load_float_data(float_ids= x$float, format = 'dataframe')
  y<-filter(y, CYCLE_NUMBER %in% x$profile)
  z<-list(x,y)
  names(z)<-c("tc_info", "profiles")
  return(z)
}

f.pressure_of_omz<-function(z){
  
  
  
  
  y<-z$tc_info
  x<-z$profiles
  c<-unique(x$CYCLE_NUMBER)
  q<-vector(mode="list", length=length(c))  
  
  for(i in 1:length(c)){
    f<-filter(x, CYCLE_NUMBER==c[i])
    p<-f.find_omz(f)
    q[[i]]<-p
    }
  
  y$omz_p<-unlist(q)
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

f.pressure_of_omz<-function(z){
  
  
  
  
  y<-z$tc_info
  x<-z$profiles
  c<-unique(x$CYCLE_NUMBER)
  z<-vector(mode="list", length=length(c))  
  
  for(i in 1:length(c)){
    f<-filter(x, CYCLE_NUMBER==c[i])
    p<-f.find_omz(f)
    z[[i]]<-p}
  
  y$omz_p<-unlist(z)
  return(y)
}


f.final_match<-function(x){
z<-f.get_profiles(x)
y<-f.pressure_of_omz(z)
y<-f.omz_change_stats(y)
return(y)    
}  




#4 change

