# functions for ch2 scripts
# need to proff read an clen this up. make function mor generic

#### general functions ####

f.ipak <- function(pkg){
  # loads packages, quietly, given by a vector of package names e.g., pkg<-c("ggplot", "tidyverse")
  # will install  packages listed , and their dependencies, if needed.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, quiet=T, verbose = F)
  sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts=F)
}


f.select_tc<-function(x, key_id){
  x<-x$tc
  x<-x[which(names(x) == key_id)]
  return(x)
}


f.ipak <- function(pkg){
  # loads packages, quietly, given by a vector of package names e.g., pkg<-c("ggplot", "tidyverse")
  # will install  packages listed , and their dependencies, if needed.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, quiet=T, verbose = F)
  sapply(pkg, require, character.only = TRUE, quietly = FALSE, warn.conflicts=F)
}


### script specific functions ####

f.single_storm_h.pts<-function(h.pts, key_id){
  tc.pts<-filter(h.pts, Key == key_id )
  conv<-1.94384 # convert to m/s
  tc.pts$max_sustained_wind_m_s<-tc.pts$Wind/conv
  return(tc.pts)
}

f.select_tc_specific_profiles<-function(list.x, key_id){
  # list. x comes from data$tc_search_radi_profiles$km_x wehr x is 100 or 200 or 500
  tc.p<-f.select_tc(x=list.x, key_id) %>%.[[1]] %>%st_as_sf()
  e<-extent(tc.p)%>%as.vector(.)
  q<-select_profiles(lon_lim=c(e[1],e[2]), lat_lim=c(e[3], e[4]), start_date = as.character(min(tc.p$time)-1), end_date =as.character(max(tc.p$time)+1),  outside="none", sensor=NULL)
}


f.reformat_profiles_for_ONE_Argo<-function(tc.p){
  # tc.p -s one lement in 
  s<-vector(mode="list")
  s$float_ids<-unique(tc.p$float)
  s$float_profs<-vector(mode="list", length=length(s$float_ids))
  for(i in 1:length(s$float_profs)){
    f<-s$float_ids[i] 
    x<-filter(tc.p, float ==f)
    s$float_profs[[i]]<-as.integer(x$profile)
  }
  return(s)}


# f.one_argo_to_oce<-function(x){
#   
#   data.argo<-load_float_data(float_ids =x$float_ids, float_profs = x$float_profs)
#   
# }