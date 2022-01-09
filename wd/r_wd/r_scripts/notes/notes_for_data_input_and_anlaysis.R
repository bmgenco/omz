### to edit to include functions
#1 line 245  rm(list=ls()[! ls() %in% c( "robj", "fig", "gis_data", "wd", "data_d", "ctdata")])
#2 line 336 (list=ls()[! ls() %in% c("x", "cent",  "robj", "fig", "gis_data", "wd", "stations", "data_d")])


## header for flourmter


files<-list.files(directory)
data<-vector(mode="list", length=length(files))

for(i in 1:length(files)){
  data[[i]]<-read.delim(file.path(directory, files[i]), sep =",", header = F)  
}

return(data)}

# here

files<-list.files(fdata)
test<-read_lines(file.path(fdata, files[1]), skip=0, n_max=1)  
test<-read.delim(file.path(fdata, files[1]), sep =",", header = F)  



### testing and compariosn, note fixed gps..... #####

ol<-data

rm(raw, data)
nl<-data
dim(c.names)
length(c.names)
dim(nl[[1]])



#here
o<-x
rm(x)

data<-f.read_cruise(corrected_fdata)

f.t<-function(x){
  
  #names(x)<-c.names
  options("digits.secs"=3)
  x$time<-mdy_hms(paste(x[,1], x[,2], tz="UTC"))
  x$time<-with_tz(x$time, tz="US/Mountain")
  x$date<-date(x$time)
  x<-subset(x, select= c(date, time,  V4))
  x<-as_tibble(x)
  names(x)<-c("date", "time", "chlorophyll")
  return(x)
}

test<-lapply(data, f.t)
test<-do.call(rbind.data.frame, test)


head(x)
summary(x$V4)

#4 v4 is concnetrated

ch<-(x$V5-0.12)
ch<-ifelse(ch < 0, 0, ch)
summary(ch)
ch<-ch*15.87







### station fix for pl;ottimng in future steps

#older option


f1<-ggplot(x, aes(x=time)) +
  geom_line(aes(y=chlorophyll * scale.chla), color=fcolor) +
  geom_line(aes(y=celsius), color=tcolor) +
  xlab("Date") +
  ggtitle("Flowthrough: SST and Chlorophyll-a") +
  scale_x_datetime(date_breaks=('2 days'))+
  scale_y_continuous(name ="SST (°C)", sec.axis = sec_axis(~./scale.chla,name="Chlorophyll-a [μg/L]")) + 
  #scale_y_continuous(name ="Chlorophyll-a [μg/L]", sec.axis = sec_axis(~./scale.chla,name="°C")) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.ticks.y.left = element_line(color = tcolor),
        axis.text.y.left = element_text(color = tcolor), 
        axis.title.y.left = element_text(color = tcolor)
  ) +
  
  theme(axis.ticks.y.right = element_line(color = fcolor),
        axis.text.y.right = element_text(color = fcolor), 
        axis.title.y.right = element_text(color = fcolor)
  )


##### Version 2 add stations bars# ####

# creat geom_rect see code block prior

rect_color<-"lightgrey"
#alpha.v<-0.005 
alpha.v <- 0.5
f2<-f1 +  geom_rect(data =x ,  aes(xmin =rects$x1[1], xmax = rects$x2[1] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha = alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[2], xmax = rects$x2[2] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[3], xmax = rects$x2[3] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[4], xmax = rects$x2[4] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[5], xmax = rects$x2[5] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[6], xmax = rects$x2[6] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_text(data=rects, x=centered, label=rects[,2], y =16, angle =90)






f2


rect_color<-"gray20"
alpha.v<-0.005 
alpha.v <- 1
f1 <- ggplot(x, aes(x = time)) + geom_rect(data =x ,  aes(xmin =rects$x1[1], xmax = rects$x2[1] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha = alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[2], xmax = rects$x2[2] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[3], xmax = rects$x2[3] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[4], xmax = rects$x2[4] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[5], xmax = rects$x2[5] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_rect(data =x ,  aes(xmin =rects$x1[6], xmax = rects$x2[6] , ymin = -Inf, ymax = Inf), fill=rect_color, alpha =alpha.v) +
  geom_text(data=rects, x=centered, label=rects[,2], y =16, angle =90)

# no alpha option....

f1 <- ggplot(x, aes(x = time)) + geom_rect(data =x ,  aes(xmin =rects$x1[1], xmax = rects$x2[1] , ymin = -Inf, ymax = Inf), fill=rect_color) +
  geom_rect(data =x ,  aes(xmin =rects$x1[2], xmax = rects$x2[2] , ymin = -Inf, ymax = Inf), fill=rect_color) +
  geom_rect(data =x ,  aes(xmin =rects$x1[3], xmax = rects$x2[3] , ymin = -Inf, ymax = Inf), fill=rect_color) +
  geom_rect(data =x ,  aes(xmin =rects$x1[4], xmax = rects$x2[4] , ymin = -Inf, ymax = Inf), fill=rect_color) +
  geom_rect(data =x ,  aes(xmin =rects$x1[5], xmax = rects$x2[5] , ymin = -Inf, ymax = Inf), fill=rect_color) +
  geom_rect(data =x ,  aes(xmin =rects$x1[6], xmax = rects$x2[6] , ymin = -Inf, ymax = Inf), fill=rect_color) +
  geom_text(data=rects, x=centered, label=rects[,2], y =16, angle =90) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.ticks.y.left = element_line(color = tcolor),
        axis.text.y.left = element_text(color = tcolor), 
        axis.title.y.left = element_text(color = tcolor)
  ) +
  
  theme(axis.ticks.y.right = element_line(color = fcolor),
        axis.text.y.right = element_text(color = fcolor), 
        axis.title.y.right = element_text(color = fcolor)
  )+

  theme(
    panel.background = element_rect(fill = NA),
    panel.ontop = TRUE
  )
  





f2<-f1 +
  geom_line(aes(y=chlorophyll * scale.chla), color=fcolor) +
  geom_line(aes(y=celsius), color=tcolor) +
  xlab("Date") +
  ggtitle("Flowthrough: SST and Chlorophyll-a") +
  scale_x_datetime(date_breaks=('2 days'))+
  scale_y_continuous(name ="SST (°C)", sec.axis = sec_axis(~./scale.chla,name="Chlorophyll-a [μg/L]"))+
  #scale_y_continuous(name ="Chlorophyll-a [μg/L]", sec.axis = sec_axis(~./scale.chla,name="°C")) + 
  theme(
    panel.background = element_rect(fill = NA),
    panel.ontop = FALSE
  )

f2


base_plot +  
  scale_x_date(breaks = function(x) seq.Date(from = min(x), 
                                             to = max(x), 
                                             by = "12 years"),
               minor_breaks = function(x) seq.Date(from = min(x), 
                                                   to = max(x), 
                                                   by = "2 years")) +
  ggtitle("(minor_)breaks = custom function")
