quartz()
CentralCoast<- get_map(location=c(lon=-122.42725,lat=37.2425),zoom=7,maptype='satellite')

ggmap(CentralCoast, fullpage = TRUE)



violent_crimes$offense <- factor(violent_crimes$offense,+   levels = c("robbery", "aggravated assault", "rape", "murder"))



ggplot(R28meansd, aes(treat, meanbiov, shape=factor(disp))) + 
  geom_point(size=5,position=pd) + 
  geom_errorbar(aes(ymin=meanbiov-sdbiov,ymax=meanbiov+sdbiov), 
                colour="black",width=.1,position=pd,size=.7) 
opts(aspect.ratio=1.4) + 
  theme_bw() + 
  opts(axis.title.y=theme_text(size=16, face="bold",
                               colour="black",angle=90,vjust=4.3)) +
  ylab(expression(paste(log~(Biovolume~µm^3~L^-1))))

get_googlemap(urlonly = TRUE)

# get_googlemap has several argument checks
get_googlemap(zoom = 13.5)
get_googlemap(scale = 3)
get_googlemap(center = c(-30,-110))

# markers and paths are easy to access
d <- function(x=-95.36, y=29.76, n,r,a){
  round(data.frame(
    lon = jitter(rep(x,n), amount = a),
    lat = jitter(rep(y,n), amount = a)
  ), digits = r)
}
df <- d(n=50,r=3,a=.3)
map <- get_googlemap(markers = df, path = df,, scale = 2)
ggmap(map)
ggmap(map, fullpage = TRUE) +
  geom_point(aes(x = lon, y = lat), data = df, size = 3, colour = 'black') +
  geom_path(aes(x = lon, y = lat), data = df)

gc <- geocode('waco, texas')
center <- as.numeric(gc)
ggmap(get_googlemap(center = center, color = 'bw', scale = 2), fullpage = T)

# the scale argument can be seen in the following
# (make your graphics device as large as possible)
ggmap(get_googlemap(center, scale = 1), fullpage = TRUE) # pixelated
ggmap(get_googlemap(center, scale = 2), fullpage = TRUE) # fine

