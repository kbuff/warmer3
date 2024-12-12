#retrieve SLR scenario. AR6 50th percentile

world=st_read("input/TM_WORLD_BORDERS-0.3.shp")

globalSLR= readRDS('input/ar6_global_SLR_50th_percentile.rds')

N=lat+0.5
S=lat-0.5
E=long+0.5
W=long-0.5

ind=which(globalSLR$lat<=N & globalSLR$lat>=S & globalSLR$long>=W & globalSLR$long<=E)

theLat=globalSLR$lat[ind]
theLong=globalSLR$long[ind]
pts = data.frame(lon=theLong, lat=theLat)
pts = st_as_sf(pts, coords=c('lon','lat'), crs=4326)
sitePt=st_as_sf(data.frame(lon=long, lat=lat), coords=c('lon', 'lat'), crs=4326)

bbox_extent = st_bbox(pts)

theDist= sqrt((lat-theLat)^2 + (long-theLong)^2)
sel=which.min(theDist)

plotSLR=ggplot()+
  geom_sf(data = world) +
  geom_sf(data=pts)+
  geom_sf(data=sitePt, col='red' )+
  geom_sf(data=pts[sel,], col='cyan' )+
  coord_sf(xlim = c(bbox_extent["xmin"], bbox_extent["xmax"]),
           ylim = c(bbox_extent["ymin"], bbox_extent["ymax"])) +
  theme_minimal()

print(plotSLR)

slr = data.frame(Year=seq(2020,2150,10),
                 Low=globalSLR$low[ind[sel],]/1000,
                 IntLow=globalSLR$intlow[ind[sel],]/1000,
                 Int=globalSLR$int[ind[sel],]/1000,
                 IntHigh=  globalSLR$inthigh[ind[sel],]/1000,
                 High=globalSLR$high[ind[sel],]/1000)

#Use a spline to interpolate the SLR projections to generate annual increments of sea level change
slr.low= gam(Low~s(Year), data=slr)
slr.intlow = gam(IntLow~s(Year), data=slr)
slr.int = gam(Int~s(Year), data=slr)
slr.inthigh = gam(IntHigh~s(Year), data=slr)
slr.high = gam(High~s(Year), data=slr)

p1 =predict(slr.low, newdata=data.frame(Year=2020:2150))
p2 =predict(slr.intlow, newdata=data.frame(Year=2020:2150))
p3 =predict(slr.int, newdata=data.frame(Year=2020:2150))
p4 =predict(slr.inthigh, newdata=data.frame(Year=2020:2150))
p5 =predict(slr.high, newdata=data.frame(Year=2020:2150))

slrDF = data.frame(Year=2020:2150, Low=p1, IntLow=p2, Int=p3, IntHigh=p4, High=p5)


