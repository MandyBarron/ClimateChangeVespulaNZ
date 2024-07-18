setwd("C:/.../Climate Change Pests/")
library(ncdf4)
library(terra)
library(RColorBrewer)

# # ncin <-nc_open("NIWAClimateChangeRCPScenarios_netcdf/Mean temp/NIWA_MeanTemp_RCPpast_1986-2005_MAM.nc")
# # ncin<-nc_open("NIWAClimateChangeRCPScenarios_netcdf/Mean temp/NIWA_TempMean_2040_RCP4_5_MAM.nc")
# # print(ncin)
# # dname<-"air_temperature"
# # dname<-"tmean"
# # lon<-ncvar_get(ncin,"longitude")
# # #dim(lon)
# # lat<-ncvar_get(ncin,"latitude")
# # #dim(lat)
# # time<- ncvar_get(ncin,"time")
# # # refDate<-as.Date("1971-01-01 00:00:00")
# # # tDates<-refDate+time
# # tunits <- ncatt_get(ncin,"time","units")
# # #dim(time)
# # 
# # # get temperature
# # tmp_array <- ncvar_get(ncin,dname)
# # dlname <- ncatt_get(ncin,dname,"long_name")
# # dunits <- ncatt_get(ncin,dname,"units")
# # fillvalue <- ncatt_get(ncin,dname,"_FillValue")
# # dim(tmp_array)
# # #tmp_array[tmp_array==fillvalue$value] <- NA
# # #print(length(which(is.na(tmp_array)==F))/ncin$dim$time$len)
# # print(length(which(is.na(tmp_array)==F)))
# # 
# # 
# # # image(lon,lat,tmp_array, col=rev(brewer.pal(10,"RdBu")))  #doesn't work cos dec lats
# # # m <- 1
# # # tmp_slice <- tmp_array[,,m]
# # # r<-rast(tmp_slice)
# # # image(lon,lat,tmp_slice, col=rev(brewer.pal(10,"RdBu")))  #doesn't work cos dec lats
# # # head(lon)
# # # head(lat)
# # 
# # 
# # #some you can just import direct as raster i.e. RCPPast, others (projections) don't work
# # r<-rast("NIWAClimateChangeRCPScenarios_netcdf/Mean temp/NIWA_MeanTemp_RCPpast_1986-2005_MAM.nc")
# # plot(r)
# # r1<-rast("NIWAClimateChangeRCPScenarios_netcdf/Mean temp/NIWA_TempMean_2040_RCP4_5_MAM.nc")
# # plot(r1)
# # 
# # r1 = rast(t(tmp_array))
# # ext(r1)=c(166.4, 178.55, -47.35, -34.35)
# # crs(r1) = "epsg:4326"
# # plot(r1)
# # #first need to change extent for first raster
# # ext(r)<-ext(r1)
# # test<-c(r,r1)
# 
# #first need to change extent for first raster and resave netcdf file
# ncin <-nc_open("NIWAClimateChangeRCPScenarios_netcdf/Mean temp/NIWA_MeanTemp_RCPpast_1986-2005_MAM.nc")
# print(ncin)
# dname<-"air_temperature"
# lon<-ncvar_get(ncin,"longitude")
# lat<-ncvar_get(ncin,"latitude")
# time<- ncvar_get(ncin,"time")
# tunits <- ncatt_get(ncin,"time","units")
# tmp_array <- ncvar_get(ncin,dname)
# tmp_array[143:155,24:35]
# # tempf = rast(t(tmp_array))
# # plot(tempf)
# dlname <- ncatt_get(ncin,dname,"long_name")
# dunits <- ncatt_get(ncin,dname,"units")
# fillvalue <- ncatt_get(ncin,dname,"_FillValue")
# dim(tmp_array)
##this is to align RCPast with future projections
# lon<-lon+0.05
# x <- ncdim_def( "longitude", "degrees_east", lon)
# y <- ncdim_def( "latitude", "degrees_north", lat)
# t <- ncdim_def( "time", "days since 1986-1-1", time)
# # Make a variable with those dimensions.  Note order: time is LAST
# meanT <- ncvar_def(dname,    dunits$value,  list(x,y,t), missval = fillvalue$value )
# # Create a netCDF file with this variable
# ncnew <- nc_create( "NIWAClimateChangeRCPScenarios_netcdf/Mean temp/NIWA_MeanTemp_RCPpast_1986-2005_2_MAM.nc", meanT)
# ncvar_put(ncnew, dname, tmp_array)
# nc_close(ncnew)
# 
# 
# folderPath<-"NIWAClimateChangeRCPScenarios_netcdf/Mean temp/NIWA_"
# fileExt<-"MAM.nc"
# varList<-c("MeanTemp_RCPpast_1986-2005_2_","TempMean_2040_RCP4_5_","TempMean_2090_RCP4_5_","TempMean_2040_RCP8_5_","TempMean_2090_RCP8_5_")
# 
# for( i in 1:length(varList)) {
#   fileNom<-paste(folderPath,varList[i],fileExt,sep="")
#   ncin <-nc_open(fileNom)
#   lon<-ncvar_get(ncin,"longitude")
#   lat<-ncvar_get(ncin,"latitude")
#   dname<-names(ncin$var)
#   tmp_array <- ncvar_get(ncin,dname)
#   nc_close(ncin)
#   tempf = rast(t(tmp_array))
#   ext(tempf)=round(c(min(lon), max(lon), min(lat), max(lat)),digits=3)
#   crs(tempf) = "epsg:4326"
#   names(tempf)<-varList[i]
#   if (i==1){
#     r<-tempf
#   } else {
#     r<-c(r,tempf)
#   }
# 
# }
# writeRaster(r,"meanMAMtemps.tif", filetype="GTiff")

# r2<-rast("meanMAMtemps.tif")
# plot(subset(r2, 2:5))
# 
# #want to add differences from projected climates to RCPpast
# r2[24:35,143:144,1:2]
# r3<-r2
# 
# for (i in 2:5) {
#   r3[[i]]<-r3[[1]]+r2[[i]]
# }
# #sanity chk
# r3[24:35,143:144,1:2]
# plot(subset(r3, 2:5))
# writeRaster(r3,"meanMAMtempsSum.tif", filetype="GTiff",overwrite=TRUE)

# ppnFn<-function(tempDev,interc=-2,slope=0.5) {
#   lp<-interc + (tempDev*slope)
#   ppn<-1/(1+exp(-lp))
#   return(ppn)
# }
# tempDev<-seq(0.5,3.8,by=0.05)
# plot(tempDev,ppnFn(tempDev))

r3<-rast("meanMAMtempsSum.tif")

ppnFn<-function(temp,interc=-5.98,slope=0.264) {
  lp<-interc + (temp*slope)
  ppn<-1/(1+exp(-lp))
  return(ppn)
}
temp<-seq(0,20,by=0.5)
plot(temp,ppnFn(temp))

newNames<-c("Historic (1986-2005)", "RCP 4.5, mid-century (2040)","RCP 4.5, end century (2090)",
            "RCP 8.5, mid-century (2040)","RCP 8.5, end century (2090)")
propPeren<-r3
for( i in 1:5) {
  propPeren[[i]]<-ppnFn(r3[[i]])
  names(propPeren[[i]])<-newNames[i]
}
plot(propPeren)
writeRaster(propPeren,"VespGermPropPeren.tif", filetype="GTiff",overwrite=TRUE)

#colfunc <- colorRampPalette(c("#31a354","#ffeda0","#feb24c","#f03b20"))

plg = list(
  title = expression(atop("Propn. nests","perennial")),
  title.cex = 0.6,
  cex = 0.6,
  shrink=0
)
#plot(subset(propPeren,2:5),range=c(0,0.35),col=colfunc(20),plg=plg)
colfunc <- colorRampPalette(c("#ffeda0","#feb24c","#f03b20"))
b <- seq(0, 0.36, 0.01)
cols<-c("white","aliceblue","antiquewhite",colfunc(34))

plot(subset(propPeren,1),range=c(0,0.36),col=cols,
     type = 'continuous',plg=plg,main=names(propPeren)[1],legend=T)
plot(subset(propPeren,2:5),range=c(0,0.36),col=cols,
     type = 'continuous',plg=plg,legend=F)


#redo process for precipitation
#first need to change extent for first raster and resave netcdf file
# ncin <-nc_open("NIWAClimateChangeRCPScenarios_netcdf/Precipitation/NIWA_Rain_RCPpast_1986-2005_SON.nc")
# print(ncin)
# dname<-"rain"
# lon<-ncvar_get(ncin,"longitude")
# lat<-ncvar_get(ncin,"latitude")
# time<- ncvar_get(ncin,"time")
# tunits <- ncatt_get(ncin,"time","units")
# tmp_array <- ncvar_get(ncin,dname)
# #tmp_array[143:155,24:35]
# dlname <- ncatt_get(ncin,dname,"long_name")
# dunits <- ncatt_get(ncin,dname,"units")
# fillvalue <- ncatt_get(ncin,dname,"_FillValue")
# dim(tmp_array)
# lon<-lon+0.05
# x <- ncdim_def( "longitude", "degrees_east", lon)
# y <- ncdim_def( "latitude", "degrees_north", lat)
# t <- ncdim_def( "time", "days since 1986-1-1", time)
# # Make a variable with those dimensions.  Note order: time is LAST
# sumRain <- ncvar_def(dname,    dunits$value,  list(x,y,t), missval = fillvalue$value )
# # Create a netCDF file with this variable
# ncnew <- nc_create( "NIWAClimateChangeRCPScenarios_netcdf/Precipitation/NIWA_Rain_RCPpast_1986-2005_2_SON.nc", sumRain)
# ncvar_put(ncnew, dname, tmp_array)
# nc_close(ncnew)
# 
# 
# folderPath<-"NIWAClimateChangeRCPScenarios_netcdf/Precipitation/NIWA_"
# fileExt<-"SON.nc"
# varList<-c("Rain_RCPpast_1986-2005_2_","RainPercChange_2040_RCP4_5_","RainPercChange_2090_RCP4_5_","RainPercChange_2040_RCP8_5_","RainPercChange_2090_RCP8_5_")
# 
# for( i in 1:length(varList)) {
#   fileNom<-paste(folderPath,varList[i],fileExt,sep="")
#   ncin <-nc_open(fileNom)
#   lon<-ncvar_get(ncin,"longitude")
#   lat<-ncvar_get(ncin,"latitude")
#   dname<-names(ncin$var)[2]
#   tmp_array <- ncvar_get(ncin,dname)
#   nc_close(ncin)
#   tempf = rast(t(tmp_array))
#   ext(tempf)=round(c(min(lon), max(lon), min(lat), max(lat)),digits=3)
#   crs(tempf) = "epsg:4326"
#   names(tempf)<-varList[i]
#   if (i==1){
#     r<-tempf
#   } else {
#     r<-c(r,tempf)
#   }
# }
# writeRaster(r,"meanSONrainPer.tif", filetype="GTiff")

# r<-rast("meanSONrainPer.tif")
# plot(r)
# r1<-r
# for( i in 2:5) {
#   r1[[i]]<-r1[[1]]*(r1[[i]]/100)
# #   names(r1[[i]])<-names(r[[i]])
# }
# plot(r)
# plot(r1)
# names(r1)<-c("Rain_RCPpast_1986-2005_2_","RainMMchange_2040_RCP4_5_",
#              "RainMMchange_2090_RCP4_5_", "RainMMchange_2040_RCP8_5_",
#              "RainMMchange_2090_RCP8_5_")
# writeRaster(r1,"meanSONrainDiff.tif", filetype="GTiff",overwrite=TRUE)

r1<-rast("meanSONrainDiff.tif")
plot(subset(r1,2:5))

# rRT<-r1
# for (i in 2:5) {
#      rRT[[i]]<-r1[[1]]+r1[[i]]
# }
# plot(rRT)
# names(rRT)<-c("Rain_SON_RCPpast_1986-2005_2_","Rain_SON_2040_RCP4_5_",
#              "Rain_SON_2090_RCP4_5_", "Rain_SON_2040_RCP8_5_",
#              "Rain_SON_2090_RCP8_5_")
# writeRaster(rRT,"meanSONrain.tif", filetype="GTiff",overwrite=TRUE)

rRT<-rast("meanSONrain.tif")
plot(rRT)

# rainFn<-function(rainDev,interc=2,slope=-0.008) {
#   lp<-interc + (rainDev*slope)
#   surv<-1/(1+exp(-lp))
#   return(surv)
# }
# rainDev<-seq(-900,900,by=50)
# plot(rainDev,rainFn(rainDev))

rainFn<-function(rain,slope=0.0006518) {
  lp<-(rain*slope)
  surv<-exp(-lp)
  return(surv)
}
rain<-seq(0,4000,by=50)
predSurv<-rainFn(rain)
plot(rain,predSurv,type="l")
#cbind(rain,predSurv)


#need to import K map based on basic ecosystems, has been recast to same extent
#and pixel size as RCP projs - this has been done in RasterFiddling2.R
rK<-rast("VespVulgKmap.tif")
plot(rK)

rm<-1.22 #1.22 me, 1.31 BBB, ln(1.5-2.5) Lester etal 2017 yet also have lambda560*0.02=2.416
newNames<-c("Historic (1986-2005)", "RCP 4.5, mid-century (2040)","RCP 4.5, end century (2090)",
            "RCP 8.5, mid-century (2040)","RCP 8.5, end century (2090)")
projK<-rRT
# projK[[1]]<-rK*1.0
# names(projK)[1]<-"Historic"
for( i in 1:5) {
  projK[[i]]<-rK*max((1+(log(rainFn(rRT[[i]]))/rm)),0)
  #projK[[i]]<-max((1+(log(rainFn(rRT[[i]]))/rm)),0)
  names(projK[[i]])<-newNames[i]
}
plot(projK)

#writeRaster(projK,"VespVulgProjKmaps.tif", filetype="GTiff",overwrite=TRUE)

projK<-rast("VespVulgProjKmaps.tif")

# breakpoints <- c(0, seq(0.1, 13, 0.1))
# colors <- c("white",brewer.pal(130, "YlOrRd"))
# plot(r, breaks = breakpoints, col = colors)
#pal <- leaflet::colorNumeric(palette = "RdBu", domain=c(-ceil, ceil), reverse = T)
colfunc <- colorRampPalette(c("blanchedalmond","#ffeda0","#feb24c","#f03b20"))
b <- seq(0, 13, 0.5)
#cols<-c("white",colfunc(13))
cols<-c("darkolivegreen3",colfunc(26))

plg = list(
  title = expression(atop("Autumn nest","density (/ha)")),
  title.cex = 0.6,
  cex = 0.6#,
  #shrink=0
)
# plot(subset(projK,1),range=c(0,13),col=cols, breaks=b,
#      type = 'interval',plg=plg,main=names(projK)[1],legend=T)
plot(subset(projK,1),range=c(0,13),col=cols, 
     type = 'continuous',plg=plg,main=names(projK)[1],legend=T)
plot(subset(projK,2:5),range=c(0,13),col=cols,
     type = 'continuous',plg=plg,legend=F)


colfunc <- colorRampPalette(c("#31a354","#ffeda0","#feb24c","#f03b20"))

plot(subset(projK,2),range=c(0,13),col=colfunc(20),
     legend.only=TRUE,
     legend.args=list(text='value', side=3, font=1, line=2.5, cex=0.8))




plg = list(loc="right",
          #ext=c(180, 182, -47, -37),
           #title= expression(atop("Autumn nest", "density (/ha)")),
           title = expression("Autumn nest density (/ha)"),
           title.cex = 0.8, cex = 0.6, shrink=0 )
#loc="right",

#plot(projK,type = 'continuous',range=c(0,13.5),col=colfunc(20),plg=plg)
# 
ext(projK)
plotext<-round(ext(projK),2)
par(mfrow=c(2,2))
for (i in 2:5){
  if (i %% 2 ==0){
    plot(projK[[i]],range=c(0,13), mar=c(2, 2, 2, 2),col=colfunc(20),
         type = 'continuous', legend=FALSE)
  } else {
    plot(projK[[i]],range=c(0,13), mar=c(2, 2, 2, 2),col=colfunc(20),
         type = 'continuous', plg = plg)
  }
}


mar1=c(4,3,1,1)
par(mfrow=c(2,2),mar=mar1,oma=c(0,0,1.5,0))

plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "")
plot(r,legend=F,mar=mar1,ext=c(5.7,6.5,49.4,50.2), axes = FALSE, add = TRUE)
fbb(xlim, ylim, xat, yat)
mtext("Map one",adj=0)

plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "")
plot(r,legend=F,mar=mar1,ext=c(6.1,6.5,49.4,50.2), axes = FALSE, add = TRUE)
fbb(xlim, ylim, xat, yat)
mtext("Map two",adj=0)


plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "")
plot(r,legend=F,mar=mar1,ext=c(5.7,6.3,49.9,50), axes = FALSE, add = TRUE)
fbb(xlim, ylim, xat, yat)
mtext("Map three",adj=0)

plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "")
plot(r,legend=F,mar=mar1,ext=c(5.7,6.5,49.4,50.2), axes = FALSE, add = TRUE)
fbb(xlim, ylim, xat, yat)
mtext("Map four",adj=0)

plot(subset(projK,2:5),range=c(0,13),col=colfunc(20),
     type = 'continuous',legend=F)
plot(subset(projK,2),range=c(0,13),col=colfunc(20),
     type = 'continuous',legend.only=TRUE, plg = plg)
     
plot(subset(projK,2:5),range=c(0,13),col=colfunc(20),
     type = 'continuous',plg=plg,legend=c(F,T,F,T))


plot(subset(projK,2),range=c(0,13),col=colfunc(20),
     legend.only=TRUE,
     legend.args=list(text='value', side=3, font=1, line=2.5, cex=0.8))
plot(subset(projK,2:5),range=c(0,13),col=colfunc(20),plg=plg,legend=c(F,T,F,T))
add_legend("bottom", legend = "centroids", pch = 20, col="red")
