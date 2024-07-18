setwd("C:/.../Climate Change Pests/")

library(dismo)
library(terra)
library(geodata)
library(rJava)
#library(rgbif)
library(readr)
library(raster)
library(maps)
# library(compGeometeR) 
# library(MASS) #

#GBIF.org (10 April 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.cprbur
# POLYGON((-12.32288 22.15805,89.61571 22.15805,89.61571 71.681,-12.32288 71.681,-12.32288 22.15805))
# locs<-read.csv("VespGerm_locs/VespGermPresNativeRange.csv",header=T)
# rgbif::occ_download_import(key=gbif_download_key,path=path_to_download)
# locs<-read_tsv("VespGerm_locs/VespGermPresNativeRange.csv",quote="")
# head(locs)
# dim(locs)
# length(which(locs$occurrenceStatus=="PRESENT"))
# min(locs$decimalLongitude)
# max(locs$decimalLongitude)
# min(locs$decimalLatitude)
# max(locs$decimalLatitude)

#clip to native range
xmin <- -10
xmax <- 50
ymin <- 27.5
ymax <- 67.3

#filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))
# VgGeo<-subset(locs, select=c("species","decimalLongitude","decimalLatitude"))
# VgGeo<-VgGeo[which(complete.cases(VgGeo)==T),]
# rm(locs)
# fIn<-which(VgGeo$decimalLongitude>xmin&VgGeo$decimalLongitude<xmax&VgGeo$decimalLatitude>ymin&VgGeo$decimalLatitude<ymax)
# VgGeo<-VgGeo[fIn,]
# #removing dupes takes out 40% records
# VgGeo <- VgGeo[!duplicated(VgGeo),]
# nrow(VgGeo)
# head(VgGeo)
# vgg<-as.data.frame(VgGeo)
# colnames(vgg)<-c("species","longitude","latitude")
# rm(VgGeo)


# ##worldclim data processing
# predictor.files <- list.files(path=paste(getwd(), '/Env_Data/wc2-5', sep=''),
#                               pattern='bil', full.names=TRUE)
# predictor.files<-predictor.files[-1]
# bioclim_data<-rast(predictor.files)
# plot(bioclim_data[[10]])
# points(cbind(vgg$longitude, vgg$latitude), col="red", cex=0.6, pch=20)
# 
# 
# e <- extent(xmin, xmax, ymin, ymax)
# envcrop <- crop(bioclim_data, e)
# plot(envcrop)
# plot(envcrop[[17:19]])
# names(envcrop)
# 
# tempLyrs<-c("bio1","bio2","bio5","bio6","bio7","bio8","bio9","bio10", "bio11")
# tempLyrsI<-which(names(envcrop) %in% tempLyrs)
# for (lyr in tempLyrsI) {
#   envcrop[[lyr]] <- envcrop[[lyr]] / 10
# }
# 
# par(mfrow=c(2,2))
# for (lyr in 1:19) {
#   plot(envcrop[[lyr]],main=names(envcrop)[lyr])
#   #points(cbind(vgg$longitude, vgg$latitude), col="red", cex=0.2, pch=20)
#   points(cbind(presence$longitude, presence$latitude), col="red", cex=0.2, pch=20)
# }
# 
# wtdLyrs<-c("bio1","bio2","bio3","bio5","bio6","bio8","bio11","bio12","bio16", "bio18")
# wtdLyrsI<-which(names(envcrop) %in% wtdLyrs)
# envcrop<-subset(envcrop,wtdLyrsI)
# writeRaster(envcrop,"C:/.../Climate Change Pests/Env_Data/EuropeData/EuropeBioClims.tif", filetype = 'GTiff',overwrite = T)
# 

####if already done clim processing and occurence cleaning
envcrop<-rast("C:/.../Climate Change Pests/Env_Data/EuropeData/EuropeBioClims.tif")
names(envcrop)
#plot(envcrop)

#set up training and testing data
# presence<-vgg[,-1]
# dim(presence)
# rm(vgg)
# ##thin occurences so only one per cell
# presence<-gridSample(presence,bioclim_data[[1]],n=1)
## fIn<-which(presence$longitude>xmin&presence$longitude<xmax&presence$latitude>ymin&presence$latitude<ymax)
## presence<-presence[fIn,]
#write.csv(presence, "VespGerm_Locs/VespGermPresNativeRange.csv", row.names = FALSE)

presence <- read.csv(file = "VespGerm_Locs/VespGermPresNativeRange.csv")

dim(presence)
nobs<-nrow(presence)

#check for co-linearity in pred varss
# all_points <- rbind(pres_train,pres_test,backg_train,backg_test)
# bioclim_extract <- extract(x = envcrop,
#                            y = all_points,
#                            ID = FALSE) # No need for an ID column
# panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
# {
#   par(usr = c(0, 1, 0, 1))
#   r <- abs(cor(x, y,use="pairwise.complete.obs"))
#   txt <- format(c(r, 0.123456789), digits = digits)[1]
#   txt <- paste0(prefix, txt)
#   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#   text(0.5, 0.5, txt, cex = cex.cor * r)
# }
# pairs(bioclim_extract,  lower.panel = panel.smooth, upper.panel = panel.cor,
#       gap=0, row1attop=FALSE)


# exclLyr<-which(names(envcrop)=="bio11")
# envcrop<-subset(envcrop,exclLyr,negate=TRUE)

#select fewer pred vars
#exclLyrs<-c("bio11","bio3","bio8","bio18")
#exclLyrs<-c("bio11","bio3","bio8","bio18","bio16")
#exclLyrs<-c("bio11","bio3","bio8","bio18","bio16")
#exclLyrs<-c("bio11","bio3","bio8","bio18","bio1")
exclLyrs<-c("bio11","bio3","bio8","bio18","bio1","bio16")
exclLyrsI<-which(names(envcrop) %in% exclLyrs)
envcrop<-subset(envcrop,exclLyrsI,negate=T)
names(envcrop)

#can't cope with SpatRaster need to go back to raster
predictors <- stack(envcrop)
predLabs<-names(predictors)


#another diversion - maybe removing locs without climate data 
#may help mahab distance algorithm, presumed removed NAs but ??
# Extract climate data for all those points
bioclim_extract <- extract(x = predictors,
                           y = presence,
                           ID = FALSE) # No need for an ID column
dim(bioclim_extract)
length(which(complete.cases(bioclim_extract)==F))
presence<-presence[which(complete.cases(bioclim_extract)==T),]
dim(presence)
nobs<-nrow(presence)


set.seed(0)
group <- kfold(presence, 5)
pres_train <- presence[group != 1, ]
pres_test <- presence[group == 1, ]

backg <-spatSample(x = envcrop, size = nobs,
                   values = FALSE, # don't need values
                   na.rm = TRUE,   # don't sample from ocean
                   xy = TRUE)      # just need coordinates
colnames(backg) = c("longitude", "latitude")
dim(backg)
# bioclim_extract <- extract(x = predictors,
#                            y = backg,
#                            ID = FALSE) # No need for an ID column
# dim(bioclim_extract)
# length(which(complete.cases(bioclim_extract)==F))
# backg<-backg[which(complete.cases(bioclim_extract)==T),]
# dim(backg)
#don't need to worry for backg cos only samples where values avail


group <- kfold(backg, 5)
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]

# plot(!is.na(envcrop[[1]]), col=c('white', 'light grey'), legend=FALSE)
# points(backg_train, pch='-', cex=0.5, col='yellow')
# points(backg_test, pch='-',  cex=0.5, col='black')
# points(pres_train, pch= '+', col='green')
# points(pres_test, pch='+', col='blue')


#Mahalanobis distance method
mm <- mahal(predictors, pres_train)
#mm <- mahal(predictors, pres_train, factors="landcover")
e <- evaluate(pres_test, backg_test, mm, predictors)
e
#w/ spatial subsampling no landcover, bioclim 1,12,2,5,6
# class          : ModelEvaluation 
# n presences    : 2089 
# n absences     : 2089 
# AUC            : 0.862778 
# cor            : 0.271807 
# max TPR+TNR at : 0.9267904 
#w/ spatial subsampling no landcover, bioclim 16,12,2,5,6
# class          : ModelEvaluation 
# n presences    : 2089 
# n absences     : 2089 
# AUC            : 0.8607445 
# cor            : 0.3481036 
# max TPR+TNR at : 0.9412704 
#w/ spatial subsampling no landcover, bioclim 12,2,5,6
# class          : ModelEvaluation 
# n presences    : 2089 
# n absences     : 2089 
# AUC            : 0.8572433 
# cor            : 0.3443403 
# max TPR+TNR at : 0.9635267 
pm = predict(predictors, mm, progress='')
par(mfrow=c(1,2))
pm[pm < -10] <- -10
map('world',xlim=c(xmin,xmax),ylim=c(ymin,ymax))
plot(pm, main='Mahalanobis distance')
points(pres_train, pch=3, cex=0.3)
threshold(e)
trmh <- threshold(e, 'spec_sens')
trmh <- threshold(e, 'prevalence')
trmh <- threshold(e, 'equal_sens_spec') 
#trmh <- threshold(e, 'no_omission')
plot(pm > trmh, main='presence/absence')
#plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch=3, cex=0.3)

#Maxent method
#maxent()
# xm <- maxent(predictors, pres_train,args=c('removeduplicates=FALSE'))
# xm <- maxent(predictors, pres_train, path=paste0(getwd(),"/maxent_out/"),
#              args=c('responsecurves=TRUE','jackknife=TRUE','writeplotdata=TRUE','removeduplicates=FALSE'))
# xm <- maxent(predictors, pres_train, factors="landcover",
#              args=c('responsecurves=TRUE','jackknife=TRUE','writeplotdata=TRUE','removeduplicates=FALSE'))
xm <- maxent(predictors, pres_train, backg_train, path=paste0(getwd(),"/maxent_out/"),
             args=c('responsecurves=TRUE','jackknife=TRUE','writeplotdata=TRUE',
                    'outputformat=cloglog','removeduplicates=FALSE'))
# xm2 <- maxent(predictors, presence, backg, path=paste0(getwd(),"/maxent_out/"),
#              args=c('responsecurves=TRUE','jackknife=TRUE','writeplotdata=TRUE',
#                     'outputformat=cloglog','removeduplicates=FALSE'))
# xm3 <- maxent(predictors, pres_train, backg_train, path=paste0(getwd(),"/maxent_out/"),
#              args=c('responsecurves=TRUE','jackknife=TRUE','writeplotdata=TRUE',
#                     'outputformat=cloglog','randomtestpoints=20','removeduplicates=FALSE'))
# xm3
par(mfrow=c(1,1))
plot(xm)
response(xm)

em <- evaluate(pres_test, backg_test, xm, predictors)
em
#w/ spatial subsampling no landcover, bioclim 1,12,2,5,6
# class          : ModelEvaluation 
# n presences    : 2089 
# n absences     : 2089 
# AUC            : 0.8875949 
# cor            : 0.682609 
# max TPR+TNR at : 0.4729085 
#w/ spatial subsampling no landcover, bioclim 16,12,2,5,6
# class          : ModelEvaluation 
# n presences    : 2089 
# n absences     : 2089 
# AUC            : 0.8884849 
# cor            : 0.6829653 
# max TPR+TNR at : 0.4866778 
#w/ spatial subsampling no landcover, bioclim 12,2,5,6
# class          : ModelEvaluation 
# n presences    : 2089 
# n absences     : 2089 
# AUC            : 0.8882717 
# cor            : 0.6822893 
# max TPR+TNR at : 0.4856429 
#why is AUC in maxent html output and evaluate here so different?
xm
evaluate(pres_train, backg_train, xm, predictors)

px <- predict(predictors, xm, progress='')
par(mfrow=c(1,2))
plot(px, main='Maxent, raw values')
plot(europe,add=TRUE)
points(pres_train, pch='+', cex=0.2)
threshold(em)
trme <- threshold(em, 'spec_sens')
# trme <- threshold(em, 'prevalence')
# trme <- threshold(em, 'equal_sens_spec')
#trme <- threshold(e, 'no_omission')
# trme<-0.128
plot(px > trme, main='presence/absence')
plot(europe,add=TRUE)
#plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+', cex=0.2)

#now predict for NZ using this model
envcrop<-rast("C:/.../Climate Change Pests/Env_Data/NZdata/NZBioClims.tif")
# exclLyr<-which(names(envcrop)=="bio11")
# envcrop<-subset(envcrop,exclLyr,negate=TRUE)
#exclLyrs<-c("bio11","bio3","bio8","bio18")
#exclLyrs<-c("bio11","bio3","bio8","bio18","bio16")
#exclLyrs<-c("bio11","bio3","bio8","bio18","bio1")
exclLyrs<-c("bio11","bio3","bio8","bio18","bio1","bio16")
exclLyrsI<-which(names(envcrop) %in% exclLyrs)
envcrop<-subset(envcrop,exclLyrsI,negate=T)

#can't cope with SpatRaster need to go back to raster
newPreds <- stack(envcrop)

obs_data <- read.csv(file = "VespGerm_Locs/VespGermPresNZ.csv")
head(obs_data)
# Pull out coordinate columns, x (longitude) first, then y (latitude) from 
presenceNZ <- obs_data[, c("lon", "lat")]
colnames(presenceNZ) <- c("longitude", "latitude")

# projNZmm<- dismo::predict(mm,newPreds, progress='')
# par(mfrow=c(1,2))
# #projNZmm[projNZmm < -10] <- -10
# # trmh<-threshold(e,'no_omission')
# # trmh<--8
# plot(projNZmm, main='Mahalanobis distance')
# points(presenceNZ, pch='+', cex=0.3)
# plot(projNZmm> trmh, main='presence/absence')
# points(presenceNZ, pch='+', cex=0.3)


projNZme <- predict(xm, newPreds, progress='')
par(mfrow=c(1,2))
# trme<-threshold(em,'equal_sens_spec')
# trme<-0.47
plot(projNZme, main='Maxent, NZ')
plot(NZc,add=T)
points(presenceNZ, pch='+', cex=0.3)
plot(projNZme > trme, main='presence/absence')
plot(NZc,add=T)
points(presenceNZ, pch='+', cex=0.3)

projNZmeBin<-projNZme>trme 
nzvals<-extract(projNZmeBin,presenceNZ)
# length(nzvals)
# length(which(nzvals==0))
# length(which(nzvals==1))
length(which(nzvals==1))/(length(which(nzvals==1))+length(which(nzvals==0)))
# 727/(727+17)
# [1] 0.9771505

writeRaster(projNZme,"projNZVgermDistHist.tif", filetype="GTiff",overwrite=TRUE)
writeRaster(projNZmeBin,"projNZVgermDistHistBin.tif", filetype="GTiff",overwrite=TRUE)

#now predict for future NZ using this model
inclLyrs<-names(newPreds)
folderStr<-c("ssp245_2041_2060","ssp245_2081_2100",
             "ssp585_2041_2060","ssp585_2081_2100")

for (k in 1:length(folderStr)){
  scenario<-folderStr[k]
  bioclim.files <- list.files(path=paste0("C:/Climate/cropToNZ/", scenario),
                              pattern='tif', full.names=TRUE)

  for (i in 1:length(bioclim.files)){
    filenom<-bioclim.files[i]
    newfilenom<-paste0(getwd(),"/ProjVgermDist/projNZ_",scenario,".tif")
    envcrop<-rast(filenom)
    inclLyrsI<-which(names(envcrop) %in% inclLyrs)
    envcrop<-subset(envcrop,inclLyrsI)
    #can't cope with SpatRaster need to go back to raster
    newPreds <- stack(envcrop)
    projNZme <- predict(xm, newPreds, progress='')
    #names(projNZme)<-substr(filenom,start=54, stop=nchar(filenom)-24)
    #newRast[[i]]<-projNZme
    if (i==1){
      newRast<-projNZme
    } else {
      newRast<-stack(newRast,projNZme)
    }
  }  
  #newRast
  writeRaster(newRast,newfilenom, filetype="GTiff",overwrite=TRUE)
} #end k/scenario loop  


# newin<-rast("ProjVgermDist/projNZ_ssp245_2041_2060.tif")
# plot(newin)
# LyrNames<-c("ACCESS.CM2", "BCC.CSM2.MR", "CMCC.ESM2", "EC.Earth3.Veg", "GISS.E2.1.G", "UKESM1.0.LL")
# rmean <- app(newin, mean)

for (k in 1:length(folderStr)){
  scenario<-folderStr[k]
  filenom<-paste0(getwd(),"/ProjVgermDist/projNZ_",scenario,".tif")
  newfilenom<-paste0(substr(filenom,start=1,stop=nchar(filenom)-4),"mean.tif")
  rIn<-rast(filenom)
  rmean <- app(rIn, mean)
  # plot(rmean)
  # plot(rmean > trme, main=scenario)
  writeRaster(rmean,newfilenom, filetype="GTiff",overwrite=TRUE)
}



VgermDist.files <- list.files(path=paste0(getwd(), "/ProjVgermDist"),
                            pattern='mean.tif', full.names=TRUE)
VgermDist.files <- c(paste0(getwd(),"/projNZVgermDistHist.tif"),VgermDist.files)
projVgermDist <-rast(VgermDist.files)
trme<-0.4856429
projVgermDist
plot(projVgermDist)

newnames<-c("Historic (1970-2000)","SSP2-4.5 2050","SSP2-4.5 2090",
            "SSP5-8.5 2050","SSP5-8.5 2090")
names(projVgermDist)<-newnames
#writeRaster(projVgermDist,"projVgermDistHistnProj.tif", filetype="GTiff",overwrite=TRUE)



cols<-c("white","indianred")
plot(subset(projVgermDist,1)> trme,col=cols, 
     type = 'classes',main=names(projVgermDist)[1],legend=F)
plot(NZc,add=T)
plot(subset(projVgermDist,2:5)> trme,col=cols,
     type = 'classes',legend=F)
plot(NZc,add=T)

projVgermDistBin<-projVgermDist
for (i in 1:5){
  projVgermDistBin[[i]]<-projVgermDist[[i]]>trme  
}
#writeRaster(projVgermDistBin,"projVgermDistHistnProjBin.tif", filetype="GTiff",overwrite=TRUE)

projVgermDistBin<-rast("projVgermDistHistnProjBin.tif")

presCounts <- data.frame(scenario=names(projVgermDistBin),
                        numOnes=rep(0,5),
                        numZeros=rep(0,5),
                        numNAs=rep(0,5))
for (i in 1:5){
  nzvals<-extract(projVgermDistBin[[i]],presenceNZ)
  presCounts[i,2]<-length(which(nzvals==1))
  presCounts[i,3]<-length(which(nzvals==0))
  presCounts[i,4]<-length(which(is.na(nzvals)))
  #print(length(which(nzvals==1))/(length(which(nzvals==1))+length(which(nzvals==0))))
}
presCounts$numVals<-presCounts$numOnes+presCounts$numZeros
presCounts$numPts<-presCounts$numVals+presCounts$numNAs
presCounts$sens<-presCounts$numOnes/presCounts$numVals
# > presCounts
# scenario numOnes numZeros numNAs numVals numPts      sens
# 1 Historic (1970-2000)     728       17     12     745    757 0.9771812
# 2        SSP2-4.5 2050     731       20      6     751    757 0.9733688
# 3        SSP2-4.5 2090     728       23      6     751    757 0.9693742
# 4        SSP5-8.5 2050     728       23      6     751    757 0.9693742
# 5        SSP5-8.5 2090     702       49      6     751    757 0.9347537

distCounts <- data.frame(scenario=names(projVgermDistBin),
                         numOnes=rep(0,5),
                         numZeros=rep(0,5),
                         numNAs=rep(0,5))
for (i in 1:5){
  cellvals<-values(projVgermDistBin[[i]])
  distCounts[i,2]<-length(which(cellvals==1))
  distCounts[i,3]<-length(which(cellvals==0))
  distCounts[i,4]<-length(which(is.na(cellvals)))
  #print(length(which(nzvals==1))/(length(which(nzvals==1))+length(which(nzvals==0))))
}
distCounts$numVals<-distCounts$numOnes+distCounts$numZeros
distCounts$numPts<-distCounts$numVals+distCounts$numNAs
distCounts$propDist<-distCounts$numOnes/distCounts$numVals
# > distCounts
# scenario                 numOnes numZeros numNAs numValsnumPts  propDist
# 1 Historic (1970-2000)   14219     3592  87116   17811 104927 0.7983269
# 2        SSP2-4.5 2050   14051     3797  87079   17848 104927 0.7872591
# 3        SSP2-4.5 2090   13631     4217  87079   17848 104927 0.7637270
# 4        SSP5-8.5 2050   13980     3868  87079   17848 104927 0.7832810
# 5        SSP5-8.5 2090   11340     6508  87079   17848 104927 0.6353653

# test<-mask(projVgermDistBin,NZc)
# plot(test[[1]])

cols<-c("white","indianred")
par(mfrow=c(1,1))
plot(subset(projVgermDistBin,1),col=cols, 
     type = 'classes',main=names(projVgermDistBin)[1],legend=F)
plot(NZc,add=T)
points(presenceNZ, pch='+', cex=0.5)

par(mfrow=c(2,2))
for (i in 2:5){
  plot(subset(projVgermDistBin,i),col=cols, 
       type = 'classes',main=names(projVgermDistBin)[i],legend=F)
  plot(NZc,add=T)
}
# y <- mask(projVgermDist[[1]], NZc,inverse=T)
# plot(y)

library(rnaturalearth)
# europe <- ne_countries(scale = 50, returnclass = "sf", continent = "Europe") 
# europe<-vect(europe)  

worldc <- ne_coastline(scale = 50, returnclass = "sf") 
worldc<-vect(worldc)  
e <- extent(-10, 50, 27.5, 67.3)
europe<- crop(worldc, e)
plot(europe)

# ne<-extent(newPreds)
# ne<-ext(projVgermDist)
ne <- extent(165.8, 179.0, -47.8, -34.0)
NZc<-crop(worldc, ne)
plot(NZc)

