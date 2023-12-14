rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
#####
# Load the North American Pollen Data Base
NPDB <- read.csv("./Data/Pollen/Raw/whitmoreetal2005_v1-72.csv")
# Make the locations a SpatVector
NPDBVect <- vect(NPDB[,c("LONDD", "LATDD")],
                 geom=c("LONDD", "LATDD"),
                 crs = "EPSG:4326")
#####

#####
# Load Present Monty Temperature from Lorentz downscale and Avg over the Baseline period 1960-1990
## 0akBP (band 2201) is 1950AD and it goes until 1990AD (band 2204)
Temp.Seson.Bsl.List <- lapply(which((-2200:3)%in%c(1:3)),
                              function(j){
                                out <- lapply(1:12,
                                              function(x){
                                                In1 <- raster("./Data/PaleoClimate/Raw/ccsm3_22-0k_temp.nc",
                                                              varname= "tmax",
                                                              band=j,
                                                              level=x)
                                                In2 <- raster("./Data/PaleoClimate/Raw/ccsm3_22-0k_temp.nc",
                                                              varname= "tmin",
                                                              band=j,
                                                              level=x)
                                                Out <-rast(mean(stack(In1 + In2)))
                                              })
                                #out <- rast(mean(do.call("stack",out)))
                                #return(out)
                                TempSummPres <- do.call("c",out)
                                Temp.SeasonalPres <- c(app(TempSummPres[[3:5]],mean),
                                                       app(TempSummPres[[6:8]],mean),
                                                       app(TempSummPres[[9:11]],mean),
                                                       app(TempSummPres[[c(12,1,2)]],mean))
                                return(Temp.SeasonalPres)
                              })
# Estimate the Mean Seasonal Temperature for the 1960,1970, & 1990 downscaled values
Temp.MnSeson.Bsl <- c(app(do.call("c",lapply(Temp.Seson.Bsl.List,function(x){x[[1]]})),mean),
                      app(do.call("c",lapply(Temp.Seson.Bsl.List,function(x){x[[2]]})),mean),
                      app(do.call("c",lapply(Temp.Seson.Bsl.List,function(x){x[[3]]})),mean),
                      app(do.call("c",lapply(Temp.Seson.Bsl.List,function(x){x[[4]]})),mean))
names(Temp.MnSeson.Bsl) <- c("Mn.Sprg","Mn.Summ","Mn.Fall","Mn.Wint")
Temp.MnSeson.Bsl
# Estimate the SD Seasonal Temperture for the 1960,1970, & 1990 downscaled values
Temp.Sd.Seson.Bsl <- c(app(do.call("c",lapply(Temp.Seson.Bsl.List,function(x){x[[1]]})),sd),
                       app(do.call("c",lapply(Temp.Seson.Bsl.List,function(x){x[[2]]})),sd),
                       app(do.call("c",lapply(Temp.Seson.Bsl.List,function(x){x[[3]]})),sd),
                       app(do.call("c",lapply(Temp.Seson.Bsl.List,function(x){x[[4]]})),sd))
names(Temp.Sd.Seson.Bsl) <- c("Sd.Sprg","Sd.Summ","Sd.Fall","Sd.Wint")
Temp.Sd.Seson.Bsl

# Get the Mean Values for NPDs sites
TempNPD <- extract(Temp.MnSeson.Bsl,
                   NPDBVect)
names(TempNPD) <- c("ID",paste0("Mn.",c("Sprg","Summ","Fall","Wint")))

# Get the SD Values for NPDs sites
SDTempNPD <- extract(Temp.Sd.Seson.Bsl,
                     NPDBVect)
names(SDTempNPD) <- c("ID",paste0("Sd.",c("Sprg","Summ","Fall","Wint")))
#####

#####
# Load Present Seasonal Precipitation from Lorentz downscale and Avg over the Baseline period 1960-1990
## 0akBP (band 2201) is 1950AD and it goes until 1990AD (band 2204)
Prec.Seson.Bsl.List <- lapply(which((-2200:3)%in%c(1:3)),
                              function(j){
                                out <- lapply(1:12,
                                              function(x){
                                                Out <- raster("./Data/PaleoClimate/Raw/ccsm3_22-0k_prcp.nc",
                                                              varname= "prcp",
                                                              band=j,
                                                              level=x)
                                                Out <-rast(mean(Out))
                                                return(Out)
                                              })
                                #out <- rast(mean(do.call("stack",out)))
                                #return(out)
                                PrecSummPres <- do.call("c",out)
                                Prec.SeasonalPres <- c(app(PrecSummPres[[3:5]],sum),
                                                       app(PrecSummPres[[6:8]],sum),
                                                       app(PrecSummPres[[9:11]],sum),
                                                       app(PrecSummPres[[c(12,1,2)]],sum))
                                return(Prec.SeasonalPres)
                              })

# Estimate the Mean Seasonal Precipitation for the 1960,1970, & 1990 downscaled values
Prec.MnSeson.Bsl <- c(app(do.call("c",lapply(Prec.Seson.Bsl.List,function(x){x[[1]]})),mean),
                      app(do.call("c",lapply(Prec.Seson.Bsl.List,function(x){x[[2]]})),mean),
                      app(do.call("c",lapply(Prec.Seson.Bsl.List,function(x){x[[3]]})),mean),
                      app(do.call("c",lapply(Prec.Seson.Bsl.List,function(x){x[[4]]})),mean))
names(Prec.MnSeson.Bsl) <- c("Mn.Sprg","Mn.Summ","Mn.Fall","Mn.Wint")
Prec.MnSeson.Bsl
# Estimate the SD Seasonal Temperture for the 1960,1970, & 1990 downscaled values
Prec.Sd.Seson.Bsl <- c(app(do.call("c",lapply(Prec.Seson.Bsl.List,function(x){x[[1]]})),sd),
                       app(do.call("c",lapply(Prec.Seson.Bsl.List,function(x){x[[2]]})),sd),
                       app(do.call("c",lapply(Prec.Seson.Bsl.List,function(x){x[[3]]})),sd),
                       app(do.call("c",lapply(Prec.Seson.Bsl.List,function(x){x[[4]]})),sd))
names(Prec.Sd.Seson.Bsl) <- c("Sd.Sprg","Sd.Summ","Sd.Fall","Sd.Wint")
Prec.Sd.Seson.Bsl

# Get the Mean Values for NPDs sites
PrecNPD <- extract(Prec.MnSeson.Bsl,
                   NPDBVect)
names(PrecNPD) <- c("ID",paste0("Mn.",c("Sprg","Summ","Fall","Wint")))

# Get the SD Values for NPDs sites
SDPrecNPD <- extract(Prec.Sd.Seson.Bsl,
                     NPDBVect)
names(SDPrecNPD) <- c("ID",paste0("Sd.",c("Sprg","Summ","Fall","Wint")))
#####

#####
# Get summary values for the seasonal temperature and precipitation values for NPDS sites
NPDB.Clim <- data.frame(TempNPD,
                        SDTempNPD[,-1],
                        PrecNPD[,-1],
                        SDPrecNPD[,-1]) 

names(NPDB.Clim)[-1] <- c(paste0("BslTempMn.",c("Sprg","Summ","Fall","Wint")),
                          paste0("BslTempSd.",c("Sprg","Summ","Fall","Wint")),
                          paste0("BslPrecMn.",c("Sprg","Summ","Fall","Wint")),
                          paste0("BslPrecSd.",c("Sprg","Summ","Fall","Wint")))
NPDB.Clim <- data.frame(NPDB[,c("ID1","SITE", "DBCODE", "SITECODE", "SITENAME", "LONDD", "LATDD","ELEVATION")],
                        NPDB.Clim[,-1])

NPDB.Clim <- NPDB.Clim[complete.cases(NPDB.Clim),]
# Save the climate (mean and SD) for the North American Pollen Data
write.csv(NPDB.Clim,
          "./Data/PaleoClimate/Processed/NPDB_Clim.csv")