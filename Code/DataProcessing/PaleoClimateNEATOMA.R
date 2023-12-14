rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")

#####
# Load the Pollen data
Old <- read.csv("./Data/Pollen/Raw/abund_spp_env_21to1bp.csv")
Old2 <- data.frame(Old[, c("Time","Clim.ID","sites","Latitude","Longitude")])
#####

#####
# Load the Paleo Seasonal Temperature and summarize as 500yrs averages
TmpPaleoSeasonalList <- lapply(1:21,
                               function(TimeBP){#(TimeBP <- 21)
                                 Summ500YrsList <- lapply(1:12,
                                                          function(x){
                                                            out <- lapply(which((-2200:3)%in%c(c(-TimeBP*100)+c(-20,-10,00,10,20))),
                                                                          function(j){
                                                                            In1 <- raster("./Data/PaleoClimate/Raw/ccsm3_22-0k_temp.nc",
                                                                                          varname= "tmax",
                                                                                          band=j,
                                                                                          level=x)
                                                                            In2 <- raster("./Data/PaleoClimate/Raw/ccsm3_22-0k_temp.nc",
                                                                                          varname= "tmin",
                                                                                          band=j,
                                                                                          level=x)
                                                                            Out <- (In1 + In2)/2
                                                                            return(Out)
                                                                            
                                                                          })
                                                            out <- rast(mean(do.call("stack",out)))
                                                            return(out)
                                                          })
                                 Summ500Yrs <- do.call("c",Summ500YrsList)
                                 SeasonalSumm500Yrs <- c(app(Summ500Yrs[[3:5]],mean),
                                                         app(Summ500Yrs[[6:8]],mean),
                                                         app(Summ500Yrs[[9:11]],mean),
                                                         app(Summ500Yrs[[c(12,1,2)]],mean))
                                 names(SeasonalSumm500Yrs) <- c("Mn.Sprg","Mn.Summ","Mn.Fall","Mn.Wint")
                                 return(SeasonalSumm500Yrs)
                               })

# Extract the Seasonal Values
NEATOMA.Seasonal.Temp <- lapply(1:dim(Old2)[1],
                                function(x){
                                  # Make a spatial vector for the Observation
                                  Old3 <- vect(Old2[x, c("Longitude","Latitude")],
                                               geom = c("Longitude","Latitude"))
                                  # Distance in the temperature dimension
                                  
                                  SesTempPaleo <- terra::extract(TmpPaleoSeasonalList[[Old2[x,"Time"]]],
                                                                 Old3)
                                  return(SesTempPaleo)
                                })
# Merge into a table
NEATOMA.Seasonal.Temp <- do.call("rbind",NEATOMA.Seasonal.Temp)
names(NEATOMA.Seasonal.Temp)[-1] <- paste0("Temp.",names(NEATOMA.Seasonal.Temp)[-1])
#####

#####
# Load the Paleo Seasonal Precipitation and summarize as 500yrs averages
PrcPaleoSeasonalList <- lapply(1:21,
                               function(TimeBP){#(TimeBP <- 21)
                                 Summ500YrsList <- lapply(1:12,
                                                          function(x){
                                                            out <- lapply(which((-2200:3)%in%c(c(-TimeBP*100)+c(-20,-10,00,10,20))),
                                                                          function(j){
                                                                            raster("./Data/PaleoClimate/Raw/ccsm3_22-0k_prcp.nc",
                                                                                   varname= "prcp",
                                                                                   band=j,
                                                                                   level=x)
                                                                          })
                                                            out <- rast(mean(do.call("stack",out)))
                                                            return(out)
                                                          })
                                 Summ500Yrs <- do.call("c",Summ500YrsList)
                                 SeasonalSumm500Yrs <- c(app(Summ500Yrs[[3:5]],sum),
                                                         app(Summ500Yrs[[6:8]],sum),
                                                         app(Summ500Yrs[[9:11]],sum),
                                                         app(Summ500Yrs[[c(12,1,2)]],sum))
                                 names(SeasonalSumm500Yrs) <- c("Mn.Sprg","Mn.Summ","Mn.Fall","Mn.Wint")
                                 return(SeasonalSumm500Yrs)
                               })

# Extract the Seasonal Values
NEATOMA.Seasonal.Prec <- lapply(1:dim(Old2)[1],
                                function(x){
                                  # Make a spatial vector for the Observation
                                  Old3 <- vect(Old2[x, c("Longitude","Latitude")],
                                               geom = c("Longitude","Latitude"))
                                  # Distance in the temperature dimension
                                  
                                  SesTempPaleo <- terra::extract(PrcPaleoSeasonalList[[Old2[x,"Time"]]],
                                                                 Old3)
                                  return(SesTempPaleo)
                                })
# Merge into a table
NEATOMA.Seasonal.Prec <- do.call("rbind",NEATOMA.Seasonal.Prec)
names(NEATOMA.Seasonal.Prec)[-1] <- paste0("Prec.",names(NEATOMA.Seasonal.Prec)[-1])

#####
# Get summary values for the seasonal temperature and precipitation values for NEATOMA sites
NEATOMA.Clim <- data.frame(NEATOMA.Seasonal.Temp,
                           NEATOMA.Seasonal.Prec[,-1]) 

NEATOMA.Clim <- data.frame(ID = NEATOMA.Clim[,1],
                        Old2,
                        NEATOMA.Clim[,-1])

# Save the CRU.TS climate for the North American Pollen Data Base
write.csv(NEATOMA.Clim,
          "./Data/PaleoClimate/Processed/NEATOMA_Clim.csv")