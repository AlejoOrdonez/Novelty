rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")

#####
# Taxa translate between NPDB and Paleo Data
TaxtTrans <- read.csv("./Data/Pollen/Raw/Taxa transation.csv")
# Which are the summarizing taxa - note the TaxtTrans$With.Traits!="" filter -> only species with traits are used
TaxaUse <- sort(unique(TaxtTrans$Translated[TaxtTrans$With.Traits!=""]))

# Load the Pollen data
NEATOMA.ALL <- read.csv("./Data/Pollen/Raw/abund_spp_env_21to1bp.csv")
# Aggregate to taxa with data 
NEATOMA1 <-  as.data.frame(t(apply(NEATOMA.ALL[,TaxaUse],1,function(x){x/sum(x,na.rm=T)})))
# Make 0 the pollen count for taxa with percentages below 5%
NEATOMA1 <-  as.data.frame(t(apply(NEATOMA1,1,function(x){ifelse(x>0.05,x,0)})))
# Aggregate to taxa with data 
NEATOMA1 <-  as.data.frame(t(apply(NEATOMA1,1,function(x){x/sum(x,na.rm=T)})))
# Remove taxa with no counts in in all periods
NEATOMA1 <- NEATOMA1[,apply(NEATOMA1,2,sum)>0]
NEATOMA1[is.na(NEATOMA1)] <- 0


# Create a new aggregated NEATOMO Data Base
NEATOMA2 <- data.frame(NEATOMA.ALL[, c("Time","Clim.ID","sites","Latitude","Longitude")],
                       NEATOMA1)

# Sate the aggregated and curated file
write.csv(NEATOMA2,"./Data/Pollen/Processed/NEATOMA_Agg.csv")
