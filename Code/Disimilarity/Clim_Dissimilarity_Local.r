rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
#####
#####
# Estimate the dissimilarity to present (TRaCE21 defined) of paleo-climates of cored sites.
# Load the climate of the NEATOMA Data Base
NEATOMA.Clim <- read.csv("./Data/PaleoClimate/Processed/NEATOMA_Clim.csv")
# Filter sites that not take 1 ka BP and have less than 6 dates
SummSyteXTime <- table(NEATOMA.Clim$sites,NEATOMA.Clim$Time)
SummSyteXTime <- SummSyteXTime[c(apply(SummSyteXTime,1,sum)>6 & SummSyteXTime[,1]==1),]


# Estimate Local climate dissimilarity - Standrized eucliden distance
LocalClimDisList <- lapply(rownames(SummSyteXTime),
                           function(Site){#(Site <- "3PINES"unique(NEATOMA.Clim$sites)[2]
                             # Estimate the standardized Euclidean distance
                             SiteStdEucDisList <- lapply(c("Sprg","Summ","Fall","Wint"),#c("Summ","Wint")
                                                         function(x){#x <- "Sprg"
                                                           # Distance for Temp
                                                           ClimMtx <- NEATOMA.Clim[NEATOMA.Clim$sites==Site,paste0("Temp.Mn.",x)]
                                                           names(ClimMtx) <- NEATOMA.Clim$Time[NEATOMA.Clim$sites==Site]
                                                           DistTemp <- (as.matrix(dist(ClimMtx)^2)/sd(ClimMtx))["1",]
                                                           # Distance for Prec
                                                           ClimMtx <- NEATOMA.Clim[NEATOMA.Clim$sites==Site,paste0("Prec.Mn.",x)]
                                                           names(ClimMtx) <- NEATOMA.Clim$Time[NEATOMA.Clim$sites==Site]
                                                           DistPrc <- (as.matrix(dist(ClimMtx)^2)/sd(ClimMtx))["1",]
                                                           # Distance Temp+Prec - > sesonal 
                                                           DistSeason <- DistTemp+DistPrc
                                                           return(DistSeason)
                                                         })
                             # Summary of standardized Euclidean distance
                             SiteStdEucDis <- data.frame(Site=Site,
                                                         Time = NEATOMA.Clim$Time[NEATOMA.Clim$sites == Site],
                                                         SESClim= (as.numeric(apply(do.call("cbind",SiteStdEucDisList),1,sum)))^0.5)
                             return(SiteStdEucDis)
                           })
# Turn list into a table
  LocalClimDis <- do.call("rbind",LocalClimDisList)
  saveRDS(LocalClimDis,"./Results/StdEuDis_ClimDis_Local.rds")

  
# Plot (taking a even sub sample of sites across periods)
rm(list=ls()); gc()
LocalClimDis <- readRDS("./Results/StdEuDis_ClimDis_Local.rds")

# Estimate the mean dissimilarity per-time period based on 10 samples
ClimDisBoot <- lapply(1:1000,
                      function(i){
                        SamplTmp <- do.call("c",lapply(1:21,
                                                       function(x){
                                                         sample(which(LocalClimDis$Time==x), 10)
                                                       }))
                        tapply(LocalClimDis$SESClim[SamplTmp],
                               LocalClimDis$Time[SamplTmp],
                               median)
                      })

ClimDisBoot2 <- do.call("rbind",ClimDisBoot)
ClimDisBootQuant <- apply(ClimDisBoot2,2,quantile,c(0.0275,0.5,0.975))

require(tidyverse)
require(ggbreak)
ClimDistCI<-data.frame(Time = -c(1:21,21:1),
                       Phylo = c(ClimDisBootQuant["2.75%",],
                                 rev(ClimDisBootQuant["97.5%",])))
ClimDistMean <- data.frame(Time = -(1:21),
                           Phylo = ClimDisBootQuant["50%",])

NEATOMA.Clim <- readRDS("./Results/StdEuDis_ClimDis.rds")
ggplot(ClimDistCI, aes(x = Time, y = Phylo)) +
  geom_polygon(fill="lightblue") + 
  geom_line(data=ClimDistMean,color="black") + 
  geom_point(data=ClimDistMean,color="black") +
  xlab("Time (kaBP)") + # for the x axis label
  ylab("Dif to present \n[Stz Euclidean Distance] ") +# for the y axis label
  ggtitle("Local Climate difference") +
  geom_hline(yintercept=NEATOMA.Clim$ROC.Cutoff$roc$Combined$optimal,
             linetype="dashed", 
             color = "red")




plot(y = rev(ClimDisBootQuant[2,]),
     x = -21:-1,
     type = "b",
     main = "Local Dissimilarity Env Change", 
     xlab ="Time (kyrBP)",
     ylab = "Disimularity (stndz Euc Dist)",
     ylim=range(ClimDisBootQuant))
lines(y = rev(ClimDisBootQuant[1,]),
      x = -21:-1)
lines(y = rev(ClimDisBootQuant[3,]),
      x = -21:-1)
NEATOMA.ClimDis <- readRDS("./Results/StdEuDis_ClimDis.rds")
abline(h=NEATOMA.ClimDis$ROC.Cutoff$roc$ Combined$ optima)



