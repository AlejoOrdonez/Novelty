rm(list=ls());gc()
require(analogue)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
#####
#####
# Estimate the dissimilarity to present of paleo-assemblages of cored sites.
# Load the Aggregated NEATOMA Data Base
NEATOMA.Comm <- read.csv("./Data/Pollen/Processed/NEATOMA_Agg.csv")
# Filter sites that not take 1 ka BP and have less than 6 dates
SummSyteXTime <- table(NEATOMA.Comm$sites,NEATOMA.Comm$Time)
SummSyteXTime <- SummSyteXTime[c(apply(SummSyteXTime,1,sum)>6 & SummSyteXTime[,1]==1),]


# Estimate Local composition dissimilarity - Squared Coord distance
LocalCompDisList <- lapply(rownames(SummSyteXTime),
                           function(Site){#(Site <- rownames(SummSyteXTime)[2])
                             CommDataTmp <- NEATOMA.Comm[NEATOMA.Comm$sites==Site,-c(1:6)]
                             rownames(CommDataTmp) <- NEATOMA.Comm$Time[NEATOMA.Comm$sites == Site]
                             DistoNPDS <- analogue::distance(CommDataTmp,
                                                             method = "SQchord")
                             # Summary of Squared Coord distance
                             SiteSQchordDis <- data.frame(Site=Site,
                                                         Time = NEATOMA.Comm$Time[NEATOMA.Comm$sites == Site],
                                                         SQchordComp= DistoNPDS[dim(DistoNPDS)[1],])
                             return(SiteSQchordDis)
                           })
# Turn list into a table
LocalCompDis <- do.call("rbind",LocalCompDisList)
saveRDS(LocalCompDis,"./Results/SQchord_CompDis_Local.rds")

# Plot (taking a even sub sample of sites across periods)
rm(list=ls()); gc()
LocalCompDis <- readRDS("./Results/SQchord_CompDis_Local.rds")

# Estimate the mean dissimilarity per-time period based on 10 samples
CompDisBoot <- lapply(1:1000,
                      function(i){
                        SamplTmp <- do.call("c",lapply(1:21,
                                                       function(x){
                                                         sample(which(LocalCompDis$Time==x), 10)
                                                       }))
                        tapply(LocalCompDis$SQchordComp[SamplTmp],
                               LocalCompDis$Time[SamplTmp],
                               median)
                      })

CompDisBoot2 <- do.call("rbind",CompDisBoot)
CompDisBootQuant <- apply(CompDisBoot2,2,quantile,c(0.0275,0.5,0.975))

require(tidyverse)
require(ggbreak)
TaxaDistCI<-data.frame(Time = -c(1:21,21:1),
                       Phylo = c(CompDisBootQuant["2.75%",],
                                 rev(CompDisBootQuant["97.5%",])))
TaxaDistMean <- data.frame(Time = -(1:21),
                           Phylo = CompDisBootQuant["50%",])

NEATOMA.CommDis <- readRDS("./Results/SQchord_TaxoDis.rds")
ggplot(TaxaDistCI, aes(x = Time, y = Phylo)) +
  geom_polygon(fill="lightblue") + 
  geom_line(data=TaxaDistMean,color="black") + 
  geom_point(data=TaxaDistMean,color="black") +
  xlab("Time (kaBP)") + # for the x axis label
  ylab("Dif to present \n[Sqr-Chord Distance] ") +# for the y axis label
  ggtitle("Local Taxonomic difference") +
  geom_hline(yintercept=NEATOMA.CommDis$ROC.Cutoff$roc$Combined$optimal,
             linetype="dashed", 
             color = "red")



plot(y = rev(CompDisBootQuant[2,]),
     x = -21:-1,
     type = "b",
     main = "Local Dissimilarity compositional change", 
     xlab ="Time (kyrBP)",
     ylab = "Disimularity (Squared-Cord Dist)",
     ylim=range(CompDisBootQuant))
lines(y = rev(CompDisBootQuant[1,]),
      x = -21:-1)
lines(y = rev(CompDisBootQuant[3,]),
      x = -21:-1)
NEATOMA.CompDis <- readRDS("./Results/SQchord_TaxoDis 2.rds")
abline(h=NEATOMA.CompDis$ROC.Cutoff$roc$ Combined$ optima)



