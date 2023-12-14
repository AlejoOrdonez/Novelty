rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
#setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
setwd("/Users/alejandroordonez/Documents/3. Ecography PaleoTraits Paper/PaleoNovelty")
#####
# Load the Aggregated North American Pollen Data Base
NPDB.Pollen <- read.csv("./Data/Pollen/Processed/NPDB_Agg.csv")
# Load the climate for the North American Pollen Data Base
NPDB.Clim <- read.csv("./Data/PaleoClimate/Processed/NPDB_Clim.csv")

# Select only sites with Pollen data and Climate data
NPDB.SiteUse <- table(c(NPDB.Clim$ID1,NPDB.Pollen$ID1))
NPDB.SiteUse <- as.numeric(names(NPDB.SiteUse)[(NPDB.SiteUse==2)])
# Climate
NPDB.Clim <- NPDB.Clim[c(NPDB.Clim$ID1%in%NPDB.SiteUse),]
# Aggregated composition
NPDB.Pollen <- NPDB.Pollen[c(NPDB.Pollen$ID1%in%NPDB.SiteUse),]
# Load the Aggregated NEATOMA Data Base
NEATOMA.Agg <- read.csv("./Data/Pollen/Processed/NEATOMA_Agg.csv")

# Which taxa to use - ensure match between NPDS and NEATOMA
TaxaUse <- table(c(names(NEATOMA.Agg[,-c(1:6)]),names(NPDB.Pollen[,-c(1:10)])))
TaxaUse <- names(TaxaUse)[TaxaUse==2]
#####


#####
# Define the analogue threshold
# Estimate the distance between NPDB sites
Comp.CordDist.NPDS <- analogue::distance(NPDB.Pollen[,TaxaUse],
                                     method = "SQchord")
# estimate the cutoff value
NPDS.ROC.Cutoff <- roc(Comp.CordDist.NPDS, # current time Taxon data.frame 
                       groups = NPDB.Pollen$BIOME # vector of group memberships
                       )
NPDS.ROC.Cutoff
#####


#####
# Estimate the SQchord to present of each Paelo site
DistToCurrent <- lapply(1:dim(NEATOMA.Agg)[1],
                        function(i){
                          DistoNPDS <- analogue::distance(NPDB.Pollen[,TaxaUse],
                                                          NEATOMA.Agg[i,TaxaUse],
                                                          method = "SQchord")
                          DistoNPDS <- as.numeric(DistoNPDS)
                          return(DistoNPDS)
                        })
# Make it a distance matrix
DistToCurrent2 <- do.call("cbind",DistToCurrent)
colnames(DistToCurrent2) <- (NPDB.Pollen$sites)
rownames(DistToCurrent2) <- NPDB.Clim$SITECODE

# Estimate the Min Distance to present for each Neatoma site
NEATOMA.Agg$minSQchord <- as.numeric(apply(DistToCurrent2,2,min))


#####

# Save the climate dissimilarity output
NEATOMA.CommDis <- list(FullDisMatrix = DistToCurrent2,
                        NEATOMA.CommDis = NEATOMA.Agg[,-c(1,6:45)],
                        NEATOMA.Agg = NEATOMA.Agg,
                        NPDB.Agg = NPDB.Pollen,
                        ROC.Cutoff = NPDS.ROC.Cutoff)
saveRDS(NEATOMA.CommDis,"./Results/SQchord_TaxoDis.rds")

# Dummy plot (taking a even sub sample of sites across periods)
rm(list=ls());gc()
NEATOMA.CommDis <- readRDS("./Results/SQchord_TaxoDis.rds")
NEATOMA.Agg <- NEATOMA.CommDis[["NEATOMA.Agg"]]
CompDisBoot <- lapply(1:1000,
                      function(i){
                        SamplTmp <- do.call("c",lapply(1:21,
                                                       function(x){
                                                         sample(which(NEATOMA.Agg$Time==x), 10)
                                                       }))
                        tapply(NEATOMA.Agg$minSQchord[SamplTmp],
                               NEATOMA.Agg$Time[SamplTmp],
                               median)
                      })

CompDisBoot2 <- do.call("rbind",CompDisBoot)
CompDisBootQuant <- apply(CompDisBoot2,2,quantile,c(0.0275,0.5,0.975))
plot(y = rev(CompDisBootQuant[2,]),
     x = -21:-1,
     type = "b",
     main = "Regional Dissimilarity Env Change", 
     xlab ="Time (kyrBP)",
     ylab = "Disimularity (stndz Euc Dist)",
     ylim=range(CompDisBootQuant))
lines(y = rev(CompDisBootQuant[1,]),
      x = -21:-1,)
lines(y = rev(CompDisBootQuant[3,]),
      x = -21:-1,)
abline(h=NEATOMA.CommDis[['ROC.Cutoff']]$roc$ Combined$ optima)

a <- table(NEATOMA.CommDis[['NEATOMA.Agg']]$minSQchord>NEATOMA.CommDis[['ROC.Cutoff']]$roc$ Combined$ optima,
           NEATOMA.CommDis[['NEATOMA.Agg']]$Time)

apply(a,2,function(x){x/sum(x)})


rm(list=ls());gc()
setwd("/Users/alejandroordonez/Documents/3. Ecography PaleoTraits Paper/PaleoNovelty/")
NEATOMA.CommDis <- readRDS("./Results/SQchord_TaxoDis.rds")
NEATOMA.Agg <- NEATOMA.CommDis[["NEATOMA.Agg"]]
CompDisBoot <- lapply(1:1000,
                      function(i){
                        SamplTmp <- do.call("c",lapply(1:21,
                                                       function(x){
                                                         sample(which(NEATOMA.Agg$Time==x), 10)
                                                       }))
                        tapply(NEATOMA.Agg$minSQchord[SamplTmp],
                               NEATOMA.Agg$Time[SamplTmp],
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

ggplot(TaxaDistCI, aes(x = Time, y = Phylo)) +
  geom_polygon(fill="lightblue") + 
  geom_line(data=TaxaDistMean,color="black") + 
  geom_point(data=TaxaDistMean,color="black") +
  xlab("Time (kaBP)") + # for the x axis label
  ylab("Dif to present \n[Sqr-Chord Distance] ") +# for the y axis label
  ggtitle("Regional Taxonomic difference") +
  geom_hline(yintercept=NEATOMA.CommDis$ROC.Cutoff$roc$Combined$optimal,
             linetype="dashed", 
             color = "red")

#"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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"


