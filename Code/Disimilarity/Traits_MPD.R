rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
require(snowfall)
require(picante)
setwd("/Users/alejandroordonez/Documents/3. Ecography PaleoTraits Paper/PaleoNovelty/")
#setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
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
TaxaUseObs <- table(c(names(NEATOMA.Agg[,-c(1:6)]),names(NPDB.Pollen[,-c(1:10)])))
TaxaUseObs <- names(TaxaUseObs)[TaxaUseObs==2]
#####

#######################################################################################################################

SrtTime1 <- Sys.time()
BootTraits <- lapply(1:100,
                     function(RunIt){
                       SrtTime <- Sys.time()
                       # Generate a possible trait realization based on the know distribution of values for a evaluated taxa
                       # Load the trait data
                       Traits <- read.csv2("./Data/Traits/Raw/TraitsSummary.csv")
                       Traits <- Traits[c(Traits$Name%in%TaxaUseObs),]
                       # Seeds 
                       SWT <- apply(Traits[,c("Min..of.Seed.Size..mg.","Max..of.Seed.Size..mg.")],
                                    1,
                                    function(x){
                                      ifelse(is.na(x[1]),
                                             NA,
                                             runif(1,log10(x[1]),log10(x[2])))
                                    })
                       SWT <- scale(SWT)
                       
                       # Height
                       Hmax <- apply(Traits[,c("Min..of.Height..m.","Max..of.Height..m.")],
                                     1,
                                     function(x){
                                       ifelse(is.na(x[1]),
                                              NA,
                                              runif(1,x[1],x[2]))
                                     })
                       Hmax <- scale(Hmax)
                       # LMA
                       LMA <- apply(Traits[,c("Min..of.LMA..log...g.m2.","Max..of.LMA..log...g.m2.")],
                                    1,
                                    function(x){
                                      ifelse(is.na(x[1]),
                                             NA,
                                             runif(1,x[1],x[2]))
                                    })
                       LMA <- scale(LMA)
                       
                       # Build a trait dataset
                       TraitsDBS <- data.frame(SWT = SWT,
                                               Hmax = Hmax,
                                               LMA = LMA)
                       row.names(TraitsDBS) <- Traits$Name
                       TraitsDBS <- TraitsDBS[complete.cases(TraitsDBS),]
                       # Which taxa to use - ensure match between NPDS and NEATOMA and TRAITS
                       TaxaUse <- table(c(TaxaUseObs,row.names(TraitsDBS)))
                       TaxaUse <- names(TaxaUse)[TaxaUse==2]
                       # Estimate the Trait distance in RAW Space                       
                       TraitsDist <- as.matrix(dist(TraitsDBS[TaxaUse,]))
                       #####
                       #####
                       # Trait Mean pairwise distance contrast
                       NPDB.MPD.ALL <- mpd(NPDB.Pollen[,colnames(TraitsDist)],
                                          TraitsDist,
                                          abundance.weighted=T)
                       names(NPDB.MPD.ALL) <- paste0(NPDB.Pollen$DBCODE,"_",NPDB.Pollen$ID1)
                       NPDB.MPD <- as.matrix(dist(NPDB.MPD.ALL[which(!is.na(NPDB.MPD.ALL))]))
                       #Sys.time()-SrtTime
                       #####
                       
                       #####
                       # Define the analogue threshold
                       # Define the Biome for the loaction to use
                       BiomeAll <- NPDB.Pollen$BIOME[which(!is.na(NPDB.MPD.ALL))]
                       # Biomes with more than two sites
                       BiomeUse <- names(table(BiomeAll))[table(BiomeAll)>2]
                       

                       # Estimate the distance between NPDB sites Removing non biome sites and biomes with only one observation
                       TraitMPDDist.ROC <- roc(NPDB.MPD[c(BiomeAll%in%BiomeUse),c(BiomeAll%in%BiomeUse)], # current time Taxon data.frame 
                                               groups = BiomeAll[c(BiomeAll%in%BiomeUse)] # vector of group memberships
                                               )
                       #####

                       #####
                       # Estimate distance of paleo assemblage to the NPDB
                       #SrtTime<-Sys.time()
                       NEATOMA.MPD <- mpd(NEATOMA.Agg[,TaxaUse],
                                                 TraitsDist,
                                                 abundance.weighted=T)
                       names(NEATOMA.MPD) <- paste0(NEATOMA.Agg$sites,"_",NEATOMA.Agg$Time)
                       NEATOMA.minMPD <- as.matrix(dist(c(NEATOMA.MPD[!is.na(NEATOMA.MPD)],NPDB.MPD.ALL[!is.na(NPDB.MPD.ALL)])))
                       NEATOMA.minMPD <- NEATOMA.minMPD[colnames(NEATOMA.minMPD)%in%names(NEATOMA.MPD),
                                                        colnames(NEATOMA.minMPD)%in%names(NPDB.MPD.ALL)]
                       NEATOMA.minMPD <- apply(NEATOMA.minMPD,1,min)
                       NEATOMA.Agg$minSCDToPres <- NA # Add the Min distance to the Agg NEATOMA DATA
                       NEATOMA.Agg$minSCDToPres[paste0(NEATOMA.Agg$sites,"_",NEATOMA.Agg$Time)%in%names(NEATOMA.minMPD)] <- NEATOMA.minMPD
                       #Sys.time()-SrtTime
                       ####
                       ####
                       # Summary fo distances per Time 
                       Out.list <- list(NEATOMA.CommDis = NEATOMA.Agg[,-c(1,7:45)],
                                        ROC.Cutoff = TraitMPDDist.ROC)
                       Sys.time()-SrtTime
                       return(Out.list)
                     })

saveRDS(BootTraits,"./Results/Euclidean_TraitMPDis.rds")



# Load and plot the MPD for traits
rm(list=ls());gc()
require(snowfall)
setwd("/Users/alejandroordonez/Documents/3. Ecography PaleoTraits Paper/PaleoNovelty/")
#setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
BootTraits <- readRDS("./Results/Euclidean_TraitMPDis.rds")

# Estimate the Variability in estimates per Trait space iteration 

TraitsMPDVar <- lapply(BootTraits,
                       function(x){#x<-BootTraits[[3]]
                         TraitsMPD <- x$NEATOMA.CommDis
                         MeanMPD <- lapply(1:1000,
                                           function(j){
                                           tapply(TraitsMPD$minSCDToPres[do.call("c",lapply(21:1,function(i){sample(which(TraitsMPD$Time==i),10)}))],
                                                  TraitsMPD$Time[do.call("c",lapply(21:1,function(i){sample(which(TraitsMPD$Time==i),10)}))],
                                                  median,na.rm=T)})
                         MeanMPD <- do.call("cbind",
                                            MeanMPD)
                         Out <- data.frame(Time = 1:21,
                                           t(apply(MeanMPD,1,quantile, c(0.025,0.5,0.975))))
                         return(Out)
                       })
TraitsMPD <- data.frame(Time = 1:21,
                        X50 = apply(sapply(TraitsMPDVar,function(x){x$X50.}),1,median),
                        X2.5 = apply(sapply(TraitsMPDVar,function(x){x$X2.5.}),1,median),
                        X97.5 = apply(sapply(TraitsMPDVar,function(x){x$X97.5.}),1,median))

plot(y = rev(TraitsMPD$X50),
     x = -21:-1,
     type = "b",
     main = "Regional Dissimilarity Trait Change", 
     xlab ="Time (kyrBP)",
     ylab = "MPD (Euclidean Dist)",
     ylim = range(TraitsMPD[-1]))
lines(y = rev(TraitsMPD$X2.5),
      x = -21:-1)
lines(y = rev(TraitsMPD$X97.5),
      x = -21:-1)

BootTraits[[1]][[2]]$roc$ Combined$ optima
abline(h=median(sapply(BootTraits,
                       function(x){
                         x[[2]]$roc$ Combined$ optima
                       })))


################################################################################
################################################################################
# Estimate the Variability in estimates per Trait space iteration 
rm(list=ls());gc()
setwd("/Users/alejandroordonez/Documents/3. Ecography PaleoTraits Paper/PaleoNovelty/")
BootTraits <- readRDS("./Results/Euclidean_TraitMPDis.rds")

# Estimate the Variability in estimates per Trait space iteration 

TraitsMPDVar <- lapply(BootTraits,
                       function(x){#x<-BootTraits[[3]]
                         TraitsMPD <- x$NEATOMA.CommDis
                         MeanMPD <- lapply(1:1000,
                                           function(j){
                                             tapply(TraitsMPD$minSCDToPres[do.call("c",lapply(21:1,function(i){sample(which(TraitsMPD$Time==i),10)}))],
                                                    TraitsMPD$Time[do.call("c",lapply(21:1,function(i){sample(which(TraitsMPD$Time==i),10)}))],
                                                    median,na.rm=T)})
                         MeanMPD <- do.call("cbind",
                                            MeanMPD)
                         Out <- data.frame(Time = 1:21,
                                           t(apply(MeanMPD,1,quantile, c(0.025,0.5,0.975))))
                         return(Out)
                       })
require(tidyverse)
require(ggbreak)
TraitsMPD2<-data.frame(Time = -c(1:21,21:1),
                       Phylo = c(apply(sapply(TraitsMPDVar,function(x){x$X2.5.}),1,median),
                                 rev(apply(sapply(TraitsMPDVar,function(x){x$X97.5.}),1,median))))
TraitsMPD3 <- data.frame(Time = -(1:21),
                         Phylo = apply(sapply(TraitsMPDVar,function(x){x$X50.}),1,median,na.rm=T))

ggplot(TraitsMPD2, aes(x = Time, y = Phylo)) +
  geom_polygon(fill="lightblue") + 
  geom_line(data=TraitsMPD3,color="black") + 
  geom_point(data=TraitsMPD3,color="black") +
  ylim(0, 0.0025) +
  xlab("Time (kaBP)") + # for the x axis label
  ylab(" MPD Dif to present MPD\n[Euclidean Distance] ") +# for the y axis label
  ggtitle("Regional Functional difference") +
  geom_hline(yintercept=median(sapply(BootTraits,
                                      function(x){
                                        x[[2]]$roc$ Combined$ optima
                                      })),
             linetype="dashed", 
             color = "red") +
  scale_y_break(c(0.0005, 0.002),scales=0.25,expand=F)

#"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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"

