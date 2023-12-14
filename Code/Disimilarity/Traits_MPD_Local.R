rm(list=ls());gc()
require(analogue)
require(snowfall)
setwd("/Users/alejandroordonez/Documents/3. Ecography PaleoTraits Paper/PaleoNovelty/")
#####
# Load the Aggregated NEATOMA Data Base
NEATOMA.Comm <- read.csv("./Data/Pollen/Processed/NEATOMA_Agg.csv")
# Filter sites that not take 1 ka BP and have less than 6 dates
SummSyteXTime <- table(NEATOMA.Comm$sites,NEATOMA.Comm$Time)
SummSyteXTime <- SummSyteXTime[c(apply(SummSyteXTime,1,sum)>6 & SummSyteXTime[,1]==1),]

# Estimate functional distance for each assemblage to conditions 1 kaBP
sfInit( parallel=TRUE, cpus=10)
sfExport("NEATOMA.Comm")
sfExport("SummSyteXTime")
SrtTime1 <- Sys.time()
BootTraits <- sfLapply(1:100,
                     function(RunIt){
                       #####
                       # Generate a possible trait realization based on the know distribution of values for a evaluated taxa
                       # Load the trait data
                       Traits <- read.csv2("./Data/Traits/Raw/TraitsSummary.csv")
                       Traits <- Traits[c(Traits$Name%in%names(NEATOMA.Comm)[-c(1:7)]),]
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
                       
                       #####
                       # Build a trait dataset
                       TraitsDBS <- data.frame(SWT = SWT,
                                               Hmax = Hmax,
                                               LMA = LMA)
                       row.names(TraitsDBS) <- Traits$Name
                       TraitsDBS <- TraitsDBS[complete.cases(TraitsDBS),]
                       
                       #####
                       # Estimate the Trait distance in RAW Space                       
                       TraitsDist <- as.matrix(dist(TraitsDBS))
                       
                       #####
                       # Remove spp with no trait data from the occurrence dataset
                       SppUse <- names(NEATOMA.Comm)[-c(1:7)][names(NEATOMA.Comm)[-c(1:7)]%in%row.names(TraitsDBS)]
                       
                       #####
                       # Estimate functional dissimilarity - Squared Coord distance
                       LocalTraitDisList <- lapply(rownames(SummSyteXTime),
                                                   function(Site){#(Site <- rownames(SummSyteXTime)[2])
                                                     CommDataTmp <- NEATOMA.Comm[NEATOMA.Comm$sites==Site,][,SppUse]
                                                     SppList <- apply(CommDataTmp,1,function(x){SppUse[x!=0]})
                                                     names(SppList) <- NEATOMA.Comm$Time[NEATOMA.Comm$sites == Site]
                                                     # Estimate the Euclidean distance of Scaled traits
                                                     SiteEucDist <- sapply(SppList,
                                                                           function(x){#(x<-SppList[[1]])
                                                                             out <- sapply(SppList,
                                                                                           function(y){#(y<-SppList[[2]])
                                                                                             a <- abs(mean(TraitsDist[x,y])-mean(TraitsDist[x,x]))
                                                                                             #a<-TraitsDist[x,y]
                                                                                             #a <- mean(a)
                                                                                             ifelse(is.na(a),0,a)
                                                                                           })
                                                                             return(out)
                                                                           })
                                                     # Summary of Euclidean distance of Scaled traits
                                                     SiteEucDis <- data.frame(Site=Site,
                                                                              Time = NEATOMA.Comm$Time[NEATOMA.Comm$sites == Site],
                                                                              EucDTrait= SiteEucDist[dim(SiteEucDist)[1],])
                                                     return(SiteEucDis)
                                                   })
                       # Turn list into a table
                       
                       LocalTraitDis <- do.call("rbind",LocalTraitDisList)
                       return(LocalTraitDis)
                     })

saveRDS(BootTraits,"./Results/Euclidean_TraitMPDis_Local.rds")
Sys.time() - SrtTime1

# Plot (taking a even sub sample of sites across periods)
rm(list=ls()); gc()
require(snowfall)
BootTraits <- readRDS("./Results/Euclidean_TraitMPDis_Local.rds")

# Estimate the Variability in estimates per Trait space iteration 
sfInit( parallel=TRUE, cpus=10)
sfExport("BootTraits")
TraitsMPDVar <- sfLapply(BootTraits,
                         function(TraitsMPD){#TraitsMPD<-BootTraits[[3]]
                           # Estimate the mean dissimilarity per-time period based on 10 samples
                           MeanMPD <- lapply(1:1000,
                                             function(j){
                                               tapply(TraitsMPD$EucDTrait[do.call("c",lapply(21:1,function(i){sample(which(TraitsMPD$Time==i),10)}))],
                                                      TraitsMPD$Time[do.call("c",lapply(21:1,function(i){sample(which(TraitsMPD$Time==i),10)}))],
                                                      median)})
                           MeanMPD <- do.call("cbind",
                                              MeanMPD)
                           Out <- data.frame(Time = 1:21,
                                             t(apply(MeanMPD,1,quantile, c(0.025,0.5,0.975),na.rm=T)))
                           return(Out)
                         })
sfStop()

# Summary of the dissimilarity per-time period based on 10 samples across 100 Trait iterations
require(tidyverse)
require(ggbreak)
TraitsMPD2<-data.frame(Time = -c(1:21,21:1),
                       Phylo = c(apply(sapply(TraitsMPDVar,function(x){x$X2.5.}),1,median),
                                 rev(apply(sapply(TraitsMPDVar,function(x){x$X97.5.}),1,median))))
TraitsMPD3 <- data.frame(Time = -(1:21),
                         Phylo = apply(sapply(TraitsMPDVar,function(x){x$X50.}),1,median,na.rm=T))
BootTraits <- readRDS("./Results/Euclidean_TraitMPDis.rds")
ggplot(TraitsMPD2, aes(x = Time, y = Phylo)) +
  geom_polygon(fill="lightblue") + 
  geom_line(data=TraitsMPD3,color="black") + 
  geom_point(data=TraitsMPD3,color="black") +
  #ylim(0, 0.0025) +
  xlab("Time (kaBP)") + # for the x axis label
  ylab(" MPD Dif to present MPD\n[Euclidean Distance] ") +# for the y axis label
  ggtitle("Local Functional difference") +
  geom_hline(yintercept=median(sapply(BootTraits,
                                      function(x){
                                        x[[2]]$roc$ Combined$ optima
                                      })),
             linetype="dashed", 
             color = "red")

TraitDisBoot <- data.frame(Time = 1:21,
                           X50 = apply(sapply(TraitsMPDVar,function(x){x$X50.}),1,median),
                           X2.5 = apply(sapply(TraitsMPDVar,function(x){x$X2.5.}),1,median),
                           X97.5 = apply(sapply(TraitsMPDVar,function(x){x$X97.5.}),1,median))


plot(y = rev(TraitDisBoot[,"X50"]),
     x = -21:-1,
     type = "b",
     main = "Local Dissimilarity Trait change", 
     xlab ="Time (kyrBP)",
     ylab = "Disimularity (Euclidean Dist)",
     ylim=range(TraitDisBoot[,-1]))
lines(y = rev(TraitDisBoot[,"X2.5"]),
      x = -21:-1)
lines(y = rev(TraitDisBoot[,"X97.5"]),
      x = -21:-1)

NEATOMA.TriatDis <- readRDS("./Results/BootTraitsMPD.rds")
abline(h=NEATOMA.TriatDis[[1]]$ROC.Cutoff$roc$ Combined$ optima)


