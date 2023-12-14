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

#####
# Super tree data
#Load a super tree
MegaTree <- read.tree("./Data/Phylogeny/ALLMB.tre")
# Spp Translation
SppTrans <- read.csv2("./Data/Phylogeny/SppUse_Pylo_PLANTS.csv")
# Remove taxa not in the Magatree
SppTrans <- SppTrans[c(SppTrans$name.for.phylo%in%MegaTree$tip.label),]
## Prune the Phylogeny
SppPhylo <- keep.tip(MegaTree,MegaTree$tip.label[c(MegaTree$tip.label%in%SppTrans$name.for.phylo)])
# Estimate the Pairwise distance between species
PhyloDist <- cophenetic(SppPhylo)
#####
# Define the taxa to use
TaxaUse <- table(c(unique(SppTrans$Translation),names(NEATOMA.Comm)[-c(1:6)]))
TaxaUse <- names(TaxaUse)[TaxaUse==2]
sfInit( parallel=TRUE, cpus=10)
sfExport("NEATOMA.Comm")
sfExport("SummSyteXTime")
sfExport("PhyloDist")
sfExport("SppTrans")
sfExport("TaxaUse")
sfLibrary(picante)
SrtTime1 <- Sys.time()
BootPhylo <- sfLapply(1:100,
                    function(RunIt){
                      SrtTime1 <- Sys.time()
                      #####
                      # Random selection of pylogeny tips
                      SppTrasTbl<- SppTrans[SppTrans$Translation%in%TaxaUse,] # Ensure the matching spp are used
                      SppTrasTbl <- table(SppTrasTbl$name.for.phylo,SppTrasTbl$Translation)
                      SppTrasTbl <- apply(SppTrasTbl,2,function(x){sample(rownames(SppTrasTbl)[x==1],1)})
                      # Generate a sample soecific phylogenetic distance matrix
                      PhyloDistTmp <- PhyloDist[as.character(SppTrasTbl),as.character(SppTrasTbl)]
                      dimnames(PhyloDistTmp)[[1]] <- dimnames(PhyloDistTmp)[[2]] <- names(SppTrasTbl)
                      #####
                      # Estimate MPD distance of each paleo assemblage X 1
                      NEATOMA.MPD <- NEATOMA.Comm[,c("sites","Time")]
                      NEATOMA.MPD <- NEATOMA.MPD[NEATOMA.MPD$sites%in%rownames(SummSyteXTime),]
                      NEATOMA.MPD$minMPD <- mpd(NEATOMA.Comm[NEATOMA.Comm$sites%in%rownames(SummSyteXTime),TaxaUse],
                                                PhyloDistTmp,
                                                abundance.weighted=T)
                      #####
                      # # Estimate MPD distance between times in each paleo assemblage to 
                      SiteCophDistList <- lapply(unique(NEATOMA.MPD$sites),
                                            function(Site){#(Site<-unique(NEATOMA.MPD$sites)[1])
                                              DistSumm <- as.matrix(dist(NEATOMA.MPD$minMPD[NEATOMA.MPD$sites==Site]))
                                              dimnames(DistSumm) <- list(NEATOMA.MPD$Time[NEATOMA.MPD$sites==Site],
                                                                         NEATOMA.MPD$Time[NEATOMA.MPD$sites==Site])
                                              out <- data.frame(sites = NEATOMA.MPD$sites[NEATOMA.MPD$sites==Site],
                                                                Time = NEATOMA.MPD$Time[NEATOMA.MPD$sites==Site],
                                                                CophDist = DistSumm[,"1"])
                                              return(out)
                                              })
                      SiteCophDist <- do.call("rbind",SiteCophDistList)
                      return(SiteCophDist)
                      Sys.time() -SrtTime1
                      return(Out.list)
                    })
sfStop()
saveRDS(BootPhylo,"./Results/Cophe_PhyloMPDis_Local.rds")

# Estimate the Variability in estimates per Trait space iteration 
rm(list=ls());gc()
BootPhylo <- readRDS("./Results/Cophe_PhyloMPDis_Local.rds")

sfInit( parallel=TRUE, cpus=10)
sfExport("BootPhylo")

PhyloMPDVar <- sfLapply(BootPhylo,
                        function(x){#(x<-BootPhylo[[3]])
                          PhyloMPD <- x
                          MeanMPD <- lapply(1:1000,
                                            function(j){
                                              tapply(PhyloMPD$CophDist [do.call("c",lapply(21:1,function(i){sample(which(PhyloMPD$Time==i),10)}))],
                                                     PhyloMPD$Time[do.call("c",lapply(21:1,function(i){sample(which(PhyloMPD$Time==i),10)}))],
                                                     median,na.rm=T)})
                          MeanMPD <- do.call("cbind",
                                             MeanMPD)
                          Out <- data.frame(Time = 1:21,
                                            t(apply(MeanMPD,1,quantile, c(0.025,0.5,0.975))))
                          return(Out)
                        })
sfStop()

require(tidyverse)
require(ggbreak)
PhyloMPD2<-data.frame(Time = -c(1:21,21:1),
                      Phylo = c(apply(sapply(PhyloMPDVar,function(x){x$X2.5.}),1,median),
                                rev(apply(sapply(PhyloMPDVar,function(x){x$X97.5.}),1,median))))
PhyloMPD3 <- data.frame(Time = -(1:21),
                        Phylo = apply(sapply(PhyloMPDVar,function(x){x$X50.}),1,median,na.rm=T))

BootPhylo <- readRDS("./Results/Cophe_PhyloMPDis.rds")
ggplot(PhyloMPD2, aes(x = Time, y = Phylo)) +
  geom_polygon(fill="lightblue") + 
  geom_line(data=PhyloMPD3,color="black") + 
  geom_point(data=PhyloMPD3,color="black") +
  xlab("Time (kaBP)") + # for the x axis label
  ylab(" MPD Dif to present MPD\n[Cophenetic Distance] ") +# for the y axis label
  ggtitle("Local Phylogentic difference") +
  geom_hline(yintercept=median(sapply(BootPhylo,
                                      function(x){
                                        x[[2]]$roc$ Combined$ optima
                                      })),
             linetype="dashed", 
             color = "red")



PhyloMPD <- data.frame(Time = 1:21,
                       X50 = apply(sapply(PhyloMPDVar,function(x){x$X50.}),1,median,na.rm=T),
                       X2.5 = apply(sapply(PhyloMPDVar,function(x){x$X2.5.}),1,median),
                       X97.5 = apply(sapply(PhyloMPDVar,function(x){x$X97.5.}),1,median))

plot(y = rev(PhyloMPD$X50),
     x = -21:-1,
     type = "b",
     main = "Local Dissimilarity Phylogentic Change", 
     xlab ="Time (kyrBP)",
     ylab = "MPD (Cophenetic Dist)",
     ylim = range(PhyloMPD[-1]))
BootPhyloReg <- readRDS("./Results/Cophe_PhyloMPDis.rds")

abline(h=median(sapply(BootPhyloReg,
                       function(x){
                         x[[2]]$roc$ Combined$ optima
                         })))

#"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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"

