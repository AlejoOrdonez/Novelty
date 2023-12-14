rm(list=ls());gc()
require(ape)
require(snowfall)
require(analogue)
require(picante)
#setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
setwd("/Users/alejandroordonez/Documents/3. Ecography PaleoTraits Paper/PaleoNovelty/")
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

#####
# Community data NPDB and NEATOMA 
# Load the Aggregated North American Pollen Data Base
NPDB.Pollen <- read.csv("./Data/Pollen/Processed/NPDB_Agg.csv")
# Load the climate for the North American Pollen Data Base
NPDB.Clim <- read.csv("./Data/PaleoClimate/Processed/NPDB_Clim.csv")
# Select only sites with Pollen data and Climate data
NPDB.SiteUse <- table(c(NPDB.Clim$ID1,NPDB.Pollen$ID1))
NPDB.SiteUse <- as.numeric(names(NPDB.SiteUse)[(NPDB.SiteUse==2)])
# Aggregated composition
NPDB.Pollen <- NPDB.Pollen[c(NPDB.Pollen$ID1%in%NPDB.SiteUse),]
# Load the Aggregated NEATOMA Data Base
NEATOMA.Agg <- read.csv("./Data/Pollen/Processed/NEATOMA_Agg.csv")
#####

#####
# List of Spp to use in all three data sets
TaxaUse <- table(c(names(NEATOMA.Agg)[-c(1:7)],
                   names(NPDB.Pollen)[-c(1:10)],
                   unique(SppTrans$Translation)))
TaxaUse <- names(TaxaUse)[c(TaxaUse==3)]
#####


BootPhylo <- lapply(1:10,
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
                #####
                # Phylogenetic Mean pairwise distance contrast
                NPDB.MPD.ALL <- mpd(NPDB.Pollen[,TaxaUse],
                                    PhyloDistTmp,
                                    abundance.weighted=T)
                names(NPDB.MPD.ALL) <- paste0(NPDB.Pollen$DBCODE,"_",NPDB.Pollen$ID1)
                NPDB.MPD <- as.matrix(dist(NPDB.MPD.ALL[which(!is.na(NPDB.MPD.ALL))]))
                #####
                # Define the analogue threshold
                # Define the Biome for the loaction to use
                BiomeAll <- NPDB.Pollen$BIOME[which(!is.na(NPDB.MPD.ALL))]
                # Biomes with more than two sites
                BiomeUse <- names(table(BiomeAll))[table(BiomeAll)>2]
                
                
                # Estimate the distance between NPDB sites Removing non biome sites and biomes with only one observation
                PhyloMPDDist.ROC <- roc(NPDB.MPD[c(BiomeAll%in%BiomeUse),c(BiomeAll%in%BiomeUse)], # current time Taxon data.frame 
                                        groups = BiomeAll[c(BiomeAll%in%BiomeUse)] # vector of group memberships
                )
                #####
                #####
                # Estimate distance of paleo assemblage to the NPDB
                NEATOMA.MPD <- mpd(NEATOMA.Agg[,TaxaUse],
                                   PhyloDistTmp,
                                   abundance.weighted=T)
                names(NEATOMA.MPD) <- paste0(NEATOMA.Agg$sites,"_",NEATOMA.Agg$Time)
                NEATOMA.minMPD <- as.matrix(dist(c(NEATOMA.MPD[!is.na(NEATOMA.MPD)],NPDB.MPD.ALL[!is.na(NPDB.MPD.ALL)])))
                NEATOMA.minMPD <- NEATOMA.minMPD[colnames(NEATOMA.minMPD)%in%names(NEATOMA.MPD),
                                                 colnames(NEATOMA.minMPD)%in%names(NPDB.MPD.ALL)]
                NEATOMA.minMPD <- apply(NEATOMA.minMPD,1,min,na.rm=T)
                NEATOMA.Agg2<-NEATOMA.Agg
                NEATOMA.Agg2$minSCDToPres <- NA # Add the Min distance to the Agg NEATOMA DATA
                NEATOMA.Agg2$minSCDToPres[paste0(NEATOMA.Agg2$sites,"_",NEATOMA.Agg2$Time)%in%names(NEATOMA.minMPD)] <- NEATOMA.minMPD
                ####
                # Summary fo distances per Time 
                Out.list <- list(NEATOMA.CommDis = NEATOMA.Agg2[,-c(1,7:45)],
                                 ROC.Cutoff = PhyloMPDDist.ROC)
                Sys.time() -SrtTime1
                return(Out.list)
              })
saveRDS(BootPhylo,"./Results/Cophe_PhyloMPDis.rds")

################################################################################
################################################################################
# Estimate the Variability in estimates per Phylo space iteration 
rm(list=ls());gc()
setwd("/Users/alejandroordonez/Documents/3. Ecography PaleoTraits Paper/PaleoNovelty/")
BootPhylo <- readRDS("./Results/Cophe_PhyloMPDis.rds")
require(snowfall)
sfInit( parallel=TRUE, cpus=10)
sfExport("BootPhylo")

PhyloMPDVar <- sfLapply(BootPhylo,
                        function(x){#x<-BootPhylo[[3]]
                          PhyloMPD <- x$NEATOMA.CommDis
                          MeanMPD <- lapply(1:1000,
                                            function(j){
                                              tapply(PhyloMPD$minSCDToPres[do.call("c",lapply(21:1,function(i){sample(which(PhyloMPD$Time==i),10)}))],
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

ggplot(PhyloMPD2, aes(x = Time, y = Phylo)) +
  geom_polygon(fill="lightblue") + 
  geom_line(data=PhyloMPD3,color="black") + 
  geom_point(data=PhyloMPD3,color="black") +
  ylim(0, 0.385) +
  xlab("Time (kaBP)") + # for the x axis label
  ylab(" MPD Dif to present MPD\n[Cophenetic Distance] ") +# for the y axis label
  ggtitle("Regional Phylogentic difference") +
  geom_hline(yintercept=median(sapply(BootPhylo,
                                      function(x){
                                        x[[2]]$roc$ Combined$ optima
                                      })),
             linetype="dashed", 
             color = "red") +
  scale_y_break(c(0.12, 0.375 ),scales=0.25,expand=F)

#"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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"
