rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
require(snowfall)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
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
# Load the trait data
Traits <- read.csv2("./Data/Traits/Raw/TraitsSummary.csv")
Traits <- Traits[c(Traits$Name%in%TaxaUse),]
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

#turn traits into PCoA axes
GowdisDfrm <- FD::gowdis(TraitsDBS)
PCoATraits <- ape::pcoa(GowdisDfrm)
PCoATraitsEgVal <- PCoATraits $ vectors

# Which taxa to use - ensure match between NPDS and NEATOMA and TRAITS
TaxaUse <- row.names(PCoATraitsEgVal)
# Estimate the Trait distance in PCoA Space                       
TraitsDist <- as.matrix(dist(PCoATraitsEgVal))
#####


#######################################################################################################################
#####
# Load the North American Pollen Data Base and its translation to the Paleo data
# raw NPDB
SrtTime1 <- Sys.time()
BootTraits <- lapply(1:100,
                     function(RunIt){
                       #####
                       # Build a species names per site for the NPDS
                       SppPerSiteNPDB <- apply(NPDB.Pollen[,TaxaUse],
                                               1,
                                               function(x){TaxaUse[x>0]})
                       #####
                       #####
                       # Trait Mean pairwise distance contrast
                       sfInit( parallel=TRUE, cpus=10)
                       sfExport("SppPerSiteNPDB")
                       sfExport("TraitsDist")
                       SrtTime<-Sys.time()
                       MPDBaseline <- lapply(SppPerSiteNPDB,
                                             function(i){
                                               sapply(SppPerSiteNPDB,
                                                      function(x){
                                                        mean(TraitsDist[i,x])
                                                      })
                                             })
                       MPDBaseline <- do.call("rbind",MPDBaseline)
                       sfStop()
                       Sys.time()-SrtTime
                       #####
                       
                       #####
                       # Define the analogue threshold
                       # Estimate the distance between NPDB sites Removing non biome sites and biomes with only one observation
                       TraitMPDDist.ROC <- roc(MPDBaseline, # current time Taxon data.frame 
                                               groups = NPDB.Pollen$BIOME # vector of group memberships
                       )
                       TraitMPDDist.ROC
                       #####
                       
                       #####
                       # Build a species per site Data for the NEATOMA dataset
                       PollenSppSiteNames <- apply(NEATOMA.Agg[,TaxaUse],
                                                   1,
                                                   function(x){TaxaUse[x>0]
                                                   })
                       #####
                       
                       #####
                       # Estimate distance to each paleo assemblage
                       sfInit( parallel=TRUE, cpus=10)
                       sfExport("SppPerSiteNPDB")
                       sfExport("PollenSppSiteNames")
                       sfExport("TraitsDist")
                       SrtTime<-Sys.time()
                       PaleoToPresentTraitDist <- sfLapply(PollenSppSiteNames,
                                                           function(i){
                                                             out <- sapply(SppPerSiteNPDB,
                                                                           function(x){
                                                                             mean(TraitsDist[i,x],na.rm=T)
                                                                           })
                                                             min(out,na.rm=T)
                                                           })
                       NEATOMA.Agg$minSCDToPres <- do.call("c",PaleoToPresentTraitDist)
                       sfStop()
                       Sys.time()-SrtTime
                       ####
                       ####
                       # Summary fo distances per Time 
                       Out.list <- list(NEATOMA.CommDis = NEATOMA.Agg[,-c(1,7:45)],
                                        ROC.Cutoff = TraitMPDDist.ROC)
                       return(Out.list)
                     })
Sys.time()-SrtTime1
saveRDS(BootTraits,"./Results/BootPCoATraitsMPD.rds")

BootTraits <- readRDS("./Results/BootPCoATraitsMPD.rds")


TraitsMPD <- lapply(BootTraits,
                    function(x){
                      x[[1]][,"minSCDToPres"]
                    })

TraitsMPD <- data.frame(BootTraits[[1]][[1]][,1:5],
                        t(apply(do.call("cbind",TraitsMPD),
                                1,
                                quantile,
                                c(0.0275,0.5,0.975)))
)



# Dummy plot (taking a even sub sample of sites across periods)
CompDisBoot <- lapply(1:1000,
                      function(i){
                        SamplTmp <- do.call("c",lapply(1:21,
                                                       function(x){
                                                         sample(which(TraitsMPD$Time==x), 10)
                                                       }))
                        tapply(TraitsMPD$X50.[SamplTmp],
                               TraitsMPD$Time[SamplTmp],
                               median)
                      })

CompDisBoot2 <- do.call("rbind",CompDisBoot)
CompDisBootQuant <- apply(CompDisBoot2,2,quantile,c(0.0275,0.5,0.975))
plot(y = rev(CompDisBootQuant[2,]),
     x = -21:-1,
     type = "b",
     main = "Regional Dissimilarity Trait Change", 
     xlab ="Time (kyrBP)",
     ylab = "Disimilarity (Sqr Cord Dist)",
     ylim=range(CompDisBootQuant))
lines(y = rev(CompDisBootQuant[1,]),
      x = -21:-1,)
lines(y = rev(CompDisBootQuant[3,]),
      x = -21:-1,)


abline(h=median(sapply(BootTraits,
                       function(x){
                         x[[2]]$roc$ Combined$ optima
                       })))