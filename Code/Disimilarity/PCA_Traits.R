rm(list=ls());gc()
require(analogue)
require(FD)
require(terra)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/Data")
#######################################################################################################################
#####
# Load the North American Pollen Data Base and its translation to the Paleo data
# raw NPDB
BootTraits <- lapply(1:100,
                     function(RunIt){
                       NPDB <- read.csv("./Pollen/NorthAmerican Pollen DBS/whitmoreetal2005_v1-72.csv")
                       # Taxa translate between NPDB and Paleo Data
                       TaxtTrans <- read.csv("./Pollen/NorthAmerican Pollen DBS/Taxa transation.csv")
                       # Which are the summarizing taxa - note the TaxtTrans$With.Traits!="" filter -> only species with traits are used
                       TaxaUse <- sort(unique(TaxtTrans$Translated[TaxtTrans$With.Traits!=""]))
                       #TaxaUse <- sort(unique(TaxtTrans$Translated))[-1]
                       #####
                       
                       #####
                       # Build a translated NPDB dataset
                       SummNPDB <- lapply(TaxaUse, # only use species with trait data (Woody species)
                                          function(x){
                                            out <- NPDB[,TaxtTrans$Original[TaxtTrans$Translated==x]]
                                            if(class(out)=="data.frame"){out <- apply(out,1,sum)}
                                            out <- as.numeric(out)
                                            return(out)
                                          })
                       names(SummNPDB) <- TaxaUse
                       SummNPDB2 <- as.data.frame(do.call("cbind",SummNPDB))
                       SummNPDB3 <- as.data.frame(t(apply(SummNPDB2,1,function(x){x/sum(x,na.rm=T)})))
                       SummNPDB3[is.na(SummNPDB3)] <- 0
                       #####
                       
                       #####
                       # Load the trait data
                       Traits <- read.csv("./Traits/TraitsSummary.csv")
                       Traits$TaxtTrans <- TaxtTrans$Translated[match(Traits$Genus,TaxtTrans$PLANTS)]
                       # Seeds 
                       SWT <- apply(Traits[,c("Min..of.Seed.Size..mg.","Max..of.Seed.Size..mg.")],
                                    1,
                                    function(x){
                                      runif(1,x[1],x[2])})
                       SWT <- scale(log10(SWT))
                       
                       # Height
                       Hmax <- apply(Traits[,c("Min..of.Height..m.","Max..of.Height..m.")],
                                     1,
                                     function(x){
                                       runif(1,x[1],x[2])})
                       Hmax <- scale(Hmax)
                       # LMA
                       LMA <- apply(Traits[,c("Min..of.LMA..log...g.m2.","Max..of.LMA..log...g.m2.")],
                                    1,
                                    function(x){
                                      runif(1,x[1],x[2])})
                       LMA <- scale(LMA)
                       
                       # Build a trait dataset
                       TraitsDBS <- data.frame(SWT = SWT,
                                               Hmax = Hmax,
                                               LMA = LMA)
                       row.names(TraitsDBS) <- Traits$TaxtTrans
                       TraitsDBS <- TraitsDBS[complete.cases(TraitsDBS),]

                       #trun traits into PCoA axes
                       GowdisDfrm <- gowdis(TraitsDBS)
                       PCoATraits <- ape::pcoa(GowdisDfrm)
                       PCoATraits <- PCoATraits $ vectors
                       # b <- vegan::wcmdscale(GowdisDfrm,
                       #                       k=15,
                       #                       eig = TRUE)                       
                       
                       # Build a occurence dataset
                       SummNPDB4 <- SummNPDB3[,rownames(TraitsDBS)]
                       SummNPDB4 <- SummNPDB4[apply(SummNPDB4,1,sum)!=0,]
                       SummNPDB4 <- SummNPDB4[c(apply(SummNPDB4,1,function(x){sum(x!=0)})>3),]
                       SummNPDB4 <- as.data.frame(t(apply(SummNPDB4,1,function(x){x/sum(x,na.rm=T)})))
                       
                       # Trait space per NPDB Site
                       NPDBTrait <- lapply(1:dim(SummNPDB4)[1],
                                           function(i){
                                             PCoATraits[names(SummNPDB4)[SummNPDB4[i,]!=0],]
                                           })
                       
                       # Community mean per NBP
                       # Raw mean
                       # NPDBTraitMean <- lapply(NPDBTrait,
                       #                         function(x){
                       #                           apply(x,2,mean)
                       #                         })
                       # NPDBTraitMean <- do.call("rbind",NPDBTraitMean)
                       
                       # Community Weighted Mean
                       NPDBTraitMean <- lapply(1:length(NPDBTrait),
                                               function(x){#x<-1
                                                 apply(NPDBTrait[[x]],
                                                       2,
                                                       function(y){
                                                         weighted.mean(x = y,
                                                                       SummNPDB4[x,rownames(NPDBTrait[[x]])])
                                                       })
                                               })
                       NPDBTraitMean <- do.call("rbind",NPDBTraitMean)
                       #####
                       # Define the analogue threshold
                       BIOMES <- vect("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/Data/WWF-Biomes/wwf_terr_ecos.shp")
                       BIomesNames <- data.frame(code = 1:14,
                                                 Name = c("Tropical and Subtropical Moist Broadleaf Forests",
                                                          "Tropical and Subtropical Dry Broadleaf Forests",
                                                          "Tropical and Subtropical Coniferous Forests",
                                                          "Temperate Broadleaf and Mixed Forests",
                                                          "Temperate Coniferous Forests",
                                                          "Boreal forests/Taiga",
                                                          "Tropical and Subtropical Grasslands, Savannas, and Shrublands",
                                                          "Temperate Grasslands, Savannas, and Shrublands",
                                                          "Flooded Grasslands and Savannas",
                                                          "Montane Grasslands and Shrublands",
                                                          "Tundra",
                                                          "Mediterranean Forests, Woodlands, and Scrub",
                                                          "Deserts and Xeric Shrublands",
                                                          "Mangroves"))
                       # Get the BIOME type for each NPDB location
                       NPDBLoc <- NPDB[rownames(SummNPDB4),c("LONDD","LATDD")]
                       names(NPDBLoc) <- c('x',"y")
                       NPDBBiome <- extract(BIOMES[,"BIOME"],
                                            NPDBLoc)
                       # Make NA values a category
                       NPDBBiome[is.na(NPDBBiome[,2]),2] <- 100
                       # Add biome to Traits
                       NPDBTraitMean <- data.frame(NPDBTraitMean,
                                                   Biome = NPDBBiome[,2])
                       # keep only complete cases
                       NPDBTraitMean <- NPDBTraitMean[complete.cases(NPDBTraitMean),]
                       # Remove non Biomes
                       NPDBTraitMean <- NPDBTraitMean[NPDBTraitMean$Biome<15,]
                       # Remove Biomes with only one observation
                       NPDBTraitMean<- NPDBTraitMean[c(!NPDBTraitMean$Biome==names(which(table(NPDBTraitMean$Biome)<3))),]
                       
                       # Estimate the distance between NPDB sites
                       TraitEuclDist <- as.matrix(dist(NPDBTraitMean[,-4]))
                       TraitEuclDist.ROC <- roc(TraitEuclDist, # current time Taxon data.frame 
                                                groups = BIomesNames$Name[NPDBTraitMean$Biome] # vector of group memberships
                       )
                       TraitEuclDist.ROC
                       #####
                       
                       #####
                       # Load the Paleo-Pollen DBS
                       Old <- read.csv("./Pollen/Ordonez & Williams 2013/abund_spp_env_21to1bp.csv")
                       Old2 <- Old[,rownames(PCoATraits)]
                       Old2 <- as.data.frame(t(apply(Old2,1,function(x){x/sum(x,na.rm=T)})))
                       Old2 <- data.frame(Old[, c("Time","Clim.ID","sites","Latitude","Longitude")],
                                          Old2)
                       Old2 <- Old2[c(apply(Old2[,rownames(PCoATraits)],1,function(x){sum(x!=0)})>3),]
                       
                       
                       
                       #####
                       # Estimate distance to each paleo assemablage
                       
                       ## Trait summary per paleo site 
                       TraitSummPaleoSites <- lapply(1:dim(Old2)[1],
                                                     function(i){
                                                       PCoATraits[row.names(PCoATraits)%in%names(SummNPDB4)[Old2[i,names(SummNPDB4)]!=0],]
                                                     })
                       # Distance of each paleo site to all sites in the present - raw mean
                       # PaleoToPresentTraitDist <- lapply(1:21,
                       #                                   function(x){
                       #                                       DistTmpList <- lapply(TraitSummPaleoSites[which(Old2$Time%in%x)],
                       #                                                             function(y){
                       #                                                               apply(y,2,mean,na.rm=T)})
                       #                                       
                       #                                       DistTmp1 <- do.call("rbind",DistTmpList)
                       #                                       DistTmp1 <- na.omit(DistTmp1)
                       #                                       DistTmp <- analogue::distance(DistTmp1,
                       #                                                                     NPDBTraitMean[,-4])
                       #                                       median(apply(DistTmp,1,min))
                       #                                     })
                       # Distance of each paleo site to all sites in the present - weighted mean
                       PaleoToPresentTraitDist <- lapply(1:21,
                                                         function(x){
                                                           DistTmpList <- lapply(which(Old2$Time%in%x),
                                                                                 function(y){#y<-3935
                                                                                   apply(TraitSummPaleoSites[[y]],
                                                                                         2,
                                                                                         function(i){
                                                                                           weighted.mean(as.numeric(i),
                                                                                                         as.numeric(Old2[y,rownames(TraitSummPaleoSites[[y]])]))
                                                                                           })
                                                                                   })
                                                           DistTmp1 <- do.call("rbind",DistTmpList)
                                                           DistTmp1 <- na.omit(DistTmp1)
                                                           DistTmp <- analogue::distance(DistTmp1,
                                                                                         NPDBTraitMean[,-dim(NPDBTraitMean)[2]])
                                                           median(apply(DistTmp,1,min))
                                                         })
                       Out.list <- list(DataFrm = data.frame(Time=-21:-1,
                                                             TraitDist = rev(do.call("c",PaleoToPresentTraitDist))),
                                        RocTresh = TraitEuclDist.ROC)
                       return(Out.list)
                     })
saveRDS(BootTraits,"./Traits/BootPCoATraitsCWM.rds")

BootTraitsPCoA <- readRDS("./Traits/BootPCoATraitsCWM.rds")

BootTraitsPCoAONLy <- lapply(BootTraitsPCoA,
                            function(x){x[[1]][,2]})


plot(y=apply(do.call("rbind",BootTraitsPCoAONLy),2,median),
     x=-21:-1,
     type="b",
     ylab = "Eucludean Distance raw-CWM",
     xlab = "Time (kyrBP)")

BootTraitsCutoff <- sapply(BootTraitsPCoA,
                           function(x){
                             x[[2]]$roc$Combined$optimal
                           })
abline(h=median(BootTraitsCutoff))


saveRDS(BootTraits,"./Traits/BootPCoATraitsCWM.rds")

BootTraitsCWM <- readRDS("./Traits/BootTraitsCWM.rds")

BootTraitsCWMONLy <- lapply(BootTraitsCWM,
                            function(x){x[[1]][,2]})


plot(y=apply(do.call("rbind",BootTraitsCWMONLy),2,median),
     x=-21:-1,
     type="b",
     ylab = "Eucludean Distance raw-CWM",
     xlab = "Time (kyrBP)")

BootTraitsCutoff <- sapply(BootTraitsCWM,
                           function(x){
                             x[[2]]$roc$Combined$optimal
                           })
abline(h=median(BootTraitsCutoff))






