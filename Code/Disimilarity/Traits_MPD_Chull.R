rm(list=ls());gc()
require(geometry)
require(snowfall)
require(analogue)
#require(raster)
#require(terra)


setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
#####
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
# Which taxa to use - ensure match between NPDS and NEATOMA and TRAITS
TaxaUse <- table(c(TaxaUse,row.names(TraitsDBS)))
TaxaUse <- names(TaxaUse)[TaxaUse==2]
#####

#####
# Define a 3D functional space for NPDB sites with more than 4 spp
# Which sites in the NPDB have more than 4 spp?
SiteNPDBUse <- which(apply(NPDB.Pollen[,TaxaUse],
                           1,
                           function(x){sum(x!=0)})>4)

# Determine the functional space convex hull vertexesNPDB
ChullNPDBVrtx <- apply(NPDB.Pollen[SiteNPDBUse,TaxaUse], 
                   1,
                   function(x){#x <- NPDB.Pollen[1,TaxaUse]
                     # Species in a site
                     SppSiteUse <- TaxaUse[x!=0]
                     # Traits of those species
                     TraitsSite <- TraitsDBS[rownames(TraitsDBS)%in%SppSiteUse,]
                     # Estimate the convex hull
                     #TraitsSite[convhulln(TraitsSite),]# Using convhulln of the geometry package (the space is set based on triangular facets)
                     #TraitsSite[chull(TraitsSite),] # Using the chull method
                     return(TraitsSite)
                   })
#####
#####
# Estimate the intersection between NPDB convex-hulls pairs
sfInit( parallel=TRUE, cpus=10)
sfExport("ChullNPDBVrtx")
SrtTime<-Sys.time()
ChullDistList <- sfLapply(1:length(ChullNPDBVrtx),
                        function(x){#x<-1
                          DistOut <- lapply(ChullNPDBVrtx,
                                            function(y){#y<-ChullNPDBVrtx[[2]]
                                              try(
                                                # Get intersection convex hulls - this function needs as input the full set of points
                                                out <- geometry::intersectn(ps1 = ChullNPDBVrtx[[x]],
                                                                          ps2 = y) 
                                                                ,silent=T)
                                              if("out"%in%ls()){
                                                ProOverLp <- out$ch$vol/out$ch2$vol #define the proportion of the base convex hull the second convex hull overlaps
                                                PairDistOut <- 1 - ProOverLp # Make the proportion overlap a dissimilarity index
                                              }
                                              else{
                                                PairDistOut<-NA
                                              }
                                              return(round(PairDistOut,5))
                                            })
                          DistFullContrast <- do.call("c",DistOut) # Turn the dissimilarity list into a vector
                          #write.csv(DistFullContrast,paste0("./temp/",x,".csv"))
                          #return(x)
                          return(DistFullContrast)
                        })
Sys.time()-SrtTime
sfStop()
# ChullDistList <- lapply(dir("./temp"),
#                         function(i){
#                           as.numeric(read.csv(paste0("./temp/",i))[,2])
#                         })
ChullDist <- do.call("cbind",ChullDistList)
ChullDist[is.na(ChullDist)] <-1 
ChullDist <- as.data.frame(ChullDist)
dimnamesChullDist <- NPDB.Pollen$SITECODE[SiteNPDBUse]

#####

#####
# Define the analogue threshold
# Estimate the distance between NPDB sites Removing non biome sites and biomes with only one observation
TraitChullMPDDist.ROC <- roc(ChullDist, # current time Taxon data.frame 
                             groups = NPDB.Pollen$BIOME[SiteNPDBUse] # vector of group memberships
                             )
TraitChullMPDDist.ROC
#####


#####
# Define a 3D functional space for sites with more than 4 spp
# Which sites in the Neatoma sites have more than 4 spp?
SiteNEATOMAUse <- which(apply(NEATOMA.Agg[,TaxaUse],
                           1,
                           function(x){sum(x!=0)})>4)

# Determine the functional space convex hull vertexesNPDB
NEATOMAChullVrtx <- apply(NEATOMA.Agg[SiteNEATOMAUse,TaxaUse], 
                          1,
                          function(x){
                            # Species in a site
                            SppSiteUse <- TaxaUse[x!=0]
                            # Traits of those species
                            TraitsSite <- TraitsDBS[rownames(TraitsDBS)%in%SppSiteUse,]
                            # Estimate the convex hull
                            #TraitsSite[convhulln(TraitsSite),]# Using convhulln of the geometry package (the space is set based on triangular facets)
                            #TraitsSite[chull(TraitsSite),] # Using the chull methos
                            return(TraitsSite)
                          })
#####
#####
# Estimate the min intersection between each NEATOMA assemblages and all NPDB convex-hulls pairs

SrtTime<-Sys.time()
sfInit( parallel=TRUE, cpus=10)
sfExport("ChullNPDBVrtx")
sfExport("NEATOMAChullVrtx")

NEATOMAChullDistList <- sfLapply(1:length(NEATOMAChullVrtx),
                          function(x){#x<-1
                            DistOut <- lapply(ChullNPDBVrtx,
                                              function(y){#y<-ChullNPDBVrtx[[2]]
                                                try(
                                                  # Get intersection convex hulls - this function needs as input the full set of points
                                                  out <- geometry::intersectn(ps1 = NEATOMAChullVrtx[[x]],
                                                                              ps2 = y) # Get intersection convex hulls  
                                                                              ,silent=T)
                                                if("out"%in%ls()){
                                                  ProOverLp <- out$ch$vol/out$ch2$vol #define the proportion of the base convex hull the second convex hull overlaps
                                                  PairDistOut <- 1 - ProOverLp # Make the proportion overlap a dissimilarity index
                                                }
                                                else{
                                                  PairDistOut<-NA
                                                }
                                                return(round(PairDistOut,5))
                                              })
                            DistFullContrast <- do.call("c",DistOut) # Turn the dissimilarity list into a vector
                            hist(DistFullContrast)
                            MnDistFullContrast <- min(DistFullContrast,na.rm=T)
                            return(MnDistFullContrast)
                          })
Sys.time()-SrtTime
sfStop()
NEATOMAChullDist <- do.call("c",NEATOMAChullDistList)


a <- NEATOMA.Agg[SiteNEATOMAUse,]
b<- data.frame(Dist =  NEATOMAChullDist,
               Time = a$Time)




plot(y = rev(tapply(b$Dist,b$Time,mean)),
     x = -21:-1,
     type = "b",
     main = "Regional Dissimilarity Trait Change", 
     xlab ="Time (kyrBP)",
     ylab = "Disimilarity (Sqr Cord Dist)")
