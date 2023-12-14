rm(list=ls());gc()
require(analogue)
require(raster)
require(terra)
setwd("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/PaleoNovelty")
#####
# Load the North American Pollen Data Base
NPDB <- read.csv("./Data/Pollen/Raw/whitmoreetal2005_v1-72.csv")
# Taxa translate between NPDB and Paleo Data
TaxtTrans <- read.csv("./Data/Pollen/Raw/Taxa transation.csv")

# Which are the summarizing taxa - note the TaxtTrans$With.Traits!="" filter -> only species with traits are used
TaxaUse <- sort(unique(TaxtTrans$Translated[TaxtTrans$With.Traits!=""]))

# Build a translated NPDB dataset for taxa with trait values
SummNPDB <- lapply(TaxaUse, # only use species with trait data (Woody species)
                   function(x){
                     out <- NPDB[,TaxtTrans$Original[TaxtTrans$Translated==x]]
                     if(class(out)=="data.frame"){out <- apply(out,1,sum)}
                     out <- as.numeric(out)
                     return(out)
                   })
names(SummNPDB) <- TaxaUse
SummNPDB2 <- as.data.frame(do.call("cbind",SummNPDB))
SummNPDB2 <- SummNPDB2[which(!is.na(apply(SummNPDB2,1,sum))),]

# Re-scale NPDS so pollen counts add up to 1
SummNPDB3 <- as.data.frame(t(apply(SummNPDB2,1,function(x){x/sum(x,na.rm=T)})))
# Make 0 the pollen count for taxa with percentages below 5%
SummNPDB3 <-  as.data.frame(t(apply(SummNPDB3,1,function(x){ifelse(x>0.05,x,0)})))
# Aggregate to taxa with data 
SummNPDB3 <-  as.data.frame(t(apply(SummNPDB3,1,function(x){x/sum(x,na.rm=T)})))
# Remove taxa not sampled in any site
SummNPDB3 <- SummNPDB3[,apply(SummNPDB3,2,sum)!=0]

# Add the Biome for each observation
# Make the locations a SpatVector
NPDBVect <- vect(NPDB[which(!is.na(apply(SummNPDB2,1,sum))),c("LONDD", "LATDD")],
                 geom=c("LONDD", "LATDD"),
                 crs = "EPSG:4326")

# Get the Biome for each observation
BIOMES <- vect("~/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/3. Ecography PaleoTraits Paper/Data/WWF-Biomes/wwf_terr_ecos.shp")
BiomesNames <- data.frame(code = c(1:14,98,99,100),
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
                                   "Mangroves",
                                   "Rock/ice",
                                   "Lake/River",
                                   "Ocean"))
BiomeNPD <- extract(BIOMES[,"BIOME"],
                    NPDBVect)
BiomeNPD$BIOME[c(is.na(BiomeNPD$BIOME))] <- 100

# Create a new aggregated North American Pollen Data Base
NPDBAgg <- data.frame(NPDB[which(!is.na(apply(SummNPDB2,1,sum))),c("ID1","SITE", "DBCODE", "SITECODE", "SITENAME", "LONDD", "LATDD","ELEVATION")],
                      BIOME  = BiomesNames$Name[c(match(BiomeNPD[,2],BiomesNames$code))],
                      SummNPDB3)
# Remove the sites located in Rock/ice, Lake/River, or Ocean.
NPDBAgg <- NPDBAgg[BiomeNPD$BIOME<15,]
# Save the Location data
write.csv(NPDBAgg,"./Data/Pollen/Processed/NPDB_Agg.csv")
