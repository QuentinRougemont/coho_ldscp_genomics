
library(magrittr)
library(dplyr)
library(corrplot)

###### DOWNLOAD ENV DATA
metadata <- read.table("01-info/metadata", h = T, sep = "\t")
metadata<- metadata %>%filter(POP !="SAI" & POP !="BNV")

geol <- read.table("02-data/env/geology", T)
#remove BNV and SAI:
geol <- geol %>%filter(POP !="SAI" & POP !="BNV")
#remove BNV and SAI:

#replace altitude of zero by 1
metadata$elevation[metadata$elevation == 0.00000 ] <- 1

#now takes into accont the standardized elevation * distance interaction
#standardiztation function:
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
metadata$normalized_distance =  metadata$elevation * metadata$dist_max_km
metadata$normalized_distance <- range01(metadata$normalized_distance)

####### PERFORM PCA ON ENVT VAR #################################################
temp <- read.table("02-data/env/temperature_coordinates_pca.txt", T)
prec <- read.table("02-data/env/precipitation_coordinates_pca.txt", T)

env1 <- dplyr::select(metadata, POP, Latitude, normalized_distance, Region)
env1 <- merge(env1, temp)
env1 <- merge(env1, prec)
env1 <- merge(env1, geol)

#remplacer avec sed directement dans metadata, geology et bioclim
env1 <- env1[order(env1$POP),]

#choose the variable we want to work with:
env <- select(env1, -POP, -Region)

#rename for readiblity of the plot:
colnames(env) <- c("Latitude","normalized_distance", 
    "TemperaturePC1", "TemperaturePC2", "TemperaturePC3", "TemperaturePC4",
    "PrecipitationPC1", "PrecipitationPC2", "PrecipitationPC3", 
    "geology")

#again write the correlation 
mat <- round(cor(env),2)

pdf(file="corrplot.pdf")
corrplot(mat, method = "number", type = "lower")
dev.off()
