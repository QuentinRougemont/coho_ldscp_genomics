#!/usr/bin/env/ Rscript

#Purpose: RDA script to test significance of the axis and environemental variable
#Author : QR
#Date: October 2021

###Verify if libs are installed:
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("magrittr" %in% rownames(installed.packages()) == FALSE)
{install.packages("magrittr", repos="https://cloud.r-project.org") }
if("reshape2" %in% rownames(installed.packages()) == FALSE)
{install.packages("reshape", repos="https://cloud.r-project.org") }
if("vegan" %in% rownames(installed.packages()) == FALSE)
{install.packages("vegan", repos="https://cloud.r-project.org") }
if("data.table" %in% rownames(installed.packages()) == FALSE)
{install.packages("data.table", repos="https://cloud.r-project.org") }
if("corrplot" %in% rownames(installed.packages()) == FALSE)
{install.packages("corrplot", repos="https://cloud.r-project.org") }

## load libs
libs <- c('dplyr','reshape2','data.table', 'magrittr','vegan','corrplot')
invisible(lapply(libs, library, character.only = TRUE))

################## DOWNLOAD SNP frequency data      ############################
freq <- fread("zcat plink.frq.strat.gz") #obtained from plink
freq2 <- dplyr::select(freq,SNP,CLST,MAF)
freq3 <- reshape2::dcast(freq2,SNP~CLST)
freq3 <- freq3[,-1]
freq4 <- t(freq3)
pop <- unique(freq2$CLST)
pop <- data.frame(pop)
snps <- unique(freq2$SNP)
colnames(freq4) <- snps

########## DOWNLOAD ENVIRONMENTAL DATA #########################################
metadata <- read.table("01-info/metadata", h = T, sep = "\t")
#remove BNV and SAI:
metadata<- metadata %>%filter(POP !="SAI" & POP !="BNV")

#remove BNV and SAI:
geol <- read.table("02-data/env/geology", T)
geol <- geol %>%filter(POP !="SAI" & POP !="BNV")

#replace altitude of zero by 1
metadata$elevation[metadata$elevation == 0.00000 ] <- 1

#now takes into accont the standardized elevation * distance interaction
#standardiztation function:
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
metadata$normalized_distance =  metadata$elevation * metadata$dist_max_km
metadata$normalized_distance <- range01(metadata$normalized_distance)

temp <- read.table("02-data/env/temperature_coordinates_pca.txt", T)
prec <- read.table("02-data/env/precipitation_coordinates_pca.txt", T)

#merge all:
env1 <- dplyr::select(metadata, POP, Latitude, normalized_distance, Region)
env1 <- merge(env1, temp)
env1 <- merge(env1, prec)
env1 <- merge(env1, geol)

#remplacer avec sed directement dans metadata, geology et bioclim
env1 <- env1[order(env1$POP),]

#ensuire that environmental variable will match the genotype data order:
colnames(pop) <- "POP"
env1 <- left_join(pop, env1) 

#choose the variable we want to work with:
env <- select(env1, -POP, -Region)

#rename for readiblity of the plot:
colnames(env) <- c("Latitude","normalized_distance", 
    "TemperaturePC1", "TemperaturePC2", "TemperaturePC3", "TemperaturePC4",
    "PrecipitationPC1", "PrecipitationPC2", "PrecipitationPC3", 
    "geology")

#write the correlation 
mat <- round(cor(env),2)

pdf(file="corrplot.pdf")
corrplot(mat, method = "number", type = "lower")
dev.off()


####################### Now perform the RDA ####################################
rda1 <- rda(freq4 ~ geology +
    TemperaturePC1 + TemperaturePC2 + TemperaturePC3 + TemperaturePC4 + 
    PrecipitationPC1 + PrecipitationPC2 + PrecipitationPC3 + normalized_distance +
    Condition(Latitude) , 
    data=env, scale=T)

######  significance testing  ##################################################

all.rda<- rda(freq4 ~ 1 , data = env, scale =T)

ord <- ordiR2step(all.rda, scope= formula(rda1), 
    direction = "forward", 
    permutations = 1000,
    Pin = 0.1, R2scope = FALSE, verbose=FALSE)
write.table(ord$anova,"ordistep_anova.txt",quote=F,sep="\t")
sink("ordistep")
print(ord)
sink()

#rda1
signif.axis <- anova.cca(rda1, by="axis", parallel=5 ) 
signif.axis

print("signif.axis.is.done")
signif.marg <- anova.cca(rda1, by="margin", parallel=5)
signif.marg
print("signif.marg.is.done")

signif.full <- anova.cca(rda1, parallel=5 ) # default is permutation=999, could increase to 9999
signif.full
print("signif.full.is.done")

sink("signif.axis")
print(signif.axis)
sink()

sink("signif.marg")
print(signif.marg)
sink()
