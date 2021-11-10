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


###### DOWLOAD SNP frequency data ###########################################################
freq <- fread("zcat 02-data/genotype/plink.frq.strat.gz") #obtained from plink
freq2 <- dplyr::select(freq,SNP,CLST,MAF)
freq3 <- reshape2::dcast(freq2,CLST ~ SNP)
freq4 <- freq3[,-1]
pop <- freq3$CLST
pop <- data.frame(pop)
snp = colnames(freq3)
##############################################################################################
###### DOWNLOAD ENVIRONMENTAL  DATA  ########################################################
pop_lat <- read.table("coho_dist_v2.txt",T)
pop_lat <- dplyr::select(pop_lat,POP_ID, elevation, dist_max_km, Latitude, Region)
#replace altitude of zero by 1
pop_lat$elevation[pop_lat$elevation == 0.00000 ] <- 1
enviro <-read.table("climat_epic4_17_07_2021_v2.txt",T)
geol <- read.table("era_rocktype_quanti_v2.txt",T)
enviro <- merge(enviro, geol, by="SITE")
 
colnames(pop) <- "SITE"
enviro <- merge(pop, enviro, by = "SITE",  sort=F)
enviro <- filter(enviro, SITE %in% pop$SITE)
####### PERFORM PCA ON ENVT VAR #################################################
X.temp <- dudi.pca(df = enviro[, 3:57], center = T, scale = T, scannf = FALSE, nf = 4)
X.prec  <- dudi.pca(df = enviro[, 58:97], center = T, scale = T, scannf = FALSE, nf = 4)

prec.ind <- get_pca_ind(X.prec)
temp.ind <- get_pca_ind(X.temp)
colnames(temp.ind$coord) <- c("Temperature1","Temperature2","Temperature3","Temperature4")
colnames(prec.ind$coord) <- c("Precipitation1","Precipitation2","Precipitation3")
#keep significant axis (4 in temperaturure and 3 in precipitation and bind population name :
temp <- cbind(enviro$SITE, temp.ind$coord[,c(1:4)])
prec <- cbind(enviro$SITE, prec.ind$coord[,c(1:3)])
colnames(temp)[1] <- "SITE"
colnames(prec)[1] <- "SITE"

env1 <- dplyr::select(enviro, SITE,ROCK)
env1 <- left_join(env1, temp  ) #$coord[,c(1:3)])
env1 <- left_join(env1, prec )  #$coord[,c(1:3)])
colnames(pop_lat)[1]<- "SITE"
env2 <- merge(env1, pop_lat, by="SITE")

#now takes into accont the standardized elevation * distance interaction
#standardiztation function:
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
env2$normalized_distance =  env2$elevation * env2$dist_max_km
env2$normalized_distance <- range01(env2$normalized_distance)
#Keep only wanted variable for RDA:
env2 <- left_join(pop, env2) #, c("pop"="SITE")) 

#remove unwanted variable :
env <- dplyr::select(env2, -SITE, -elevation, -Region, -dist_max_km)
colnames(env)[1] <- "Geology"

#write the correlation 
mat <- round(cor(env),2)

pdf(file="corrplot.pdf")
corrplot(mat, method = "number", type = "lower")
dev.off()

####################### Now perform the RDA ####################################
rda1 <- rda(freq4 ~ Geology +
    Temperature1 + Temperature2 + Temperature3 + Temperature4 + 
    Precipitation1 + 
    Precipitation2 + 
    Precipitation3 + normalized_distance +
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
