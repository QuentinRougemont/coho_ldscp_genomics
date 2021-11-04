#!/usr/bin/env/ Rscript

#Purpose: RDA script to extract outliers:
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

## load libs
libs <- c('dplyr','reshape2','data.table', 'magrittr','vegan')
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


####################### Now perform the RDA ####################################
rda1 <- rda(freq4 ~ geology +
    TemperaturePC1 + TemperaturePC2 + TemperaturePC3 + TemperaturePC4 + 
    PrecipitationPC1 + PrecipitationPC2 + PrecipitationPC3 + normalized_distance +
    Condition(Latitude) , 
    data=env, scale=T)

################################################################################
#R²
RsquareAdj(rda1)
summary(eigenvals(rda1, model= "constrained"))

#check the vif:
vif.cca(rda1)

#screeplot
pdf(file="screeplot.pdf")
screeplot(rda1)
dev.off()

###############   extract the scores and find outliers #########################
#7 RDA axis are significant (see script for signficance)
#and we keep them
signif <- 7
load.rda <- scores(rda1, choices = c(1:signif), display = "species") 

pdf(file="hist2.pdf")
for (i in 1:signif){ 
hist(load.rda[,i], main = paste("Loadings on RDA", i ,sep = ""))
}
dev.off()

##### find outliers 
##### TO DO: add Robust mahalanobis distance to extract score and -log10(pval) for each variable

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}
thresh = 3

for(i in 1:signif){
   nam <- paste("cand", i, sep = "")
   assign(nam, outliers(load.rda[,i],thresh) )
}

#get total number of outliers 
tot = ls(pattern="cand") 
total = NULL
tmp = NULL
for(i  in 1:signif){
tmp <- length(get(tot[i]))
  total <-rbind(total, tmp)
}
ncand <- sum(total[,1])

cand1 <- cbind.data.frame(rep(1, times = length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2, times = length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3, times = length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4, times = length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5, times = length(cand5)), names(cand5), unname(cand5))
cand6 <- cbind.data.frame(rep(6, times = length(cand6)), names(cand6), unname(cand6))
cand7 <- cbind.data.frame(rep(7, times = length(cand7)), names(cand7), unname(cand7))

colnames(cand1) <- colnames(cand2)<- colnames(cand3) <- colnames(cand4) <- c("axis","snp","loading")
colnames(cand5) <- colnames(cand6)<- colnames(cand7) <- colnames(cand1)

cand <- rbind(cand1,cand2,cand3,cand4,cand5,cand6, cand7) 
cand$snp <- as.character(cand$snp)
foo <- matrix(nrow=(ncand), ncol=ncol(env))  # 8 columns for 8 predictors
foo <- matrix(nrow=(ncand), ncol=ncol(env[-1]))  # 8 columns for 8 predictors

#colnames(foo) <- colnames(env)
colnames(foo) <- colnames(env[-1]) #to rm lat

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- freq4[,nam]
  #foo[i,] <- apply(env,2,function(x) cor(x,snp.gen))
  foo[i,] <- apply(env[,-1],2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)
#head(cand)

length(cand$snp[duplicated(cand$snp)])  # 7 duplicate detections
foo <- cbind(cand$axis, duplicated(cand$snp))

cand <- cand[!duplicated(cand$snp),]  #remove duplicate
col<-ncol(cand)

#correlation:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,col+1] <- names(which.max(abs(bar[4:col]))) # gives the variable
  cand[i,col+2] <- max(abs(bar[4:col]))              # gives the correlation
}

colnames(cand)[col+1] <- "predictor"
colnames(cand)[col+2] <- "correlation"
table(cand$predictor)

write.table(cand,
    paste("candidate_outliers_with_var",thresh,".txt",sep=""),
    quote=F, row.names=F, col.names=T)


snp <- "population_for_GEA.renamed.vcfsnp"
snp <- read.table(snp) %>% 
    select(V1, V2, V3) %>% 
    set_colnames(., c("CHR","POS","snp"))


snp_cand <- merge(snp, cand)

write.table(snp_cand,
    paste("candidate_outliers_with_var",thresh,"_with_chr.txt",sep=""),
    quote=F, 
    row.names=F, col.names=T)


