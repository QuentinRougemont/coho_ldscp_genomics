#!/usr/bin/env Rscript

#### load libs
libs <- c('dplyr','reshape2','ade4','data.table', 'magrittr', 'factoextra','vegan')
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
pop_lat <- read.table("02-data/env/coho_dist_v2.txt",T)
pop_lat <- dplyr::select(pop_lat,POP_ID, elevation, dist_max_km, Latitude, Region)
#replace altitude of zero by 1
pop_lat$elevation[pop_lat$elevation == 0.00000 ] <- 1
enviro <-read.table("02-data/env/climat_epic4_wanted_pop.txt",T)
geol <- read.table("02-data/env/era_rocktype_quanti_v2.txt",T)
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

#verify that the order of pop in env match the genotype matrix:
setdiff(env2$SITE, pop$SITE)  #there should be no difference 
setdiff(pop$SITE,env2$SITE)   #idem
table(pop$SITE==env2$SITE)    #if false occur then there's difference

#remove unwanted variable :
env <- dplyr::select(env2, -SITE, -elevation, -Region, -dist_max_km)
colnames(env)[1] <- "Geology"


#perform RDA
rda1 <- rda(freq4 ~ Geology +
    Temperature1 + Temperature2 + Temperature3 + Temperature4 + 
    Precipitation1 + 
    Precipitation2 + 
    Precipitation3 + normalized_distance +
    Condition(Latitude) , 
    data=env, scale=T)

#R²
sink("r2_eigenval_vif.txt")
RsquareAdj(rda1)
summary(eigenvals(rda1, model= "constrained"))
#check the vif:
vif.cca(rda1)
sink()

#screeplot
pdf(file="screeplot.pdf")
screeplot(rda1)
dev.off()

####extract the scores ###########
#the 6 RDA axis are significant (see script for signficance)
#and we keep them
signif <- 7
load.rda <- scores(rda1, choices=c(1:signif), display="species") 
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}
thresh = 3

for(i in 1:signif){
   nam <- paste("cand",i,sep="")
   assign(nam, outliers(load.rda[,i],thresh) )
}

#get total number of outliers 
tot = ls(pattern="cand") 
total=NULL
tmp=NULL
for(i  in 1:signif){
tmp <- length(get(tot[i]))
  total <-rbind(total, tmp)
}
ncand <- sum(total[,1])

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))
cand6 <- cbind.data.frame(rep(6,times=length(cand6)), names(cand6), unname(cand6))
cand7 <- cbind.data.frame(rep(7,times=length(cand7)), names(cand7), unname(cand7))

colnames(cand1) <- colnames(cand2)<- colnames(cand3) <- c("axis","snp","loading")
colnames(cand4) <- colnames(cand5) <- colnames(cand6)<-colnames(cand7) <- colnames(cand1)
cand <- rbind(cand1,cand2,cand3,cand4,cand5,cand6,cand7)
cand$snp <- as.character(cand$snp)
foo <- matrix(nrow=(ncand), ncol=ncol(env))  # 8 columns for 8 predictors
colnames(foo) <- colnames(env)

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- freq4[,nam]
  foo[i,] <- apply(env,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)
length(cand$snp[duplicated(cand$snp)])  # 7 duplicate detections

foo <- cbind(cand$axis, duplicated(cand$snp))
table(foo[foo[,1]==1,2]) #
table(foo[foo[,1]==2,2]) # 

cand <- cand[!duplicated(cand$snp),]  #remove duplicate
col<-ncol(cand)

#correlation extraction:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,col+1] <- names(which.max(abs(bar[4:col]))) # gives the variable
  cand[i,col+2] <- max(abs(bar[4:col]))              # gives the correlation
}
colnames(cand)[col+1] <- "predictor"
colnames(cand)[col+2] <- "correlation"
d = data.frame(table(cand$predictor))
write.table(d, "outlier_per_predictor.txt", quote =F , row.names = F, col.names = F)

write.table(cand,
    paste("candidate_outliers_with_var",thresh,".txt",sep=""),
    quote=F, row.names=F, col.names=T)

snp <- ("02-data/genotype/population_for_GEA.vcfsnp")
snp <- read.table(snp) %>% 
    select(V1, V2, V3) %>% 
    set_colnames(., c("CHR","POS","snp"))

snp_cand <- merge(snp, cand)

write.table(snp_cand,
    paste("candidate_outliers_with_var",thresh,"_with_chr.txt",sep=""),
    quote=F, 
    row.names=F, col.names=T)


#distribution of outliers
dat <-  data.frame(table(snp_cand$CHR))
write.table(dat,"distribution_outliers_along_chr.txt",quote=F,row.names=F,col.names=F)

d = read.table("distribution_outliers_along_chr.txt")
chr = read.table("../length.txt.gz")
#lfmm = read.table("LFMM/significant_outlier_control_forLatitudeK20.txt",T)
#rda = read.table("01.RDA/candidate_outliers_with_var3_with_chr.txt",T)

chr_d <- merge(chr,d, by="V1")
summary(lm(chr_d$V2.x ~ chr_d$V2.y))

print("searching epas1")
filter(snp_cand, snp == "4244396:56:+")
