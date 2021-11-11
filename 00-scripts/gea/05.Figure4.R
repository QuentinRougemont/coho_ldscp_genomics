#Script to test for parallelism
#Author: Q.R
#date: 20-04-20
###Verify if libs are installed:
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("magrittr" %in% rownames(installed.packages()) == FALSE)
{install.packages("magrittr", repos="https://cloud.r-project.org") }
if("ade4" %in% rownames(installed.packages()) == FALSE)
{install.packages("ade4", repos="https://cloud.r-project.org") }
if("reshape" %in% rownames(installed.packages()) == FALSE)
{install.packages("reshape", repos="https://cloud.r-project.org") }
if("factoextra" %in% rownames(installed.packages()) == FALSE)
{install.packages("factoextra", repos="https://cloud.r-project.org") }
if("vegan" %in% rownames(installed.packages()) == FALSE)
{install.packages("vegan", repos="https://cloud.r-project.org") }
if("data.table" %in% rownames(installed.packages()) == FALSE)
{install.packages("data.table", repos="https://cloud.r-project.org") }

## load libs
libs <- c('dplyr','reshape2','ade4','data.table', 'magrittr', 'factoextra','vegan','cowplot','ggsci')
invisible(lapply(libs, library, character.only = TRUE))

freq <- fread("zcat 02-data/genotype/plink.frq.strat.gz")
freq2 <- dplyr::select(freq,SNP,CLST,MAF)
pop <- unique(freq2$CLST)
pop <- data.frame(pop)

########## DOWNLOAD ENVIRONMENTAL DATA #########################################
#metadata <- read.table("01-info/metadata", h = T, sep = "\t")
pop_lat <- read.table("02-data/env/coho_dist_v2.txt",T)
pop_lat$elevation[pop_lat$elevation == 0.00000 ] <- 1

metadata <- dplyr::select(pop_lat,POP_ID, elevation, dist_max_km, Latitude, Region)
colnames(metadata)[1] <-c("POP")
#remove BNV and SAI:
metadata<- metadata %>%filter(POP !="BNV")
#standardiztation function:
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
metadata$normalized_distance =  metadata$elevation * metadata$dist_max_km
metadata$normalized_distance <- range01(metadata$normalized_distance)

#merge all:
env1 <- dplyr::select(metadata, POP, Latitude, normalized_distance, Region)
env1 <- env1[order(env1$POP),]

#ensure that environmental variable will match the genotype data order:
colnames(pop) <- "POP"
env1 <- left_join(pop, env1) 

### read rda results ###########################################################
rdares <- read.table("candidate_outliers_with_var3_with_chr.txt",T) #%>% select(-CHR, -POS)
colnames(rdares)[1] <- 'SNP'

lfmm <- read.table("significant_outlier_control_forLatitudeK20.txt",T)
lfmm_dist <- filter(lfmm, var=="Dist")

all <- full_join(lfmm, rdares) #, c("SNP"="snp")) #combine data to exclude for random pca

dist_rda <- dplyr::filter(rdares, predictor=="normalized_distance")
dist_rda_snp <- dplyr::select(dist_rda, SNP) #keep only SNPs
rdafull_lfmm <- merge(lfmm_dist,dist_rda_snp) # by.x="SNP", by.y="V3")

rda_lfmm <- merge(lfmm, rdares)
tab <- as.data.frame(table(rda_lfmm$var, rda_lfmm$predictor))
write.table(rda_lfmm, "shared_outlier.txt", quote = F, row.names = F)
write.table(rdafull_lfmm, "shared_outlier_distance.txt", quote = F, row.names = F)
write.table(tab, "table_of_shared_outliers.txt", quote = F ) 
sink("table_of_shared_outliers.txt")
table(rda_lfmm$var, rda_lfmm$predictor)
sink()


#keep strong outlier
rda_strong <- filter(dist_rda, correlation > 0.64)
lfmm_strong <- filter(lfmm_dist, BH <0.009)
rda_strong_snp <- select(rda_strong , SNP)
#strong <- merge(lfmm_dist, rda_strong)
strong <- merge(lfmm_strong, rda_strong)

write.table(strong,"shared_outlier_distance_strong.txt",quote=F,row.names=F)

rdafull_lfmm <- strong

#rda_lfmm <- rda_lfmm %>% filter(predictor == "normalized_distance" | var == "Dist")
################################################################################
################################################################################
## subset frequencies on random samples:
wRandom <- dplyr::filter(freq2, ! SNP %in% all$SNP)
freq3 <- reshape2::dcast(wRandom,CLST~SNP)
rownames(freq3) = freq3[,1]
freqRandom <- freq3[,-1]
freqRandom <- sample(freqRandom,500) #choose random individuals
pcaRandom <- dudi.pca(freqRandom,scale=FALSE,scannf=FALSE)

## working on LFMM -RDA strong outlier
shared_snp <- select(rdafull_lfmm,SNP)
wLFMM_RDA <- dplyr::filter(freq2, SNP %in% shared_snp$SNP)
freq3 <- reshape2::dcast(wLFMM_RDA,CLST~SNP)
rownames(freq3) = freq3[,1]
freqLFMM_RDA <- freq3[,-1]
pcaLFMM_RDA <- dudi.pca(freqLFMM_RDA,scale=FALSE,scannf=FALSE)
eig.val <- get_eigenvalue(pcaLFMM_RDA) #first we get the eigen value in an easy form        
eig.val$eig <- eig.val$variance.percent/100                                          

#get individual coordinates for all PCAs
res.indRandom <- get_pca_ind(pcaRandom)
res.indLFMM_RDA <- get_pca_ind(pcaLFMM_RDA)

### Will use broken sticks from vegan package:                                       
eig.val <- get_eigenvalue(pcaLFMM_RDA) #first we get the eigen value in an easy form        
eig.val$eig <- eig.val$variance.percent/100                                          
expected <- bstick(length(eig.val$eig) ) #get expected value given number of axis    
signif <- eig.val$eig > expected                                                     
signif[1:10] ##only the first axis should be significant if we successfully capture the effect of distance 
write.table(eig.val, "eigenvalue", quote=F, row.names=F)
write.table(signif, "siginificativity_of_pca.txt", quote = F, row.names = F)
######## Prepare data for PLOT #########################################
library(ggsci)
strata <- env1
pRandom <- fviz_pca_ind(pcaRandom, habillage=strata$Region,geom="text",labelsize=3,
                      addEllipses=F, ellipse.level=0.95) + 
  theme(legend.position="none") + scale_color_igv() +
  theme_minimal()

#shared and correlated:
pLFMM_RDA <- fviz_pca_ind(pcaLFMM_RDA, habillage=strata$Region,geom="text",
                      addEllipses=F, ellipse.level=0.95) + 
  theme(legend.position="none") + scale_color_igv() +
  theme_minimal()

pdf(file="pca_based_on_random_snp.pdf")
pRandom
dev.off()
pdf(file="pca_distance_outliers.pdf")
pLFMM_RDA
dev.off()
########################################################################################
#  now test the effect of Distance*elevation and Latitude on our inferences
#  Latitude, if well controlled for should be non-significant when using outliers
#  On the contrary Latitude should be significant when using random SNPs
#  A significant interaction reveals both parallelism & non-parallelism (i.e. location specific
# outliers)
#######################################################################################

metadata <- select(metadata, POP, Region, Latitude, elevation, normalized_distance, dist_max_km)

#test on random SNPs:
tmp <- data.frame(row.names(
  res.indRandom$contrib), 
  res.indRandom$contrib[,1], 
  res.indRandom$coord[,1]) %>% set_colnames(.,c("POP_ID","contrib","coord"))
tmpdist <- merge(tmp, metadata,  by.x="POP_ID", by.y="POP")
tmpdist$grp = ifelse(tmpdist$dist_max_km >438 | tmpdist$elevation>403, "long", "short")
tmpdist %>% group_by(grp) %>% summarise(m_con = mean(contrib), m_coord= mean(coord), n=n())

#test de comparaison des modèles:
res0 = (lm(tmpdist$coord ~ tmpdist$normalized_distance ))
res1 = (lm(tmpdist$coord ~ tmpdist$normalized_distance + tmpdist$Latitude))
res2 = (lm(tmpdist$coord ~ tmpdist$Latitude))
res3 = (lm(tmpdist$coord ~ tmpdist$normalized_distance + tmpdist$normalized_distance * tmpdist$Latitude + tmpdist$Latitude))
res4 = (lm(tmpdist$coord ~ 1))
AIC(res0,res1,res2,res3, res4)

summary(res3)
sink("resultats_lm_random.txt")
print(summary(res3))
sink()
tmpdist$grp <- gsub("short","short-distance",tmpdist$grp)
tmpdist$grp <- gsub("long","long-distance",tmpdist$grp)
tmpdist$Region[tmpdist$Region=="Washington&Oregon"] <- "Cascadia"
tmpdistRandom <- tmpdist

####################  LFMM-RDA shared strong outliers:    ################
tmp <- data.frame(row.names(
  res.indLFMM_RDA$contrib), 
  res.indLFMM_RDA$contrib[,1], 
  res.indLFMM_RDA$coord[,1]) %>% set_colnames(.,c("POP_ID","contrib","coord"))
tmpdist <- merge(tmp, metadata,  by.x="POP_ID", by.y="POP")
tmpdist$grp = ifelse(tmpdist$dist_max_km >438 | tmpdist$elevation>403, "long", "short")
tmpdist %>% group_by(grp) %>% summarise(m_con = mean(contrib), m_coord= mean(coord), n=n())

#test de comparaison des modèles:
res0 = (lm(tmpdist$coord ~ tmpdist$normalized_distance ))
res1 = (lm(tmpdist$coord ~ tmpdist$normalized_distance + tmpdist$Latitude))
res2 = (lm(tmpdist$coord ~ tmpdist$Latitude))
res3 = (lm(tmpdist$coord ~ tmpdist$normalized_distance + tmpdist$normalized_distance * tmpdist$Latitude + tmpdist$Latitude))
res4 = (lm(tmpdist$coord ~ 1))
AIC(res0,res1,res2,res3, res4)

summary(res3)
sink("resultats_lm_outliers.txt")
print(summary(res3))
sink()

tmpdist$grp <- gsub("short","short-distance",tmpdist$grp)
tmpdist$grp <- gsub("long","long-distance",tmpdist$grp)
tmpdist$Region[tmpdist$Region=="Washington&Oregon"] <- "Cascadia"
tmpdistLFMM_RDA <- tmpdist

###################################################################################
#now do the plots:
#################################################################################
#usual colors:
myColors <- c("blue","orange","red","green","darkviolet","springgreen4")
names(myColors) <- levels(tmpdistLFMM_RDA$Region)
colScale <- scale_colour_manual(name = "REGION",values = myColors)

#plot of randomSNP:
plot_Random <- ggplot(tmpdistRandom, aes(x=coord, y=normalized_distance, label=POP_ID) ) +
  geom_point()+ geom_label(aes(colour=factor(Region)))+colScale +
  theme(legend.position="none") +
  theme(strip.background = element_blank(),strip.text.x = element_blank() )+
  ylab("Normalized migratory distance") + xlab("PC1 coordinates") +
  theme_classic() + 
  theme(axis.title.x=element_text(size=12, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=12,family="Helvetica",face="bold", angle=0, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=12, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=12,family="Helvetica",face="bold"),
        strip.text.x = element_text(size=12)) +
  theme(panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())

#plot of outliers
plot_LFMM_RDA <- ggplot(tmpdistLFMM_RDA, aes(x=coord, y=normalized_distance, label=POP_ID) ) +
  geom_point()+ geom_label(aes(colour=factor(Region)))+colScale +
  theme(legend.position="none") +
  theme(strip.background = element_blank(),strip.text.x = element_blank() )+
  ylab("Normalized migratory distance") + xlab("PC1 coordinates") +
  theme_classic() + 
  theme(axis.title.x=element_text(size=12, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=12,family="Helvetica",face="bold", angle=0, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=12, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=12,family="Helvetica",face="bold"),
        strip.text.x = element_text(size=12)) +
  theme(panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())

first_row = plot_grid(plot_Random, labels = c('A'))
second_row = plot_grid(pRandom, labels = c('B'), nrow = 1)
gg_all = plot_grid(first_row, second_row, labels=c('', ''), ncol=1)

###########################################################################################
## 	look at the MAF and frequency shift
############################################################################################
freqshift <- cbind(tmpdistLFMM_RDA,freqLFMM_RDA) %>%
   select( -contrib, -coord, -Latitude,-elevation, -dist_max_km, -normalized_distance,-Region)  %>%
   reshape2::melt() %>% group_by(grp,variable) %>% summarise(mean=mean(value))
freqshiftgrp <- gsub("long-distance","1-long"  , freqshift$grp)
freqshift$grp <- gsub("short-distance","0-short", freqshift$grp )
freqshiftLFMM_RDA <- freqshift

## random SNP
freqshift <- cbind(tmpdistRandom,freqRandom) %>%
   select( -contrib, -coord, -Latitude,-elevation, -dist_max_km, -normalized_distance,-Region)  %>%
   reshape2::melt() %>% group_by(grp,variable) %>% summarise(mean=mean(value))

freqshiftgrp <- gsub("long-distance","1-long"  , freqshift$grp)
freqshift$grp <- gsub("short-distance","0-short", freqshift$grp )
freqshiftRandom <- freqshift

delta <-  dcast(freqshiftLFMM_RDA, variable ~ grp)
delta$diff = delta$`long-distance` - delta$`0-short`
print("looking at shift in allele frequency")
summary(delta$diff)

sink("allele_frequency_shit.txt")
summary(delta$diff)
sink()

####################################################################################################
#plot all shared_strong_LFMM SNPs:
###################################################################################################
shift_LFMM_RDA <- ggplot(freqshiftLFMM_RDA, aes(x=grp, y=mean, group=variable)) + 
  geom_point(aes(color=variable)) + 
  geom_line(aes(color=variable)) +
  scale_color_manual(values=c(rep(c("#999999", "#E69F00"),605), "black") ) +
  ylab("Allele frequency shift") +  xlab("Normalized Migratory distance") +
  theme_classic() +   theme(legend.position="none")  + 
  theme(axis.title.x=element_text(size=12, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=12,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=12, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=12,family="Helvetica",face="bold"),
        strip.text.x = element_text(size=12)) + 
  theme(panel.grid.major = element_blank() )

shift_Random <- ggplot(freqshiftRandom, aes(x=grp, y=mean, group=variable)) + 
  geom_point(aes(color=variable)) + 
  geom_line(aes(color=variable)) +
  scale_color_manual(values=c(rep(c("#999999", "#E69F00"),605), "black") ) +
  ylab("Allele frequency shift") +  xlab("Normalized Migratory distance") +
  theme_classic() +   theme(legend.position="none")  + 
  theme(axis.title.x=element_text(size=12, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=12,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=12, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=12,family="Helvetica",face="bold"),
        strip.text.x = element_text(size=12)) + 
  theme(panel.grid.major = element_blank() )

pdf(file="Figure4ABstrong.pdf",12,7)
plot_grid(
  plot_LFMM_RDA,
  shift_LFMM_RDA, 
  labels = c("A", "B"), ncol = 2,
  label_size = 12)
dev.off()

first_row = plot_grid(plot_Random, shift_Random, labels = c('A','B'))
second_row = plot_grid(pRandom, labels = c('C'), nrow = 1)
gg_all = plot_grid(first_row, 
    second_row, 
    labels = c('', ''), 
    ncol = 1,
    rel_heights = c(1,1.25) )


pdf(file="FigureS11.pdf",10,12)
gg_all

dev.off()

