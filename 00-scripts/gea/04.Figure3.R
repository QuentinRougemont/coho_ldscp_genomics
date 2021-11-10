#purpose: script to reproduce the Fig3 of our paper
#author: QR
#date: October 2021

## load libs
libs <- c('dplyr','reshape2','ade4','data.table', 'magrittr', 'factoextra','vegan','ggplot2', 'cowplot','CMplot', 'gridGraphics')
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

#verify that the order of pop in env match the genotype matrix:
setdiff(env2$SITE, pop$SITE)  #there should be no difference 
setdiff(pop$SITE,env2$SITE)   #idem
table(pop$SITE==env2$SITE)    #if false occur then there's difference

#remove unwanted variable :
env <- dplyr::select(env2, -SITE, -elevation, -Region, -dist_max_km)
colnames(env) <- c( "geology", 
    "TemperaturePC1", "TemperaturePC2", "TemperaturePC3", "TemperaturePC4",
    "PrecipitationPC1", "PrecipitationPC2", "PrecipitationPC3",
    "Latitude","normalized_distance" )

#all.rda<- rda(freq4 ~ 1 , data = env, scale =T)
rda1 <- rda(freq4 ~ geology +
    TemperaturePC1 + TemperaturePC2 + TemperaturePC3 + TemperaturePC4 + 
    PrecipitationPC1 + PrecipitationPC2 + PrecipitationPC3 + normalized_distance +
    Condition(Latitude) , 
    data=env, scale=T)

########################################################################
#### now do the plot of individuals/SNPs ###############################
########################################################################
env2$Region[env2$Region=="Washington&Oregon"] <- "Cascadia"

eco <- factor(env2$Region)
bg <- c("Darkblue","Orange","Red","green","Darkviolet","Springgreen4")

tmp <- round(summary(eigenvals(rda1, model= "constrained"))[2]*100,2)   
RDA1var <- paste("RDA1 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[5]*100,2)   
RDA2var <- paste("RDA2 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[8]*100,2)   
RDA3var <- paste("RDA3 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[11]*100,2)   
RDA4var <- paste("RDA4 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[14]*100,2)   
RDA5var <- paste("RDA5 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[17]*100,2)   
RDA6var <- paste("RDA6 ",tmp,"%",sep="")                                

###############   extract the scores and find outliers #########################
#the 7 RDA axis are significant (see script for signficance)
#and we keep them
signif <- 7
load.rda <- scores(rda1, choices = c(1:signif), display = "species") 

pdf(file="hist2.pdf")
for (i in 1:signif){ 
hist(load.rda[,i], main = paste("Loadings on RDA", i ,sep = ""))
}
dev.off()

##### find outliers 
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
foo <- matrix(nrow=(ncand), ncol=ncol(env))  
#foo <- matrix(nrow=(ncand), ncol=ncol(env[-1]))  
colnames(foo) <- colnames(env)
#colnames(foo) <- colnames(env[-1]) #to rm lat

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- freq4[,nam]
  foo[i,] <- apply(env,2,function(x) cor(x,snp.gen))
  #foo[i,] <- apply(env[,-1],2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)

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

cand12 <- dplyr::filter(cand, axis==1 |axis==2)
cand12 <- dplyr::filter(cand12,  predictor != "Latitude")
sel  <- data.frame(cbind(cand12$snp, cand12$predictor ))

colnames(sel) <- c("SNP","predictor")
# color by predictor:
col.pred <- rownames(rda1$CCA$v)  %>% 
  as.data.frame() %>% 
  set_colnames(., "SNP")# pull the SNP names

###############################################################################
col.pred <- left_join(col.pred , sel) #, by=c("SNP"="V1"))
col.pred$color <- 
                ifelse(col.pred$predictor=='geology',            '#1f78b4',
                ifelse(col.pred$predictor=='TemperaturePC1',     '#a6cee3' ,
                ifelse(col.pred$predictor=='TemperaturePC2',     '#33a02c',
                ifelse(col.pred$predictor=='TemperaturePC3',     '#ffff33',
                ifelse(col.pred$predictor=='TemperaturePC4',      '#ff7f00',
                ifelse(col.pred$predictor=='PrecipitationPC1',   '#fb9a99',
                ifelse(col.pred$predictor=='PrecipitationPC2',   '#6a3d9a',
                ifelse(col.pred$predictor=='PrecipitationPC3',   '#e31a1c',
                ifelse(col.pred$predictor=='normalized_distance','#b2df8a',
                  '#f1eef6') ))))))))
col.pred[is.na(col.pred)] <-'#f1eef6'
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")

#only a few predictor are associate to axes 1 & 2 so we don't need all colors actually
bgsnp <- c(
    '#1f78b4','#a6cee3',
    '#33a02c','#ffff33',
    '#ff7f00', '#fb9a99',
    '#6a3d9a','#e31a1c', 
    '#b2df8a')
    
bgsnp <- c(
    '#ffff33',
    '#ff7f00', 
    '#fb9a99',
    '#e31a1c', 
    '#b2df8a')

#table(col.pred$predictor)
#          #f1eef6 normalized_distance    PrecipitationPC1    PrecipitationPC3 
#            58174                1110                  17                   9 
#   TemperaturePC3      TemperaturePC4 
#               29                 114 

legsnp <- colnames(env)[c(-1,-2,-3,-7,-9)] #remove latitude and geology
################################################################################
## Figure 3D
#distribution of outliers
#snp <- "02-data/genotype/population_for_GEA.vcfsnp"
#without thompson
snp <- "02-data/genotype/population_for_GEA.renamed.nothompson.recode.vcfsnp"
snp <- read.table(snp) %>% 
    select(V1, V2, V3) %>% 
    set_colnames(., c("CHR","POS","snp"))

snp_cand <- merge(snp, cand)

outlier_distribution <- data.frame(table(snp_cand$CHR)) %>% 
    set_colnames(., c("CHR","COUNT"))

#correlation with chromosome length:
chromo_length <- read.table("length.txt.gz") %>% set_colnames(., c("CHR","LEN"))
outlier_distribution <- merge(outlier_distribution, chromo_length)

summary(lm(outlier_distribution$COUNT ~ outlier_distribution$LEN))

##    faire le plot ############################################################
myColors <- rep(c("darkviolet"),nrow(outlier_distribution))

outlier_distribution <- outlier_distribution  %>% filter(!(grepl("scaff", CHR)))
outlier_distribution$CHR<- droplevels(outlier_distribution$CHR)

outlier_distribution$CHR<- factor(outlier_distribution$CHR, 
    levels = levels(outlier_distribution$CHR))

names(myColors) <- levels(outlier_distribution$CHR)
colScale <- scale_colour_manual(name = "CHR",values = myColors)

yleg <- "RDA outliers"

library(ggplot2)
p <- ggplot(outlier_distribution, aes(x = LEN, y = COUNT) ) +
  stat_smooth(method="lm") +
  geom_point(aes(color=factor(CHR) ) )+
  colScale +
  ylab(yleg) +  xlab("chromosome length - bp") +
  theme(legend.position="none") +
  theme( strip.background = element_blank(),
    strip.text.x = element_blank()  )
p <- p + theme_bw()
p  <- p + theme(legend.position="none")
p <- p + theme(axis.title.x=element_text(size=16, family="Helvetica",face="bold"),
               axis.text.x=element_text(size=12,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
               axis.title.y=element_text(size=16, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
               axis.text.y=element_text(size=12,family="Helvetica",face="bold"),
               strip.text.x = element_text(size=12))
p <- p + theme(panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())
p <- p + annotate(geom="text", x=31000000, y=120, label= "paste(italic(R) ^ 2, \" = .90 p < 0.001***\")", parse = TRUE,
                  color="black", size=8)

### Figure 3C : manahatan plot ###
#### LFMM results ##############################################################
pval_loc <- fread("zcat LFMM/adjust_pvaluesBH_lfkmm_K.20.txt.gz")

#we will remove scaffold to plot only the 30 chromosome + the 8 dups.
pvalsub <- pval_loc %>% filter(!(grepl("scaff", V1)))
pvalsub$V1<- gsub("Okis","",pvalsub$V1)
pvalsub$V1<- gsub("LG","",pvalsub$V1)
pvalsub$V1<- gsub("_hom","",pvalsub$V1)

distplot = select(pvalsub,V3,V1,V2,normalized_distance) %>% set_colnames(.,c("SNP","Chromosome","Position","distance"))
distplot <- as.data.frame(distplot)

library(gridGraphics)
library(CMplot)

CMplot(distplot, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=NULL, mar = c(5,8,3,3),
       chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),signal.cex=c(1.5,1.5),
       signal.pch=c(19,19),file="pdf",memo="",dpi=300,
       file.output=F,
       verbose=TRUE,
       width=14,height=6, chr.labels.angle=45)

p5 <- recordPlot()

#########   now put everything together ########################################
plot(rda1, scaling = 3,  choices = c(1,2), type  = "n", cex.lab = 1,
     cex.axis = 1.2,
     xlab = RDA1var, ylab = RDA2var)
points(rda1, display = "sites",
       pch = 21, cex = 1.2,
       col = "gray32",
       scaling = 3, bg = bg[eco], choices = c(1,2))              
text(rda1, scaling = 3, display = "bp", col = "#0868ac", cex = 1.1)   
legend("bottomright", legend = levels(eco) ,bty = "n",
       col = "gray32", pch = 21, cex = 1.25, pt.bg = bg)

p1 <-  recordPlot()
#par(mfrow = c(1,2))
plot(rda1, type = "n", scaling = 3,
     cex.lab = 1.2,
     xlab = RDA1var, ylab = RDA2var,
     xlim = c(-1,1), ylim = c(-1,1))
points(rda1, display = "species",
       pch = 21, cex = 1, col = "gray32", bg = col.pred$color, scaling = 3)
text(rda1, scaling = 3, display = "bp", col = "#0868ac", cex=1.1)
legend("topleft", legend = legsnp[c(2,3,4,5,1)],
       bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bgsnp)
p2 <-  recordPlot()

top_row <- plot_grid(p1, p2,
  labels = c("A", "B"),
  #rel_widths = c(1.5, 1.5 ),
  #rel_heigths = c(1.5, 1.5),
  nrow = 1
)

midlle_row <- plot_grid(p5, labels=c("C"))

bottom_row <- plot_grid(p, ep,
  labels = c("D", "E"),
  rel_widths = c(1, 1),
  nrow = 1
)


#export final figure
pdf(file="Fig3.pdf",15,20)
plot_grid(top_row, midlle_row, bottom_row, ncol = 1, align = "v")
dev.off()



################################################################################
#### Axis 34 et 56
cand34 <- dplyr::filter(cand, axis==3 |axis==4)
cand34 <- dplyr::filter(cand34,  predictor != "Latitude")
sel34  <- data.frame(cbind(cand34$snp, cand34$predictor ))

cand56 <- dplyr::filter(cand, axis==5 |axis==6)
cand56 <- dplyr::filter(cand56,  predictor != "Latitude")
sel56  <- data.frame(cbind(cand56$snp, cand56$predictor ))
colnames(sel56) <- colnames(sel34) <- c("SNP","predictor")

# color by predictor:
col.pred34 <- rownames(rda1$CCA$v) %>% as.data.frame() %>% set_colnames(., "SNP") 
col.pred56 <- rownames(rda1$CCA$v) %>% as.data.frame() %>% set_colnames(., "SNP") 
col.pred34 <- left_join(col.pred34 , sel34) #, by=c("SNP"="V1"))
col.pred56 <- left_join(col.pred56 , sel56) #, by=c("SNP"="V1"))

col.pred34$color <- 
                ifelse(col.pred34$predictor=='geology',            '#1f78b4',
                ifelse(col.pred34$predictor=='TemperaturePC1',     '#a6cee3' ,
                ifelse(col.pred34$predictor=='TemperaturePC2',     '#33a02c',
                ifelse(col.pred34$predictor=='TemperaturePC3',     '#ffff33',
                ifelse(col.pred34$predictor=='TemperaturePC4',      '#ff7f00',
                ifelse(col.pred34$predictor=='PrecipitationPC1',   '#fb9a99',
                ifelse(col.pred34$predictor=='PrecipitationPC2',   '#6a3d9a',
                ifelse(col.pred34$predictor=='PrecipitationPC3',   '#e31a1c',
                ifelse(col.pred34$predictor=='normalized_distance','#b2df8a',
                  '#f1eef6') ))))))))
col.pred34[is.na(col.pred34)] <-'#f1eef6'
empty <- col.pred34
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
#only a few predictor are associate to axes 1 & 2 so we don't need all colors actually
table(col.pred34$predictor)    
bgsnp34 <- c(
    '#1f78b4','#a6cee3',
    '#33a02c','#ffff33',
    '#ff7f00', '#fb9a99',
    '#6a3d9a','#e31a1c', 
    '#b2df8a')
    

legsnp34 <- colnames(env)[c(-1,-7,-9)] #remove latitude and geology
col.pred56$color <- 
                ifelse(col.pred56$predictor=='geology',            '#1f78b4',
                ifelse(col.pred56$predictor=='TemperaturePC1',     '#a6cee3' ,
                ifelse(col.pred56$predictor=='TemperaturePC2',     '#33a02c',
                ifelse(col.pred56$predictor=='TemperaturePC3',     '#ffff33',
                ifelse(col.pred56$predictor=='TemperaturePC4',      '#ff7f00',
                ifelse(col.pred56$predictor=='PrecipitationPC1',   '#fb9a99',
                ifelse(col.pred56$predictor=='PrecipitationPC2',   '#6a3d9a',
                ifelse(col.pred56$predictor=='PrecipitationPC3',   '#e31a1c',
                ifelse(col.pred56$predictor=='normalized_distance','#b2df8a',
                  '#f1eef6') ))))))))
col.pred56[is.na(col.pred56)] <-'#f1eef6'
empty56 <- col.pred56
empty56[grep("#f1eef6",empty56)] <- rgb(0,1,0, alpha=0) # transparent
empty56.outline <- ifelse(empty56=="#00FF0000","#00FF0000","gray32")
#only a few predictor are associate to axes 1 & 2 so we don't need all colors actually
table(col.pred56$predictor)

bgsnp56 <- c(
    '#1f78b4','#a6cee3',
    '#33a02c','#ffff33',
    '#ff7f00', '#fb9a99',
    '#6a3d9a','#e31a1c', 
    '#b2df8a')
    
legsnp56 <- colnames(env)[c(-9)] #remove latitude and geology
#########   now plot axis34 and axis56 ########################################

#export final figure
#pdf(file="FigS11.pdf", 15, 15)
#par(mfrow = c(3,2)) #no thompson

pdf(file="FigS08.pdf", 15, 15)
par(mfrow = c(2,2)) #full

#âxes12: #only without thompson
#plot(rda1, scaling = 3,  choices = c(1,2), type  = "n", cex.lab = 1,
#     cex.axis = 1.2,
#     xlab = RDA1var, ylab = RDA2var)
#points(rda1, display = "sites",
#       pch = 21, cex = 1.2,
#       col = "gray32",
#       scaling = 3, bg = bg[eco], choices = c(1,2))              
#text(rda1, scaling = 3, display = "bp", col = "#0868ac", cex = 1.1)   
#legend("bottomright", legend = levels(eco) ,bty = "n",
#       col = "gray32", pch = 21, cex = 1.25, pt.bg = bg)
#plot(rda1, type = "n", scaling = 3,
#     cex.lab = 1.2,
#     xlab = RDA1var, ylab = RDA2var,
#     xlim = c(-1,1), ylim = c(-1,1))
#points(rda1, display = "species",
#       pch = 21, cex = 1, col = "gray32", bg = col.pred$color, scaling = 3)
#text(rda1, scaling = 3, display = "bp", col = "#0868ac", cex=1.1)
#legend("topleft", legend = legsnp[c(2,3,4,5,1)],
       bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bgsnp)

#axes 34
plot(rda1, scaling = 3,  choices = c(3,4), type  = "n", cex.lab = 1,
     cex.axis = 1.2,
     xlab = RDA3var, ylab = RDA4var)
points(rda1, display = "sites",
       pch = 21, cex = 1.2,
       col = "gray32",
       scaling = 3, bg = bg[eco], choices = c(3,4))
text(rda1, scaling = 3, display = "bp", col = "#0868ac", cex = 1.1)   
legend("topleft", legend = levels(eco) ,bty = "n",
       col = "gray32", pch = 21, cex = 1.25, pt.bg = bg)

plot(rda1, type = "n", scaling = 3,
     cex.lab = 1.2, choices = c(3,4),
     xlab = RDA3var, ylab = RDA4var,
     xlim = c(-1,1), ylim = c(-1,1) )
     
points(rda1, display = "species", choices = c(3,4), 
       pch = 21, cex = 1, col = "gray32", bg = col.pred34$color, scaling = 3)
text(rda1, scaling = 3, display = "bp", col = "#0868ac", cex=1.1)
legend("topleft", legend = legsnp34, #[c(2,3,4,5,1)],
       bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bgsnp34)

## axes 56
plot(rda1, scaling = 3,  choices = c(5,6), type  = "n", cex.lab = 1,
     cex.axis = 1.2,
     xlab = RDA5var, ylab = RDA6var)
points(rda1, display = "sites",
       pch = 21, cex = 1.2,
       col = "gray32",
       scaling = 3, bg = bg[eco], choices = c(5,6))
text(rda1, scaling = 3, display = "bp", col = "#0868ac", cex = 1.1)   
legend("topleft", legend = levels(eco) ,bty = "n",
       col = "gray32", pch = 21, cex = 1.25, pt.bg = bg)

plot(rda1, type = "n", scaling = 3,
     cex.lab = 1.2, choices = c(5,6),
     xlab = RDA5var, ylab = RDA6var,
     xlim = c(-1,1), ylim = c(-1,1) )
     
points(rda1, display = "species", choices = c(5,6), 
       pch = 21, cex = 1, col = "gray32", bg = col.pred56$color, scaling = 3)
text(rda1, scaling = 3, display = "bp", col = "#0868ac", cex=1.1)
legend("topleft", legend = legsnp56, #[c(2,5,6,5,1)],
       bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bgsnp56)

dev.off()


### then we repeat all the script in the folder without thompson samples
