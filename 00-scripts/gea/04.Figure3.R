#purpose: script to reproduce the Fig3 of our paper
#author: QR
#date: October 2021

## load libs
libs <- c('dplyr','reshape2','ade4','data.table', 'magrittr', 
    'factoextra','vegan','ggplot2', 'cowplot','CMplot',
    'gridGraphics')

invisible(lapply(libs, library, character.only = TRUE))

###### DOWLOAD SNP frequency data
freq <- fread("zcat plink.frq.strat.gz") #obtained from plink
freq2 <- dplyr::select(freq,SNP,CLST,MAF)
freq3 <- reshape2::dcast(freq2,SNP~CLST)
freq3 <- freq3[,-1]
freq4 <- t(freq3)
pop <- unique(freq2$CLST)
pop <- data.frame(pop)
snps <- unique(freq2$SNP)
colnames(freq4) <- snps

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

#ensuire that environmental variable will match the genotype data order:
colnames(pop) <- "POP"
env1 <- left_join(pop, env1) 

#choose if I want to work with elevation as well interaction dist*elev
env <- select(env1, -POP, -Region)

colnames(env) <- c("Latitude","normalized_distance", 
    "TemperaturePC1", "TemperaturePC2", "TemperaturePC3", "TemperaturePC4",
    "PrecipitationPC1", "PrecipitationPC2", "PrecipitationPC3", 
    "geology")

#all.rda<- rda(freq4 ~ 1 , data = env, scale =T)
rda1 <- rda(freq4 ~ geology +
    TemperaturePC1 + TemperaturePC2 + TemperaturePC3 + TemperaturePC4 + 
    PrecipitationPC1 + PrecipitationPC2 + PrecipitationPC3 + normalized_distance +
    Condition(Latitude) , 
    data=env, scale=T)

########################################################################
#### now do the plot of individuals/SNPs ###############################
########################################################################
env1$Region[env1$Region=="Washington&Oregon"] <- "Cascadia"

eco <- factor(env1$Region)

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
#col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
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
    '#a6cee3',
    '#ffff33',
    '#ff7f00', 
    '#fb9a99',
    '#e31a1c', 
    '#b2df8a')

table(col.pred$predictor)
#normalized_distance    PrecipitationPC1    PrecipitationPC3      TemperaturePC3 
#               1170                  19                  10                  35 
#     TemperaturePC4 
#                 88 

legsnp <- colnames(env)[c(-1,-4,-8,-10)] #remove latitude and geology


## Figure 3D
#distribution of outliers
snp <- "population_for_GEA.renamed.vcfsnp"
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

##    faire le plot :
myColors <- rep(c("darkviolet"),nrow(outlier_distribution))

outlier_distribution$CHR<- factor(outlier_distribution$CHR, 
    levels = levels(outlier_distribution$CHR))

names(myColors) <- levels(outlier_distribution$CHR)
colScale <- scale_colour_manual(name = "CHR",values = myColors)

yleg <- "RDA outliers"

library(ggplot2)
p <- ggplot(outlier_distribution, aes(x = LEN, y = COUNT) ) +
  stat_smooth(method="lm") +
  geom_point(aes(color=factor(CHR) ) )+
  #geom_line(aes(x=LENGTH,y=OUTLIERS)) + 
  colScale +
  ylab(yleg) +
  xlab("chromosome length - bp") +
  theme(legend.position="none") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
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

p

## Figure 3E: Epas1
#subseter le freq4 pour garder le candidat:
##Epas1:OkisLG25 28154607 4244396:56:+ A T . PASS NS=8867;AF=0.117 
epas <- filter(freq2, SNP == "4244396:56:+")
epas <- merge(epas, env1, by.x = "CLST", by.y = "POP")

summary(lm(epas$MAF ~ epas$normalized_distance))

epas$Region[epas$Region=="Washington&Oregon"] <- "Cascadia"
#div_dist$REGION=droplevels(div_dist$REGION)   

myColors <- c("blue","orange","red","green",
              "darkviolet",
              "springgreen4")
names(myColors) <- levels(epas$Region)
colScale <- scale_colour_manual(name = "Region",values = myColors)

ep <- ggplot(epas, aes(x = normalized_distance, y = MAF) ) +
  stat_smooth(method="lm") +
  geom_point(aes(color=factor(Region) ), size = 2 )+
  #geom_line(aes(x=LENGTH,y=OUTLIERS)) + 
  colScale +
  ylab("Derived Allele Frequency") +
  xlab("normalized distance") +
  theme(legend.position="none") +
  #theme(axis.title.x=element_blank(),
  #    axis.text.x=element_blank(),
  #    axis.ticks.x=element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
ep  <- ep + theme_bw()
ep  <- ep + theme(legend.position="none") + ylim(0,0.9)
ep  <- ep + theme(axis.title.x=element_text(size=16, family="Helvetica",face="bold"),
               axis.text.x=element_text(size=12,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
               axis.title.y=element_text(size=16, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
               axis.text.y=element_text(size=12,family="Helvetica",face="bold"),
               strip.text.x = element_text(size=12))
ep <- ep + theme(panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())
ep <- ep + annotate(geom="text", x=0.3, y=0.6, label= "paste(italic(R) ^ 2, \" = .41 p < 0.001***\")", parse = TRUE,
                  color="black", size=8)
ep

### Figure 4C : manahatan plot ###
#### données LFMM
pval_loc <- fread("zcat 04-RDA_LFMM/adjust_pvaluesBH_lfkmm_K.5.txt.gz")

#we will remove scaffold to plot only the 30 chromosome + the 8 dups.
pvalsub <- pval_loc %>% filter(!(grepl("scaff", V1)))
pvalsub$V1<- gsub("Okis","",pvalsub$V1)
pvalsub$V1<- gsub("LG","",pvalsub$V1)
pvalsub$V1<- gsub("_hom","",pvalsub$V1)

distplot = select(pvalsub,V3,V1,V2,normalized_distance) %>% set_colnames(.,c("SNP","Chromosome","Position","distance"))
distplot <- as.data.frame(distplot)

#CMplot(distplot,type="p",plot.type="m",LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,mar = c(5,6,3,3),
#       file.output=F,verbose=TRUE,width=16,height=6,chr.labels.angle=45)
#p5 <- recordPlot()

library(gridGraphics)
library(CMplot)

#CMplot(distplot,type="p",plot.type="m",LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,mar = c(5,6,3,3),
#       file.output=TRUE,verbose=TRUE,width=16,height=6,chr.labels.angle=45)

CMplot(distplot, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=NULL, mar = c(5,8,3,3),
       chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),signal.cex=c(1.5,1.5),
       signal.pch=c(19,19),file="pdf",memo="",dpi=300,
       file.output=F,
       verbose=TRUE,
       width=14,height=6, chr.labels.angle=45)

p5 <- recordPlot()


#########   now put everything together ########################################
#library(gridBase)
#library(grid)
#ex: https://stackoverflow.com/questions/14124373/combine-base-and-ggplot-graphics-in-r-figure-window

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

#plot.new()                                  # Create empty plot
#p3 <- recordPlot()
#plot.new()                                  # Create empty plot
#p4 <- recordPlot()

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
plot_grid(top_row, midlle_row, bottom_row, ncol = 1, aligned = "v")
dev.off()
