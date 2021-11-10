
library(dplyr);library(data.table);library(lme4);library(MuMIn)

####################################################################################################
#### loading recombination data ##############
recomb <-read.table("ALL.recomb250kb.gz")
colnames(recomb) <- c('POP','CHR','START','MID','END','N','RECOMB')
m_rec <- recomb %>% group_by(CHR,START,END) %>% summarise(mean=mean(RECOMB), n =mean(N))
#length(which(m_rec$n <100)) #56 avec moins de 100 SNPs par windows
m_rec <- m_rec %>% filter(n >100) 
pop_rec = recomb %>% group_by(POP,CHR) %>% summarise(mean=mean(RECOMB), n =mean(N)) 
pop_rec_chr = recomb %>% group_by(POP) %>% summarise(mean=mean(RECOMB), n =mean(N)) 

#on voit que la recomb est significativement réduite dans la pop SALMO et l'inverse et vrai pour TSO
#on peut garder les pop qui sont dans la moyenne +/1 SD de distrib pour faire un LD "moyen"
mean(pop_rec_chr$mean)+sd(pop_rec_chr$mean)
mean(pop_rec_chr$mean)-sd(pop_rec_chr$mean)

#sur cette base on exclue donc TSO et SALMO .....
m_rec2 = recomb %>% 
    filter(POP !="TSO" & POP !="SAL") %>%
    group_by(CHR,START,END) %>% 
    summarise(mean=mean(RECOMB), n =mean(N))

m_rec2$START = m_rec2$START*1000
m_rec2$END =   m_rec2$END*1000
m_rec2$CHR <- gsub("Okis","",m_rec2$CHR) #reshape name to have matching name with the RDA

rho_chro <- m_rec2 %>%
    group_by(CHR) %>%
    summarise(rhomean = mean(mean),
        rho_m2sd = mean(mean) - 5*(sd(mean)/sqrt(length(mean))),
        rho_p2sd = mean(mean) + 5*(sd(mean)/sqrt(length(mean))) ,
    n = n() )
#       rho_m2sd = mean(mean) - 8*(sd(mean)/sqrt(length(mean))),
#        rho_p2sd = mean(mean) + 8*(sd(mean)/sqrt(length(mean))) ,

recomb <- left_join(m_rec2, rho_chro, by="CHR")
recomb$state <- ifelse(recomb$mean < recomb$rho_m2sd,"coldspot" ,
        ifelse(recomb$mean > recomb$rho_p2sd,"hotspot","normal"))
table(recomb$state)

###################### WORK ON OUTLIERS ###########################################################
##### Now we will test the influence of recombination on outliers #############
#RDA results 
RDA <- read.table("01.RDA_LFMM/candidate_outliers_with_var3_with_chr.txt",T)
colnames(RDA)[1] <- "SNP"
colnames(RDA)[2] <- "CHR"
colnames(RDA)[3] <- "POS"
RDA$CHR<- gsub("Okis","",RDA$CHR)  

## read lfmm results:
lfmm <- read.table("01.RDA_LFMM/significant_outlier_control_forLatitudeK20.txt" ,T)
lfmm$CHR<- gsub("Okis","",lfmm$CHR)  
rdalfmm <- merge(RDA,lfmm) 
RDA <- rdalfmm
########################## COMBINE OUTLIER AND RECOMBINATION ######################################
x <- data.table(recomb[,c(1:3)])
x <- setkey(x)

RDA <- filter (RDA, !grepl("scaf",CHR))
y <- data.table(select(RDA, CHR, POS))
z  <- data.table(select(lfmm, CHR, POS))
w <- data.table(select(rdalfmm, CHR, POS))
z[, POS2 := POS]
y[, POS2 := POS]
w[, POS2 := POS]
xy <- foverlaps(y, x, by.x = names(y),type = "any", mult = "all", nomatch = 0L)
xz <- foverlaps(z, x, by.x = names(y),type = "any", mult = "all", nomatch = 0L)
xw <- foverlaps(w, x, by.x = names(y),type = "any", mult = "all", nomatch = 0L)

RDA1 = left_join(RDA, xy ,by=c("CHR","POS")) #merge to have the position of the recomb data 
RDA1 = left_join(RDA1,recomb,by=c("CHR","START","END")) #merge RDA with recomb data
RDA_rec = RDA1 %>% group_by(CHR) %>% summarise(mean=mean(mean))
#LFMM = 

#perform the chi-square to test for Randomness of distribution:
dat_state <- rbind(table(recomb$state), table(RDA1$state))
(Xsq <- chisq.test(dat_state))


#X-squared = 21.311, df = 2, p-value = 2.357e-05
snp <- read.table("../01.RDA_LFMM/20.DEFINITIF/02-data/genotype/population_for_GEA.vcfsnp")[,c(1:3)]

colnames(snp) <- c("CHR","POS","SNP")
snp$CHR<- gsub("Okis","",snp$CHR)  
RDA_ID=RDA[,c(1:3)]
RDA_ID=unique(RDA_ID)
RDA_ID$outlier <- "outlier"

snp2 = full_join(snp,RDA_ID)

y <- data.table(snp2[,c(1,2)])
y[, POS2 := POS]
xy <- foverlaps(y, x,
    by.x = names(y),
    type = "any",
    mult = "all",
    nomatch = 0L)

snp_out <- left_join(snp2, xy)
snp_out <- left_join(snp_out,recomb,by=c("CHR","START","END")) 
snp_out <- filter (snp_out, !grepl("scaf",CHR))
snp_out$outlier[is.na(snp_out$outlier)] <- "neutral"

wilcox.test(snp_out[snp_out[,4]=="outlier",8], snp_out[snp_out[,4]=="neutral",8])
#W = 5068012, p-value = 8.59e-06
#W = 4249648, p-value = 0.0003289

wilcox.test(snp_out[,8] ~ snp_out[,4])
#W = 5854058, p-value = 0.0003289

snp_out %>%group_by(outlier) %>%summarise(mean=mean(mean))
setwd("~/postdoc/coho/10K_coho/EXPLO/FULL/")
m_out <- filter(snp_out, outlier=="outlier") %>% summarise(mean=mean(mean), n=n())
m_neut <- filter(snp_out, outlier=="neutral") %>% summarise(mean=mean(mean), n=n())

#faire un plot des moyennes
library(ggplot2)
#violin plot:
p <- ggplot(snp_out, aes(x=outlier,y=mean, fill=outlier)) + 
  geom_violin() + geom_boxplot(width=0.1)+  #scale_fill_manual(values=c("#999999", "#E69F00")) +
  scale_fill_brewer(palette="RdBu") + theme_minimal() +
  theme(legend.position="none") +
  #theme(strip.background = element_blank(),strip.text.x = element_blank() )+
  ylab("mean recombination rate") + xlab("SNPs status") 
p <- p + stat_summary(fun.data=mean_sdl, mult=1, 
             geom="pointrange", color="red")
pdf(file="differenceinrecombination.pdf")
p
dev.off()

  #geom_jitter(width = 0.4)
RDAout <- filter(snp_out, outlier=="outlier")
RDAneut <- filter(snp_out, outlier=="neutral")
t.test(RDAout$mean,RDAneut$mean)
wilcox.test(RDAout$mean,RDAneut$mean, pairded=F,alternative="less", conf.int=T)

library(lme4)
library(MuMIn)
snp_out$out <- ifelse(snp_out$outlier=="neutral",0,1)

res0 <- glmer(snp_out$out ~  1+ (1 | snp_out$CHR), family=binomial(link="logit"))
res1 <- glmer(snp_out$out ~ scale(snp_out$mean) +
     (1 |snp_out$CHR), family=binomial(link="logit"))
anova(res0,res1)
r.squaredGLMM(res1)

mchr_snp_out <- snp_out %>% group_by(CHR,outlier) %>% summarise(mea_rho=mean(mean), n=n())

