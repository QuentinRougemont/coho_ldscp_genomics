library(data.table)
library(dplyr)
library(magrittr)
library(CMplot)
rm(list=ls())
ALB_LCA <- fread("zcat all_ALB_LCA_SAI.OkisLG01.pbs.50.txt.gz") #1 Thompson
POR_SHK <- fread("zcat all_POR_SHK_SAI.OkisLG01.pbs.50.txt.gz") #2 alaska
SFK_GOO <- fread("zcat all_SFK_GOO_SAI.OkisLG01.pbs.50.txt.gz") #3 alkaska #2
SUS_OON <- fread("zcat all_SUS_OON_SAI.OkisLG01.pbs.50.txt.gz") #4 BC North

ALB_LCA  <- dplyr::select(ALB_LCA, -Fst02, -Fst12, -nSite, -V3, -PBS1, -PBS2)
POR_SHK  <- dplyr::select(POR_SHK, -Fst02, -Fst12, -nSite, -V3, -PBS1, -PBS2)
SFK_GOO  <- dplyr::select(SFK_GOO, -Fst02, -Fst12, -nSite, -V3, -PBS1, -PBS2)
SUS_OON  <- dplyr::select(SUS_OON, -Fst02, -Fst12, -nSite, -V3, -PBS1, -PBS2)

### Effectuer un manathan plots:
ALB_LCA2 <- dplyr::select(ALB_LCA, chr, midPos, PBS0)
POR_SHK2 <- dplyr::select(POR_SHK, chr, midPos, PBS0)
SFK_GOO2 <- dplyr::select(SFK_GOO, chr, midPos, PBS0)
SUS_OON2 <- dplyr::select(SUS_OON, chr, midPos, PBS0)

ALB_POR <- full_join(ALB_LCA2, POR_SHK2, by=c("chr","midPos"))
ALB_POR_SFK <- full_join(ALB_POR, SFK_GOO2, by=c("chr","midPos"))
ALB_POR_SFK_SUS <- full_join(ALB_POR_SFK, SUS_OON2, by=c("chr","midPos"))

distplot <- ALB_POR_SFK_SUS
distplot$SNP <- seq(1:nrow(distplot) )

distplot <- dplyr::select(distplot, SNP,chr,midPos,
     PBS0.x, PBS0.y, PBS0.x.x, PBS0.y.y) %>% 
     set_colnames(.,c("SNP","Chromosome","Position","ALB","POR","SFK","SUS"))

distplot <- filter(distplot, !grepl("scaf", distplot$Chromosome))
#distplot <- filter(distplot, !grepl("hom", distplot$Chromosome))
distplot$Chromosome <- gsub("Okis","",distplot$Chromosome)
distplot[distplot<0] <- 0

CMplot(distplot,type="p",plot.type="m", multracks=TRUE,
    LOG10=F,threshold=NULL,memo="",
    ylab="PBS",cex=0.6,cex.lab=0.8,
    file.output=TRUE,file="pdf",dpi=300,
    verbose=TRUE,width=20,height=4,chr.labels.angle=45)

################## FST CMPLOT ##################################################
ALB_LCA2 <- dplyr::select(ALB_LCA, chr, midPos, Fst01)
POR_SHK2 <- dplyr::select(POR_SHK, chr, midPos, Fst01)
SFK_GOO2 <- dplyr::select(SFK_GOO, chr, midPos, Fst01)
SUS_OON2 <- dplyr::select(SUS_OON, chr, midPos, Fst01)

ALB_POR <- full_join(ALB_LCA2, POR_SHK2, by=c("chr","midPos"))
ALB_POR_SFK <- full_join(ALB_POR, SFK_GOO2, by=c("chr","midPos"))
ALB_POR_SFK_SUS <- full_join(ALB_POR_SFK, SUS_OON2, by=c("chr","midPos"))

distplot <- ALB_POR_SFK_SUS
distplot$SNP <- seq(1:nrow(distplot) )

distplot <- dplyr::select(distplot, SNP,chr,midPos,
            Fst01.x, Fst01.y, Fst01.x.x, Fst01.y.y ) %>% 
	set_colnames(.,c("SNP","Chromosome","Position" ,"ALB","POR","SFK","SUS"))

distplot <- filter(distplot, !grepl("scaf", distplot$Chromosome))
distplot$Chromosome <- gsub("Okis","",distplot$Chromosome)
distplot[distplot<0] <- 0

CMplot(distplot,type="p",plot.type="m", multracks=TRUE,
       LOG10=F,threshold=NULL,memo="",
       ylab="Fst",cex=0.6,cex.lab=0.8,
       file.output=TRUE,file="pdf",dpi=300,
       verbose=TRUE,width=20,height=4,chr.labels.angle=45)
