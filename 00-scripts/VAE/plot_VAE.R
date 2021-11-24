#!/usr/bin/env Rscript

#Author: QR
#Date:   October 2021
#Purpose:Plot VAE results, look at correlation with latitude 

library(dplyr);library(cowplot);library(ggplot2);library(data.table);library(tidyr);library(ggsci)

pd <- read.table("fulltest_latent_coords.txt",T)
names(pd)[1:2] <- c("LD1","LD2")
meta <- read.table("metadata",T, sep="\t" ) %>% select(POP, Latitude, Longitude, Region)
pd1 <- separate(pd,sampleID,c("POP","ind"),"_")
pd <- merge(pd1, meta ) #.x="POP",by.y="POP_ID")

############################################################################################### 
# test correlation  with geography:
summary(lm(pd$LD1 ~ pd$Latitude))
summary(lm(pd$LD2 ~ pd$Latitude))
summary(lm(pd$LD2 ~ pd$Longitude))
summary(lm(pd$LD1 ~ pd$Longitude))
summary(lm(pd$LD1 ~ pd$Latitude + pd$Longitude ))
summary(lm(pd$LD1 ~ pd$Latitude + pd$Longitude + pd$Latitude * pd$Longitude ))

#on prepare les couleurs
myColors <- c("blue","orange","red","darkviolet","chocolate4", "springgreen4","green")
names(myColors) <- levels(pd$Region)
colScale <- scale_colour_manual(name = "Region",values = myColors)

#on fait un plot des 2 ld
p <- pd %>% 
  dplyr::rename_at("POP",~"river") %>%
  #separate(REGION2,c("Rivers",NA), sep ="-" ) #)
  #ggplot(aes(x=LD1,y=LD2, col= river))+
  ggplot(aes(x = LD1, y = LD2, col = Region))+
  colScale +
  geom_point(size=1) + theme_classic()

haida <- dplyr::filter(pd, Region=="HaidaGwaii") %>% 
  dplyr::rename_at("POP",~"river") %>%
  ggplot(aes(x=LD1,y=LD2, col=river))+ #labs(title = "B") +
  geom_text(aes ( label= river), size=3.5) + theme_classic() +
  theme(legend.position = "none")  
h <- haida + scale_color_igv()

all <- ggdraw(p + theme_half_open(12)) +
    draw_plot(h , .62, .62, .35, .35) +
        draw_plot_label(
         c("A", "B"),
         c(0, 0.45),
         c(1, 0.95),
         size = 12)

pdf(file="vae_LD1_LD2.pdf")
all 
dev.off()
