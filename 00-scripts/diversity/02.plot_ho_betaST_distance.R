

library(ggplot2)
library(dplyr)
library(magrittr)
library(data.table)
library(cowplot)
rm(list=ls())


hierfs <- fread("zcat hierfstats.hs.gz")

hs  <-  data.frame((colMeans(hierfs[,-1], na.rm=T))) %>% set_colnames(.,c("hs"))
hs$POP <- rownames(hs) 
hs[which(hs$POP == "DRA" ), 2] <- "WAC"

beta <- read.table("beta_iovl_ci",T) %>% select(beta) %>% set_colnames(.,c("beta"))

beta$POP <- rownames(beta)
beta[which(beta$POP == "DRA" ), 2] <- "WAC"

meta <- read.table("metadata.txt", T, sep = "\t") %>% 
    select(POP, Region, dist_max_km)
mouth_dist <- read.table("embouchure_coho_dist.txt.gz", T) %>% 
    select(-LON_emb, -LAT_emb)
meta <- merge(meta, mouth_dist, by.x = "POP", by.y = "SITE")
meta$dist_total <- meta$dist_max_km + meta$dist_SCO_3k

#combine with hs and beta:
divall <- merge(merge(hs, meta) , beta)
summary(lm(divall$dist_total ~ divall$beta))

#remove the sample from Russia as it does not follow the pattern:
d <- divall %>% filter(POP !="SAI")
summary(lm(d$dist_total ~ d$beta))

## prepare data for plotting:
divall$Region[divall$Region=="Washington&Oregon"] <- "Cascadia"
myColors <- c("blue","orange","red","green","darkviolet","brown","springgreen4")
names(myColors) <- levels(divall$Region)
colScale <- scale_colour_manual(name = "Region",values = myColors)


############## now plot the data ######################################################
p <- ggplot(data = divall,
    aes(x = dist_total, y = hs, label = POP))
p <- p + geom_vline(colour = "darkgrey",linetype = 5, size = 3,
    data = NULL, xintercept = 1050)
p <- p + stat_smooth(method = "lm")
p <- p + geom_point(aes(colour = factor(Region)),size = 2)
#p <- p + geom_text(aes(label = POP), size = 2)
p <- p + colScale + theme_bw()
p <- p + labs(x ="Distance from southernmost site (km)", y = "Hs")  + 
     expand_limits(y = 0.005, x = 8300) +   scale_x_continuous(limits = c(0, 8300))
p <- p + theme(axis.title.x = element_text(size = 12, family = "Helvetica", face ="bold"),
     axis.text.x = element_text(size = 12, family = "Helvetica",
     face = "bold", angle = 90, hjust = 0, vjust = 0.5),
     axis.title.y = element_text(size = 12, family = "Helvetica",
     face = "bold",angle = 0, hjust = 0, vjust = 0.5),
     axis.text.y = element_text(size = 12,family = "Helvetica",face ="bold"),
     strip.text.x = element_text(size = 12))
p  <- p + theme(panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())
p <- p + annotate(geom = "text", x = 3800, y = 0.01,
     label = "paste(italic(R) ^ 2, \" = .4  p < 2e-16, slope = -1.65e-06\")",
     parse = TRUE, color = "black", size = 3)
p <- p + theme(legend.text = element_text(size = 8, colour = "black", family="Helvetica"))
p <- p + theme(legend.title = element_text(size = 8, face = "bold"))
p <- p + theme(legend.position = c(.01, .002),
    legend.justification = c("left", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(4,4,4,4)
    )

p


pdf(file = "Ho_corrected_distances_with_labels.pdf", 20, 6)
p
dev.off()

#p <- p + geom_text(aes(label=POP_ID),hjust=0, vjust=0,size=3)
pdf(file="Ho_corrected_distance_withlabels.pdf",18,10)
p + geom_text(aes(label=POP_ID),hjust=0, vjust=0,size=1)
dev.off()


d <- ggplot(data = divall,
    aes(x = dist_total, y = beta))
d <- d + geom_vline(colour = "darkgrey",linetype = 5, size = 3,
    data = NULL, xintercept = 1050)
d <- d + stat_smooth(method = "lm")
d <- d + geom_point(aes(colour = factor(Region)) , size = 2)
d <- d + colScale + theme_bw()
d <- d + labs(x = "Distance from southernmost site (km)", y = expression(beta["ST"]) )  + 
         expand_limits(y = 0.005, x = 8300) +  scale_x_continuous(limits = c(0, 8300))
d <- d + theme(axis.title.x = element_text(size = 12, family = "Helvetica",face = "bold"),
              axis.text.x = element_text(size = 12,family = "Helvetica",
              face = "bold", angle = 90, hjust = 0, vjust = 0.5),
              axis.title.y = element_text(size = 12, family = "Helvetica",
              face = "bold",angle = 0, hjust = 0, vjust = 0.5),
              axis.text.y = element_text(size=12,
              family = "Helvetica",face="bold"),
              strip.text.x = element_text(size=12))
d <- d + theme(panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())
d <- d + annotate(geom = "text",x = 4000, y = 0.01,
     label = "paste(italic(R) ^ 2, \" = .39, p < 2e-16, slope = 3.36e-05\")",
     parse = TRUE,
     color = "black", size = 3)
d  <- d + theme(legend.position = "none")

pdf(file="figure2.pdf", width = 20, height = 8)
plot_grid(
    p,
    d, 
    labels = c("A", "B"), ncol = 2,
    label_size = 12)
dev.off()
