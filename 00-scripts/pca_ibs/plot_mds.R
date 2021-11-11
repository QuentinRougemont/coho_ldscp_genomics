

libs <- c('ggplot2','data.table', 'magrittr','dplyr')
invisible(lapply(libs, library, character.only = TRUE))

strata <- read.table("01-info/pop_ind_region_latitute_corrected3.txt") 
strata <- select(strata, V2,V6) %>% 
    set_colnames(.,c("FID","REGION"))
strata <- unique(strata)
strata <- strata[order(strata$FID),]  #reorder because population DRA has been renamed into WAC

mds = fread("zcat mds_to_plot.mds.gz")
#for IBM
#mds = fread("zcat ibm_mds_to_plot.mds.gz")

strata = left_join(mds, strata)

#personalizing the colors:
myColors <- c("blue","orange","red","darkviolet","chocolate4", "springgreen4","green")
names(myColors) <- levels(strata$REGION)
colScale <- scale_colour_manual(name = "REGION",values = myColors)


#Full data
p <- ggplot(mds) + 
  geom_point(aes(x=C1, y=C2 ,  colour=factor(strata$REGION)), alpha =0.2) + #plot les points noirs
  colScale + #add la couleur suivant le dataset
  theme_bw() 
p <- p +theme(legend.text = element_text(size = 10, face = "bold")) +
    theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
    colScale + 
    ggtitle("MDS on IBS") + 
     theme(plot.title = element_text(lineheight=.8, face="bold" , hjust = 0.5))

pdf(file="mds_point.pdf",12,12)
p
dev.off()

