
## load libs
libs <- c('dplyr','reshape2','ade4','data.table', 'magrittr', 'factoextra','vegan','ggplot2', 'cowplot','CMplot', 'gridGraphics')
invisible(lapply(libs, library, character.only = TRUE))

pval_loc <- fread("zcat LFMM/adjust_pvaluesBH_lfkmm_K.20.txt.gz")

#we will remove scaffold to plot only the 30 chromosome + the 8 dups.
pvalsub <- pval_loc %>% filter(!(grepl("scaff", V1)))
pvalsub$V1<- gsub("Okis","",pvalsub$V1)
pvalsub$V1<- gsub("LG","",pvalsub$V1)
pvalsub$V1<- gsub("_hom","",pvalsub$V1)

#temperature
tempplot <- select(pvalsub,V3,V1,V2,Temperature1, Temperature2, Temperature3, Temperature4) %>% 
    set_colnames(.,c("SNP","Chromosome","Position",
    "TemperaturePCaxis1","TemperaturePCaxis2","TemperaturePCaxis3","TemperaturePCaxis4"))

templot <- as.data.frame(tempplot)

#precipitation
precplot <- select(pvalsub,V3,V1,V2,Precipitation1, Precipitation2, Precipitation3) %>% 
    set_colnames(.,c("SNP","Chromosome","Position",
    "PrecipitationPCaxis1","PrecipitationPCaxis2","PrecipitationPCaxis3"))

precplot <- as.data.frame(precplot)

#geology
geolplot <- select(pvalsub,V3,V1,V2,Geology) %>% 
    set_colnames(.,c("SNP","Chromosome","Position","Geology"))

geolplot <- as.data.frame(geolplot)


CMplot(geolplot, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=NULL ,
       chr.den.col = c("darkgreen", "yellow", "red"),
       signal.pch = c(19,19),file="pdf",memo="",dpi=300,
       file.output = F,
       verbose = TRUE,
       width=14,height=6, chr.labels.angle=45)

p3 <- recordPlot()
CMplot(templot, plot.type="m", multracks = TRUE, LOG10=TRUE, ylim=NULL, threshold=NULL, mar = c(5,8,3,3),
       chr.den.col = c("darkgreen", "yellow", "red"),
       signal.pch = c(19,19),file="pdf",memo="",dpi=300,
       file.output = F,
       verbose = TRUE,
       width=14,height=6, chr.labels.angle=45)
p1 <- recordPlot()


CMplot(precplot, plot.type="m", multracks = TRUE, LOG10=TRUE, ylim=NULL, threshold=NULL, mar = c(5,8,3,3),
       chr.den.col = c("darkgreen", "yellow", "red"),
       signal.pch = c(19,19),file="pdf",memo="",dpi=300,
       file.output = F,
       verbose = TRUE,
       width=14,height=6, chr.labels.angle=45)
p2 <- recordPlot()

#export final figure
pdf(file="FigS09LFMM.pdf",15,20)
plot_grid(p1, p2,p3, 
    ncol = 1, 
    align = "v", 
    labels=c("A - Temperature", "B - Precipitation", "C - Geology") )
dev.off()
