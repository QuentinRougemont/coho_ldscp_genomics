file="plink.frq.strat.renamed.gz"

#load libs
libs <- c('dplyr','reshape2','ade4','data.table', 'factoextra', 'magrittr','ggsci')
invisible(lapply(libs, library, character.only = TRUE))

file <- paste0("zcat ",file) 
freq <- fread(file)
freq2 <- dplyr::select(freq,SNP,CLST,MAF) 
freq3 <- reshape2::dcast(freq2,SNP~CLST)
freq3 <- freq3[,-1] 
freq4 <- t(freq3)
pop <- unique(freq$CLST)   
pop <- data.frame(pop)

#perform PCA
pca1 <- dudi.pca(freq4,scale=FALSE,scannf=FALSE)

#load strata
strata <- read.table("pop_ind_region_latitute_corrected3.txt") 
strata <- select(strata, V2,V6) %>% 
    set_colnames(.,c("POP","REGION"))
strata <- unique(strata)
strata <- strata[order(strata$POP),]  #reorder because population DRA has been renamed into WAC

#personalizing the colors:
myColors <- c("blue","orange","red","darkviolet","springgreen4","green")
names(myColors) <- levels(strata$REGION)
colScale <- scale_colour_manual(name = "REGION",values = myColors)

#plot with river label:
p <- fviz_pca_ind(pca1, label="none", pointsize = 0.0) +
    geom_text(aes(label=strata$POP, 
        colour=factor(strata$REGION)),
        size = 2 )
#p <- p + theme_minimal() + theme(legend.position = "none")  + colScale
#p <- p + scale_color_igv() + theme_minimal() + theme(legend.position = "none")
p <- p +theme(legend.text = element_text(size = 10, face = "bold")) +
	theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
	colScale
pdf(file="pca_on_freq_population.pdf")
p
dev.off()

