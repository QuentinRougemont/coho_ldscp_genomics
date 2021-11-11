#!/usr/bin/env/ Rscript

#Purpose: RDA script to test significance of the axis and environemental variable
#Author : QR
#Date: October 2021

## Verify if libs are installed:
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("ade4" %in% rownames(installed.packages()) == FALSE)
{install.packages("ade4", repos="https://cloud.r-project.org") }
if("factoextra" %in% rownames(installed.packages()) == FALSE)
{install.packages("factoextra", repos="https://cloud.r-project.org") }
if("vegan" %in% rownames(installed.packages()) == FALSE)
{install.packages("vegan", repos="https://cloud.r-project.org") }
if("cowplot" %in% rownames(installed.packages()) == FALSE)
{install.packages("cowplot", repos="https://cloud.r-project.org") }


## load libs
libs <- c('dplyr','ade4','factoextra','vegan','cowplot')
invisible(lapply(libs, library, character.only = TRUE))

###### DOWNLOAD ENV DATA
pop <- read.table("01-info/population.txt", h = T)

bioclim <- read.table("02-data/env/climat_epic4_wanted_pop.txt", h = T, sep = "\t")

#preserve only the wanted population for the analyses using a joint:
bioclim <- merge(pop, bioclim, by = "SITE", sort = F)

####### PERFORM PCA ON ENVT VAR #################################################
#will do a PCA on the dataaset that contains only climatic variables:
#separate temperature and precipitation and keep only significant axis or only first axis
X.temp <- dudi.pca(df = bioclim[, 3:57], center = T, scale = T, scannf = FALSE) #nf = 3)
X.prec  <- dudi.pca(df = bioclim[, 58:97], center = T, scale = T, scannf = FALSE) # nf = 3)

## singificance of axis
eig.val <- get_eigenvalue(X.temp) #first we get the eigen value in an easy form
eig.val$eig <- eig.val$variance.percent/100 #percent
expected <- bstick(length(eig.val$eig) )
signif <- eig.val$eig > expected #get signifcicant axis #3axes sont signif
signif
#4 axis are significant

eig.val.p <- get_eigenvalue(X.prec) #first we get the eigen value in an easy form
eig.val.p$eig <- eig.val.p$variance.percent/100 #percent
expected <- bstick(length(eig.val.p$eig) )
signif.p <- eig.val.p$eig > expected #get signifcicant axis #3axes sont signif
signif.p
#3 axis are signficant

#visualisation of individuals cos^2:
X.temp <- dudi.pca(df = bioclim[, 2:53], center = T, scale = T, scannf = FALSE, nf = 4)
X.prec  <- dudi.pca(df = bioclim[, 54:95], center = T, scale = T, scannf = FALSE,  nf = 3)

fviz_pca_ind(X.prec, col.ind="cos2", geom = "point") +
    scale_color_gradient2(low="white", mid="blue",
    high="red", midpoint=0.6)+ theme_minimal()

p12 <- fviz_pca_var(X.prec, col.var="steelblue")+
    theme_minimal()  +  ggtitle("Precipitation PC axis 12\nvariables") + 
     theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))

p23 <- fviz_pca_var(X.prec, axes = c(2,3) , col.var="steelblue")+
    theme_minimal() +  ggtitle("Precipitation PC axis 23\nvariables") + 
     theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))

t12 <- fviz_pca_var(X.temp, col.var="steelblue")+
    theme_minimal() +  ggtitle("Temperature PC axis 12\nvariables") + 
     theme(plot.title = element_text(lineheight=.8, face="bold" , hjust = 0.5))

t34 <- fviz_pca_var(X.temp, axes = c(3,4) , col.var="steelblue")+
    theme_minimal() +  ggtitle("Temperature PC axis 34\nvariables") + 
     theme(plot.title = element_text(lineheight=.8, face="bold" , hjust = 0.5))


t_eig  <- fviz_eig(X.temp, addlabels = TRUE, ylim = c(0, 50)) + 
    ggtitle("Temperature eigenvalues") + 
     theme(plot.title = element_text(lineheight=.8, face="bold" , hjust = 0.5))


p_eig  <- fviz_eig(X.prec, addlabels = TRUE, ylim = c(0, 50)) + 
    ggtitle("Precipitation eigenvalues") + 
     theme(plot.title = element_text(lineheight=.8, face="bold" , hjust = 0.5))



pdf(file = "temperature_and_precipiration.pca.pdf", 18,12)
plot_grid(t12, t34, t_eig, p12, p23,p_eig, 
    labels = "auto", ncol = 3)
dev.off()

#We will keep all three axis (significant, but could also work on the first only:
prec.ind <- get_pca_ind(X.prec)
temp.ind <- get_pca_ind(X.temp)
colnames(temp.ind$coord) <- c("TemperaturePC1","TemperaturePC2","TemperaturePC3","TemperaturePC4")
colnames(prec.ind$coord) <- c("PrecipitationPC1","PrecipitationPC2","PrecipitationPC3")

prec = cbind(bioclim$POP, prec.ind$coord)
temp = cbind(bioclim$POP, temp.ind$coord)


write.table(temp, "02.data/env/temperature_coordinates_pca.txt" , quote = F, row.names = F )
write.table(prec, "02.data/env/precipitation_coordinates_pca.txt" , quote = F, row.names = F )

