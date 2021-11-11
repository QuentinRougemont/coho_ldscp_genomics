#!/usr/bin/env Rscript

#script to run LFMM and plot the results.
#needs to be clean
#author: QR
###Verify if libs are installed:
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("magrittr" %in% rownames(installed.packages()) == FALSE)
{install.packages("magrittr", repos="https://cloud.r-project.org") }

## load libs
libs <- c('dplyr','data.table','magrittr', 'ade4','vegan','factoextra', 'lfmm')
invisible(lapply(libs, library, character.only = TRUE))

## read lfmm file:
argv <- commandArgs(T)

file <- argv[1]
snp  <- argv[2]
k <- as.numeric(argv[3])

input = paste0("zcat ", file)
dat <- fread(input)
loci <- read.table(snp) %>% 
    select(V1,V2,V3) #%>% 
    #set_colnames(.,c("CHR","POS","SNP")) 

##################  DOWNLOAD ENV DATA ########################  
pop_lat <- read.table("02-data/env/coho_dist_v2.txt",T)
pop_lat <- dplyr::select(pop_lat,POP_ID, elevation, dist_max_km, Latitude, Region)
#replace altitude of zero by 1
pop_lat$elevation[pop_lat$elevation == 0.00000 ] <- 1
enviro <-read.table("02-data/env/climat_epic4_wanted_pop.txt",T)
geol   <- read.table("02-data/env/era_rocktype_quanti_v2.txt",T)
enviro <- merge(enviro, geol, by="SITE")
pop    <- read.table("01-info/population.txt", h = T)
enviro <- merge(pop, enviro, sort=F)

####### PERFORM PCA ON ENVT VAR #################################################
X.temp <- dudi.pca(df = enviro[, 3:57], center = T, scale = T, scannf = FALSE, nf = 4)
X.prec  <- dudi.pca(df = enviro[, 58:97], center = T, scale = T, scannf = FALSE, nf = 4)
prec.ind <- get_pca_ind(X.prec)
temp.ind <- get_pca_ind(X.temp)
colnames(temp.ind$coord) <- c("Temperature1","Temperature2","Temperature3","Temperature4")
colnames(prec.ind$coord) <- c("Precipitation1","Precipitation2","Precipitation3")
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
env2$normalized_distance <-  env2$elevation * env2$dist_max_km
env2$normalized_distance <- range01(env2$normalized_distance)
#Keep only wanted variable for RDA:
env2 <- left_join(pop, env2) #, c("pop"="SITE")) 
env  <- dplyr::select(env2,  -elevation, -Region, -dist_max_km)
colnames(env)[2] <- "Geology"
#set them to an individual level now:
ind <- read.table("01-info/strata_gea.txt") %>% set_colnames(., c("SITE","IND"))
env <- filter(env, SITE %in% pop$SITE)
#ici il faut réordooner env1 pour qu'il matche l'ordre des ind?
X <- merge(ind, env, sort = F) # by.x="POP", by.y="SITE", sort = F)
X <-  dplyr::select(X, -SITE, -IND)                                                                                                                                                         

### data imputation for LFMM1:
#dat[dat == 9 ] <- NA   
Y <- apply(dat, 
    2, 
    function(x) {
    replace(x, is.na(x), as.numeric(names(which.max(table(x)))))
    })


#pc <- prcomp(Y)

#pdf(file="pca_variance.pdf")
#plot(pc$sdev[1:100]^2, xlab = 'PC', ylab = "Variance explained")
#dev.off()

#run lfmm
for (K in k){
mod.lfmm <- lfmm_ridge(Y = Y,
                        X = X,
                        K = K) #we could increase to 10 groups at most

    pv <- lfmm_test(Y = Y,
                X = X,
                lfmm = mod.lfmm,
                calibrate = "gif")
    pvalues <- pv$calibrated.pvalue

    #look at qqplot
    pdf(file = paste0("qplot_ridgeK",K,".pdf")) 
    qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
        pch = 19, cex = .4)
    abline(0,1)
    dev.off()

    #associat pvalue to loci
    pval_loc <- cbind(loci[,c(1:3)],pvalues)                              
    pvadj <- matrix(ncol = ncol(pvalues),nrow = nrow(pvalues), NA) 

    #correct value with BH method
    for(i in 1:ncol(pvalues)){
    pvadj[,i] <- p.adjust(pvalues[,i], method = 'BH')  
    }

    colnames(pvadj) <- colnames(pvalues)
    padj_loc <- cbind(loci[,c(1:3)],pvadj)                              

    #€xport value
    write.table(pval_loc,
    paste("pvalues_lfkmm_K",K,"txt", sep = "."),
    sep = "\t", quote = F, row.names = F)
    
    #corrected value
    write.table(padj_loc,
    paste("adjust_pvaluesBH_lfkmm_K",K,"txt", sep = "."),
    sep = "\t", quote = F, row.names = F)

    #Extract significiant variable using a stringeant p-value cutoff!
    geology <- filter(padj_loc, geology < 0.01) %>% select(V1,V2,V3,geology) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Temp1 <- filter(padj_loc,Temp1 < 0.01) %>% select(V1,V2,V3,Temp1) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH")) 
    Temp2 <- filter(padj_loc,Temp2 < 0.01) %>% select(V1,V2,V3,Temp2) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Temp3 <- filter(padj_loc,Temp3 < 0.01) %>% select(V1,V2,V3,Temp3) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Temp4 <- filter(padj_loc,Temp4 < 0.01) %>% select(V1,V2,V3,Temp4) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Prec1 <- filter(padj_loc,Prec1 < 0.01) %>% select(V1,V2,V3,Prec1) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Prec2 <- filter(padj_loc,Prec2 < 0.01) %>% select(V1,V2,V3,Prec2) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Prec3 <- filter(padj_loc,Prec3 < 0.01) %>% select(V1,V2,V3,Prec3) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Latitude <- filter(padj_loc,Latitude < 0.01) %>% select(V1,V2,V3,Latitude) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Dist <- filter(padj_loc,normalized_distance < 0.01) %>% 
        select(V1,V2,V3,normalized_distance) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH")) 
    
    geology$var = "geology"
    Temp1$var = "Temp1"
    Temp2$var = "Temp2"
    Temp3$var = "Temp3"
    Temp4$var = "Temp4"
    Prec1$var = "Prec1"
    Prec2$var = "Prec2"
    Prec3$var = "Prec3"
    Latitude$var = "Latitude"
    Dist$var = "Dist"
    
    all <- rbind(geology, Temp1, Temp2, Temp3, 
        Prec1, Prec2, Prec3, Dist)
    
    #remove variable that covary with latitude:
    all_corrected <- anti_join(all, Latitude, c("SNP" = "SNP"))   
    
    write.table(all_corrected,
    paste0("significant_outlier_control_forLatitudeK", K, ".txt"),
    quote = F, row.names = F, sep = "\t")

    write.table(all,
        paste0("significant_outlierK",K,".txt")
        ,quote = F,row.names = F, sep = "\t")

}
