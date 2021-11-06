#script to run LFMM and plot the results.
#needs to be clean
#author: QR
###Verify if libs are installed:
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("magrittr" %in% rownames(installed.packages()) == FALSE)
{install.packages("magrittr", repos="https://cloud.r-project.org") }
if("reshape" %in% rownames(installed.packages()) == FALSE)
{install.packages("reshape", repos="https://cloud.r-project.org") }
if("data.table" %in% rownames(installed.packages()) == FALSE)
{install.packages("data.table", repos="https://cloud.r-project.org") }
##for LFMM1:
if("lfmm" %in% rownames(installed.packages()) == FALSE)
{install.packages("lfmm", repos="https://cloud.r-project.org") }
##for LFMM2:
if("LEA" %in% rownames(installed.packages()) == FALSE)
{BiocManager::install("LEA") }

## load libs
libs <- c('dplyr','reshape2','data.table', 'magrittr',  'lfmm')
invisible(lapply(libs, library, character.only = TRUE))

## read lfmm file:
argv <- commandArgs(T)

file <- argv[1]
snp  <- argv[2]

input = paste0("zcat ", file)
dat <- fread(input)
loci <- read.table(snp) %>% 
    select(V1,V2,V3) #%>% 
    #set_colnames(.,c("CHR","POS","SNP")) 

##################  DOWNLOAD ENV DATA ########################  
metadata <- read.table("01-info/metadata", h = T, sep = "\t")
metadata <- metadata %>% filter(POP !="SAI" & POP !="BNV")

geol <- read.table("02-data/env/geology", T)
#remove BNV and SAI:
geol <- geol %>% filter(POP !="SAI" & POP !="BNV")
#remove BNV and SAI:

#replace altitude of zero by 1
metadata$elevation[metadata$elevation == 0.00000 ] <- 1

#now takes into accont the standardized elevation * distance interaction
#standardiztation function:
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
metadata$normalized_distance =  metadata$elevation * metadata$dist_max_km
metadata$normalized_distance <- range01(metadata$normalized_distance)

####### READ PCA axis #################################################
temp <- read.table("02-data/env/temperature_coordinates_pca.txt", T)
prec <- read.table("02-data/env/precipitation_coordinates_pca.txt", T)

#merge all environmental variable:
env1 <- dplyr::select(metadata, POP, Latitude, normalized_distance)
env1 <- merge(env1, temp)
env1 <- merge(env1, prec)
env1 <- merge(env1, geol)

#set them to an individual level now:
ind <- read.table("01-info/strata.txt") %>% set_colnames(., c("POP","IND"))

#ici il faut réordooner env1 pour qu'il matche l'ordre des ind?

X <- left_join(ind,env1)
X <- merge(ind,env1, sort = F)

#WAC: 1593-1638
#check the position of WAC in env data:
#X[1593:1638,1:5]

X <-  dplyr::select(X, -POP, -IND)                                                                                                                                                         
#X <-  dplyr::select(X, -POP, -IND, -Temp4)  

### data imputation for LFMM1:
#dat[dat == 9 ] <- NA   
Y <- apply(dat, 
    2, 
    function(x) {
    replace(x, is.na(x), as.numeric(names(which.max(table(x)))))
    })


pc <- prcomp(Y)

pdf(file="pca_variance.pdf")
plot(pc$sdev[1:100]^2, xlab = 'PC', ylab = "Variance explained")
dev.off()

#run lfmm
for (K in c(5, 6, 12, 20, 40, 60, 100, 200)){
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
    geology <- filter(padj_loc, geology <= 0.02) %>% select(V1,V2,V3,geology) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Temp1 <- filter(padj_loc,Temp1 <= 0.02) %>% select(V1,V2,V3,Temp1) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH")) 
    Temp2 <- filter(padj_loc,Temp2 <= 0.02) %>% select(V1,V2,V3,Temp2) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Temp3 <- filter(padj_loc,Temp3 <= 0.02) %>% select(V1,V2,V3,Temp3) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Temp4 <- filter(padj_loc,Temp4 <= 0.02) %>% select(V1,V2,V3,Temp4) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Prec1 <- filter(padj_loc,Prec1 <= 0.02) %>% select(V1,V2,V3,Prec1) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Prec2 <- filter(padj_loc,Prec2 <= 0.02) %>% select(V1,V2,V3,Prec2) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Prec3 <- filter(padj_loc,Prec3 <= 0.02) %>% select(V1,V2,V3,Prec3) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Latitude <- filter(padj_loc,Latitude <= 0.02) %>% select(V1,V2,V3,Latitude) %>% 
        set_colnames(.,c("CHR","POS","SNP","BH"))
    Dist <- filter(padj_loc,normalized_distance <= 0.02) %>% 
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
    
    all <- rbind(geology, Temp1, Temp2, Temp3, Temp4,
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
