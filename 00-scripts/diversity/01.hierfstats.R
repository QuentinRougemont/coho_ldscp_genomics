
#Purpose:script to compute basic genetic diversity statistics
#date: October 2021
#author: QR

if("hierfstat" %in% rownames(installed.packages()) == FALSE)
{install.packages("hierfstat", repos="https://cloud.r-project.org") }


library(hierfstat)

#read data
data <-read.table("hierfstat.data.txt", header=T ) #, sep="\t")

stats <- basic.stats(data, diploid=TRUE)
write.table(stats$perloc,"hierfstat.diversity.index.",quote=F,row.names=T,col.names=T)
write.table(stats$Ho, "hierfstat.ho",quote=F,row.names=T,col.names = T)
write.table(stats$Hs, "hierfstat.hs",quote=F,row.names=T,col.names = T)
write.table(stats$Fis,"hierfstats.fis",quote=F, row.names=T, col.names=T)
write.table(stats$overall,"hierfstats.overall",quote=F,row.names=T,col.names=T)


beta1 <- betas(data, nboot=100)
write.table(t(rbind(beta1$betaiovl,beta1$ci)),"beta_iovl_ci", col.names=T, row.names=T) 
write.table(t(rbind(beta1$betaiovl,beta1$ci)),"beta_iovl_ci2", col.names=c("beta","2.5%","97.5%"), row.names=T)

fis <- boot.ppfis(data, nboot = 1000) 

write.table(b$fis.ci, "hierfstat.confidence.interval.fis.txt",
    quote = F, row.names = T, col.names = T)


ho <- read.table("hierfstat.ho")
hs <- read.table("hierfstat.hs")

ho2 <- data.frame(colMeans(ho[,-1], na.rm=T))
colnames(ho2) <- "ho"

hs2 <- data.frame(colMeans(hs[,-1], na.rm=T))
colnames(hs2) <- "hs"

write.table(ho2,"mean_Ho_hierfstat.txt",quote=F)
write.table(hs2,"mean_Hs_hierfstat.txt",quote=F)

