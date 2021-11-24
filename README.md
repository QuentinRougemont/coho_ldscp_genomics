# landscape_genomics_coho
scripts to reproduce our landscape genomics results in coho salmon

See paper: Long distance migration is a major factor driving local adaptation at continental scale in a Pacific Salmon by Rougemont et al.

all necessary steps are highlighted below:

# Dependancies: 

Alignment and SNP calling: 
**bwa** software avialble [here](https://sourceforge.net/projects/bio-bwa/files/)

**fastp** software availalbe [here](https://github.com/OpenGene/fastp)  

**Samtools** software available [here](http://www.htslib.org/)  

**[htslib](http://www.htslib.org/download/)**  

**cutadapt** [cutadapt](https://github.com/marcelm/cutadapt)

**stacks** available [here](http://catchenlab.life.illinois.edu/stacks/) (see details of installation on the website)

for filtration:
**bcftools** software available [here](http://www.htslib.org/download/) 

**vcftools** available [here](http://vcftools.sourceforge.net/)

ANGSD, R & python

#I like to use the following fonction in my .bashrc :

```bash
function strata () { cat <(grep "CHR" "$1" ) |\
    cut -f 10- |\                                  
    perl -pe 's/\t/\n/g' | \
              sed 's/_/\t/g' | awk '{print $1"\t"$1"_"$2}' > "$2"; }
```

# DATA:

part of the RAW data are deposited on NCBI [PRJNA647050](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA647050). The remaining samples should be available soon  

**RAW vcf** and **filtered vcf** have been depostied on dryad (link to come soon)  

The **allele frequency file** used in the RDA will be available on dryad  

Other intermediary file (e.g. LFMM or hierfstat input) can be easily constructed with commands below   

metadata are available for GEA in the folder `02-data/env`  

additional metadata (individual name, river abbreviation, Latitude, Longitude, region of sampling, etc) can be fond in `01-info`   


# Steps :

## 1. Aligned Read and Call SNPs
  * Use [stacks_workflow](https://github.com/QuentinRougemont/stacks_v2_workflow) to perform the following:
    *   demultiplexing and trimming 
    *   align read,
    *   run stacks
    *   (filter vcf)
   
   note: I don't update stacks workflow any more and recommand to use the one of Eric Normandeau available [here](https://github.com/enormandeau/stacks_workflow) 
   
  * Use ANSGD to produce likelihood file (used to compute Fst and PBS latter)

## 2. Filter vcf and perform quality checks

To reproduce the filtering you can access the raw vcf from dryad

### 2.1 Filtering vcf


```bash

############################# Filtering RAW VCF ###########################

# in stacks populations was runn requiring:
# SNP Genotyped in : >70% of the samples
#                    >70% of the population


vcftools \
      --gzvcf populations.snps.vcf.gz
        --remove blacklisted_individuals.txt
        --recode-INFO-all
        --mac 1
        --max-meanDP 120
        --min-meanDP 10
        --out populations.DP10_120
        --recode


vcftools \
        --vcf populations.10_120.recode.vcf
        --recode-INFO-all
        --max-missing 0.85
        --out populations.10_120.miss85
        --recode  

# compute allele frequency and keep one SNP per loci:  

input="populations.10_120.miss85.recode.vcf"  

vcftools \
        --vcf $input\
        --freq --out out  

sed 1d out.frq | 's/:/\t/g' |cut -f 1,2,7,8 out.frq > frq2  

#then I use a horrible bash script to create a whitelist of independant SNPs. 
#A script in Eric Normandeau's [workflow](https://github.com/enormandeau/stacks_workflow/blob/master/00-scripts/utility_scripts/extract_snp_with_max_maf.py) should do the same more efficiently.

#compute missing rate per individuals  
vcftools --vcf $input --missing-indv --out individuals

#we remove individuals with too many missing data with awk
awk awk '$5>0.20 {print $1}'  individuals.imiss > blacklisted.individuals.txt

#remove individuals with too much missing data and keep one SNP per locus:  
vcftools \
     --vcf $input\
        --remove blacklisted.individuals.txt\
        --recode-INFO-all\
	--mac 1 \
        --out populations.mac1\
        --recode\
        --snps snp_id.txt\

#then perform a last quality filter:
vcftools \
    --vcf populations.mac1.recode.vcf \
    --max-missing 0.95 \
    --out populations.mac1miss0.95\
    --min-meanDP 10 --max-meanDP 120\
    --recode \
    --recode-INFO-all

#We don't want to filter on HWE since we may be interested in finding markers under selection 
#Moreover a global HWE filter is meaningless in highly structured populations

strata populations.mac1miss0.95.recode.vcf
bgzip populations.mac1miss0.95.recode.vcf; tabix -p vcf populations.mac1miss0.95.recode.vcf.gz

#then we need to rename part of the population HOP into HOD:
cut -f 1 strata.txt |sed '2576,2607s/HOP/HOD/g' > new.ind.tmp
cit -f 1 strata.txt > ind.tmp

paste ind.tmp new.ind.tmp > new_sample_name.txt
rm *tmp 

input=populations.mac1miss0.95.recode.vcf.gz

#then we rename individuals in the vcf 
bcftools reheader --samples new_sample_name.txt -o ${input%.recode.vcf.gz}.renamed.vcf.gz $input

#this is our final vcf ready for estimating genetic diversity and structure 
#see below for filtering for GEA

```

### 2.2. Quality check

working on the same cleaned file (populations.mac1miss0.95.renamed.vcf)   


	* compute depth of sequencing, missing rate, genotyping rate  

    use:
		* vcftools site-mean-depth 
		* vcftools missing-indv  


	* verifiy the absence of a pattern in missing data :


```bash
plink --file populations.mac1miss0.95.renamed.vcf \
    --allow-extra-chr \
    --mds-plot 4  --cluster missing  \
    --out ibm_mds_to_plot
```

then plot it using:  
```R
Rscript ./00-scripts/quality_check/plot_mds.R
```

we obtained the following plot:
![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/ibm_mds_point.git.png) 


we see that missing data are low and that these missing data do not show a particular pattern of structure.  


## 3. Compute genetic diversity and plot it

first convert vcf into hierfstat and other usefull input 

```bash
vcf=$1 #vcf here

#check compression
if file --mime-type "$vcf" | grep -q gzip$; then
   echo "$vcf is gzipped"
   gunzip "$vcf"
   INPUT=${vcf%}.gz
else
   echo "$vcf is not gzipped"
   INPUT=$vcf
fi

#extract population and individuals
strata $INPUT strata.txt ; cut -f 1 strata.txt > pop.tmp

vcf2geno $INPUT
geno2lfmm ${INPUT%.vcf}.geno
snp=${INPUT%.vcf}.vcfsnp
INPUT=${INPUT%.vcf}.geno

sed -e 's/2/22/g' $INPUT |\
    sed -e 's/9/NA/g' -e 's/1/12/g' -e 's/0/11/g' > out.tmp
paste pop out.tmp > hierfstat.data.tmp
echo "loc" > loc.tmp

cut -d " " -f 3 $snp |perl -pe "s/\n/\t/g" > snp_id.tmp
paste loc.tmp snp_id.tmp > loc.tmp
cat loc.tmp hierfstat.data.tmp > hierfstat.data.txt
rm \*tmp

```

then use R to compute basic statistics (Hs, Ho, Fis, Bst, in hierfstat) 

the scripts:  

```
00-scripts/diversity/00.vcf2hierfstats.sh
00-scripts/diversity/01.hierfstats.R
```

should procude all the results. 

The second script run on a cluster with >50Gb of RAM  :
```R
Rscript ./00-scripts/diversity/01.hierfstats.R
```

Then we will perform plots of the correlation between the distance to the southernmost site and Bst and Hs statistics. 

### Plotting diversity

 
Run:
```R
Rscript ./00-scripts/diversity/02.plot_ho_betaST_distance.R
```

on the command line or alternatively in Rstudio,  
this will test for correlation among diversity and distance and produce the following graph:
![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/figureS07.git.png) 



## 4. Perform PCA and VAE analyses  

### 4a. PCA

PCA was performed on individuals using [this script](https://github.com/QuentinRougemont/utility_scripts/blob/master/02.PCA/pca_on_vcf.R) and populations allele frequency using [this script](https://github.com/QuentinRougemont/utility_scripts/blob/master/02.PCA/pca_on_freq.R). 
I also made an IBS plot using plink with this [simple script](https://github.com/QuentinRougemont/utility_scripts/blob/master/07.random_scripts/plink_cluster_IBS.sh). 

To plot and fit to the Coho salmon data here's a simple modification of the script.
```R

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
pop <- data.frame(pop) %>% set_colnames(.,c("POP"))


#perform PCA
pca1 <- dudi.pca(freq4,scale=FALSE,scannf=FALSE)

#load strata
strata <- read.table("pop_ind_region_latitute_corrected3.txt") 
strata <- select(strata, V2,V6) %>% 
    set_colnames(.,c("POP","REGION"))
strata <- unique(strata)
strata <- strata[order(strata$POP),]  #reorder because population DRA has been renamed into WAC

#ensure that the order in strata will match order of the freq file:
strata = left_join(pop, strata)

#personalizing the colors:
myColors <- c("blue","orange","red","darkviolet","chocolate4", "springgreen4","green")
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
    colScale + 
    ggtitle("PCA on allele frequencies") + 
     theme(plot.title = element_text(lineheight=.8, face="bold" , hjust = 0.5))

    
pdf(file="pca_on_freq_population.pdf")
p
dev.off()


```

you should obtained an image as follows 
![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/pca_on_freq_population.git.png)  



in addition it is possible to test for population structure using a PCA at the individual level or an IBS plot using a similar command to the one use for IBM:


```bash
plink --file  populations.mac1miss0.95.renamed.vcf \
    --allow-extra-chr \
    --mds-plot 4  --cluster\
    --out mds_to_plot
```

Run this very simple script (similar to the one used for IBM):

```R
Rscript 00-scripts/pca_ibs/plot_mds.R
```

then you'll obtain the following graph:
![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/mds_point.git.png) 


This highlights particularly well the continuity of our data, where the axis of a PCA are highly correlated with latitude and longitude (see our paper and our previous one).  
Interestingly, we see that the Russian samples falls within Alaska and is not much divergent.  



### 4b. Running VAE:

Install vae available (here)[https://github.com/kr-colab/popvae]  
Due to its stochasticity, you'll obtained slightly different results than I did.
For instance depanding on the parameters used, you may obtained a distinction of nearly all river in separate cluster (Figure 2 of our paper)
In other case, you'll obtained a continuous distribution of the genotype entierely correlated with latitude. 
These two results reflect the existence of fine scale structure (recent postglacial homing) and of strong IBD following the postglacial expansion  



```bash

#install VAE and activate
#install tensorflow as well!

#create a working directory:
mkdir VAE
cd VAE
#copy the vcf 
ln -s ../02-data/populations.mac1miss0.95.renamed.vcf.gz .
#create folder to store the results

mkdir results
#run popvae with default parameters
popvae.py --infile populations.mac1miss0.95.renamed.vcf.gz --out results/coho --seed 12345

##more advanced use:
popvae.py --infile  populations.mac1miss0.95.renamed.vcf.gz\
    --out results2/coho \ 
    --seed 12345 \
    --search_network_sizes \ 
    --train_prop 0.90


```

#advice: try different run to see how results can differ

Then plot the results in R 

here's an example with many cluster:

![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/Fig2.git.png)

this is the graph I obtained at first run of the VAE and which is in our manuscript. Almost all sampling localities have been separated in discrete cluster). 
 
Populations further away in space are also further away on the graph.  

In the panel B/C I zoomed in the location of HaidaGwaii because there are few samples there. We can see that we are able to finely discriminate all samples and that those very closed in space are very closed in the Graph.   

The same results could be obtained by running many separate PCA. Here only one run was necessary.

At the same time therese strong IBD in our daya.


here's an example with pure geography (IBD): 

![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/vae_LD1_LD2.png)

accordingly we obtain a strong correlation with latitude:

![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/vae_LD1_Latitude.png)



## 5. GEA Analyses:


## 5.1. Extract environmental data  

The environmental data are located in the ```01-info``` folder with two text file:
* metadata 
* bioclimatic_var

Data included the following:
* Distance to the spawning site (km) 
* Elevation (m) 
* Watershed areas 
* Climate data (precipitation and temperature from WorldClim) 
* Geological data (extracted from the USGS database)  (included in ```02-data/env/geology```

Climate data include the mean, min, max range, sdt of 19 variables associated to precipitation and temperature.
These are obviously redundant and we summarized them and decorelate the variable using a PCA. 

The following simple code available at ```00-scripts/gea/00.pca_of_env.R```  

```R
## load libs
libs <- c('dplyr','ade4','factoextra','vegan','cowplot')
invisible(lapply(libs, library, character.only = TRUE))

###### DOWNLOAD ENV DATA
metadata <- read.table("01-info/metadata", h = T, sep = "\t")

bioclim <- read.table("01-info/bioclimatic_var", h = T, sep = "\t" )
#remove BNV and SAI:
bioclim <- bioclim %>%filter(POP !="SAI" & POP !="BNV")

####### PERFORM PCA ON ENVT VAR #################################################
#will do a PCA on the dataaset that contains only climatic variables:
#separate temperature and precipitation and keep only significant axis or only first axis
X.temp <- dudi.pca(df = bioclim[, 2:56], center = T, scale = T, scannf = FALSE) #nf = 3)
X.prec  <- dudi.pca(df = bioclim[, 57:96], center = T, scale = T, scannf = FALSE) # nf = 3)

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
colnames(temp.ind$coord) <- c("Temp1","Temp2","Temp3","Temp4")
colnames(prec.ind$coord) <- c("Prec1","Prec2","Prec3")

prec = cbind(bioclim$POP, prec.ind$coord)
temp = cbind(bioclim$POP, temp.ind$coord)

write.table(temp, "02.data/env/temperature_coordinates_pca.txt" , quote = F, row.names = F )
write.table(prec, "02.data/env/precipitation_coordinates_pca.txt" , quote = F, row.names = F )

```

this will create the two input-files   
* temperature_coordinates_pca.txt  
* precipitation_coordinates_pca.txt  

 located in 02.data/env/ 

This will be used in the GEA  

we can check that the data are not correlated using a simple correlation plot in R for instance.
These few lines of code will do the job:

```R

library(magrittr)
library(dplyr)
library(corrplot)

###### DOWNLOAD ENV DATA
metadata <- read.table("01-info/metadata", h = T, sep = "\t")
metadata<- metadata %>%filter(POP !="SAI" & POP !="BNV")

geol <- read.table("02-data/env/geology", T)
#remove BNV and SAI:
geol <- geol %>%filter(POP !="SAI" & POP !="BNV")
#remove BNV and SAI:

#replace altitude of zero by 1
metadata$elevation[metadata$elevation == 0.00000 ] <- 1

#now takes into accont the standardized elevation * distance interaction
#standardiztation function:
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
metadata$normalized_distance =  metadata$elevation * metadata$dist_max_km
metadata$normalized_distance <- range01(metadata$normalized_distance)

####### PERFORM PCA ON ENVT VAR #################################################
temp <- read.table("02-data/env/temperature_coordinates_pca.txt", T)
prec <- read.table("02-data/env/precipitation_coordinates_pca.txt", T)

env1 <- dplyr::select(metadata, POP, Latitude, normalized_distance, Region)
env1 <- merge(env1, temp)
env1 <- merge(env1, prec)
env1 <- merge(env1, geol)

#remplacer avec sed directement dans metadata, geology et bioclim
env1 <- env1[order(env1$POP),]

#choose the variable we want to work with:
env <- select(env1, -POP, -Region)

#rename for readiblity of the plot:
colnames(env) <- c("Latitude","normalized_distance", 
    "TemperaturePC1", "TemperaturePC2", "TemperaturePC3", "TemperaturePC4",
    "PrecipitationPC1", "PrecipitationPC2", "PrecipitationPC3", 
    "geology")

#again write the correlation 
mat <- round(cor(env),2)

pdf(file="corrplot.pdf")
corrplot(mat, method = "number", type = "lower")
dev.off()

```

and produce the following graph:
![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/corrplot.png) 




## 5.2. Prepare genomics data  

 We first need to filter the data to exclude 2 populations, namely BNV and SAI for which environmental data are not available. 

```bash
 #################################################################################
#                           preparing data for RDA and GEA
#################################################################################

mkdir GEA
cd GEA

#we create a file containing population SAI and BNV that we want to exclude
grep "SAI\|BNV" ../strata > SAI_BNV.txt 

vcftools --gzvcf ../populations.24_09.missing0.92mac15.indep.clean2.recode.vcf.gz \
    --remove SAI_BNV.txt \
    --mac 15  \
    --max-missing 0.95 \
    --min-meanDP 10 --max-meanDP 120 \
    --out population_for_GEA \
    --recode \
    --recode-INFO-all

input=population_for_GEA.recode.vcf
strata $input strata.txt
bgzip $input

input=$input.gz

cut -f 2 strata.txt > ind.tmp; cp ind.tmp new.ind.tmp 

#in this file we need to rename the 32 first individuals from HOP into HOD as the "HOP" population contains individual sample from both "HOP" and "HOD"
sed -i '2536,2567s/HOP/HOD/g' new.ind.tmp
sed -i 's/DRA/WAC/g' new.ind.tmp 

paste ind.tmp new.ind.tmp > new_sample_name.txt 

#then we rename individuals in the vcf 
bcftools reheader --samples new_sample_name.txt -o ${input%.recode.vcf.gz}.renamed.vcf.gz $input 

gunzip ${input%.recode.vcf.gz}.renamed.vcf.gz

INPUT="${input%.recode.vcf.gz}.renamed.vcf"

#create usefull intput for LFMM:
vcf2geno $INPUT 
geno2lfmm ${INPUT%.vcf}.geno 
sed -i 's/9/NA/g' ${INPUT%.vcf}.lfmm

#save space
gzip ${INPUT%.vcf}.lfmm
gzip ${INPUT%.vcf}.geno


plink --vcf $INPUT \
    --allow-extra-chr \
    --out ${INPUT%.recode.vcf} \
    --recode
    --recode

#prepare a cluster file
rm strata.txt
strata $INPUT strata.txt
cut -d "_" -f2 new.ind.tmp > col2
cut -d "_" -f1 new.ind.tmp > col1
paste col1 col2 col1 > cluster.dat

#now compute allele frequency for RDA:
plink --file  ${INPUT%.recode.vcf}  \
    --allow-extra-chr --freq \
    --within cluster.dat 

gzip plink.frq.strat

rm *tmp col*
#now we can run the RDA and LFMM


################################################################################
#           prepare the data by excluding the Thompson samples                 #
################################################################################
# Since the Thompson samples populations are highly divergent they may be driving part of the signal, we will therefore replicate the analyses by removing them from the samples.

mkdir no_thompson
cd no_thompson

awk -F"\t" '$10=="Thompson" {print $0}' ../metadata > thompson
awk -F"\t" '$10=="Thompson" {print $1}' ../metadata > thompson.pop #to keep only thompson rivers

grep -Ff thompson.pop ../newnames  > individus.thompson #to extract individuals

input=population_for_GEA.renamed.vcf.gz
output=population_for_GEA.renamed.nothompson

vcftools --gzvcf ../$input \
    --remove individus.thompson \
    --mac 15  \
    --max-missing 0.95 \
    --min-meanDP 10 --max-meanDP 120 \
    --out $output \
    --recode \
    --recode-INFO-all


strata $output.recode.vcf strata.txt

INPUT="${output}.recode.vcf"

vcf2geno $INPUT 
geno2lfmm ${INPUT%.vcf}.geno 
sed -i 's/9/NA/g' ${INPUT%.vcf}.lfmm

# save space

gzip ${INPUT%.vcf}.lfmm
gzip ${INPUT%.vcf}.geno

plink --vcf $INPUT \
    --allow-extra-chr \
    --out ${INPUT%.recode.vcf} \
    --recode
    --recode

# prepare a cluster file

strata $INPUT strata.txt
cut -f1 strata.txt > col1
cut -d "\_" -f2 strata.txt > col2
paste col1 col2 col1 > cluster.dat

# now compute allele frequency for RDA

plink --file  ${INPUT%.recode.vcf} \
    --allow-extra-chr --freq \
    --within cluster.dat 

gzip plink.frq.strat

rm \*tmp col\* 
#now we can run the RDA and LFMM


```

## 5.3. Running the GEAs


## 5.3.1 LFMM

LFMM will use a mixed model analyses to identify outliers

simply use the script in 

```R
./00-scripts/gea/01.lfmm_step1.R 
```

this will run rather rapidly with approximately 20-30Gb of RAM  



## 5.3.2 RDA


Next we run the RDA in 2 steps

* 1. Significance testing 

First we need to test the significance of environmental data and of the RDA axes.
This is done through an ANOVA like permutation test of the RDA using the ```anova.cca``` function in R 
The ```ordiR2step``` function also allows to test the best model and choose significant variable of the RDA through a set of permutation test.
 
Run the script : 
```R
./00-scripts/gea/02.RDA_significance_testing.R
```

it takes a few hours (~one night) on a cluster

* 2. identify outliers:

run the script : 

```R
Rscript 00-scripts/gea/03.RDA_identify_outliers.R
```

this will write the number of outlier in a file as well as other useful data such as the number of outliers in each chromosomes.  

We used the approach of Forester et al. (2018) to identify outliers.  



## 5.3.3. Create a RDA Figure and LFMM plot:

To reproduce the Figure 3 simply run:

```R
Rscript 00-scripts/gea/Figure3.R
``` 

This should automatically produce the Figure3 from our manuscript


 
## 6. Looking for parallelism

Run the Rscript 
```R
00-scripts/gea/05.Figure4.R 
```

This should automatically produce the Figure4 from our manuscript as well as some statisticall tests
See explanation in our methods of the manuscript  
I'll document all of this later.  


## 7. plotting LFMM results:

this can be done with CMplot 
see script 00-scripts/gea/LFMM_plot_FigureS10.R or simple code below:

```R
## load libs
libs <- c('dplyr','data.table', 'magrittr', 'cowplot','CMplot', 'gridGraphics')
invisible(lapply(libs, library, character.only = TRUE))

pval_loc <- fread("zcat adjust_pvaluesBH_lfkmm_K.20.txt.gz")

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
pdf(file="FigureS10.pdf",15,20)
plot_grid(p1, p2,p3, 
    ncol = 1, 
    align = "v", 
    labels=c("A - Temperature", "B - Precipitation", "C - Geology") )
dev.off()

```

## 8. PBS analyses 

*  **ANGSD analyses:**

The major reason for using ANGSD was to compute the PBS score from Yi et al. 2010 (although it can be computed from the vcf directly)

All steps to run ANGSD are found [here](https://github.com/QuentinRougemont/utility_scripts/tree/master/00.ANGSD). They necessitate to start from the bam files.  

* **important note**

before running Fst or PBS I used a three-species outgroup sequence for proper folding made of :
 * 1. [Atlantic salmon](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA287458) n = 5 individuals 
 * 2. [Rainbow trout](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386519) n = 5 individuals and 
 * 3. [Sockeye salmon](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA530256) n = 3 individuals

Whole genome sequences were:
 * processed with [fastp](https://github.com/QuentinRougemont/gatk_haplotype/blob/master/gatk4/01_scripts/01_fastp.sh) 
 * mapped to the coho salmon genome using [bwa-mem](https://github.com/QuentinRougemont/gatk_haplotype/blob/master/gatk4/01_scripts/02_bwa_mem2_align_reads_PE.sh). 
 * marked for duplicates with [picard](https://github.com/QuentinRougemont/gatk_haplotype/blob/master/gatk4/01_scripts/03_remove_duplicates.sh)  

The resulting bam were processed to reconstruct an ancestral sequence with [ancestral_seq.sh](https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/ancestral_seq.sh) 

Then Simply follow the script in the order of their number:  

 * 1. [00_angsd_bam_to_saf.sh](https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/00_angsd_bam_to_saf.sh) 
 * 2. [01_1dsfs.sh](https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/01_1dsfs.sh)
 * 3. [02_2dsfs.sh](https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/02_2dsfs.sh)
 * 4. [04.pairwise_fst_sliding_window.sh](https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/04.pairwise_fst_sliding_window.sh) #to get only Fst
 * 5. [05_PBS_sliding_windows.sh](https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/05_PBS_sliding_windows.sh) #to get Fst and PBS



* Produce graph: I used [CMplot](https://github.com/YinLiLin/CMplot) to make all manhattan plot. 

The used is straightforward:

```R
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
ALB\_LCA2 <- dplyr::select(ALB\_LCA, chr, midPos, Fst01)
POR\_SHK2 <- dplyr::select(POR\_SHK, chr, midPos, Fst01)
SFK\_GOO2 <- dplyr::select(SFK\_GOO, chr, midPos, Fst01)
SUS\_OON2 <- dplyr::select(SUS\_OON, chr, midPos, Fst01)

ALB\_POR <- full\_join(ALB\_LCA2, POR\_SHK2, by=c("chr","midPos"))
ALB\_POR\_SFK <- full\_join(ALB\_POR, SFK\_GOO2, by=c("chr","midPos"))
ALB\_POR\_SFK\_SUS <- full\_join(ALB\_POR\_SFK, SUS\_OON2, by=c("chr","midPos"))

distplot <- ALB\_POR\_SFK\_SUS
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
```



## 9. Association between recombination and outliers

1. **estimate recombination**

Recombination was estimated using 30X whole Genome data from this [preprint](https://www.biorxiv.org/content/10.1101/732750v2).  
(note that the WGS data were removed from the published article and will be published later)   
LDhat was used with this [pipeline](https://github.com/QuentinRougemont/LDhat_workflow). To reproduce these estimate simply use the pipeline.  

Rho was estimated for each population separately  

Then results were summarized into 250 kb windows and each populations were concatenated together  

The scripts [here](https://github.com/QuentinRougemont/LDhat_workflow/blob/master/02-scripts/06.run_slidingwindows.sh) should enables to summarize the data with a bit of modifications depending on your data


2. **identify outliers**  

outliers were identified in the steps above for the RDA

3. **test for how outlier (GEA and Fst/PBS) are influenced by recombination:**

Look at the script 

```
00-scripts/recombination/association_recombination_GEA.R
```

This will allow to perform statistical test: 

	* wilcoxon test of difference in mean recombination rate between outliers vs neutral region
	* mixed models test 

A plot like the one in FigS15 will be produced. 
![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/FigS15.git.png)

Here we only have the population scale recombination rate (rho = 4*Ne*r)  

In addition it is possible to explore this by considering either:  
	* shared outliers between LFMM & RDA  
	* RDA only outliers  
	* LFMM only outliers  
	* outliers that fall into areas of residual tetraploidy can be investigated separately 
	since there's less variance in recombination in these region  

	There was no difference in recombination when considering the few outliers on region of residual tetraploidy.  
        Obviously, removing these region increased the strength of the signal   


## 10. looking for candidate

  * No GO enrichment instead I only use SNPeff and look for meaningfull outliers.

  SNPeff can be found [here](https://pcingola.github.io/SnpEff/) and was run with standard parameters



## 11. Other stuff

 * below is a zoomable map for BC/Thompson/HaidaGwaii that include all sample sites used here. 
 * may be usefull for users  
 * This map was obtained with very simple code in R

![example_graph](https://github.com/QuentinRougemont/coho_ldscp_genomics/blob/main/pictures/sampling_map.git.png) 

