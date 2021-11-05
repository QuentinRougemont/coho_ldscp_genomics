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

### 2.1. Filtering vcf


    (to fill)


### 2.2. Quality check

working on the same cleaned file (populations.mac1miss0.95.renamed.vcf)   


	* compute depth of sequencing, missing rate, genotyping rate  

    use:
		* vcftools site-mean-depth 
		* vcftools missing-indv  


	* verifiy the absence of a pattern in missing data


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

we see an absence that missing data are low and that these missing data do not show a particular pattern of structure.  


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
rm *tmp
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


### 4b. Running VAE:

Install vae available (here)[https://github.com/kr-colab/popvae]  
Due to its stochasticity, you'll obtained slightly different results than I did.
The analysis is really straighforward and I used default parameters here.


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



## 5.3.2. RDA  

Next we run the RDA in 2 steps

* 1. Significance testing 

First we need to test the significance of environmental data and of the RDA axes.
This is done through an ANOVA like permutation test of the RDA using the ```anova.cca``` function in R 
The ```ordiR2step``` function also allows to test the best model and choose significant variable of the RDA through a set of permutation test.
 
Run the script : 
```R
./00-scripts/gea/02.RDA_significance_testing.R
```

it takes a few hours on a cluster

* 2. identify outliers:

run the script : 

```R
Rscript 00-scripts/gea/03.RDA_identify_outliers.R
```

this will write the number of outlier in a file


## 5.3.3. Create a RDA Figure and LFMM plot:

To reproduce the Figure 3 simply run:

```R
Rscript 00-scripts/gea/Figure3.R
``` 

This should automatically produce the Figure3 from our manuscrit

 
## 6. Looking for parallelism
(To fill)


## 7. PBS analyses 
(To fill)


## 8. Association between recombination and outliers
(To fill)


## 9. looking for candidate
  * No GO enrichment instead I only use SNPeff and look for meaningfull outliers.
