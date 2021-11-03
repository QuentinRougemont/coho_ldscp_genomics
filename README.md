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
    (to fill)

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
``

should procude all the results. 

The second script run on a cluster with >50Gb of RAM  :
```
Rscript ./00-scripts/diversity/01.hierfstats.R
```


Then we will perform plots of the correlation between the distance to the southernmost site and Bst and Hs statistics.  


## 4. Perform PCA and VAE analyses  
     (to fill)


## 5. GEA Analyses:

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
#Since the Thompson samples populations are highly divergent they may be driving part of the signal, we will therefore replicate the analyses by removing them from the samples.

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
#save space
gzip ${INPUT%.vcf}.lfmm
gzip ${INPUT%.vcf}.geno

plink --vcf $INPUT \
    --allow-extra-chr \
    --out ${INPUT%.recode.vcf} \
    --recode
    --recode

#prepare a cluster file
strata $INPUT strata.txt
cut -f1 strata.txt > col1
cut -d "_" -f2 strata.txt > col2
paste col1 col2 col1 > cluster.dat

#now compute allele frequency for RDA:
plink --file  ${INPUT%.recode.vcf}  \
    --allow-extra-chr --freq \
    --within cluster.dat 

gzip plink.frq.strat

rm *tmp col*
#now we can run the RDA and LFMM


```
 
 
## 6. Looking for parallelism
(To fill)


## 7. PBS analyses 
(To fill)

## 8. Association between recombination and outliers
(To fill)

## 9. looking for candidate
  * No GO enrichment instead I only use SNPeff and look for meaningfull outliers.
