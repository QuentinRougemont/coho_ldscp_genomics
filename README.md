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
Rscript ./00-scripts/01.hierfstats.R

     
## 4. Perform PCA and VAE analyses
     (to fill)
## 5. GEA Analyses:
  * RDA analyses:
    * use the Rscript to perform the GEA
    * Test significance with the following script
    * Plot de data
  *  LFMM analyses
    *  use the following script:
* To complete
## 6. Looking for parallelism
(To fill)
## 7. PBS analyses 
(To fill)

## 8. Association between recombination and outliers
(To fill)

## 9. looking for candidate
  * No GO enrichment instead I only use SNPeff and look for meaningfull outliers.
