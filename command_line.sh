#!/bin/bash


#################################################################################
#                           Filter the whole RAW vcf                            #
#################################################################################






#################################################################################
#                           preparing data for RDA and GEA
#################################################################################

vcftools --gzvcf ../populations.24_09.missing0.92mac15.indep.clean2.recode.vcf.gz \
    --remove ../Russia_BNV \
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
sed -i '2536,2567s/HOP/HOD/g' new.ind.tmp
paste ind.tmp new.ind.tmp > new_sample_name.txt 

bcftools reheader --samples new_sample_name.txt -o ${input%.recode.vcf.gz}.renamed.vcf.gz $input 

gunzip ${input%.recode.vcf.gz}.renamed.vcf.gz

INPUT="${input%.recode.vcf.gz}.renamed.vcf"

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
#

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


