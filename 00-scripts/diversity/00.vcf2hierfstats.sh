#!/bin/bash

source ~/.bashrc

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
INPUT=${INPUT%.vcf}.lfmm

sed -e 's/2/22/g' $INPUT |\
    sed -e 's/9/NA/g' -e 's/1/12/g' -e 's/0/11/g' > out.tmp

paste pop.tmp  out.tmp > hierfstat.data.tmp
echo "loc" > loc.tmp

cut -d " " -f 3 $snp |perl -pe "s/\n/\t/g" > snp_id.tmp
paste loc.tmp snp_id.tmp > loc1.tmp
cat loc1.tmp hierfstat.data.tmp > hierfstat.data.txt

rm *tmp

