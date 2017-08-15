#!bin/bash

### Define project directories
proj="BgeVars"

gh_dir="${GIT_PATH}/${proj}"
local_dir="${GIT_DATA}/${proj}"

# extract bi-allelic SNPs and filter
# After filtering, kept 12968754 out of a possible 18957644 Sites
# vcftools --vcf "${local_dir}/BGE.VB.vcf" --remove-indels --recode --recode-INFO-all --out "${local_dir}/BGE.VB.snp.vcf"
# # convert to TSV
# vcf2tsv  "${local_dir}BGE.VB.snp.vcf" > "${local_dir}/BGE.VB.snp.tsv"
# # get AO and RO columns
awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.snp.tsv" > "${local_dir}/BGE.VB.snp.AORO.tsv"

### continue filtering

# filter out variants that have all alternative observations coming from one pair of the reads
#  3394272
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s "${local_dir}/BGE.VB.snp.vcf" > "${local_dir}/BGE.VB.snp.fil1.vcf"
# convert to TSV
vcf2tsv  "${local_dir}BGE.VB.snp.fil1.vcf" > "${local_dir}/BGE.VB.snp.fil1.tsv"
# get AO and RO columns
awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.snp.fil1.tsv" > "${local_dir}/BGE.VB.snp.fil1.AORO.tsv"

# filter out variants that have quality less than 1/4 of the depth
# 3392758
vcffilter -f "QUAL / DP > 0.25" "${local_dir}/BGE.VB.snp.fil1.vcf" > "${local_dir}/BGE.VB.snp.fil2.vcf"
# convert to TSV
vcf2tsv  "${local_dir}/BGE.VB.snp.fil2.vcf" > "${local_dir}/BGE.VB.snp.fil2.tsv"
# get AO and RO columns
awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.snp.fil2.tsv" > "${local_dir}/BGE.VB.snp.fil2.AORO.tsv"

# recode to a new VCF
# After filtering, kept 3392758 out of a possible 3392758 Sites
vcftools --vcf  "${local_dir}/BGE.VB.snp.fil2.vcf" --recode-INFO-all --out "${local_dir}/BGE.VB.snp.fil3" --recode 
# convert to TSV
vcf2tsv  "${local_dir}BGE.VB.snp.fil3.recode.vcf" > "${local_dir}/BGE.VB.snp.fil3.tsv"
# get AO and RO columns
awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.snp.fil3.tsv" > "${local_dir}/BGE.VB.snp.fil3.AORO.tsv"

### can filter allele counts in R