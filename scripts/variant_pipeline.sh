#!bin/bash

### Define project directories
proj="BgeVars"

gh_dir="${GIT_PATH}/${proj}"
local_dir="${GIT_DATA}/${proj}/Bge"

### variant calling
# -Ou pipe as uncompressed BCF (binary)
# -f reference
# -C50 downgrade mapping quality for reads with too many mistmatches (50 is recommended value for BWA alignments)
# -A output all alternate alleles even if they don't appear in the genotype (I think this important when we don't know ploidy)
# -m alt model for calling multiple alleles
# -v output variant sites only
bcftools mpileup -Ou -C50 --threads 4 -f "${local_dir}/BglaB1.5.fa" "${local_dir}/CA301ANXX.bam" | bcftools call -mv --threads 4 --ploidy 8 -A -Oz  > "${local_dir}/bcftools.raw.bcf"

### convert VCF to TSV
# vcf2tsv  "${local_dir}/BGE.VB.vcf" > "${local_dir}/BGE.VB.tsv"
### get AO and RO columns
# awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.tsv" > "${local_dir}/BGE.VB.AORO.tsv"

### extract SNPs and filter
# After filtering, kept 12968754 out of a possible 18957644 Sites
# vcftools --vcf "${local_dir}/BGE.VB.vcf" --remove-indels --recode --recode-INFO-all --out "${local_dir}/BGE.VB.snp.vcf"
# mv  "${local_dir}/BGE.VB.snp.vcf.recode.vcf"  "${local_dir}/BGE.VB.snp.vcf"
### convert to TSV
# vcf2tsv  "${local_dir}/BGE.VB.snp.vcf" > "${local_dir}/BGE.VB.snp.tsv"
### get AO and RO columns
# awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.snp.tsv" > "${local_dir}/BGE.VB.snp.AORO.tsv"

### continue filtering

### filter out variants that have all alternative observations coming from one pair of the reads
#  3394272
# vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s "${local_dir}/BGE.VB.snp.vcf" > "${local_dir}/BGE.VB.snp.fil1.vcf"
### convert to TSV
# vcf2tsv  "${local_dir}/BGE.VB.snp.fil1.vcf" > "${local_dir}/BGE.VB.snp.fil1.tsv"
### get AO and RO columns
# awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.snp.fil1.tsv" > "${local_dir}/BGE.VB.snp.fil1.AORO.tsv"

### filter out variants that have quality less than 1/4 of the depth
# 3392758
# vcffilter -f "QUAL / DP > 0.25" "${local_dir}/BGE.VB.snp.fil1.vcf" > "${local_dir}/BGE.VB.snp.fil2.vcf"
### convert to TSV
# vcf2tsv  "${local_dir}/BGE.VB.snp.fil2.vcf" > "${local_dir}/BGE.VB.snp.fil2.tsv"
### get AO and RO columns
# awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.snp.fil2.tsv" > "${local_dir}/BGE.VB.snp.fil2.AORO.tsv"

### recode to a new VCF
# After filtering, kept 3392758 out of a possible 3392758 Sites
# vcftools --vcf  "${local_dir}/BGE.VB.snp.fil2.vcf" --recode-INFO-all --out "${local_dir}/BGE.VB.snp.fil3" --recode 
### convert to TSV
# vcf2tsv  "${local_dir}BGE.VB.snp.fil3.recode.vcf" > "${local_dir}/BGE.VB.snp.fil3.tsv"
# get AO and RO columns
# awk '{print $1 "\t" $2 "\t" $6 "\t" $13 "\t" $36}' "${local_dir}/BGE.VB.snp.fil3.tsv" > "${local_dir}/BGE.VB.snp.fil3.AORO.tsv"

### can filter allele counts in R