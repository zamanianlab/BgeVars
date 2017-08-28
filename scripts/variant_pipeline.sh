#!bin/bash

### Define project directories
proj="BgeVars"

gh_dir="${GIT_PATH}/${proj}"
local_dir="${GIT_DATA}/${proj}"

### variant calling
# -Ou pipe as uncompressed BCF (binary)
# -f reference
# -C50 downgrade mapping quality for reads with too many mistmatches (50 is recommended value for BWA alignments)
# -A output all alternate alleles even if they don't appear in the genotype (I think this important when we don't know ploidy)
# -m alt model for calling multiple alleles
# -v output variant sites only
# bcftools mpileup -Ou -C50 --threads 4 -f "${local_dir}/BglaB1.5.fa" "${local_dir}/CA301ANXX.bam" | bcftools call -mv --threads 4 -A -Oz  > "${local_dir}/bcftools.raw.bcf"
# bcftools view -Oz --threads 4 "${local_dir}/bcftools.raw.bcf" -o "${local_dir}/bcftools.vcf.gz"

### get relevant INFO tags
# bcftools query -f '%CHROM %POS %REF %ALT %DP4\n' "${local_dir}/bcftools.vcf" > "${local_dir}/bcftools.dp4.vcf"

### extract SNPs and filter
# After filtering, kept 10489715 out of a possible 12618001 Sites
# vcftools --vcf "${local_dir}/bcftools.vcf" --remove-indels --recode --recode-INFO-all --out "${local_dir}/bcftools.snp.vcf"
# mv  "${local_dir}/bcftools.snp.vcf.recode.vcf"  "${local_dir}/bcftools.snp.vcf"
### get relevant INFO tags
# bcftools query -f '%CHROM %POS %REF %ALT %DP4\n' "${local_dir}/bcftools.snp.vcf" > "${local_dir}/bcftools.snp.dp4.vcf"

### continue filtering
### only keep bi- or tri-allelic sites
# After filtering, kept 10478307 out of a possible 10489715 Sites
vcftools --vcf "${local_dir}/bcftools.snp.vcf"  --min-alleles 2 --max-alleles 3 --recode --recode-INFO-all --out "${local_dir}/bcftools.snp.fil1.vcf"
mv  "${local_dir}/bcftools.snp.fil1.vcf.recode.vcf"  "${local_dir}/bcftools.snp.fil1.vcf"
### get relevant INFO tags
bcftools query -f '%CHROM %POS %REF %ALT %DP4\n' "${local_dir}/bcftools.snp.fil1.vcf" > "${local_dir}/bcftools.snp.fil1.dp4.vcf"

### filter out variants that have quality less than 1/4 of the depth
# 10361386 out of 10478307
vcffilter -f "QUAL / DP > 0.25" "${local_dir}/bcftools.snp.fil1.vcf" > "${local_dir}/bcftools.snp.fil2.vcf"
### get relevant INFO tags
bcftools query -f '%CHROM %POS %REF %ALT %DP4\n' "${local_dir}/bcftools.snp.fil2.vcf" > "${local_dir}/bcftools.snp.fil2.dp4.vcf"

### recode to a new VCF
# After filtering, kept 10371439 out of a possible 10371439 Sites
vcftools --vcf  "${local_dir}/bcftools.snp.fil2.vcf" --recode-INFO-all --out "${local_dir}/bcftools.snp.fil3.vcf" --recode 
mv  "${local_dir}/bcftools.snp.fil3.vcf.recode.vcf"  "${local_dir}/bcftools.snp.fil3.vcf"
### get relevant INFO tags
bcftools query -f '%CHROM %POS %REF %ALT %DP4\n' "${local_dir}/bcftools.snp.fil3.vcf" > "${local_dir}/bcftools.snp.fil3.dp4.vcf"

### can filter allele counts in R