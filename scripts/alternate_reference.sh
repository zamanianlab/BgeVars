#!bin/bash

### Define project directories
boxdr=~/Box\ Sync
proj="Local_BgeVars"

gh_dir="${boxdr}/GitHub/${proj}"
local_dir="${boxdr}/GHdata/${proj}"

# Fetch and extract Bgla reference genome
wget -nc -O "${local_dir}/GCF_000457365.1_ASM45736v1_genomic.fna.gz" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_genomic.fna.gz
gunzip -dv "${local_dir}/GCF_000457365.1_ASM45736v1_genomic.fna.gz" > "${local_dir}/GCF_000457365.1_ASM45736v1_genomic.fna"

reference="${local_dir}/GCF_000457365.1_ASM45736v1_genomic.fna"
vcf="${local_dir}/BGE.vcf"

# Create genome dict and index viles
picard CreateSequenceDictionary R="${reference}" O="${reference}.dict"
rename 's/.fna//' "${reference}.dict"
samtools faidx "${reference}"
gatk -T FastaAlternateReferenceMaker -R "${reference}" -o "${local_dir}"/Bge_GCR_00457365.fa -V "${vcf}"


