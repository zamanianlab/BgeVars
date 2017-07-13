#!bin/bash

### Define project directories
proj="Local_BgeVars"

gh_dir="${GIT_PATH}/${proj}"
local_dir="${GIT_DATA}/${proj}"


# Specify a genome 'build' (Species_Verion)
build="BglaB1.5"
# Fetch genome path
genome_path="`brew info snpEff | grep '/data' | cut -f 7 -d ' '`"

# Create directories 
mkdir -p ${genome_path}/${build}
mkdir -p ${genome_path}/genomes

# Update config file
echo "${build}.genome : B_glabrata" >> $genome_path/../snpEff.config

# Download genome / Extract sequence 
wget -nc -O ${genome_path}/genomes/${build}.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_genomic.fna.gz
zcat ${genome_path}/genomes/${build}.fa.gz > ${genome_path}/genomes/${build}.fa

# Download/extract CDS
wget -nc -O ${genome_path}/${build}/cds.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_cds_from_genomic.fna.gz
zcat ${genome_path}/${build}/cds.fa.gz > ${genome_path}/${build}/cds.fa

# Download/extract proteins
wget -nc -O ${genome_path}/${build}/protein.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_protein.faa.gz
zcat ${genome_path}/${build}/protein.fa.gz > ${genome_path}/${build}/protein.fa

# Download/extract gff3
wget -nc -O ${genome_path}/${build}/genes.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_genomic.gff.gz
zcat ${genome_path}/${build}/genes.gff.gz > ${genome_path}/${build}/genes.gff


# Build genome 
snpEff build -gff3 -v ${build} -c ${genome_path}/../snpEff.config

# Move VCF file to working directory
cp ~/data/Bge/BGE.vcf ${genome_path}/${build}
# snpEff -v -csvStats -c ${genome_path}/../snpEff.config ${local_dir}/BGE.vcf ${build} > ${genome_path}/${build}/BGE.ann.vcf
