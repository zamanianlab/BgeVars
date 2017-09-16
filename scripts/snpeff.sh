#!bin/bash

### Define project directories
proj="BgeVars"

gh_dir="${GIT_PATH}/${proj}"
local_dir="${GIT_DATA}/${proj}"


# Specify a genome 'build' (Species_Verion)
build="BglaB1.5"
# Fetch genome path
genome_path="`brew info snpeff | grep '/data' | cut -f 7 -d ' '`"

# Create directories 
mkdir -p ${genome_path}/${build}
mkdir -p ${genome_path}/genomes

# Update config file
echo "${build}.genome : B_glabrata" >> $genome_path/../snpEff.config

# Download genome / Extract sequence 
wget -nc -O ${genome_path}/genomes/${build}.fa.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02scaffoldsbglab1fagz
zcat ${genome_path}/genomes/${build}.fa.gz > ${genome_path}/genomes/${build}.fa

# Download/extract CDS
# wget -nc -O ${genome_path}/${build}/cds.fa.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02transcriptsbglab15fagz
# zcat ${genome_path}/${build}/cds.fa.gz > ${genome_path}/${build}/cds.fa

# Download/extract proteins
wget -nc -O ${genome_path}/${build}/protein.fa.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02peptidesbglab15fagz
zcat ${genome_path}/${build}/protein.fa.gz > ${genome_path}/${build}/protein.fa

# Download/extract gtf
wget -nc -O ${genome_path}/${build}/genes.gff.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02basefeaturesbglab15gtfgz
zcat ${genome_path}/${build}/genes.gff.gz > ${genome_path}/${build}/genes.gtf


# Build genome 
snpEff build -gtf22 -v ${build}

# Move VCF file to working directory
# cp ~/data/Bge/BGE.vcf ${genome_path}/${build}
# snpEff -v -csvStats -c ${genome_path}/../snpEff.config ${local_dir}/BGE.vcf ${build} > ${genome_path}/${build}/BGE.ann.vcf







