#!bin/bash

# Initialize directories
boxdr=~/Box\ Sync
proj="Local_BgeVars"

gh_dir="${boxdr}/GitHub/${proj}"
local_dir="${boxdr}/GHdata/${proj}"

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
wget -O ${genome_path}/genomes/${build}.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_genomic.fna.gz
gzcat ${genome_path}/genomes/${build}.fa.gz > ${genome_path}/genomes/${build}.fa

# Download/extract CDS
wget -O ${genome_path}/${build}/cds.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_cds_from_genomic.fna.gz
gzcat ${genome_path}/${build}/cds.fa.gz > ${genome_path}/${build}/cds.fa

# Download/extract proteins
wget -O ${genome_path}/${build}/protein.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_protein.faa.gz
gzcat ${genome_path}/${build}/protein.fa.gz > ${genome_path}/${build}/protein.fa

# Download/extract gff3
wget -O ${genome_path}/${build}/genes.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_genomic.gff.gz
gzcat ${genome_path}/${build}/genes.gff.gz > ${genome_path}/${build}/genes.gff


# Build genome 
snpEff build -gff3 -v ${build} -c ${genome_path}/../snpEff.config

# Move VCF file to working directory
#cp ~/Box\ Sync/working/snail/genome/BGE.vcf ${genome_path}/${build}
snpEff -v -csvStats -c ${genome_path}/../snpEff.config ${local_dir}/BGE.vcf ${build} > ${genome_path}/${build}/BGE.ann.vcf

#measure length of each contig
#cat bgl_1.5.fa | awk '$0 ~ ">" {print c; c=0; printf $1 "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > key.txt

#sort according to length
#sort -rgk2 key.txt > sorted.key.txt