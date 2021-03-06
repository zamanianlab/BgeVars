#!bin/bash

### Define project directories
proj="BgeVars"

gh_dir="${GIT_PATH}/${proj}"
local_dir="${GIT_DATA}/${proj}"


# Specify a genome 'build' (Species_Verion)
build="BglaB1.6"
# Fetch genome path
# genome_path="/home/linuxbrew/.linuxbrew/Cellar/snpeff/4.3i/share/snpeff/data"
genome_path="/Users/nic/bin/snpEff/data"

# Create directories 
mkdir -p ${genome_path}/${build}
mkdir -p ${genome_path}/genomes

# Update config file
echo "${build}.genome : B_glabrata" >> $genome_path/../snpEff.config

# Download genome / Extract sequence 
wget -nc -O ${genome_path}/genomes/${build}.fa.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02scaffoldsbglab1fagz
zcat ${genome_path}/genomes/${build}.fa.gz > ${genome_path}/genomes/${build}.fa
cat ${genome_path}/genomes/${build}.fa | sed "s/>\(.*\) dna:.*/>\1/" > ${genome_path}/genomes/${build}2.fa
mv ${genome_path}/genomes/${build}2.fa ${genome_path}/genomes/${build}.fa

# Download/extract CDS
wget -nc -O ${genome_path}/${build}/cds.fa.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02transcriptsbglab16fagz
zcat ${genome_path}/${build}/cds.fa.gz > ${genome_path}/${build}/cds.fa
cat ${genome_path}/${build}/cds.fa | perl -pe "s/>(.*?) .*/>\1/" ${genome_path}/${build}/cds.fa > ${genome_path}/${build}/cds2.fa
mv ${genome_path}/${build}/cds2.fa ${genome_path}/${build}/cds.fa

# Download/extract proteins
wget -nc -O ${genome_path}/${build}/protein.fa.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02peptidesbglab16fagz
zcat ${genome_path}/${build}/protein.fa.gz > ${genome_path}/${build}/protein.fa
cat ${genome_path}/${build}/protein.fa | perl -pe "s/>(.*?)-P(.) .*/>\1-R\2/" ${genome_path}/${build}/protein.fa > ${genome_path}/${build}/protein2.fa
mv ${genome_path}/${build}/protein2.fa ${genome_path}/${build}/protein.fa

# Download/extract transcripts
wget -nc -O ${genome_path}/${build}/transcripts.fa.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02transcriptsbglab16fagz
zcat ${genome_path}/${build}/transcripts.fa.gz > ${genome_path}/${build}/transcripts.fa

# Download/extract gtf
wget -nc -O ${genome_path}/${build}/genes.gtf.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02basefeaturesbglab16gtfgz
zcat ${genome_path}/${build}/genes.gtf.gz > ${genome_path}/${build}/genes.gtf

# Download/extract gff
wget -nc -O ${genome_path}/${build}/genes.gff.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02basefeaturesbglab16gff3gz
zcat ${genome_path}/${build}/genes.gff.gz > ${genome_path}/${build}/genes.gff

# correct GTF (removes transcript lines, start_codon, stop_codon, and UTRs)
cd ${genome_path}/${build}/
gffread -E genes.gtf -T -o genes2.gtf
gffread -V -H  -N -J genes.gtf -T -o genes2.gtf -g ../genomes/BglaB1.6.fa
mv genes2.gtf genes.gtf

# Build genome 
cd ${genome_path}
snpEff build -gtf22 -v ${build}
# snpEff build -gff3 -v ${build}

# run snpEff
cp "${local_dir}"/bcftools.snp.fil3.vcf ${genome_path}/..
cd ${genome_path}/${build}
# snpEff -c ../../snpEff.config -v BglaB1.6 bcftools.snp.fil3.vcf > bcftools.snp.fil3.ann.vcf
cd ~/bin/snpEff
java -jar -Xmx16g snpEff.jar -c snpEff.config -v BglaB1.6 bcftools.snp.fil3.vcf > bcftools.snp.fil3.ann.vcf






