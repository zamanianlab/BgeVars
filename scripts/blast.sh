#!bin/bash

### Define project directories
boxdr=~/Box\ Sync
proj="Local_BgeVars"

gh_dir="${boxdr}/GitHub/${proj}"
local_dir="${boxdr}/GHdata/${proj}"

### BLAST
# Create blast db
makeblastdb -out GCF_000457365.1_ASM45736v1_genomic -in GCF_000457365.1_ASM45736v1_genomic.fna -dbtype nucl
makeblastdb -out GCF_000457365.1_ASM45736v1_protein -in GCF_000457365.1_ASM45736v1_protein.faa -dbtype prot


# # tblastn for Bge actin
tblastn -db GCF_000457365.1_ASM45736v1_genomic -query Cg_actin.fasta -out actin.blastout -outfmt 6
tblastn -db GCF_000457365.1_ASM45736v1_genomic -query Cg_actin.fasta -out actin.verbose.blastout

blastp -db GCF_000457365.1_ASM45736v1_protein -query Cg_actin.fasta -out actin.blastpout -outfmt 6
blastp -db GCF_000457365.1_ASM45736v1_protein -query Cg_actin.fasta -out actin.verbose.blastpout

# get locus start/stop information
cat actin.blastout | awk '{print $2 " " $9 " " $10 " " $11}' > actin.blastout.locus

# get scaffold number from Bge_GCR_00457365.fa
cat actin.blastout | awk '{print $2}' | sort -uk1 > actin.contig.id
grep -Fwf actin.contig.id Bge_GCR_00457365.fa > contigs

bedtools getfasta -fi Bge_GCR_00457365.fa -bed actin.range.bed > actin.extract.fa
bedtools getfasta -fi Bge_GCR_00457365.fa -bed promoter.range.bed >> actin.extract.fa