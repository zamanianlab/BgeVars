#!bin/bash

### Define project directories
boxdr=~/Box\ Sync
proj="Local_BgeVars"

gh_dir="${boxdr}/GitHub/${proj}"
local_dir="${boxdr}/GHdata/${proj}"

### BLAST
# Create blast db
# makeblastdb -out "${local_dir}"/auxillary/Bge_GCR_00457365 -in "${local_dir}"/auxillary/Bge_GCR_00457365.fa  -dbtype prot 

# # blastp for Bge actin
# tblastn -db "${local_dir}"/auxillary/Bge_GCR_00457365 -query "${gh_dir}"/auxillary/Cg_actin.fasta -out "${gh_dir}"/actin.blastout -outfmt 6

# get locus start/stop information
cat actin.blastout | awk '{print $2 " " $9 " " $10 " " $11}' > actin.blastout.locus

# get scaffold number from Bge_GCR_00457365.fa
cat actin.blastout | awk '{print $2}' | sort -uk1 > actin.contig.id
grep -Fwf actin.contig.id Bge_GCR_00457365.fa > contigs

bedtools getfasta -fi Bge_GCR_00457365.fa -bed actin.range.bed > actin.extract.fa
bedtools getfasta -fi Bge_GCR_00457365.fa -bed promoter.range.bed >> actin.extract.fa