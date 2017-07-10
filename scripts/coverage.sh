# BioSamples IDs for B. glabrata genome (paired end 454 reads)
# SRR_list="SRR024007 SRR024008 SRR024017 SRR024018 SRR024019 SRR024020 SRR024021 SRR024022 SRR024023 SRR024024 SRR024025 SRR024026 SRR024027 SRR024028 SRR024031 SRR024032 SRR024033 SRR024034 SRR024035 SRR024036 SRR024037 SRR024038 SRR024039 SRR024040 SRR024041"

# for id in SRR_list; do
# 	fastq-dump $id --gzip -O ~/data/SRA/Bgla 
# done

wget -O Bgla_nt.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_genomic.fna.gz
zcat Bgla_nt.fa.gz > Bgla_nt.fa

bwa index Bgla_nt.fa
# bwa mem -t 8 -P Bgla_nt SRR024007_1.fastq SRR024007_2.fastq

### Calculate coverage 
# overall coverage: 
# bedtools genomecov -ibam ~/data/Bge/GCF_000457365.1_ASM45736v1_genomic.sorted.bam -bga > ~/data/Bge/snail.bam.bedgraph
# vs (per base)
# bedtools genomecov -ibam ~/data/Bge/GCF_000457365.1_ASM45736v1_genomic.sorted.bam -d > ~/data/Bge/snail.bam.bedgraph


