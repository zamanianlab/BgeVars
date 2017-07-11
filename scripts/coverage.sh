# BioSamples IDs for B. glabrata genome (paired end 454 reads)
SRR_list="SRR024007 SRR024008 SRR024017 SRR024018 SRR024019 SRR024020 SRR024021 SRR024022 SRR024023 SRR024024 SRR024025 SRR024026 SRR024027 SRR024028 SRR024031 SRR024032 SRR024033 SRR024034 SRR024035 SRR024036 SRR024037 SRR024038 SRR024039 SRR024040 SRR024041"

# wget -nc -O ~/data/Bgla/Bgla_nt.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_genomic.fna.gz
# zcat ~/data/Bgla/Bgla_nt.fa.gz > ~/data/Bgla/Bgla_nt.fa

# bwa index ~/data/Bgla/Bgla_nt.fa

# for id in $SRR_list; do
	#fastq-dump --split-files $id -O ~/data/Bgla/SRA
	# cat ~/data/Bgla/SRA/${id}_2.fastq ~/data/Bgla/SRA/${id}_4.fastq > ~/data/Bgla/SRA/${id}.fastq
	# bwa mem -t 8 ~/data/Bgla/Bgla_nt.fa ~/data/Bgla/SRA/${id}.fastq > ~/data/Bgla/SAM/${id}.sam
	# samtools view -bS ~/data/Bgla/SAM/${id}.sam > ~/data/Bgla/BAM/${id}.unsorted.bam
	# samtools flagstat ~/data/Bgla/BAM/${id}.unsorted.bam
	# samtools sort -@ 8 -o ~/data/Bgla/BAM/${id}.bam ~/data/Bgla/BAM/${id}.unsorted.bam 
	# samtools index -b ~/data/Bgla/BAM/${id}.bam
# done
	
# samtools merge merge.bam ~/data/Bgla/BAM/*.bam 
# bedtools genomecov -ibam ~/data/Bgla/BAM/merge.bam -d > ~/data/Bgla/BAM/merge.d.bedgraph
# bedtools genomecov -ibam ~/data/Bgla/BAM/merge.bam -bga > ~/data/Bgla/BAM/merge.bga.bedgraph
bedtools genomecov -ibam ~/data/Bge/GCF_000457365.1_ASM45736v1_genomic.sorted.bam -bga > ~/data/Bgla/BAM/BGE.bga.bedgraph
# bedtools genomecov -ibam ~/data/Bge/GCF_000457365.1_ASM45736v1_genomic.sorted.bam -d > ~/data/Bgla/BAM/BGE.d.bedgraph
