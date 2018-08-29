### Project
proj="BgeVars"

# BioSamples IDs for B. glabrata genome (paired end 454 reads)
SRR_list="SRR024007 SRR024008 SRR024017 SRR024018 SRR024019 SRR024020 SRR024021 SRR024022 SRR024023 SRR024024 SRR024025 SRR024026 SRR024027 SRR024028 SRR024031 SRR024032 SRR024033 SRR024034 SRR024035 SRR024036 SRR024037 SRR024038 SRR024039 SRR024040 SRR024041"

wget -nc -O ~/data/Bgla/Bgla_nt.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/457/365/GCF_000457365.1_ASM45736v1/GCF_000457365.1_ASM45736v1_genomic.fna.gz
zcat ~/data/Bgla/Bgla_nt.fa.gz > ~/data/Bgla/Bgla_nt.fa

for id in $SRR_list; do
	zcat ${GIT_DATA}/${proj}/snail/${id}_2.fastq.gz > ${GIT_DATA}/${proj}/snail/${id}_2.fastq
	zcat ${GIT_DATA}/${proj}/snail/${id}_4.fastq.gz > ${GIT_DATA}/${proj}/snail/${id}_4.fastq
	cat ${GIT_DATA}/${proj}/snail/${id}_2.fastq ${GIT_DATA}/${proj}/snail/${id}_4.fastq > ${GIT_DATA}/${proj}/snail/${id}.fastq
	bwa mem -t 4 ${GIT_DATA}/${proj}/BglaB1.fa ${GIT_DATA}/${proj}/snail/${id}.fastq > ${GIT_DATA}/${proj}/snail/${id}.sam
	samtools view -@ 4 -bS ${GIT_DATA}/${proj}/snail/${id}.sam > ${GIT_DATA}/${proj}/snail/${id}.unsorted.bam
	samtools flagstat ${GIT_DATA}/${proj}/snail/${id}.unsorted.bam
	samtools sort -@ 4 -m 2G ${GIT_DATA}/${proj}/snail/${id}.unsorted.bam -o ${GIT_DATA}/${proj}/snail/${id}.bam
	samtools index -@ 4 -b -m ${GIT_DATA}/${proj}/snail/${id}.bam
done

samtools merge -@ 4 ${GIT_DATA}/${proj}/snail/snail.bam ${GIT_DATA}/${proj}/snail/*.bam
samtools sort -@ 4 -m 3G ${GIT_DATA}/${proj}/snail/snail.bam -o ${GIT_DATA}/${proj}/snail/snail.sorted.bam
bedtools genomecov -ibam ~/data/Bgla/BAM/merge.bam -d > ~/data/Bgla/BAM/merge.d.bedgraph
bedtools genomecov -ibam ~/data/Bgla/BAM/merge.bam -bga > ~/data/Bgla/BAM/merge.bga.bedgraph
bedtools genomecov -ibam ${GIT_DATA}/${proj}/snail/snail.sorted.bam -bga > ${GIT_DATA}/${proj}/snail/snail.bga.bedgraph
bedtools genomecov -ibam ~/data/Bge/GCF_000457365.1_ASM45736v1_genomic.sorted.bam -d > ~/data/Bgla/BAM/BGE.d.bedgraph

bedtools coverage -sorted -a scaffold_windows.bed -b snail.sorted.bam > snail.cov.bed
