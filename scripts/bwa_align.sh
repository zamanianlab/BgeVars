### Project
proj="Local_BgeVars"

# trimmomatic PE -phred33 -threads 4 ${GIT_DATA}/${proj}/CA301ANXX_1.fastq.gz ${GIT_DATA}/${proj}/CA301ANXX_2.fastq.gz ${GIT_DATA}/${proj}/CA301ANXX_P1.fq.gz ${GIT_DATA}/${proj}/CA301ANXX_U1.fq.gz ${GIT_DATA}/${proj}/CA301ANXX_P2.fq.gz ${GIT_DATA}/${proj}/CA301ANXX_U2.fq.gz ILLUMINACLIP:${GIT_PATH}/${proj}/auxillary/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

# wget -nc -O ${GIT_DATA}/${proj}/BglaB1.fa.gz https://www.vectorbase.org/download/biomphalaria-glabrata-bb02scaffoldsbglab1fagz
# zcat ${GIT_DATA}/${proj}/BglaB1.fa.gz > ${GIT_DATA}/${proj}/BglaB1.fa

# bwa index ${GIT_DATA}/${proj}/BglaB1.fa
# zcat  ${GIT_DATA}/${proj}/CA301ANXX_P1.fq.gz > ${GIT_DATA}/${proj}/CA301ANXX_P1.fastq
# zcat ${GIT_DATA}/${proj}/CA301ANXX_P2.fq.gz > ${GIT_DATA}/${proj}/CA301ANXX_P2.fastq
# bwa mem -t 4 ${GIT_DATA}/${proj}/BglaB1.fa ${GIT_DATA}/${proj}/CA301ANXX_P1.fastq ${GIT_DATA}/${proj}/CA301ANXX_P2.fastq > ${GIT_DATA}/${proj}/CA301ANXX.sam
# samtools view -@ 4 -bS ${GIT_DATA}/${proj}/CA301ANXX.sam > ${GIT_DATA}/${proj}/CA301ANXX.unsorted.bam
# samtools flagstat ${GIT_DATA}/${proj}/CA301ANXX.unsorted.bam
samtools sort -@ 4 -m 2G -o ${GIT_DATA}/${proj}/CA301ANXX.bam ${GIT_DATA}/${proj}/CA301ANXX.unsorted.bam 
# rm ${GIT_DATA}/${proj}/CA301ANXX.unsorted.bam
# zcat -c ${GIT_DATA}/${proj}/CA301ANXX.bam > ${GIT_DATA}/${proj}/CA301ANXX.bam.gz
samtools index -@ 4 -b ${GIT_DATA}/${proj}/CA301ANXX.bam
	
# samtools merge merge.bam ~/data/Bgla/BAM/*.bam 
# bedtools genomecov -ibam ~/data/Bgla/BAM/merge.bam -d > ~/data/Bgla/BAM/merge.d.bedgraph
# bedtools genomecov -ibam ~/data/Bgla/BAM/merge.bam -bga > ~/data/Bgla/BAM/merge.bga.bedgraph
bedtools genomecov -ibam ${GIT_DATA}/${proj}/CA301ANXX.bam -bga > ${GIT_DATA}/${proj}/CA301ANXX.bga.bedgraph
# bedtools genomecov -ibam ~/data/Bge/GCF_000457365.1_ASM45736v1_genomic.sorted.bam -d > ~/data/Bgla/BAM/BGE.d.bedgraph


# [zamanian@brc6 Bge]$ samtools flagstat CA301ANXX.bam
# 336645748 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 12391422 + 0 supplementary
# 0 + 0 duplicates
# 331826214 + 0 mapped (98.57% : N/A)
# 324254326 + 0 paired in sequencing
# 162127163 + 0 read1
# 162127163 + 0 read2
# 277408126 + 0 properly paired (85.55% : N/A)
# 316177310 + 0 with itself and mate mapped
# 3257482 + 0 singletons (1.00% : N/A)
# 31014804 + 0 with mate mapped to a different chr
# 18144990 + 0 with mate mapped to a different chr (mapQ>=5)