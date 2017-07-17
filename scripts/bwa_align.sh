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
