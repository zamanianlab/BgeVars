# BgeVars
Repository for analyzing variants between the Bge cell line and the *B. glabrata* reference genome

figures.R - includes all R commands in the analysis, preparation, and plotting for all figures

bwa_align.sh - bash script for trimming paired-end Illumina reads from the Bge cell line and aligning these to the *B. glabrata* BB02 reference genome and measuring read depth coverage (RDC) across the genome

coverage.sh - bash script for fetching paired-end Illumina reads from the *B. glabrata* BB02 referenence, re-aligning these to the assembled reference, and measuring RDC

snpeff.sh - bash scdript for running SnpEff analysis

variant_pipeline - bash script that inlcudes all commands and parameters for filterings variants