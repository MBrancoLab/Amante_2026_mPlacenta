#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=72:0:0


module load use.own
module load cellranger


#merge gencode annotation with expressed IAPs
#cat ~/genomes/mm10/gencode.vM25.annotation.gtf expressed_IAPs.gtf > gencode_iap.gtf

#make cell ranger reference
cellranger mkref \
    --genome=mm10 \
    --fasta=~/genomes/mm10/mm10.fa \
    --genes=gencode_iap.gtf

#count reads
cellranger count --id=placenta_iaps \
	--fastqs=fastq \
	--sample=E14 \
	--transcriptome=mm10

