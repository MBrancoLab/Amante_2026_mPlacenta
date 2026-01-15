#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=12:0:0
#$ -l h_vmem=12G
#$ -t 1-5


r1=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
r2=$(echo $r1 | sed s/_1.fq.gz/_2.fq.gz/)
base=$(basename $r1 | sed s/_1.fq.gz//)


########
#SQuIRE#
########


module load anaconda2

rlength=$(zcat $r1 | head | sed -n "2p" | tr -d '\n' | wc -m)

##activate virtual environment
source activate squire

##map
squire Map -1 $r1 -2 $r2 -f ~/genomes/mm10/squire_fetch -r $rlength -n $base -o ${base}_map -p 4 -v

##count
squire Count -m squire_map -c ~/genomes/mm10/squire_clean \
	-f ~/genomes/mm10/squire_fetch \
	-r $rlength \
	-n $base \
	-m ${base}_map \
	-o ${base}_count \
	-p 4 -s 1 -v
source deactivate squire


###############
#TEtranscripts#
###############

module unload anaconda2
module load python/3.8.5

##TEtranscripts
~/.local/bin/TEcount -b ${base}_map/$base.bam \
	--GTF ~/genomes/mm10/squire_fetch/mm10_refGene.gtf \
	--TE ~/genomes/mm10/TEtranscripts/mm10_rmsk_TE.gtf \
	--stranded reverse \
	--sortByPos \
	--project ${base}_TEtranscripts

