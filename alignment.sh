#!/bin/bash

HOME_PATH=/PATH/TO/ANALYSIS/DIRECTORY
FASTQ_PATH=$HOME_PATH/fastq
BAM_PATH=$HOME_PATH/bam
THREADS=24
BWA_INDEX=$$RESOURCES_PATH/hs37d5/hs37d5.fa

if [ -d $BAM_PATH ]
then
    mkdir -p $BAM_PATH
fi

for FILE in `ls $FASTQ_PATH/*_1.fastq.gz`
do
    BASE=`basename $FILE | sed s/_1\.fastq\.gz//`
    F1=$FASTQ_PATH/$BASE"_1.fastq.gz"
    F2=$FASTQ_PATH/$BASE"_2.fastq.gz"
    
    RG="@RG\tID:"$BASE"\tSM:"$BASE"\tLB:WES\tPL:ILLUMINA"
    
    $BWA_PATH/bwa mem -t $THREADS -R $RG $BWA_INDEX $F1 $F2 | \
        $SAMTOOLS_PATH/samtools view -bS -o $BAM_PATH/$BASE".bam" -
done
