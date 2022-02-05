#!/bin/bash

BAM_PATH=$HOME_PATH/bam
CORES=16

for SAMPLE in `ls $BAM_PATH` 
do
    echo "Processing $SAMPLE"
    $SAMTOOLS_PATH/samtools sort -n -@ $CORES -m 4G \
        $BAM_PATH/$SAMPLE".bam" | \
    $SAMTOOLS_PATH/samtools fixmate -m – 
    $BAM_PATH/$SAMPLE"_fixmate.bam"
done

for SAMPLE in `ls $BAM_PATH` 
do
    echo "Processing $SAMPLE"
    $SAMTOOLS_PATH/samtools sort -n -@ $CORES -m 4G \
        $BAM_PATH/$SAMPLE".bam" | \
    $SAMTOOLS_PATH/samtools fixmate -m – 
    $BAM_PATH/$SAMPLE"_fixmate.bam"
done

for SAMPLE in `ls $BAM_PATH`
do
    echo "Processing $SAMPLE"
    $SAMTOOLS_PATH/samtools sort -@ $CORES -m 4G \
        $BAM_PATH/$SAMPLE"_fixmate.bam" | \
        $SAMTOOLS_PATH/samtools markdup - $BAM_PATH/$SAMPLE".bam"
    echo "Indexing $SAMPLE"
    $SAMTOOLS_PATH/samtools index $BAM_PATH/$SAMPLE".bam"
done
