#!/bin/bash

#!/bin/bash

BAM_PATH=$HOME_PATH/bam
CORES=16

for FILE in `ls $BAM_PATH/*.uns` 
do
    SAMPLE=`basename $FILE | sed s/\.uns//`
    echo "Processing $SAMPLE"
    $SAMTOOLS_PATH/samtools sort -n -@ $CORES -m 4G \
        $BAM_PATH/$SAMPLE".uns" | \
    $SAMTOOLS_PATH/samtools fixmate -m - \
    $BAM_PATH/$SAMPLE"_fixmate.bam"
done
rm $BAM_PATH/*.uns

for FILE in `ls $BAM_PATH/*_fixmate.bam`
do
    SAMPLE=`basename $FILE | sed s/_fixmate\.bam//`
    echo "Processing $SAMPLE"
    $SAMTOOLS_PATH/samtools sort -@ $CORES -m 4G \
        $BAM_PATH/$SAMPLE"_fixmate.bam" | \
        $SAMTOOLS_PATH/samtools markdup - $BAM_PATH/$SAMPLE".bam"
    echo "Indexing $SAMPLE"
    $SAMTOOLS_PATH/samtools index $BAM_PATH/$SAMPLE".bam"
done

## For single-end BAM files
#for FILE in `ls $BAM_PATH/*.uns`
#do
#    SAMPLE=`basename $FILE | sed s/\.uns//`
#    echo "Processing $SAMPLE"
#    $SAMTOOLS_PATH/samtools sort -@ $CORES -m 4G \
#        $BAM_PATH/$SAMPLE".uns" | \
#    $SAMTOOLS_PATH/samtools markdup - $BAM_PATH/$SAMPLE".bam"
#    echo "Indexing $SAMPLE"
#    $SAMTOOLS_PATH/samtools index $BAM_PATH/$SAMPLE".bam"
#done
