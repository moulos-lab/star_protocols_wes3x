#!/bin/bash

BAM_PATH=$HOME_PATH/bam
TRACKS_PATH=$HOME_PATH/tracks
GENOME_SIZE=$BEDTOOLS_PATH/../genomes/human.hg19.genome

if [ -d $TRACKS_PATH ]
then
    mkdir -p $TRACKS_PATH
fi

for FILE in `ls $BAM_PATH/*_fixmate.bam`
do 
    SAMPLE=`basename $FILE | sed s/_fixmate\.bam//`
    echo "Processing $SAMPLE"
    $BEDTOOLS_PATH/bedtools genomecov -bg \
    -ibam $BAM_PATH/$SAMPLE/$SAMPLE".bam" | \
        grep -vP 'chrU|rand|hap|loc|cox|GL|NC|hs37d5' | \
        awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' | \
        sed s/chrMT/chrM/g | \
        sort -k1,1 -k2g,2 > $TRACKS_PATH/$SAMPLE".bedGraph" &
done

wait

for FILE in `ls $TRACKS_PATH/*.bedGraph`
do
    echo "Processing $FILE"
    SAMPLE=`basename $FILE | sed s/\.bedGraph//`
    $UCSCTOOLS_PATH/bedGraphToBigWig $FILE $GENOME_SIZE $TRACKS_PATH/$SAMPLE".bigWig" &
done

wait

rm $TRACKS_PATH/*.bedGraph
