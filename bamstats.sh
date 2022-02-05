#!/bin/bash

FASTQ_PATH=$HOME_PATH/fastq
CAPTURE_KIT=$HOME_PATH/resources/MGI_Exome_Capture_V4.bed
BAM_PATH=$HOME_PATH/bam
SAMTOOLS_PATH=PATH/TO/SAMTOOLS/samtools
BEDTOOLS_PATH=/PATH/TO/BEDTOOLS/bedtools
REPORT=$HOME_PATH/reports/finalbamstats.txt

printf "%s\t%s\t%s\t%s\t%s\t%s%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "name" \
    "total reads" "total reads pairs" "aligned reads" \
     "properly paired aligned pairs" "uniquely aligned reads (q>20)" \
     "properly paired uniquely aligned reads" "chimeric reads" \
      "reads overlapping targets" "total bases" "aligned bases" \
      "uniquely aligned bases" "bases overlapping targets" > $REPORT

for SAMPLE in `ls $BAM_PATH`
do
    echo "Processing $SAMPLE"
    
    BAM=$BAM_PATH/$SAMPLE/$SAMPLE"_markdup.bam"
    
    printf "%s\t" $SAMPLE >> $REPORT
    
    echo "  total reads"
    printf "%d\t" `$SAMTOOLS_PATH/samtools view -c -F2048 $BAM` >> $REPORT
    
    echo "  total read pairs"
    printf "%d\t" `$SAMTOOLS_PATH/samtools view -c -F2048 $BAM | awk '{print $1/2}'` \
        >> $REPORT
    
    echo "  aligned reads"
    printf "%d\t" `$SAMTOOLS_PATH/samtools view -c -F2052 $BAM` >> $REPORT
    
    echo "  properly paired aligned pairs"
    printf "%d\t" `$SAMTOOLS_PATH/samtools view -c -f66 -F2048 $BAM` \
        >> $REPORT
        
    echo "  uniquely aligned reads (q>20)"
    printf "%d\t" `$SAMTOOLS_PATH/samtools view -c -F2052 -q20 $BAM` >> \
        $REPORT
    
    echo "  properly paired uniquely aligned reads"
    printf "%d\t" `$SAMTOOLS_PATH/samtools view -c -f66 -F2048 -q20 $BAM` \
        >> $REPORT
    
    echo "  chimeric reads"
    printf "%d\t" `
        $SAMTOOLS_PATH/samtools flagstat $BAM | \
        perl -e 'my @in;' \
            -e 'while(<>) { chomp $_; push(@in,$_); }' \
            -e 'my @tmp = split("\\\+",pop(@in));' \
            -e '$tmp[0] =~ s/\s+$//;' \
            -e 'print STDOUT $tmp[0];'
    ` >> $REPORT
    
    echo "  reads overlapping targets"
    printf "%d\t" `
        $BEDTOOLS_PATH/bedtools intersect -a $CAPTURE_KIT -b $BAM -c | \
            awk 'BEGIN {tot=0}{tot+=$4} END {print tot}'
    ` >> $REPORT
    
    echo "  total bases"
    printf "%d\t" `
        $SAMTOOLS_PATH/samtools view $BAM | cut -f10 | \
            awk 'BEGIN {tr=0}{tr+=length($0)} END {print tr}'
    ` >> $REPORT
    
    echo "  aligned bases"
    printf "%d\t" `
        $SAMTOOLS_PATH/samtools view -F2052 $BAM | cut -f10 | \
            awk 'BEGIN {tr=0}{tr+=length($0)} END {print tr}'
    ` >> $REPORT
    
    echo "  uniquely aligned bases"
    printf "%d\t" `
        $SAMTOOLS_PATH/samtools view -F2052 -q20 $BAM | cut -f10 | \
            awk 'BEGIN {tr=0}{tr+=length($0)} END {print tr}'
    ` >> $REPORT
    
    echo "  bases overlapping targets"
    printf "%d\n" `
        $BEDTOOLS_PATH/bedtools coverage -a $CAPTURE_KIT -b $BAM -d | \
            awk 'BEGIN {tr=0} {tr+=$5} END {print tr}'
    ` >> $REPORT
    
done
