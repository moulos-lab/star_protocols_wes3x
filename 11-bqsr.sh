#!/bin/bash

BAM_PATH=$HOME_PATH/bam
CAPTURE_KIT=$RESOURCES_PATH/panel/Agilent_SureSelect_All_Exon_V2.bed
INTERVAL_LIST_PATH=$HOME_PATH/resources/interval_scatter
BWA_INDEX=$RESOURCES_PATH/hs37d5/hs37d5.fa
DBSNP=$RESOURCES_PATH/dbSNP/dbSNP151.vcf
GNOMAD=$RESOURCES_PATH/gnomAD/gnomad.exomes.r2.1.1.sites.vcf.bgz
CORES=16
PADDING=50

# Process dbSNP
$HTSLIB_PATH/bgzip $DBSNP
$HTSLIB_PATH/tabix $DBSNP”.gz”
DBSNP=$RESOURCES_PATH/dbSNP151.vcf.gz

mkdir -p $HOME_PATH/reports
META_REPORT=$HOME_PATH/reports/bsqr_current.log

echo "=== Splitting intervals" > $META_REPORT
if [ -d $INTERVAL_LIST_PATH ]
then
    echo "    Cleaning previous intervals" >> $META_REPORT
    rm -r $INTERVAL_LIST_PATH
fi
mkdir -p $INTERVAL_LIST_PATH

# Firstly split exome intervals for parallel BSQR

$GATK_PATH/gatk SplitIntervals \
    --reference $BWA_INDEX \
    --intervals $CAPTURE_KIT \
    --interval-padding $PADDING \
    --scatter-count $CORES \
    --output $INTERVAL_LIST_PATH \
    --QUIET

echo "=== Calculating BQSR tables" >> $META_REPORT
for FILE in `ls $BAM_PATH/*_fixmate.bam`
do
    SAMPLE=`basename $FILE | sed s/_fixmate\.bam//`
    echo "Processing $SAMPLE" >> $META_REPORT
    
    BAM=$BAM_PATH/$SAMPLE/$SAMPLE".bam"
    mkdir -p $BAM_PATH/$SAMPLE
    BQSR_PART_OUT=$BAM_PATH/$SAMPLE/bqsr_parts
    
    if [ -d $BQSR_PART_OUT ]
    then
        echo "    Cleaning previous tables" >> $META_REPORT
        rm -r $BQSR_PART_OUT
    fi
    mkdir -p $BQSR_PART_OUT
    
    # Calculate BQSR over intervals
    for INTERVAL in `readlink -f $INTERVAL_LIST_PATH/*`
    do      
        BQSR_NAME=`basename $INTERVAL | sed s/\-scattered\.interval_list//`
        echo "  Processing $BQSR_NAME" >> $META_REPORT
        $GATK_PATH/gatk BaseRecalibrator \
            --input $BAM \
            --reference $BWA_INDEX \
            --output $BQSR_PART_OUT/$BQSR_NAME".tab" \
            --known-sites $DBSNP \
            --known-sites $GNOMAD \
            --intervals $INTERVAL \
            --interval-padding $PADDING \
            --QUIET &
    done
    
    # Wait for individuals to complete before moving to the  next thread
    wait
done

echo "=== Gathering BQSR reports" >> $META_REPORT
for FILE in `ls $BAM_PATH/*_fixmate.bam`
do
    SAMPLE=`basename $FILE | sed s/_fixmate\.bam//`
    echo "Processing reports for $SAMPLE" >> $META_REPORT
    
    BQSR_PART_OUT=$BAM_PATH/$SAMPLE/bqsr_parts
    
    for TAB in `readlink -f $BQSR_PART_OUT/*`
    do
        echo "--input $TAB" >> $BAM_PATH/$SAMPLE/gather_bqsr.arg
    done
    
    # Gather reports
    $GATK_PATH/gatk GatherBQSRReports \
        --arguments_file $BAM_PATH/$SAMPLE/gather_bqsr.arg \
        --output $BAM_PATH/$SAMPLE/bqsr.tab \
        --QUIET &
done

# Wait for BQSR tables to be merged for each sample
wait

echo "=== Applying BQSR to BAM files" >> $META_REPORT
for FILE in `ls $BAM_PATH`
do
    SAMPLE=`basename $FILE | sed s/_fixmate\.bam//`
    echo "Processing BAM file $SAMPLE" >> $META_REPORT
    
    BAM=$BAM_PATH/$SAMPLE/$SAMPLE".bam"
    BQSR_TABLE=$BAM_PATH/$SAMPLE/bqsr.tab
    
    # Apply BQSR to BAM files
    $GATK_PATH/gatk ApplyBQSR  \
        --input $BAM \
        --reference $BWA_INDEX \
        --bqsr-recal-file $BQSR_TABLE \
        --output $BAM_PATH/$SAMPLE/$SAMPLE"_bqsr.bam" \
        --QUIET &
done

# Wait for new BAM files to be created before reporting finished
wait

echo "=== Finished!" >> $META_REPORT
