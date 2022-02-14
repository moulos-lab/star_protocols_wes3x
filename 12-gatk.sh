#!/bin/bash

export VCF_PATH=$HOME_PATH/vcf
BAM_PATH=$HOME_PATH/bam
INTERVAL_LIST_PATH=$RESOURCES_PATH/panel/interval_scatter
BWA_INDEX=$RESOURCES_PATH/hs37d5/hs37d5.fa
CORES=16
PADDING=50

META_REPORT=$HOME_PATH/reports/haca_current.log

echo "=== Calling variants" > $META_REPORT
for SAMPLE in `ls $BAM_PATH`
do
    echo "Processing $SAMPLE" >> $META_REPORT
    
    BAM=$BAM_PATH/$SAMPLE/$SAMPLE"_bqsr.bam"
    GVCF_PART_OUT=$BAM_PATH/$SAMPLE/gvcf_parts
    
    if [ -d $GVCF_PART_OUT ]
    then
        echo "    Cleaning previous gVCFs" >> $META_REPORT
        rm -r $GVCF_PART_OUT
    fi
    mkdir -p $GVCF_PART_OUT
    
    # Call variants over intervals
    for INTERVAL in `readlink -f $INTERVAL_LIST_PATH/*`
    do      
        GVCF_NAME=`basename $INTERVAL | sed s/\-scattered\.interval_list//`
        echo "  Processing $GVCF_NAME" >> $META_REPORT
        $GATK_PATH/gatk HaplotypeCaller \
            --input $BAM \
            --reference $BWA_INDEX \
            --intervals $INTERVAL \
            --interval-padding $PADDING \
            --output $GVCF_PART_OUT/$GVCF_NAME".g.vcf" \
            --emit-ref-confidence GVCF \
            --create-output-variant-index false \
            --QUIET &
    done
    
    # Wait for individuals to complete before moving to the next thread
    wait
done

# Then GVCFs must be consolidated
echo "=== Merging gVCFs" >> $META_REPORT
for SAMPLE in `ls $BAM_PATH`
do
    echo "Processing interval gVCFs for $SAMPLE"
    
    GVCF_PART_OUT=$BAM_PATH/$SAMPLE/gvcf_parts
    
    if [ -f $BAM_PATH/$SAMPLE/interval_gvcfs.txt ]
    then
        echo "    Cleaning previous gVCFs input file" >> $META_REPORT
        rm $BAM_PATH/$SAMPLE/interval_gvcfs.txt
    fi
    for GVCF in `readlink -f $GVCF_PART_OUT/*.g.vcf`
    do
        echo "$GVCF" >> $BAM_PATH/$SAMPLE/interval_gvcfs.txt
    done
    
    # Get the gVCF header and strip the GATK command 
    GVFH=`readlink -f $GVCF_PART_OUT/*.g.vcf | head -1`
    grep "^#" $GVFH | grep -v "^##GATKCommand" > $BAM_PATH/$SAMPLE/gvcf.header
    
    # Cat the gVCFs
    for GVCF in `readlink -f $GVCF_PART_OUT/*.g.vcf`
    do
        echo "  Concatenating $GVCF"
        #echo "  Concatenating $GVCF" >> $META_REPORT
        grep -v "^#" $GVCF >> $BAM_PATH/$SAMPLE/gvcf.tmp
    done
    
    # Place the header
    echo "  Creating final gVCF"
    #echo "  Creating final gVCF" >> $META_REPORT
    cat $BAM_PATH/$SAMPLE/gvcf.header $BAM_PATH/$SAMPLE/gvcf.tmp > \
        $BAM_PATH/$SAMPLE/$SAMPLE".u.g.vcf"
        
    rm $BAM_PATH/$SAMPLE/gvcf.tmp $BAM_PATH/$SAMPLE/gvcf.header
done

# Sort gVCFs
echo "=== Sorting gVCFs" >> $META_REPORT
for SAMPLE in `ls $BAM_PATH`
do
    echo "Sorting gVCF for $SAMPLE" >> $META_REPORT
        
    $GATK_PATH/gatk SortVcf \
        --INPUT $BAM_PATH/$SAMPLE/$SAMPLE".u.g.vcf" \
        --OUTPUT $BAM_PATH/$SAMPLE/$SAMPLE".g.vcf.gz" \
        --QUIET &
done

# Wait for sorting to finish before cleaning unsorted
wait 

# Some cleanup
echo "=== Deleting unsorted gVCFs" >> $META_REPORT
for FILE in `ls $BAM_PATH/*_fixmate.bam`
do
    SAMPLE=`basename $FILE | sed s/_fixmate\.bam//`
    echo "Deleting unsorted gVCF for $SAMPLE" >> $META_REPORT
    rm $BAM_PATH/$SAMPLE/$SAMPLE".u.g.vcf"
    echo "Compression gVCF parts for $SAMPLE" >> $META_REPORT
    pigz $BAM_PATH/$SAMPLE/gvcf_parts/*
    echo "Compression BQSR reports for $SAMPLE" >> $META_REPORT
    pigz $BAM_PATH/$SAMPLE/bqsr_parts/*
done

# Gather VCFs
echo "=== Combining sorted population gVCFs" >> $META_REPORT
if [ ! -d $VCF_PATH ]
then
    mkdir $VCF_PATH
fi
# Delete the .arg file as it will get multiple entries
if [ -f $VCF_PATH/combine_gvcf.arg ]
then
    rm $VCF_PATH/combine_gvcf.arg
fi

for FILE in `ls $BAM_PATH/*_fixmate.bam`
do
    SAMPLE=`basename $FILE | sed s/_fixmate\.bam//`
    GVCF=`readlink -f $BAM_PATH/$SAMPLE/$SAMPLE".g.vcf.gz"`
    echo "--variant $GVCF" >> $VCF_PATH/combine_gvcf.arg
done

# Combine gVCFs
$GATK_PATH/gatk CombineGVCFs \
 --reference $BWA_INDEX \
 --arguments_file $VCF_PATH/combine_gvcf.arg \
 --output $VCF_PATH/haplotypecaller_full.g.vcf.gz

# Genotype VCFs
echo "=== Genotyping gVCFs" >> $META_REPORT
$GATK_PATH/gatk GenotypeGVCFs \
    --reference $BWA_INDEX \
    --variant $VCF_PATH/haplotypecaller_full.g.vcf.gz \
    --output $VCF_PATH/haplotypecaller_full.vcf.gz

# Apply basic GATK hard filters
echo "=== Applying GATK hard filters" >> $META_REPORT
$BCFTOOLS_PATH/bcftools view \
    --include 'QUAL>20 & INFO/QD>2 & INFO/MQ>40 & INFO/FS<60 & INFO/SOR<3 & INFO/MQRankSum>-12.5 & INFO/ReadPosRankSum>-8 & TYPE="snp"' \
    --output-type z \
    --output-file $VCF_PATH/haplotypecaller_filtered_snp.vcf.gz \
    $VCF_PATH/haplotypecaller_full.vcf.gz &

# The normalization step is potentially not required but it is harmless
$BCFTOOLS_PATH/bcftools view \
    --include 'QUAL>20 & INFO/QD>2 & INFO/ReadPosRankSum>-20 & INFO/InbreedingCoeff>-0.8 & INFO/FS<200 & INFO/SOR<10 & TYPE~"indel"' \
    $VCF_PATH/haplotypecaller_full.vcf.gz | \
    $BCFTOOLS_PATH/bcftools norm \
    --fasta-ref $BWA_INDEX \
    --output-type z \
    --output $VCF_PATH/haplotypecaller_filtered_norm_indel.vcf.gz &
wait

echo "=== Merging GATK filtered SNPs and INDELs" >> $META_REPORT
$GATK_PATH/gatk MergeVcfs \
    --INPUT $VCF_PATH/haplotypecaller_filtered_snp.vcf.gz \
    --INPUT $VCF_PATH/haplotypecaller_filtered_norm_indel.vcf.gz \
    --OUTPUT $VCF_PATH/haplotypecaller_filtered_norm.vcf.gz \
    --QUIET

rm $VCF_PATH/haplotypecaller_filtered_snp.vcf.gz \
    $VCF_PATH/haplotypecaller_filtered_norm_indel.vcf.gz

#echo "=== Finished!"
echo "=== Finished!" >> $META_REPORT

