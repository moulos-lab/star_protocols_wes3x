#!/bin/bash

export VCF_PATH=$HOME_PATH/vcf
BAM_PATH=$HOME_PATH/bam
CAPTURE_KIT_DIR=$RESOURCES_PATH/resources
CAPTURE_KIT=$RESOURCES_PATH/panel/Agilent_SureSelect_All_Exon_V2.bed
DV_VERSION=0.9.0
BWA_INDEX_DIR=$RESOURCES_PATH/hs37d5
BWA_INDEX=$RESOURCES_PATH/hs37d5/hs37d5.fa
CORES=32

META_REPORT=/media/raid/tmp/tmp/medex/scripts/deepvariant_current.log

echo "=== Calling variants" > $META_REPORT
for SAMPLE in `ls $BAM_PATH`
do
    echo "Processing $SAMPLE" >> $META_REPORT
    
    BAM=$BAM_PATH/$SAMPLE/$SAMPLE".bam"
    
    #GVCF_NAME=$BAM_PATH/$SAMPLE/$SAMPLE"_DV.g.vcf"
    #VCF_NAME=$BAM_PATH/$SAMPLE/$SAMPLE"_DV.vcf"
    
    docker run \
        -v "$BAM_PATH":"/data" \
        -v "$BWA_INDEX_DIR":"/reference" \
        -v "$CAPTURE_KIT_DIR":"/capture_kit" \
        google/deepvariant:$DV_VERSION \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WES \
        --ref="/reference/hs37d5.fa" \
        --reads="/data/$SAMPLE/$SAMPLE.bam" \
        --regions="/capture_kit/MGI_Exome_Capture_V4.bed" \
        --output_vcf="/data/$SAMPLE/$SAMPLE'_DV.vcf'" \
        --output_gvcf="/data/$SAMPLE/$SAMPLE'_DV.g.vcf'" \
        --num_shards=$CORES
    
    ## Wait for indiviuals to complete before moving to the  next thread
    #wait
done

echo "=== Creating list of gVCF files" >> $META_REPORT
for SAMPLE in `ls $BAM_PATH`
do
    GVCF=`readlink -f $BAM_PATH/$SAMPLE/$SAMPLE"_DV.g.vcf"`
    echo "$GVCF" >> $VCF_PATH/deepvariant_gvcf_list.txt
done

echo "=== Gathering gVCFs" >> $META_REPORT
rm -r GLnexus.DB
$GLNEXUS_PATH/glnexus_cli \
    --config DeepVariantWES \
    --bed $CAPTURE_KIT \
    --list $VCF_PATH/deepvariant_gvcf_list.txt \
    --threads $CORES | \
    $BCFTOOLS_PATH/bcftools view --include 'QUAL>=20' - | \
    $BCFTOOLS_PATH/bcftools norm \
    --fasta-ref $BWA_INDEX \
    --output-type z \
    --output $VCF_PATH/deepvariant_filtered_norm.vcf.gz
$HTSLIB_PATH/tabix $VCF_PATH/deepvariant_filtered_norm.vcf.gz

echo "=== Finished!" >> $META_REPORT
