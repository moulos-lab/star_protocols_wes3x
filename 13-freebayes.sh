#!/bin/bash

export VCF_PATH=$HOME_PATH/vcf
BAM_PATH=$HOME_PATH/bam
CAPTURE_KIT=$RESOURCES_PATH/panel/Agilent_SureSelect_All_Exon_V2.bed
INTERVAL_LIST_PATH=$RESOURCES_PATH/resources/interval_scatter_bed
BWA_INDEX=$RESOURCES_PATH/hs37d5/hs37d5.fa
CORES=32
PADDING=50

META_REPORT=$HOME_PATH/reports/freebayes_current.log

echo "=== Splitting intervals" > $META_REPORT
if [ -d $INTERVAL_LIST_PATH ]
then
    echo "    Cleaning previous intervals" >> $META_REPORT
    rm -r $INTERVAL_LIST_PATH
fi
mkdir -p $INTERVAL_LIST_PATH

# Firstly split exome intervals for parallel freebayes

$GATK_PATH/gatk SplitIntervals \
    --reference $BWA_INDEX \
    --intervals $CAPTURE_KIT \
    --interval-padding $PADDING \
    --scatter-count $CORES \
    --extension .pre \
    --output $INTERVAL_LIST_PATH \
    --QUIET

echo "=== Converting intervals" >> $META_REPORT
for INTERVAL in `ls $INTERVAL_LIST_PATH`
do
    BED=`basename $INTERVAL | sed s/\.pre//`
    INTERVAL_FILE=$INTERVAL_LIST_PATH/$INTERVAL
    grep -vP "^@" $INTERVAL_FILE | awk '{print $1"\t"$2"\t"$3}' > \
        $INTERVAL_LIST_PATH/$BED".bed" &
done

# Wait and clear intermediate intervals
wait
rm $INTERVAL_LIST_PATH/*.pre

# Prepare BAM file list for freebayes
echo "=== Preparing BAM file list" >> $META_REPORT
BAMLIST=/media/raid/tmp/tmp/medex/scripts/bamlist.txt
if [ -f $BAMLIST ]
then
    rm $BAMLIST
fi
for FILE in `ls $BAM_PATH/*_fixmate\.bam`
do
    SAMPLE=`basename $FILE | sed s/_fixmate\.bam//`
    BAM=$BAM_PATH/$SAMPLE/$SAMPLE".bam"
    echo "$BAM" >> $BAMLIST
done

echo "=== Calling variants with FreeBayes" >> $META_REPORT
if [ -d $VCF_PATH/fb_parts ]
then
    rm -r $VCF_PATH/fb_parts
fi
mkdir -p $VCF_PATH/fb_parts

for TARGET in `ls $INTERVAL_LIST_PATH`
do
    NAME=`basename $TARGET | sed s/\.bed//`
    echo "Processing interval list $NAME" >> $META_REPORT
    INTERVAL=$INTERVAL_LIST_PATH/$TARGET
    $FREEBAYES_PATH/freebayes \
        --fasta-reference $BWA_INDEX \
        --bam-list $BAMLIST \
        --targets $INTERVAL \
        --vcf $VCF_PATH/fb_parts/$NAME".vcf" &
done

# Wait before gathering the results
wait
echo "=== Merging VCFs" >> $META_REPORT
cat $VCF_PATH/*.vcf | \
    $VCFLIB_PATH/scripts/vcffirstheader | \
    $VCFLIB_PATH/bin/vcfstreamsort -w 1000 | \
    $VCFLIB_PATH/bin/vcfuniq > \
    $VCF_PATH/freebayes_full.vcf

echo "=== Compressing and indexing final VCF" >> $META_REPORT
$HTSLIB_PATH/bgzip $VCF_PATH/freebayes_full.vcf
$HTSLIB_PATH/tabix $VCF_PATH/freebayes_full.vcf.gz

### Basic filtering before decomposing and normalization

# Determine a quality and depth cutoff pre-filter based on 99th percentile of
# the respective distributions
echo "=== Determining QUAL and DP hard pre-filters" >> $META_REPORT
$BCFTOOLS_PATH/bcftools query \
    --include 'QUAL>20' \
    --format '%QUAL\n' $VCF_PATH/freebayes_full.vcf.gz > quals.tmp &
$BCFTOOLS_PATH/bcftools query \
    --include 'QUAL>20' \
    --format '%INFO/DP\n' $VCF_PATH/freebayes_full.vcf.gz | \
    awk -F "," '{print $1}' > $VCF_PATH/dps.tmp &
wait

Rscript -e '
    vp <- Sys.getenv("VCF_PATH")
    dps <- as.numeric(readLines(file.path(vp,"dps.tmp")));
    quals <- as.numeric(readLines(file.path(vp,"quals.tmp")));
    qudp <- unname(round(quantile(dps,0.99)));
    ququ <- unname(quantile(quals,0.99));
    write(qudp,file.path(vp,"dpt.tmp"));
    write(ququ,file.path(vp,"qut.tmp"));
'
QUALUP=`cat $VCF_PATH/qut.tmp`
DPUP=`cat $VCF_PATH/dpt.tmp`
rm $VCF_PATH/qut.tmp $VCF_PATH/dpt.tmp $VCF_PATH/dps.tmp $VCF_PATH/quals.tmp

# Apply the filters, decompose complex variants and normalize
echo "=== Applying filters and normalizing" >> $META_REPORT
$BCFTOOLS_PATH/bcftools view \
    --include 'QUAL>20 & INFO/DP>10 & QUAL<'$QUALUP' & INFO/DP<'$DPUP' & (QUAL/(INFO/DP))>2' $VCF_PATH/freebayes_full.vcf.gz | \
    $VCFLIB_PATH/bin/vcfallelicprimitives -kg | \
    $BCFTOOLS_PATH/bcftools norm \
    --fasta-ref $BWA_INDEX \
    --output-type z \
    --output $VCF_PATH/freebayes_filtered_norm.vcf.gz
$HTSLIB_PATH/tabix $VCF_PATH/freebayes_filtered_norm.vcf.gz

echo "=== Finished!" >> $META_REPORT
