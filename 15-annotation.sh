#!/bin/bash

export VCF_PATH=$HOME_PATH/vcf
DBSNP_FILE=$RESOURCES_PATH/dbSNP/dbSNP151.vcf.gz
DBNSFP_FILE=$RESOURCES_PATH/dbNSFP/dbNSFP2.9.3.txt.gz
GNOMAD_FILE=$RESOURCES_PATH/gnomAD/gnomad.exomes.r2.1.1.sites.vcf.bgz

if [ ! -d $SNPEFF_PATH/data ]
then
    java -jar $SNPEFF_PATH/snpEff.jar download GRCh37.75
fi

## Haplotype Caller
# Variant effect annotation
java -Xmx4096m -jar $SNPEFF_PATH/snpEff.jar ann \
    -v -noLog -noStats -noLof GRCh37.75 \
    $VCF_PATH/haplotypecaller_filtered_norm.vcf.gz > $VCF_PATH/haplotypecaller_filtered_norm_eff.vcf 
$HTSLIB_PATH/bgzip $VCF_PATH/haplotypecaller_filtered_norm_eff.vcf

$HTSLIB_PATH/tabix $VCF_PATH/haplotypecaller_filtered_norm_eff.vcf.gz

# Annotation with dbSNP
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar annotate \
    -v -id $DBSNP_FILE \
    $VCF_PATH/haplotypecaller_filtered_norm_eff.vcf.gz > $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp.vcf 
    $HTSLIB_PATH/bgzip
$VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp.vcf.gz    
$HTSLIB_PATH/tabix $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp.vcf.gz

# Annotation with dbNSFP
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar dbnsfp \
    -v -m -db $DBNSFP_FILE \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp.vcf.gz > $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf
    $HTSLIB_PATH/bgzip $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf
$HTSLIB_PATH/tabix $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
    
# Annotation with gnomAD
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar annotate \
    -v $GNOMAD_FILE \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz > $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf
    $HTSLIB_PATH/bgzip $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf
$HTSLIB_PATH/tabix $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf.gz

# FreeBayes
# Variant effect annotation
java -Xmx4096m -jar $SNPEFF_PATH/snpEff.jar ann \
    -v -noLog -noStats -noLof GRCh37.75 \
    $VCF_PATH/freebayes_filtered_norm.vcf.gz > $VCF_PATH/freebayes_filtered_norm_eff.vcf
    $HTSLIB_PATH/bgzip  $VCF_PATH/freebayes_filtered_norm_eff.vcf
    $HTSLIB_PATH/tabix \
    $VCF_PATH/freebayes_filtered_norm_eff.vcf.gz

# Annotation with dbSNP
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar annotate \
    -v -id $DBSNP_FILE \
    $VCF_PATH/freebayes_filtered_norm_eff.vcf.gz > $VCF_PATH/freebayes_filtered_norm_eff_dbsnp.vcf
    $HTSLIB_PATH/bgzip\
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp.vcf.gz

    $HTSLIB_PATH/tabix \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp.vcf.gz

# Annotation with dbNSFP
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar dbnsfp \
    -v -m -db $DBNSFP_FILE \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp.vcf.gz > $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf 
    $HTSLIB_PATH/bgzip\
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf \
    $HTSLIB_PATH/tabix \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz

# Annotation with gnomAD
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar annotate \
    -v $GNOMAD_FILE \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz > $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf
    $HTSLIB_PATH/bgzip\
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf  
    $HTSLIB_PATH/tabix \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf.gz

## Deep Variant
# Variant effect annotation
java -Xmx4096m -jar $SNPEFF_PATH/snpEff.jar ann \
    -v -noLog -noStats -noLof GRCh37.75 \
    $VCF_PATH/deepvariant_filtered_norm.vcf.gz > $VCF_PATH/deepvariant_filtered_norm_eff.vcf.gz
    $HTSLIB_PATH/bgzip \
    $VCF_PATH/deepvariant_filtered_norm_eff.vcf
    $HTSLIB_PATH/tabix \
    $VCF_PATH/deepvariant_filtered_norm_eff.vcf.gz

# Annotation with dbSNP
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar annotate \
    -v -id $DBSNP_FILE \
    $VCF_PATH/deepvariant_filtered_norm_eff.vcf.gz > $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp.vcf
    $HTSLIB_PATH/bgzip\
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp.vcf
    $HTSLIB_PATH/tabix \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp.vcf.gz

# Annotation with dbNSFP
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar dbnsfp \
    -v -m -db $DBNSFP_FILE \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp.vcf.gz >     $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf
    $HTSLIB_PATH/bgzip\
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
    $HTSLIB_PATH/tabix \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
    
# Annotation with gnomAD
java -Xmx4096m -jar $SNPEFF_PATH/SnpSift.jar annotate \
    -v $GNOMAD_FILE \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz > $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf.gz  
    $HTSLIB_PATH/bgzip\
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf
    $HTSLIB_PATH/tabix \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp_gnomad.vcf.gz

# Remove intermediate files
rm $VCF_PATH/haplotypecaller_filtered_norm_eff.vcf.gz \
    $VCF_PATH/haplotypecaller_filtered_norm_eff.vcf.gz.tbi \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp.vcf.gz \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp.vcf.gz.tbi \
    $VCF_PATH/freebayes_filtered_norm_eff.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff.vcf.gz.tbi \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp.vcf.gz.tbi \
    $VCF_PATH/deepvariant_filtered_norm_eff.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff.vcf.gz.tbi \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp.vcf.gz.tbi
