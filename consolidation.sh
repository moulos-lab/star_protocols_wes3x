#!/bin/bash

export VCF_PATH=$HOME_PATH/vcf

# 1
$BCFTOOLS_PATH/bcftools isec \
    --prefix 1 \
    --output-type z \
    --nfiles ~100 \
    --collapse none \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz

# 2
$BCFTOOLS_PATH/bcftools isec \
    --prefix 2 \
    --output-type z \
    --nfiles ~010 \
    --collapse none \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
 
# 3
$BCFTOOLS_PATH/bcftools isec \
    --prefix 3 \
    --output-type z \
    --nfiles ~001 \
    --collapse none \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
    
# 4
$BCFTOOLS_PATH/bcftools isec \
    --prefix 4 \
    --output-type z \
    --nfiles ~110 \
    --collapse none \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
    
# 5
$BCFTOOLS_PATH/bcftools isec \
    --prefix 5 \
    --output-type z \
    --nfiles ~011 \
    --collapse none \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
    
# 6
$BCFTOOLS_PATH/bcftools isec \
    --prefix 6 \
    --output-type z \
    --nfiles ~101 \
    --collapse none \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
    
# 7
$BCFTOOLS_PATH/bcftools isec \
    --prefix 7 \
    --output-type z \
    --nfiles ~111 \
    --collapse none \
    $VCF_PATH/haplotypecaller_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/freebayes_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz \
    $VCF_PATH/deepvariant_filtered_norm_eff_dbsnp_dbnsfp.vcf.gz
