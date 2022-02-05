#!/bin/bash

# Set resources path
RESOURCES_PATH=/home/user/resources
mkdir -p $ RESOURCES_PATH
CWD=`pwd`

# Download reference genome
cd $RESOURCES_PATH
mkdir hs37d5
cd hs37d5
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
cd $CWD

# Download exome capture kit coordinates
cd $RESOURCES_PATH
mkdir panel
cd panel
wget --no-check-certificate https://figshare.com/ndownloader/files/33961505
mv 33961505 Agilent_SureSelect_All_Exon_V2.bed.gz
gunzip Agilent_SureSelect_All_Exon_V2.bed.gz
cd $CWD

# Download and configure dbSNP, dbSNFP and gnomAD
cd $RESOURCES_PATH

mkdir dbSNP
cd dbSNP
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
gunzip 00-All.vcf.gz
cd ..

mkdir dbNSFP
cd dbNSFP
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv2.9.3.zip
unzip dbNSFPv2.9.3.zip
cd ..

mkdir gnomAD
cd gnomAD
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi

CWD=`pwd`
