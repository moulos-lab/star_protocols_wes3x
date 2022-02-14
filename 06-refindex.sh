#!/bin/bash

cd $RESOURCES_PATH/hs37d5
$BWA_PATH/bwa index hs37d5.fa
$SAMTOOLS_PATH/samtools faidx hs37d5.fa
$SAMTOOLS_PATH/samtools dict hs37d5.fa > hs37d5.dict
cd $CWD
