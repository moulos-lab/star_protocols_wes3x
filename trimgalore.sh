#!/bin/bash

HOME_PATH=/PATH/TO/ANALYSIS/DIRECTORY
FASTQ_PATH=$HOME_PATH/fastq
TRIMGALORE_COMMAND=$TRIMGALORE_PATH/trim_galore
CUTADAPT_COMMAND=$CUTADAPT_PATH/cutadapt
TRIMGALORE_OUTPUT=$HOME_PATH/fastq_qual
CORES=4

if [ ! -d $TRIMGALORE_OUTPUT ]
then
    mkdir -p $TRIMGALORE_OUTPUT
fi

for FILE in $FASTQ_PATH/*_1.fastq.gz
do
    BASE=`basename $FILE | sed s/_1\.fastq\.gz//`
    echo "Processing $BASE"
    mkdir -p $TRIMGALORE_OUTPUT
    F1=$FASTQ_PATH/$BASE"_1.fastq.gz"
    F1=$FASTQ_PATH/$BASE"_2.fastq.gz"
    $TRIMGALORE_COMMAND \
        --quality 30 \
        --length 50 \
        --output_dir $TRIMGALORE_OUTPUT/$BASE \
        --path_to_cutadapt $CUTADAPT_COMMAND \
        --cores 4 \
        --paired \
        --fastqc \
        --trim-n $F1 $F2
            
    mv $TRIMGALORE_OUTPUT/$BASE"_1_val_1.fq.gz" \
        $TRIMGALORE_OUTPUT/$BASE"_1.fq.gz"
    mv $TRIMGALORE_OUTPUT/$BASE"_2_val_2.fq.gz" \
        $TRIMGALORE_OUTPUT/$BASE"_2.fq.gz"
    mv $TRIMGALORE_OUTPUT/$BASE"_1_val_1_fastqc.html" \
        $TRIMGALORE_OUTPUT/$BASE"_1_fastqc.html"
    mv $TRIMGALORE_OUTPUT/$BASE"_1_val_1_fastqc.zip" \
        $TRIMGALORE_OUTPUT/$BASE"_1_fastqc.zip"
    mv $TRIMGALORE_OUTPUT/$BASE"_2_val_2_fastqc.html" \
        $TRIMGALORE_OUTPUT/$BASE"_2_fastqc.html"
    mv $TRIMGALORE_OUTPUT/$BASE"_2_val_2_fastqc.zip" \
        $TRIMGALORE_OUTPUT/$BASE"_2_fastqc.zip"
done
