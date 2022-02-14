#!/bin/bash

# Setup paths
HOME_PATH=/home/user/analysis
FASTQ_PATH=$HOME_PATH/fastq
mkdir -p $ FASTQ_PATH
cd $CWD

# Download raw data from 1000 genomes project
cd $FASTQ_PATH

# HG00119
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099967/SRR099967_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099967/SRR099967_2.fastq.gz
mv SRR099967_1.fastq.gz HG00119_1.fastq.gz
mv SRR099967_2.fastq.gz HG00119_2.fastq.gz

# HG00133
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099969/SRR099969_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099969/SRR099969_2.fastq.gz
mv SRR099969_1.fastq.gz HG00133_1.fastq.gz
mv SRR099969_2.fastq.gz HG00133_2.fastq.gz

# HG00145
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_2.fastq.gz
mv SRR099957_1.fastq.gz HG00145_1.fastq.gz
mv SRR099957_2.fastq.gz HG00145_2.fastq.gz

# HG00239
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099958/SRR099958_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099958/SRR099958_2.fastq.gz
cd $DATA_PATH

# HG00119
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099967/SRR099967_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099967/SRR099967_2.fastq.gz
mv SRR099967_1.fastq.gz HG00119_1.fastq.gz
mv SRR099967_2.fastq.gz HG00119_2.fastq.gz

# HG00133
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099969/SRR099969_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099969/SRR099969_2.fastq.gz
mv SRR099969_1.fastq.gz HG00133_1.fastq.gz
mv SRR099969_2.fastq.gz HG00133_2.fastq.gz

# HG00145
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_2.fastq.gz
mv SRR099957_1.fastq.gz HG00145_1.fastq.gz
mv SRR099957_2.fastq.gz HG00145_2.fastq.gz

# HG00239
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099958/SRR099958_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099958/SRR099958_2.fastq.gz
mv SRR099958_1.fastq.gz HG00239_1.fastq.gz
mv SRR099958_2.fastq.gz HG00239_2.fastq.gz

# HG00258
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099954/SRR099954_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099954/SRR099954_2.fastq.gz
mv SRR099954_1.fastq.gz HG00258_1.fastq.gz
mv SRR099954_2.fastq.gz HG00258_2.fastq.gz

# HG00265
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099968/SRR099968_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099968/SRR099968_2.fastq.gz
mv SRR099968_1.fastq.gz HG00265_1.fastq.gz
mv SRR099968_1.fastq.gz HG00265_2.fastq.gz

cd $CWD
