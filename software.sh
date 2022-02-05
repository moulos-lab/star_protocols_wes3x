#!/bin/bash

# Main installation path
INSTALL_PATH=/home/user/tools
mkdir -p $INSTALL_PATH
CWD=`pwd`

# Download and configure FastQC
LINK=https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
unzip $ARCHIVE
export FASTQC_PATH=$INSTALL_PATH/FastQC
rm $ARCHIVE
cd $CWD

# Install MultiQC
pip install multiqc
export MULTIQC_PATH=/home/user/.local/bin

# Install cutadapt
pip install –upgrade cutadapt
export CUTADAPT_PATH=/home/user/.local/bin

# Download and configure TrimGalore
LINK=https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz -O TrimGalore_v0.6.7.tar.gz
cd $INSTALL_PATH
wget $LINK
ARCHIVE=TrimGalore_v0.6.7.tar.gz
tar -xvf $ARCHIVE
export TRIMGALORE_PATH=$INSTALL_PATH/ TrimGalore-0.6.7
rm $ARCHIVE
cd $CWD

# Download and install bwa
LINK=https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
tar -xvf $ARCHIVE
export BWA_PATH=$INSTALL_PATH/bwa-0.7.17
rm $ARCHIVE
cd $BWA_PATH
make
cd $CWD

# Download and configure GATK
LINK=https://github.com/broadinstitute/gatk/releases/download/4.2.4.1/gatk-4.2.4.1.zip
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
unzip $ARCHIVE
export GATK_PATH=$INSTALL_PATH/gatk-4.2.4.1
rm $ARCHIVE
cd $CWD

# Download and install FreeBayes
LINK=https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
mkdir freebayes-1.3.6
mv $ARCHIVE ./freebayes-1.3.6/
cd freebayes-1.3.6
gunzip $ARCHIVE
chmod +x freebayes-1.3.6-linux-amd64-static
mv freebayes-1.3.6-linux-amd64-static freebayes
export FREEBAYES_PATH=$INSTALL_PATH/freebayes-1.3.6
cd $CWD

# Download, configure and install Docker
sudo apt remove docker docker-engine docker.io containerd runc

sudo apt update

sudo apt install ca-certificates curl gnupg lsb-release

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt update

sudo apt install docker-ce docker-ce-cli containerd.io

sudo usermod -aG docker ${USER}

BIN_VERSION="1.3.0"
sudo docker pull google/deepvariant:"${BIN_VERSION}"

# Download and configure SnpEff
LINK=https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
unzip $ARCHIVE
export SNPEFF_PATH=$INSTALL_PATH/snpEff
rm $ARCHIVE
cd $SNPEFF_PATH
chmod +x snpEff.jar SnpSift.jar
cd $CWD

# Download and install samtools
LINK=https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
tar -xvf $ARCHIVE
export SAMTOOLS_PATH=$INSTALL_PATH/samtools-1.14
rm $ARCHIVE
cd $SAMTOOLS_PATH
make
cd $CWD

# Download and install bcftools
LINK=https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
tar -xvf $ARCHIVE
export BCFTOOLS_PATH=$INSTALL_PATH/bcftools-1.14
rm $ARCHIVE
cd $BCFTOOLS_PATH
make
cd $CWD

# Download and install htslib
LINK=https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
tar -xvf $ARCHIVE
export HTSLIB_PATH=$INSTALL_PATH/htslib-1.14
rm $ARCHIVE
cd $HTSLIB_PATH
make
cd $CWD

# Download and install BEDTools
LINK=https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
tar -xvf $ARCHIVE
rm $ARCHIVE
export BEDTOOLS_PATH=$INSTALL_PATH/bedtools2/bin
cd BEDTOOLS_PATH/..
make
cd $CWD

# Download and configure UCSC tools
cd $INSTALL_PATH
mkdir ucsc_tools
cd ucsc_tools
rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ./
export UCSCTOOLS_PATH=$INSTALL_PATH/ucsc_tools
cd $CWD

# Download and configure vcflib
LINK=https://github.com/vcflib/vcflib/releases/download/v1.0.1/vcflib-1.0.1-src.tar.gz
cd $INSTALL_PATH
wget $LINK
ARCHIVE=`basename $LINK`
tar -xvf $ARCHIVE
rm $ARCHIVE
mv vcflib-1.0.1-src vcflib-1.0.1
export VCFLIB_PATH=$INSTALL_PATH/vcflib-1.0.1/bin
cd VCFLIB_PATH/..
make
cd $CWD

# Download and configure GLnexus
cd $INSTALL_PATH
mkdir GLnexus
cd GLnexus
wget https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.4.1/glnexus_cli
chmod +x glnexus_cli
cd ..
export GLNEXUS_PATH=$INSTALL_PATH/GLnexus
cd $CWD
