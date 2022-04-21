#!/bin/sh

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#This is a reference to set up the environment for the pipeline.
#Only need to run necessary steps before the first time of usage.
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

###################
#Link the ocean dir
###################
ln -s /ocean/projects/mcb180074p/sdai clip_proj

###############
#Setup anaconda
###############
module load anaconda3/2020.11
mv $HOME/.conda $PROJECT/.conda
ln -s $PROJECT/.conda $HOME/.conda

###################
#Setup anaconda env
###################
#python=2.7
conda create --name py2 python=2.7
#install packages
pip install --user pysam==0.15.4

#######################################
#Download data
#######################################
#download eclip TARDBP K562 from ENCODE
cd clip_proj
xargs -L 1 curl -O -J -L < files.txt

#rename the sequencing files
mv ENCFF734UEC.fastq.gz rep1.IP.umi.r1.fq.gz
mv ENCFF147JYD.fastq.gz rep1.IP.umi.r2.fq.gz
mv ENCFF661TYX.fastq.gz rep2.IP.umi.r1.fq.gz
mv ENCFF218BOC.fastq.gz rep2.IP.umi.r2.fq.gz

#Download genome reference: https://www.encodeproject.org/files/ENCFF159KBI/
#Download male genome reference: https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/
#Download RepBase file homo_sapiens_repbase_fixed_v2.fasta from Canvas
