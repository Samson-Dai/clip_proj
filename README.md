# Functional analysis of RNA binding sites using eCLIP PE data

<!---
[eCLIP](https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/eclip.html) is a sequencing technique developed by Illumina, which sequences mRNA strands that are bound by RNA-binding proteins. It involves in-vivo cross-linking of RNA and protein. The protein is then immunoprecipitated and the attached mRNA is sequenced. This sequencing data is available as fastq files on [ENCODE](https://www.encodeproject.org/eclip/). This project aims to create a pipeline for extracting the mRNA sequences from these fastq files. The project comes with an API that aims to enable end-users to extract sequences from their own eCLIP experiments and compare it against the database of eCLIP experiments on ENCODE. 
--->

## Technology Stack

* **Language**: Shell script
* **Computing Platform**: HPC (PSC Bridge2) using SLURM workload manager

## Pipeline

<!---
The pipeline takes the fastq files from [ENCODE](https://www.encodeproject.org/eclip/) and the human hg38 reference genome from the [Ensembl genome browser](https://useast.ensembl.org/index.html). It uses [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to check the quality of the fastq files. Then it uses [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove adapters from the reads. It then aligns the trimmed reads against the reference genome using [STAR](https://github.com/alexdobin/STAR) aligner. The bam files from the aligner are indexed using [Samtools](http://www.htslib.org/) before being fed to [PureCLIP](https://github.com/skrakau/PureCLIP) for peak-calling. [Bedtools](https://bedtools.readthedocs.io/en/latest/) is used to extract raw sequences from this output. The output is stored in a BLAST database, which allows us to [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) a sequence against this database of sequences.
--->

<!---
A detailed user manual for this API can be found [here](https://docs.google.com/document/d/1oN8uVp0X6dJDNgbZuuwEoN9pRx6a43hob3GOVItfcQU/edit?usp=sharing).
--->

<div align="center"><img src="https://github.com/Samson-Dai/clip_project/blob/main/workflow.png" width="600" height="400"></div>
<div align="center"><b>Pipeline workflow</b></div>

## Input files and naming conventions
Before running the pipeline, required input files should be prepared, renamed following the naming conventions and placed in the INPUT directory. Required input files and naming conventions are listed as follows, missing input files or failure to follow the naming conventions will lead to execution fault for the pipeline:

- **eCLIP sequencing data** naming as: 
```js
	req1.r1.fq
	req1.r2.fq
	req2.r1.fq
	req2.r2.fq
```
The Pair-end eCLIP sequencing data can be downloaded from the [ENCORE](https://www.encodeproject.org/encore-matrix/?type=Experiment&status=released&internal_tags=ENCORE) database. The data should contain 2 replicates and each replicate should contain 2 reads. Hence, in total, there should be 4 different eCLIP sequencing files.

- **RepBase repetitive RNA sequences** naming as: 
```js
	homo_sapiens_repbase.fasta
```
The RepBase file contains the sets of common repeat elements for several species and it's used to remove the repetitive elements. The sample  homo sapiens RepBase file can be downloaded from the Canvas page.

- **Reference genome .fasta file** naming as: 
```js
	ref_genome.fasta
```
The reference genome is used as a reference for the non-repeat reads to map to. We use [GRCh38](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/) reference genome as a sample genome. 

- **Reference genome .gtf information file** naming as: 
```js
	genome_info.gtf
```
The reference genome .gtf file contains information about gene structures and helps to map non-repeat reads to the reference genome. We use [GRCh38 gtf](https://www.encodeproject.org/files/ENCFF159KBI/) reference genome .gtf as a sample .gtf. 

After downloading and renaming all the required input files, a input directory should be created to hold all the input files. Create a sample input directory via the following command:
```js
	mkdir sample_input_dir
```

After moving all input files into the sample input directory, the tree structure of the sample input directory is as follows:
```
├── sample_input_dir
│   ├── req1.r1.fq
│   ├── req1.r2.fq
│   ├── req2.r1.fq
│   ├── req2.r2.fq
│   ├── homo_sapiens_repbase.fasta
│   ├── ref_genome.fasta
│   └── genome_info.gtf
```

## Usage


## Summary

<!---
We developed a pipeline that could extract mRNA sequences bound by RBP from the reads data available at ENCODE. The API for this pipeline could be used by end-users to extract sequences from their experiments and compare it against the sequences available within ENCODE.
--->

## Team

* [Songcheng Dai](https://github.com/Samson-Dai/clip_project)
*  Tianyu Xia
*  Yejie Yun
*  Shuhao Zhang