
# Functional analysis of RNA binding sites using eCLIP PE data

[CLIP-seq](https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/hits-clip-clip-seq-ptb-seq.html) is a crosslinking and immunoprecipitation (CLIP) sequencing technique that is used to identify RNA binding sites of RNA binding proteins (RBPs). CLIP-seq is frequently used to understand the protein-RNA interactions as well as its functional downstream effects. [eCLIP](https://www.illumina.com/science/sequencing-method-explorer/kits-and-arrays/eclip.html) data is a recent and improved version of CLIP-seq that can identify binding sites of RBPs in vivo. Compared to earlier CLIP-seq methods such as iCLIP, eCLIP tends to have higher efficiency and quality of library production as iCLIP methods often result in high duplication rates and low library complexity [(Nostrand, 2016)](https://www.nature.com/articles/nmeth.3810) .

This pipeline aims to use paired-end eCLIP data replicates to identify peak sequence regions where RBPs are bound and perform functional analysis on the RNA regions shown to interact with the RBP of interest. Steps taken to process raw peCLIP sequence data were adopted from the existing eCLIP pipeline of the Yeo lab [(Blue, 2022)](https://pubmed.ncbi.nlm.nih.gov/35322209/). The raw eCLIP reads (available as fastq files on [ENCODE](https://www.encodeproject.org/eclip/)) were processed and aligned against the reference genome to create .bam files. The .bam files were then used to call peaks and perform functional analysis on the peak regions. This pipeline is designed to perform in an HPC environment (PSC bridges-2).

## Technology Stack
* **Language**: Shell script
* **Computing Platform**: HPC (PSC Bridge2) using SLURM workload manager

## Pipeline

The pipeline takes fastq files from [ENCODE](https://www.encodeproject.org/eclip/) , the human hg38 reference genome and .gtf from the [Ensembl genome browser](https://useast.ensembl.org/index.html).  It uses [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove adapters from the reads. Then it maps the reads to a human specific version of
[RepBase](https://www.girinst.org/repbase/) used to remove repetitive elements, helps control for spurious artifacts from rRNA (&other) repetitive reads. The unique reads are aligned against the reference genome using [STAR](https://github.com/alexdobin/STAR) aligner. It then removes PCR duplicates by barcodecollapsepe.py. The .bam files are merged using [STAR](https://github.com/alexdobin/STAR) before being fed to [PEAKachu](https://github.com/tariks/peakachu) for peak-calling. [Homer](http://homer.ucsd.edu/homer/ngs/peaks.html) will then take the list of peaks (initial_peaks.csv) to annotate the peaks with its Entrez gene ID as well as identify binding motifs. At last, [RCAS](https://academic.oup.com/nar/article/45/10/e91/3038237) takes the annotated list of peaks to perform functional analysis of the transcript and generates a .png file containing the bar graph representing the properties of the transcripts that bind to the protein of interest. A .csv file with GO annotations and gene enrichment analysis results are generated as well.


<div align="center"><img src="https://github.com/Samson-Dai/clip_project/blob/main/test_results/workflow.png"></div>
<div align="center"><b>Pipeline workflow</b></div>

## Input files and naming conventions
Before running the pipeline, required input files should be prepared, renamed following the naming conventions and placed in the INPUT directory. Required input files and naming conventions are listed as follows, missing input files or failure to follow the naming conventions will lead to execution fault for the pipeline:

- **eCLIP sequencing data** naming as: 
```
	rep1.r1.fq
	rep1.r2.fq
	rep2.r1.fq
	rep2.r2.fq
```
The paired-end eCLIP sequencing data can be downloaded from the [ENCORE](https://www.encodeproject.org/encore-matrix/?type=Experiment&status=released&internal_tags=ENCORE) database. The data should contain 2 replicates and each replicate should contain 2 reads. Hence, in total, there should be 4 different eCLIP sequencing files.

- **RepBase repetitive RNA sequences** naming as: 
```
	homo_sapiens_repbase.fasta
```
The RepBase file contains the sets of common repeat elements for several species and it's used to remove the repetitive elements. The sample  homo sapiens RepBase file can be downloaded from the Canvas page.

- **Reference genome .fasta file** naming as: 
```
	ref_genome.fasta
```
The reference genome is used as a reference for the non-repeat reads to map to. We use [GRCh38](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/) reference genome as a sample genome. 

- **Reference genome .gtf information file** naming as: 
```
	genome_info.gtf
```
The reference genome .gtf file contains information about gene structures and helps to map non-repeat reads to the reference genome. We use [GRCh38 gtf](https://www.encodeproject.org/files/ENCFF159KBI/) reference genome .gtf as a sample .gtf. 

- **Control .bam data** naming as: 
```
	ctrl.bam
```
 We use [ENCFF948OYU.bam](https://www.encodeproject.org/files/ENCFF948OYU/@@download/ENCFF948OYU.bam) as our control .bam data. 

After downloading and renaming all the required input files, an input directory should be created to hold all the input files. Create a sample input directory via the following command:
```
	mkdir sample_input_dir
```

After moving all input files into the sample input directory, the tree structure of the sample input directory is as follows:
```
├── sample_input_dir
│   ├── rep1.r1.fq
│   ├── rep1.r2.fq
│   ├── rep2.r1.fq
│   ├── rep2.r2.fq
│   ├── homo_sapiens_repbase.fasta
│   ├── ref_genome.fasta
│   └── genome_info.gtf
```

## Environment setup
There're two bash files that can help to set up the environment and download data in the PSC bridges2. 
The first one is `shell_script/pre_setup.sh`, which can be run before the initial use of the pipeline.  It contains commands to install all the necessary packages as well as its dependencies. This bash script will also configure anaconda and create conda environments for each package to run properly. It also includes a guideline to downloading and renaming datasets such as raw eCLIP reads and RepBase. 
The second one is `shell_scripts/set_up.sh`, which can be run everytime a new bridges2 session starts. It will load all necessary modules and activate the anaconda environment.

## Usage
The pipeline consists of 2 parts. 
### Part1: eClip data processing and peak calling (PSC Bridges2)
The first part of the pipeline contains all steps for eClip data processing to peak calling and motif identification. This part is supposed to be run on PSC Bridges2.
Before running the pipeline, please make sure all required input files are renamed following the naming conventions and placed in the INPUT directory.  It's recommended to provide an output directory to hold the output file and a temp direcory to hold all the intermediate files. Create a sample output and temp directories via the following command:
```
	mkdir sample_output_dir
	mkdir sample_temp_dir
```
The pipeline can be run by directly calling the execution file. Please note that all directory paths are **absolute paths**: 
```
./eclippe [Options]

Options:
  [ -i INPUT_DIR ],          Required. Absolute path to input files directory.
  [ -o OUTPUT_DIR ],         Optional. Absolute path to output files directory. Default as INPUT_DIR.
  [ -t TEMP_DIR ],           Optional. Absolute path to temporary files directory. Default as INPUT_DIR.
  [ -h ],                    Help manuals.
```
Running the sample pipeline on PSC bridges2 requires a submission script that contains the commands to run the execution. The submission script looks like this:
```
#!/bin/sh
#SBATCH -p RM-shared
#SBATCH -t 5:00:00

./eclippe -i sample_input_di -o sample_output_dir -t sample_temp_dir
```
**Note: The user must personalize the submission script before run the pipeline.**
Then we can submit the submission script to slurm via the following command:
```
sbatch sample_submission.sh
```
This part will produce 3 output files:

-`pe_clip.fq.merged.bam`  A merged .bam file for 2 replicated to call peak on.
-`initial_peaks.bed` A .bed initial peak file for [RCAS](https://academic.oup.com/nar/article/45/10/e91/3038237) to perform functional analysis on.
-`homer_peaks_annot.txt` The annotated peaks.

### Part2: Functional analysis using  RCAS(local machine)
The user can run the RCAS.R on the local machine to conduct functional analysis. The user should place the peaks .bed file, UCSC annotation .gtf file or ENSEMBL .gtf file into one directory. The directory should have the following structure::
```
├── local_dir
│   ├── RCAS.R
│   ├── initial_peaks.bed
│   ├── hg38_UCSC.gtf
│   └── hg38_ENSEMBL.gtf
```
The RCAS.R script can be run locally by using the following command:
```
Rscript RCAS.R [Options]
Options:
  [ -b INPUT_BED ],      Required. Absolute path to input .bed file directory.
  [ -g INPUT_GTF ],      Required. Absolute path to input .gtf file directory.

```

The final output, including a .png file containing the bar graph representing the properties of the transcripts that bind to the protein of interest and also a .csv file including the enriched GO terms and p values, will be generated in the output directory.

<div align="center"><img src="https://github.com/Samson-Dai/clip_project/blob/main/test_results/RCAS_result_ENSEMBL.png"></div>
<div align="center"><img src="https://github.com/Samson-Dai/clip_project/blob/main/test_results/GO_result_screenshot.png"></div>
<div align="center"><b>Sample output</b></div>

## Team
* [Songcheng Dai](https://github.com/Samson-Dai/clip_project)
*  Tianyu Xia
*  Yejie Yun
*  Shuhao Zhang

## Reference
* Van Nostrand, E. L., Pratt, G. A., Shishkin, A. A., Gelboin-Burkhart, C., Fang, M. Y., Sundararaman, B., ... & Yeo, G. W. (2016). Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP). Nature methods, 13(6), 508-514.

  
* Blue, S. M., Yee, B. A., Pratt, G. A., Mueller, J. R., Park, S. S., Shishkin, A. A., ... & Yeo, G. W. (2022). Transcriptome-wide identification of RNA-binding protein binding sites using seCLIP-seq. Nature Protocols, 1-43.
