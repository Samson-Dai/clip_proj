#!/bin/sh

#############
#load package
#############
module load cutadapt/2.10
module load STAR/2.7.6a
module load samtools/1.13.0
module load anaconda3/2020.11
module load bedtools/2.29.2
module load MEME-suite/5.4.1
module load FastQC

###############
#setup anaconda
###############
mv $HOME/.conda $PROJECT/.conda
ln -s $PROJECT/.conda $HOME/.conda

#############
#Trim adapter
#############
##Rep1
#Round 1
cutadapt \
-f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 1 \
--quality-cutoff 6 \
-m 18 \
-a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-g CTTCCGATCTACAAGTT \
-g CTTCCGATCTTGGTCCT \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o rep1.IP.umi.r1.fqTr.fq \
-p rep1.IP.umi.r2.fqTr.fq \
rep1.IP.umi.r1.fq \
rep1.IP.umi.r2.fq

fastqc -t 2 --extract -k 7 rep1.IP.umi.r1.fqTr.fq -o quality
fastqc -t 2 --extract -k 7 rep1.IP.umi.r2.fqTr.fq -o quality

#Round 2
cutadapt \
-f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 5 \
--quality-cutoff 6 \
-m 18 \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o rep1.IP.umi.r1.fqTrTr.fq \
-p rep1.IP.umi.r2.fqTrTr.fq \
rep1.IP.umi.r1.fqTr.fq \
rep1.IP.umi.r2.fqTr.fq

fastqc -t 2 --extract -k 7 rep1.IP.umi.r1.fqTrTr.fq -o quality
fastqc -t 2 --extract -k 7 rep1.IP.umi.r2.fqTrTr.fq -o quality

##Rep2
#Round 1
cutadapt \
-f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 1 \
--quality-cutoff 6 \
-m 18 \
-a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-g CTTCCGATCTACAAGTT \
-g CTTCCGATCTTGGTCCT \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o rep2.IP.umi.r1.fqTr.fq \
-p rep2.IP.umi.r2.fqTr.fq \
rep2.IP.umi.r1.fq \
rep2.IP.umi.r2.fq

#Round 2
cutadapt \
-f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 5 \
--quality-cutoff 6 \
-m 18 \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o rep2.IP.umi.r1.fqTrTr.fq \
-p rep2.IP.umi.r2.fqTrTr.fq \
rep2.IP.umi.r1.fqTr.fq \
rep2.IP.umi.r2.fqTr.fq

###########################
#Remove repetitive elements
###########################
#Generate the repeat index
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeSAindexNbases 9 \
--genomeDir homo_sapiens_repbase_v2 \
--genomeFastaFiles homo_sapiens_repbase_fixed_v2.fasta

#Rep1
#Run alignment
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir homo_sapiens_repbase_v2 \
--genomeLoad NoSharedMemory \
--alignEndsType EndToEnd \
--outSAMunmapped Within \
--outFilterMultimapNmax 30 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix rep1.IP.umi.r1.fqTrTr.sorted.STAR \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outBAMcompression 10 \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo \
--outSAMattributes All \
--outSAMmode Full \
--outStd Log \
--readFilesIn rep1.IP.umi.r1.fqTrTr.fq rep1.IP.umi.r2.fqTrTr.fq

#Rename files
mv rep1.IP.umi.r1.fqTrTr.sorted.STARAligned.out.bam rep1.IP.umi.r1.fq.repeat-mapped.bam
mv rep1.IP.umi.r1.fqTrTr.sorted.STARUnmapped.out.mate1 rep1.IP.umi.r1.fq.repeat-unmapped.fq
mv rep1.IP.umi.r1.fqTrTr.sorted.STARUnmapped.out.mate2 rep1.IP.umi.r2.fq.repeat-unmapped.fq

#Rep2
#Run alignment
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir homo_sapiens_repbase_v2 \
--genomeLoad NoSharedMemory \
--alignEndsType EndToEnd \
--outSAMunmapped Within \
--outFilterMultimapNmax 30 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix rep2.IP.umi.r1.fqTrTr.sorted.STAR \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outBAMcompression 10 \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo \
--outSAMattributes All \
--outSAMmode Full \
--outStd Log \
--readFilesIn rep2.IP.umi.r1.fqTrTr.fq rep2.IP.umi.r2.fqTrTr.fq

#Rename files
mv rep2.IP.umi.r1.fqTrTr.sorted.STARAligned.out.bam rep2.IP.umi.r1.fq.repeat-mapped.bam
mv rep2.IP.umi.r1.fqTrTr.sorted.STARUnmapped.out.mate1 rep2.IP.umi.r1.fq.repeat-unmapped.fq
mv rep2.IP.umi.r1.fqTrTr.sorted.STARUnmapped.out.mate2 rep2.IP.umi.r2.fq.repeat-unmapped.fq

##############################
#Map non-repeats to the genome
##############################
## Download genome reference: https://www.encodeproject.org/files/ENCFF159KBI/
## Download male genome reference: https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/
#Generate index 
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir genome_STARindex \
--sjdbGTFfile ENCFF159KBI.gtf \
--genomeFastaFiles GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--sjdbOverhang 99

##Rep1
#Alignment
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir genome_STARindex \
--genomeLoad NoSharedMemory \
--readFilesIn \
rep1.IP.umi.r1.fq.repeat-unmapped.fq \
rep1.IP.umi.r2.fq.repeat-unmapped.fq \
--outSAMunmapped Within \
--outFilterMultimapNmax 1 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix rep1.IP.umi.r1.fq.genome-mapped \
--outSAMattributes All \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo \
--outStd Log \
--alignEndsType EndToEnd \
--outBAMcompression 10 \
--outSAMmode Full

#Re-name BAM: rename genome-mapped outputs
mv rep1.IP.umi.r1.fq.genome-mappedAligned.out.bam rep1.IP.umi.r1.fq.genome-mapped.bam

##Rep2
#Alignment
#!/bin/sh
STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir genome_STARindex \
--genomeLoad NoSharedMemory \
--readFilesIn \
rep2.IP.umi.r1.fq.repeat-unmapped.fq \
rep2.IP.umi.r2.fq.repeat-unmapped.fq \
--outSAMunmapped Within \
--outFilterMultimapNmax 1 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix rep2.IP.umi.r1.fq.genome-mapped \
--outSAMattributes All \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--outSAMattrRGline ID:foo \
--outStd Log \
--alignEndsType EndToEnd \
--outBAMcompression 10 \
--outSAMmode Full

#Re-name BAM: rename genome-mapped outputs
mv rep2.IP.umi.r1.fq.genome-mappedAligned.out.bam rep2.IP.umi.r1.fq.genome-mapped.bam

###############
#Name sort BAM
###############
#Sort output from STAR by name to ensure read pairs are adjacent.
##Rep1
samtools \
sort \
-n \
-o rep1.IP.umi.r1.fq.genome-mappedSo.bam \
rep1.IP.umi.r1.fq.genome-mapped.bam

##Rep2
samtools \
sort \
-n \
-o rep2.IP.umi.r1.fq.genome-mappedSo.bam \
rep2.IP.umi.r1.fq.genome-mapped.bam

######################
#Remove PCR duplicates
######################
#use python/2.7
#install packages
conda create --name py2 python=2.7
conda activate py2
pip install --user pysam==0.15.4

##Rep1
python barcodecollapsepe.py \
-o rep1.IP.umi.r1.fq.genome-mappedSo.rmDup.bam \
-m rep1.IP.umi.r1.fq.genome-mappedSo.rmDup.metrics \
-b rep1.IP.umi.r1.fq.genome-mappedSo.bam

##Rep2
python barcodecollapsepe.py \
-o rep2.IP.umi.r1.fq.genome-mappedSo.rmDup.bam \
-m rep2.IP.umi.r1.fq.genome-mappedSo.rmDup.metrics \
-b rep2.IP.umi.r1.fq.genome-mappedSo.bam

########################
#Sort resulting bam file
########################
##Rep1
samtools sort -o rep1.IP.umi.r1.fq.genome-mappedSo.rmDupSo.bam rep1.IP.umi.r1.fq.genome-mappedSo.rmDup.bam
##Rep2
samtools sort -o rep2.IP.umi.r1.fq.genome-mappedSo.rmDupSo.bam rep2.IP.umi.r1.fq.genome-mappedSo.rmDup.bam

#################
#Merge replicates
#################
samtools merge pe_clip.fq.genomemappedSo.rmDupSo.merged.bam rep1.IP.umi.r1.fq.genome-mappedSo.rmDupSo.bam rep2.IP.umi.r1.fq.genome-mappedSo.rmDupSo.bam

############
#Call peaks
############
#dowload test control data
get https://www.encodeproject.org/files/ENCFF948OYU/@@download/ENCFF948OYU.bam
mv ENCFF948OYU.bam ctrl.bam

#install peakachu
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n peakachu peakachu python=3
conda activate peakachu

#run peakachu
peakachu adaptive -M 200 -m 0.0 -f 2.0 -Q 0.05 -c ctrl.bam -t pe_clip.fq.genomemappedSo.rmDupSo.merged.bam

#convert peakachu peaks to homer peak file
awk -F'\t' -v OFS='\t' 'NR ==1 {print "ID", $0; next} {print (NR-1), $0} ' peaks.txt > peaks_homer_1.txt
cut -f6,7,8,9 --complement peaks_homer_1.txt > homer_peaks.txt

#############
#Find motifs 
#############
#download hg38 genome and load module
module load homer
module load bedtools 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

#converst peak file to bed file
awk 'NR>1{
    if ($2 < $3) {
        print $1, $2, $3, "clip_peak_"NR-1,$8,$4;
        }
    else {
        print $1, $3, $2, "clip_peak_"NR-1,$8,$4;
        }
    }' initial_peaks.csv | tr ' ' '\t' > initial_peaks.bed

#extract genomic DNA of peaks
bedtools getfasta -fi hg38.fa -bed initial_peaks.bed > initial_peaks_seq.fa

#find motifs 
findMotifs.pl initial_peaks_seq.fa fasta homer_output/ 

###############
#Annotate peaks
###############
#download genome .gtf file
wget -O hg38_UCSC.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
gunzip hg38_UCSC.gtf.gz

#annotate peaks
annotatePeaks.pl homer_peaks.txt hg38.fa -gtf hg38_UCSC.gtf > homer_peaks_annot.txt

#########################
#Gene enrichment analysis
#########################
#add RCAS rscript

