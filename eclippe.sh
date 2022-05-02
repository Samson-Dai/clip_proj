#!/bin/sh
#SBATCH -N 2
#SBATCH -p RM-shared
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node=128

INPUT=""
TEMP=""
OUTPUT=""

_usage="
Usage: $0 [OPTIONS]
Options:
  [ -i INPUT_DIR ],          Required. Absolute path to input files directory.
  [ -o OUTPUT_DIR ],         Optional. Absolute path to output files directory. Defualt as INPUT_DIR.
  [ -t TEMP_DIR ],           Optional. Absolute path to temporary files directory. Defualt as INPUT_DIR.
  [ -h ],                    Help manuals.
"

usage() {                                 # Function: Print a help message.
  echo "$_usage" 1>&2 
}

exit_abnormal() {                         # Function: Exit with error.
  usage
  exit 1
}

while getopts ":i:o:t:h" options; do         # use silent error checking;
                                          
  case "${options}" in                    
    i)                                    # If the option is i,
      INPUT=${OPTARG}                     # set $INPUT to specified value.
      ;;
    o)                                    # If the option is o,
      OUTPUT=${OPTARG}                    # set $OUTPUT to specified value.
      ;;
    t)                                    # If the option is t,
      TEMP=${OPTARG}                      # set $TEMP to specified value.
      ;;
    h)                                    # If the option is h,
      usage
      exit 0                      		    # print help manuals.
      ;;
    :)                                    # If expected argument omitted:
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal                       # Exit abnormally.
      ;;
    *)									  # If unknown (any other) option:
      echo "Error: Unknow option provided."                                    
      exit_abnormal                       # Exit abnormally.
      ;;
  esac
done


if [ "$INPUT" = "" ]; then                 # If $INPUT is an empty string,
  echo "Error: must provide INPUT_DIR."
  exit_abnormal                       
fi

if [ "$TEMP" = "" ]; then                  # If $TEMP is an empty string,
  TEMP=$INPUT                        
fi

if [ "$OUTPUT" = "" ]; then                # If $OUTPUT is an empty string,
  OUTPUT=$INPUT                       
fi

################
################
#Pipeline starts
################
################

########################
#Define input file paths
########################
ECLIP_RAW_1_1="${INPUT}/rep1.r1.fq"
ECLIP_RAW_1_2="${INPUT}/rep1.r2.fq"
ECLIP_RAW_2_1="${INPUT}/rep2.r1.fq"
ECLIP_RAW_2_2="${INPUT}/rep2.r2.fq"
REF_GENOME="${INPUT}/ref_genome.fasta"
GENOME_GTF="${INPUT}/genome_info.gtf"
REPBASE="{INPUT}/homo_sapiens_repbase.fasta"

#############
#Trim adapter
#############
##Rep1
#Round 1
TRIMMED1_1_1="${TEMP}/rep1.r1.fqTr.fq"
TRIMMED1_1_2="${TEMP}/rep1.r2.fqTr.fq"

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
-o "$TRIMMED1_1_1" \
-p "$TRIMMED1_1_2" \
"$ECLIP_RAW_1_1" \
"$ECLIP_RAW_1_2"

#Round 2
TRIMMED2_1_1="${TEMP}/rep1.r1.fqTrTr.fq"
TRIMMED2_1_2="${TEMP}/rep1.r2.fqTrTr.fq"

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
-o "$TRIMMED2_1_1" \
-p "$TRIMMED2_1_2" \
"$TRIMMED1_1_1" \
"$TRIMMED1_1_2"

##Rep2
#Round 1
TRIMMED1_2_1="${TEMP}/rep2.r1.fqTr.fq"
TRIMMED1_2_2="${TEMP}/rep2.r2.fqTr.fq"

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
-o "$TRIMMED1_2_1" \
-p "$TRIMMED1_2_2" \
"$ECLIP_RAW_2_1" \
"$ECLIP_RAW_2_2"

#Round 2
TRIMMED2_2_1="${TEMP}/rep2.r1.fqTrTr.fq"
TRIMMED2_2_2="${TEMP}/rep2.r2.fqTrTr.fq"

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
-o "$TRIMMED2_2_1" \
-p "$TRIMMED2_2_2" \
"$TRIMMED1_2_1" \
"$TRIMMED1_2_2"

exit 0                                     # Exit normally.