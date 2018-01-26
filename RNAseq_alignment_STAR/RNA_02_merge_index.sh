#!/bin/bash

# Author: Naveed Ishaque (n.ishaqaue@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Naveed Ishaque

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

clear

## SOFTWARE STACK - these binaries were used for the analysis

STAR_BINARY="STAR_2.3.0e_r291"
SAMTOOLS_BINARY="samtools-0.1.17-r973:277"
MBUFFER_BINARY="mbuffer-0.9.9"

## PATHS and FILES - this has to change with your files structure

DATABASE_PATH="/path/to/references/and/indexes"
GENE_MODELS_GTF="$DATABASE_PATH/gencode.v13.annotation_nochr.gtf"
GENOME_FASTA="$DATABASE_PATH/hs37d5.fa"
STAR_INDEX="$DATABASE_PATH/STAR_2.3.0e_Gencode13"
THREADS=8
SAMPLE="SAMPLE_XYZ"
FASTQ_FILE="$SAMPLE.fastq_R1.fq.gz"
OUTPUT_DIR="/path/to/results/$SAMPLE/alignments"

## SORT all SAM files

echo
echo "######################"
echo "### SORT SAM FILES ###"
echo "######################"
echo

read -p "#INPUT REQUIRED: please define BAM file preffix (default $SAMPLE): " SAMPLE_temp
if [ ${#SAMPLE_temp} -gt 5 ]
then
  SAMPLE=$SAMPLE_temp
else
  echo "#ERROR: bam preffix should be longer thatn 5 chars '$SAMPLE_temp'... using default bam preffix $SAMPLE"
  echo
fi

read -p "#INPUT REQUIRED: please define path to output folder for alignments (default $OUTPUT_DIR): " OUTPUT_DIR_temp
if [ -d $OUTPUT_DIR_temp ]
then
  OUTPUT_DIR=$OUTPUT_DIR_temp
else
  echo "#ERROR: cannot find '$OUTPUT_DIR_temp'... using default STAR index folder $OUTPUT_DIR"
  echo
fi

SAM_FILES=(`ls $OUTPUT_DIR/*sam`)

iterator=0;

for file in ${SAM_FILES[@]}
do 
  iterator=$[$iterator+1]
  echo "# SAM file $iterator: $file"
done

BAM_FILE_LIST=""

if [ $iterator -lt 1 ]
then
  echo "#ERROR: no sam files found in $OUTPUT_DIR... exiting!"
  exit
else
  for file in ${SAM_FILES[@]}
  do 
    echo
    echo "$SAMTOOLS_BINARY view -Sb $file | $SAMTOOLS_BINARY sort - $file.bam"
    echo "$SAMTOOLS_BINARY index $file.bam"
    BAM_FILE_LIST="$BAM_FILE_LIST $file.bam"
    echo
  done
fi

if [ $iterator -gt 0 ]
then
  echo
  echo "$SAMTOOLS_BINARY merge $SAMPLE.merged.bam $BAM_FILE_LIST"
  echo "$SAMTOOLS_BINARY index $SAMPLE.merged.bam"
  echo
fi
