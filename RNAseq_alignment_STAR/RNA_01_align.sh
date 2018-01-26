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

## CREATE INDEX

echo 
echo "###############################################################"
echo "## Creating STAR index call (this only needs to be run once)..."
echo "###############################################################"
echo 

read -p "#INPUT REQUIRED: please define path to genome models GTF file (default $GENE_MODELS_GTF): " GENE_MODELS_GTF_temp
if [ -e "$GENE_MODELS_GTF_temp" ]
then
  GENE_MODELS_GTF=$GENE_MODELS_GTF_temp
else
  echo "#ERROR: cannot find \"$GENE_MODELS_GTF_temp\"... using default gene model GTF $GENE_MODELS_GTF"
  echo
fi

read -p "#INPUT REQUIRED: please define path to genome fasta file (default $GENOME_FASTA): " GENOME_FASTA_temp
if [ -e "$GENOME_FASTA_temp" ]
then
  GENOME_FASTA=$GENOME_FASTA_temp
else
  echo "#ERROR: cannot find '$GENOME_FASTA_temp'... using default genome fasta $GENOME_FASTA"
  echo
fi

read -p "#INPUT REQUIRED: please define path to STAR index folder (default $STAR_INDEX): " STAR_INDEX_temp
if [ -d "$STAR_INDEX_temp" ]
then
  STAR_INDEX=$STAR_INDEX_temp
else
  echo "#ERROR: cannot find '$STAR_INDEX_temp'... using default STAR index folder $STAR_INDEX"
  echo
fi

echo 
echo $STAR_BINARY --runMode genomeGenerate --genomeDir $STAR_INDEX --sjdbGTFfile $GENE_MODELS_GTF --sjdbOverhang 50 --genomeFastaFiles $GENOME_FASTA --runThreadN 4
echo 

## ALIGN RNA USING STAR

echo
echo "######################################"
echo "### Generating STAR alignment call ###"
echo "######################################"
echo

read -p "#INPUT REQUIRED: please define path to fastq file (default $FASTQ_FILE): " FASTQ_FILE_temp
if [ -e "$FASTQ_FILE_temp" ]
then
  FASTQ_FILE=$FASTQ_FILE_temp
else
  echo "#ERROR: cannot find '$FASTQ_FILE_temp'... using default STAR index folder $FASTQ_FILE"
  echo
fi

read -p "#INPUT REQUIRED: please define path to output folder for alignments (default $OUTPUT_DIR): " OUTPUT_DIR_temp
if [ -d "$OUTPUT_DIR_temp" ]
then
  OUTPUT_DIR=$OUTPUT_DIR_temp
else
  echo "#ERROR: cannot find '$OUTPUT_DIR'... using default STAR index folder $OUTPUT_DIR"
  echo
fi

echo
echo $STAR_BINARY --genomeDir $STAR_INDEX --sjdbFileChrStartEnd $STAR_INDEX --readFilesIn $FASTQ_FILE --runThreadN 8 --outFileNamePrefix --outFileNamePrefix $OUTPUT_DIR/$SAMPLE.aligned. --genomeLoad LoadAndRemove --alignIntronMax 500000 --alignMatesGapMax 500000 --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --outFilterMismatchNoverLmax 0.3 --sjdbOverhang 50 --chimSegmentMin 15 --chimScoreMin 1 --chimScoreJunctionNonGTAG 0 --chimJunctionOverhangMin 15 --readFilesCommand zcat --outStd SAM
echo
