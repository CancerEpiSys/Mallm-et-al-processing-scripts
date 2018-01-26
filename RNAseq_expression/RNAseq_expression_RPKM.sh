#!/bin/bash

# Author: Naveed Ishaque (n.ishaqaue@dkfz.de). Based on original work by Barbara Hutter (b.hutter@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Naveed Ishaque

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

clear

## SOFTWARE STACK - these binaries were used for the analysis

SAMTOOLS_BINARY="samtools-0.1.17-r973:277"
SAM_VIEW_OPTS="-bu -q 1"
MBUFFER_BINARY="mbuffer-0.9.9"
COVERAGEBED="coverageBed-2.16.2"
COVBED_OPTS="-split"

## PATHS and FILES - this has to change with your files structure

REFEXONS=/path/to/RefSeq_Nov15_2011_from_annovar_Exons_plain.bed.gz
REFGENES=/path/to/RefSeq_Nov15_2011_from_annovar_Genes_plain.bed.gz

bam=/path/to/MySample.bam
read -p "#INPUT REQUIRED: please define path to RNAseq BAM file (default $bam): " bam

bamstem=`echo $bam | perl -ne '{$_=~ s/\.bam$//; print $_;}'`
outfile=$bamstem.Exons_RPKM.bed 
generpkm=$bamstem.Exons_RPKM.bed

# WORKFLOW
echo
echo "$SAMTOOLS_BINARY view $SAM_VIEW_OPTS $bam | $COVERAGEBED $COVBED_OPTS -abam stdin -b $REFEXONS > ${outfile}"
echo "echo -e \"#chrom\tchromStart\tchromEnd\tname\texonNr\tstrand\tlength\treads\tbases_covered\tlength\tcoverage\tRPKM\" > ${outfile}.tmp"
echo "awk '{a+=\$(NF-3);}END{print a}' ${outfile} | perl RPKM.pl ${outfile} - | sort -k1,1d -k2,2n"
echo "mv ${outfile}.tmp ${outfile}"
echo "perl countcoverageBedRPKM.pl ${outfile} <(bgzip -c -d $REFGENES) $total | awk 'NR==1; NR > 1 {print $0 | "sort -k1,1d -k2,2n"}' > $generpkm"
echo
