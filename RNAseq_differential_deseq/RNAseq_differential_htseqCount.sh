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
HTSEQ_BINARY="htseq-count-0.6.0"
THREADS=8

GENE_ANNOTATION="/path/to/gencode.v17.annotation.gtf"

bam=/path/to/MySample.bam
read -p "#INPUT REQUIRED: please define path to RNAseq BAM file (default $bam): " bam

bamstem=`echo $bam | perl -ne '{$_=~ s/\.bam$//; print $_;}'`

#

echo
echo "$SAMTOOLS_BINARY sort -n -m 5G -@ $THREADS -l 0 -o $bam $bamstem.namesorted"
echo "$HTSEQ_BINARY -m intersection-nonempty -s reverse -a 1 -t exon -i gene_name -s yes $bamstem.namesorted.bam $GENE_ANNOTATION > $bamstem.countTable_rev_all.txt"
echo
