#

# Author: Sandra Koser (s.koser@dkfz.de).

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Naveed Ishaque

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS 
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

library( "DESeq" )
samplesheet = "/path/to/example_samplesheet.txt"
reulsts_dir="/path/to/results"
df = read.table(samplesheet, header=FALSE )
cds = newCountDataSetFromHTSeqCount( df, reulsts_dir=)
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
# cds = estimateDispersions( cds )
cds = estimateDispersions( cds, method="blind", sharingMode="fit-only" )
str( fitInfo( cds ) )
plotDispEsts( cds )
head( fData( cds ) )

# write.csv( counts( cds, normalized=TRUE ), file="expressionValues.csv" )

res_1_1 = nbinomTest( cds, "control1", "cancer" ) 

res_1_1Sig = res_1_1[ res_1_1$padj < 0.1, ] 

head( res_1_1Sig[ order( res_1_1Sig$pval), ]) #most signifcantly diferentially expressed genes
head( res_1_1Sig[ order( res_1_1Sig$foldChange, -res_1_1Sig$baseMean ), ] ) #most strongly down-regulated of the signifcant genes
head( res_1_1Sig[ order( -res_1_1Sig$foldChange, -res_1_1Sig$baseMean ), ] ) #most strongly up-regulated ones

write.csv( res_1_1, file="results_diff_c1_pAll.csv" )
write.csv( res_1_1Sig[ order(res_1_1Sig$pval), ], file="differentiallyExpr_c1_pAll.csv" )
