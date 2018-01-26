#

# Author: Naveed Ishaque (n.ishaqaue@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Naveed Ishaque

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

7library("DiffBind")
tamoxifen<-dba(sampleSheet="config.csv", 
                minOverlap=1,
                config=data.frame(RunParallel=TRUE, reportInit="DBA"),
                skipLines=0, 
                bAddCallerConsensus=FALSE, 
                bRemoveM=TRUE, 
                bRemoveRandom=TRUE, 
                bCorPlot=FALSE, 
                attributes=c(DBA_ID),
                bLowerScoreBetter=TRUE,
                peakCaller="raw") 

pdf("CES_H3K27ac_panabinostat_enhancer_panobinostat.peaks.pdf", width = 10, height = 10)
plot(tamoxifen)
dev.off()

pdf("CES_H3K27ac_panabinostat_enhancer_panobinostat.reads.DBA_SCORE_READS_FOLD.bScale_control.pdf", width = 10, height = 10)
tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_READS_FOLD, bScaleControl=TRUE)
dev.off()

save.image(file = "diffbind_config_CES_control.RData")

############################################################################################################################################

tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION, minMembers=2, block=DBA_TISSUE)

save.image(file = "diffbind_config_CES_control.RData")

load("diffbind_config_CES_control.RData")

#tamoxifen = dba.analyze(tamoxifen,bSubControl=FALSE,bParallel=TRUE,bFullLibrarySize=TRUE)
#tamoxifen = dba.analyze(tamoxifen,bSubControl=FALSE,bParallel=TRUE,bFullLibrarySize=FALSE)
#tamoxifen = dba.analyze(tamoxifen,bSubControl=TRUE,bParallel=TRUE,bFullLibrarySize=TRUE)
tamoxifen = dba.analyze(tamoxifen,bSubControl=TRUE,bParallel=TRUE,bFullLibrarySize=FALSE,method=DBA_ALL_METHODS)

############################################################################################################################################

save.image(file = "diffbind_config_CES_control.RData")

tamoxifen.DB = dba.report(tamoxifen, contrast=1, file="CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp.24T_VS_24C.normalized_F.DESEQ2_BLOCK.FDR_0.01", ext="csv", th=0.01, bUsePval=FALSE, bCounts=TRUE, bNormalized=FALSE,method=DBA_DESEQ2_BLOCK)
tamoxifen.DB = dba.report(tamoxifen, contrast=1, file="CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp.24T_VS_24C.normalized_F.DESEQ2_BLOCK.FDR_1.00", ext="csv", th=1, bUsePval=FALSE, bCounts=TRUE, bNormalized=FALSE,method=DBA_DESEQ2_BLOCK)
tamoxifen.DB = dba.report(tamoxifen, contrast=1, file="CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp.24T_VS_24C.normalized_T.DESEQ2_BLOCK.FDR_0.01", ext="csv", th=0.01, bUsePval=FALSE, bCounts=TRUE, bNormalized=TRUE,method=DBA_DESEQ2_BLOCK)
tamoxifen.DB = dba.report(tamoxifen, contrast=1, file="CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp.24T_VS_24C.normalized_T.DESEQ2_BLOCK.FDR_1.00", ext="csv", th=1, bUsePval=FALSE, bCounts=TRUE, bNormalized=TRUE,method=DBA_DESEQ2_BLOCK)

tamoxifen.DB = dba.report(tamoxifen, contrast=1, file="CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp.24T_VS_24C.normalized_F.EDGER_BLOCK.FDR_0.01", ext="csv", th=0.01, bUsePval=FALSE, bCounts=TRUE, bNormalized=FALSE,method=DBA_EDGER_BLOCK)
tamoxifen.DB = dba.report(tamoxifen, contrast=1, file="CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp.24T_VS_24C.normalized_F.EDGER_BLOCK.FDR_1.00", ext="csv", th=1, bUsePval=FALSE, bCounts=TRUE, bNormalized=FALSE,method=DBA_EDGER_BLOCK)
tamoxifen.DB = dba.report(tamoxifen, contrast=1, file="CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp.24T_VS_24C.normalized_T.EDGER_BLOCK.FDR_0.01", ext="csv", th=0.01, bUsePval=FALSE, bCounts=TRUE, bNormalized=TRUE,method=DBA_EDGER_BLOCK)
tamoxifen.DB = dba.report(tamoxifen, contrast=1, file="CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp.24T_VS_24C.normalized_T.EDGER_BLOCK.FDR_1.00", ext="csv", th=1, bUsePval=FALSE, bCounts=TRUE, bNormalized=TRUE,method=DBA_EDGER_BLOCK)


############################################################################################################################################

#pdf("CES_H3K27ac_panabinostat_enhancer_panobinostat.BOX_plot.FDR_0.01.4_fold.pdf", width = 10, height = 10)
#pvals = dba.plotBox(tamoxifen, bUsePval=FALSE, th=0.01, fold=4)
#dev.off()

#save.image(file = "diffbind_config_CES_control.RData")

############################################################################################################################################

pdf("CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp_sites.24T_VS_24C.FDR_0.01.EDGER_BLOCK.pdf", width = 10, height = 10)
dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE, bUsePval=FALSE, th=0.01, method=DBA_EDGER_BLOCK)
dev.off()
pdf("CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp_sites.24T_VS_24C.FDR_0.01.DESEQ2_BLOCK.pdf", width = 10, height = 10)
dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE, bUsePval=FALSE, th=0.01, method=DBA_DESEQ2_BLOCK)
dev.off()

pdf("CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp_sites.24T_VS_24C.FDR_0.05.EDGER_BLOCK.pdf", width = 10, height = 10)
dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE, bUsePval=FALSE, th=0.05,method=DBA_EDGER_BLOCK)
dev.off()
pdf("CES_H3K27ac_panabinostat_enhancer_panobinostat.diff_exp_sites.24T_VS_24C.FDR_0.05.DESEQ2_BLOCK.pdf", width = 10, height = 10)
dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE, bUsePval=FALSE, th=0.05,method=DBA_DESEQ2_BLOCK)
dev.off()

save.image(file = "diffbind_config_CES_control.RData")
