#

# Author: Naveed Ishaque (n.ishaqaue@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Naveed Ishaque

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

library("DiffBind")
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

png("CES_H3K27ac_enhancer_interaction.peaks.png", width = 500, height = 500)
plot(tamoxifen)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_READS.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_READS)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_READS_FOLD.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_READS_FOLD)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_READS_MINUS.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_READS_MINUS)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_RPKM.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_RPKM)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_RPKM_FOLD.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_RPKM_FOLD)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_TMM_READS_FULL.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_TMM_READS_FULL)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_TMM_READS_EFFECTIVE.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_TMM_READS_EFFECTIVE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_TMM_MINUS_FULL.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_TMM_MINUS_FULL)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_TMM_MINUS_EFFECTIVE.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_TMM_MINUS_EFFECTIVE)
dev.off()

#bScaleControl

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_READS.bScale_control.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_READS, bScaleControl=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_READS_FOLD.bScale_control.png", width = 500, height = 500)
tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_READS_FOLD, bScaleControl=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_READS_MINUS.bScale_control.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_READS_MINUS, bScaleControl=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_RPKM.bScale_control.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_RPKM, bScaleControl=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_RPKM_FOLD.bScale_control.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_RPKM_FOLD, bScaleControl=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_TMM_READS_FULL.bScale_control.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_TMM_READS_FULL, bScaleControl=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_TMM_READS_EFFECTIVE.bScale_control.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_TMM_READS_EFFECTIVE, bScaleControl=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_TMM_MINUS_FULL.bScale_control.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_TMM_MINUS_FULL, bScaleControl=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.reads.DBA_SCORE_TMM_MINUS_EFFECTIVE.bScale_control.png", width = 500, height = 500)
#tamoxifen = dba.count (tamoxifen, minOverlap=1, fragmentSize=180, bParallel=tamoxifen$config$RunParallel, bCorPlot=TRUE, score=DBA_SCORE_TMM_MINUS_EFFECTIVE, bScaleControl=TRUE)
dev.off()

save.image(file = "diffbind_config_CES_control.RData")

############################################################################################################################################

tamoxifen = dba.contrast(tamoxifen, categories=DBA_CONDITION,minMembers=2)

save.image(file = "diffbind_config_CES_control.RData")

#png("CES_H3K27ac_enhancer_interaction.DBA_SCORE_TMM_READS_FULL.diff_exp_cor.subcont_F.fulllib_T.png", width = 500, height = 500)
#tamoxifen = dba.analyze(tamoxifen,bSubControl=FALSE,bParallel=TRUE,bFullLibrarySize=TRUE)
#dev.off()

png("CES_H3K27ac_enhancer_interaction.DBA_SCORE_TMM_READS_FULL.diff_exp_cor.subcont_F.fulllib_F.png", width = 500, height = 500)
tamoxifen = dba.analyze(tamoxifen,bSubControl=FALSE,bParallel=TRUE,bFullLibrarySize=FALSE)
dev.off()

#png("CES_H3K27ac_enhancer_interaction.DBA_SCORE_TMM_READS_FULL.diff_exp_cor.subcont_T.fulllib_T.png", width = 500, height = 500)
#tamoxifen = dba.analyze(tamoxifen,bSubControl=TRUE,bParallel=TRUE,bFullLibrarySize=TRUE)
#dev.off()

#png("CES_H3K27ac_enhancer_interaction.DBA_SCORE_TMM_READS_FULL.diff_exp_cor.subtcont_T.fulllib_F.png", width = 500, height = 500)
#tamoxifen = dba.analyze(tamoxifen,bSubControl=TRUE,bParallel=TRUE,bFullLibrarySize=FALSE)
#dev.off()

############################################################################################################################################

save.image(file = "diffbind_config_CES_control.RData")
tamoxifen.DB = dba.report(tamoxifen,file="CES_H3K27ac_enhancer_interaction.diff_exp.normalized_F.FDR_1", ext="csv", th=1, bUsePval=FALSE, bCounts=TRUE, bNormalized=FALSE)
save.image(file = "diffbind_config_CES_control.RData")

tamoxifen.DB = dba.report(tamoxifen,file="CES_H3K27ac_enhancer_interaction.diff_exp.normalized_T.FDR_1", ext="csv", th=1, bUsePval=FALSE, bCounts=TRUE, bNormalized=TRUE)
save.image(file = "diffbind_config_CES_control.RData")

############################################################################################################################################

png("CES_H3K27ac_enhancer_interaction.MA_plot.FDR_0.01.4_fold.normalized_T.png", width = 500, height = 500)
dba.plotMA(tamoxifen, bUsePval=FALSE, th=0.01, fold=4, bNormalized=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.MA_plot.FDR_0.01.3_fold.normalized_T.png", width = 500, height = 500)
dba.plotMA(tamoxifen, bUsePval=FALSE, th=0.01, fold=3, bNormalized=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.MA_plot.FDR_0.01.2_fold.normalized_T.png", width = 500, height = 500)
dba.plotMA(tamoxifen, bUsePval=FALSE, th=0.01, fold=2, bNormalized=TRUE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.MA_plot.FDR_0.01.4_fold.normalized_F.png", width = 500, height = 500)
dba.plotMA(tamoxifen, bUsePval=FALSE, th=0.01, fold=4, bNormalized=FALSE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.MA_plot.FDR_0.01.3_fold.normalized_F.png", width = 500, height = 500)
dba.plotMA(tamoxifen, bUsePval=FALSE, th=0.01, fold=3, bNormalized=FALSE)
dev.off()

png("CES_H3K27ac_enhancer_interaction.MA_plot.FDR_0.01.2_fold.normalized_F.png", width = 500, height = 500)
dba.plotMA(tamoxifen, bUsePval=FALSE, th=0.01, fold=2, bNormalized=FALSE)
dev.off()

save.image(file = "diffbind_config_CES_control.RData")

############################################################################################################################################

png("CES_H3K27ac_enhancer_interaction.BOX_plot.FDR_0.01.4_fold.png", width = 500, height = 500)
pvals = dba.plotBox(tamoxifen, bUsePval=FALSE, th=0.01, fold=4)
dev.off()

png("CES_H3K27ac_enhancer_interaction.BOX_plot.FDR_0.01.3_fold.png", width = 500, height = 500)
pvals = dba.plotBox(tamoxifen, bUsePval=FALSE, th=0.01, fold=3)
dev.off()

png("CES_H3K27ac_enhancer_interaction.BOX_plot.FDR_0.01.2_fold.png", width = 500, height = 500)
pvals = dba.plotBox(tamoxifen, bUsePval=FALSE, th=0.01, fold=2)
dev.off()

save.image(file = "diffbind_config_CES_control.RData")

############################################################################################################################################

png("CES_H3K27ac_enhancer_interaction.diff_exp_sites.FDR_0.01.png", width = 500, height = 500)
dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE, bUsePval=FALSE, th=0.01)
dev.off()

save.image(file = "diffbind_config_CES_control.RData")
