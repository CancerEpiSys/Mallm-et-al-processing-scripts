###########################################################################
###          Rsript: Plot color bar for quality scores (FRiP)           ###
###              (author: Lara Klett    date: 2017-04-04)               ###
###########################################################################

########
writeLines("\n***\nCreated on 2017-07-24\n\nUpdated on 2018-02-23\n\nauthor: Lara Klett\n***\n"); writeLines("Hello User!"); writeLines("\nPlot color bar for quality scores (FRiP) ...\n")
########


##################################################
###  Plot color bar for quality scores (FRiP)  ###
##################################################
library("gplots", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")

folder="your_path"
setwd(folder)

info=read.csv("sample-names.txt", sep="\t", header=F)
colnames(info)=c("individual","ID","health_status", "mut_status")
rownames(info)=info[,1]

frip=read.csv("FRiPscores.csv")


frip$FRiP_score[frip$FRiP_score>0.15]=0.15

mycol=colorRampPalette(c("white","#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(150)
col_breaks = seq(from=0, to=0.15, by=0.001)

par(mfrow=c(1,1))

info2=info[c(1,3:37),]

ma=frip[match(rownames(info2),frip[,1]),2]


pdf("Result_ATAC_QC_FRiP_colorbar_MuratColors.pdf", height=3.5, width=32)
par(mar=c(5,5,5,5))
heatmap.2(t(cbind(ma,ma)), 
          lhei=rep(0.1,2),
          col=mycol, 
          dendrogram="none", 
          Rowv=F, Colv=F,
          labCol=info2[,1],
          labRow="", 
          density.info="none",
          trace="none", 
          breaks=col_breaks, 
          sepcol="black", 
          colsep=c(0:dim(info2)[1]), 
          rowsep=c(0,2),
          srtCol=45,
          cexCol=2.5,
          key.par = list(cex=1),keysize=0.5)
dev.off()




