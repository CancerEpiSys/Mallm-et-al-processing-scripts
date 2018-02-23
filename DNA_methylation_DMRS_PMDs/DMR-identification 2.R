options(max.print=200)
options(stringsAsFactors=FALSE)

library(GenomicRanges)
#library(Gviz)
library(rtracklayer)

library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)
library(caTools)
library(parallel)
library(exomeCopy)
library(RColorBrewer)
library(zoo)
library(ggplot2)
library(data.table)
library(bsseq)
library(DSS)
library(VennDiagram)


library(BSgenome.Hsapiens.UCSC.hg19)
setwd("~/Projects/PBS1/methylation/results-mem-1000gen/scripts")
#itrack <- IdeogramTrack(genome = gen, chromosome = chr)

load("publicly-available-data/ucsc.gap.gr.20140513.RData")
##################################################
#preparing tracks:
#chip-seq and rna-seq tracks
#loading epigenomics tools:
source("epigenomics-tools.git/classes.R")
source("epigenomics-tools.git/accessors.R")
source("epigenomics-tools.git/extractors.R")
source("epigenomics-tools.git/constructors.R")
source("epigenomics-tools.git/loaders.R")
source("epigenomics-tools.git/writers.R")
source("epigenomics-tools.git/modifiers.R")
source("epigenomics-tools.git/modifiers_merge.R")
source("epigenomics-tools.git/modifiers_bin.R")


library("BSgenome.Hsapiens.UCSC.hg19")
##########################################################################################################################################
explabel=read.table("../results/downstream-analysis/supp-fig-sample-summary/experiment-labels.txt")
#explabel=read.table("experiment-labels.txt")
exporder=unique(explabel[,2])

########################################################################################################




###################################################################################################################
#Preparing methylation data for BSseq analysis:

# sample-related info
explabel=read.table("../results/downstream-analysis/supp-fig-sample-summary/experiment-labels.txt")
#explabel=read.table("experiment-labels.txt")
exporder=unique(explabel[,2])


#getting Beta values and perform t-test for averages:
l2=list.files("../results/wgbs",pattern="merged.dupmarked.RData.Beta.RData",full.names = TRUE,recursive = TRUE)
l2=l2[grep("wgbs/CES_",l2)]
l2=l2[grep("CES_GN",l2,invert=TRUE)]
l2=l2[grep("CES_WS-2",l2,invert=TRUE)]
l2=l2[grep("CES_YF-2",l2,invert=TRUE)]
l2=l2[grep("CES_VL-2",l2,invert=TRUE)]
allBetas=list()
for(i in 1:length(l2))
{
  print(i)
  load(l2[i])
  allBetas[[i]]=Beta
}
lname=basename(l2)
lname=sub("dm.combined.","",lname)
lname=sub("_complete_filtered_merged.dupmarked.RData.Beta.RData","",lname)

#getting methylation Coverage and beta values:
l2=list.files("../results/wgbs",pattern="_complete_filtered_merged.dupmarked.CG.all.txt",full.names = TRUE,recursive = TRUE)
l2=l2[grep("wgbs/CES_",l2)]
l2=l2[grep("CES_GN",l2,invert=TRUE)]
l2=l2[grep("CES_WS-2",l2,invert=TRUE)]
l2=l2[grep("CES_YF-2",l2,invert=TRUE)]
l2=l2[grep("CES_VL-2",l2,invert=TRUE)]
allBetas=list()
allCov=list()
allM=list()
for(i in 1:length(l2))
{
  print(i)
  a=read.table(pipe(paste0("cut -f5,6,7 ",l2[i])))
  a[a[,1]>0.25,2]=NA
  a[a[,1]>0.25,3]=NA
  cov1=a[seq(from=1,to=length(a[,1]),by=2),2]+a[seq(from=2,to=length(a[,1]),by=2),2]
  cov2=a[seq(from=1,to=length(a[,1]),by=2),3]+a[seq(from=2,to=length(a[,1]),by=2),3]
  allCov[[i]]=cov1+cov2
  allM[[i]]=cov1
}
save(allCov,allM,file="../results/downstream-analysis/dmr-tf-enrichment.Rdata")
lname=basename(l2)
lname=sub("CES_","",lname)
lname=sub("_complete_filtered_merged.dupmarked.CG.all.txt","",lname)
lname=sub("-1","",lname)
lab=explabel[match(lname,explabel[,1]),3]
sel=c(1:23)[!(c(1:23)%in%c(grep("-B",lname),grep("-T",lname)))]
allCov2=allCov[sel]
allM2=allM[sel]
lname2=lname[sel]
lab=explabel[match(lname2,explabel[,1]),3]
plab=data.frame(lab)
rownames(plab)=lname2
allCov3=lapply(allCov2,function(x)data.frame(x))
allM3=lapply(allM2,function(x)data.frame(x))
allM3mat=matrix(unlist(rbindlist(allM3)),nc=length(lname2))
allCov3mat=matrix(unlist(rbindlist(allCov3)),nc=length(lname2))
save(allM3mat,file="../results/downstream-analysis/allM3mat.Rdata")
save(allCov3mat,file="../results/downstream-analysis/allCov3mat.Rdata")

i=1
a=read.table(pipe(paste0("cut -f1,2 ",l2[i])))
a=a[seq(from=1,to=length(a[,1]),by=2),]
save(a,allM3mat,allCov3mat,lname2,plab,file="../results/downstream-analysis/bsseq-input.RData")


BStmp <- BSseq(chr = a[,1], pos = a[,2], M = allM3mat, Cov = allCov3mat, sampleNames = lname2,pData=plab)
#Coverage filtering:

BS.cov <- getCoverage(BStmp)

keepLoci.ex <- sum(rowSums(allCov3mat[, plab == "CLL"] >= 5) >= 4 & rowSums(allCov3mat[, plab == "H"] >= 5) >= 4,na.rm=TRUE)
keepLoci.ex <- which(rowSums(is.na(allCov3mat[, plab == "CLL"])) == 0 & rowSums(is.na(allCov3mat[, plab == "H"])) == 0 & rowSums(allCov3mat[, plab == "CLL"] >= 5) >= 4 & rowSums(allCov3mat[, plab == "H"] >= 5) >= 4)
length(keepLoci.ex)
allCov4mat=allCov3mat[keepLoci.ex,]
allM4mat=allM3mat[keepLoci.ex,]
a2=a[keepLoci.ex,]
colnames(allCov4mat)=lname2
colnames(allM4mat)=lname2

BSfilt <- BSseq(chr = a2[,1], pos = a2[,2], M = allM4mat, Cov = allCov4mat, sampleNames = lname2,pData=plab)
save(BSfilt,file="../results/downstream-analysis/BSfilt.RData")

BSsmooth <- BSmooth(BSfilt, ns=11 ,h=500 , maxGap=2000, mc.cores = 40, verbose = TRUE)

BS.cancer.ex.tstat <- BSmooth.tstat(BSfilt, group1 = c("CX", "FG", "IO","KL","PP","QU","RC","TH","WW","XP","ZV"), group2 = c("AY", "BK", "EM","VL","WS","YF"), estimate.var = "group2", local.correct = TRUE, verbose = TRUE,mc.cores=40)


# Using DDS package:
load("../results/downstream-analysis/dmr-identification/BSfilt.RData")
dmlTest <- DMLtest(BSfilt, group1 = c("CX", "FG", "IO","KL","PP","QU","RC","TH","WW","XP","ZV"), group2 = c("AY", "BK", "EM","VL","WS","YF"),smoothing=FALSE)
dmls <- callDML(dmlTest, p.threshold=0.05)
dmrs2 <- callDMR(dmlTest, delta=0.1, p.threshold=0.05)

load("../results/downstream-analysis/dmr-identification/dmlTest.RData")

load("../results/downstream-analysis/dmr-identification/dmrs2.RData")
####################################################################################################################
####################################################################################################################
####################################################################################################################
# loading UMRs:
# get all file names with UMR.
l=list.files("../results/wgbs/LMR",pattern="UMR.bedGraph",recursive="TRUE",full.names=TRUE)
l=l[grep("old-analysis",l,invert=TRUE)]
l=l[grep("PMD-DMV-list",l,invert=TRUE)]
# outlier
l=l[grep("CES_GN",l,invert=TRUE)]
# not interested in T cells
l=l[grep("-T",l,invert=TRUE)]
# SD-B is behaving different than other B cells.
l=l[grep("SD-B",l,invert=TRUE)]


umrlist=list()
for(i in 1:length(l))
{
  print(i)
  a=read.table(l[i])
  # filter out chrY
  a=a[!a[,1]%in%"chrY",]
  umrlist[[i]]=GRanges(seqnames=a[,1],IRanges(start=as.numeric(a[,2]),end=as.numeric(a[,3])))
  # filter out misidentified ranges that contains genomic gaps.
  umrlist[[i]]=umrlist[[i]][countOverlaps(umrlist[[i]],ucsc.gap.gr)==0]
}
samplenames=basename(l)
samplenames=sub("dm.combined.","",samplenames)
samplenames=sub("-1_complete_filtered_merged.dupmarked.RData.UMR.bedGraph","",samplenames)
samplenames=sub("-1_complete_filtered_merged.dupmarked.UMR.bedGraph","",samplenames)
samplenames=sub("CES_","",samplenames)

samplecond=explabel[match(samplenames,explabel[,1]),3]

allumr=GRanges()
for(i in 1:length(umrlist))
{
  allumr=c(allumr,umrlist[[i]])
}
nonredumr=reduce(allumr)
#by checking redundant umrs, we already end up with 35K out of 51k sites.
nonredumr2=nonredumr[countOverlaps(nonredumr,allumr)>1]

##########################################################################################################################################
#get PMD information
# to order samples based on the total length of their PMDs
# to prepare a summary PMD track to plot.
l= list.files("../results/downstream-analysis/PMD-vs-chip-seq/pmds/",full.names=TRUE,recursive=TRUE)
l=l[grep("shuffle",l,invert=TRUE)]
pmdgr=GRanges()
for(i in 1:length(l))
{
  pmd=read.table(l[i])
  a=GRanges(seqnames=pmd[,1],IRanges(start=pmd[,2],end=pmd[,3]))
  a=a[countOverlaps(a,ucsc.gap.gr)==0]
  
  pmdgr=c(pmdgr,a)
}
redpmdgr=pmdgr
pmdgr2=reduce(pmdgr)
a=countOverlaps(pmdgr2,pmdgr)
#we require that at least half of the cancer samples show PMD in that region.
pmdgr=pmdgr2[a>5]
pmd=data.frame(chr=seqnames(pmdgr),start=start(pmdgr),end=end(pmdgr))


#loading K27ac
data=read.csv("../../../ChIP-Seq/diffbind/H3K27ac/DBA_CES_H3K27ac_interaction.diff_exp.normalized_T.FDR_1.00.csv")
# Conc_WT refers to wild type and if fold change is high, it means h3k27ac is higher in controls.
#data=data[data[,1]%in%paste0("chr",c(1:22)),]
data=data[data[,1]%in%c(1:22),]
k27all=GRanges(seqnames=paste("chr",as.character(data[,1]),sep=""),IRanges(start=as.numeric(data[,2]),end=as.numeric(data[,3])))
mcols(k27all)=data[,6]
data2=data[data[,9]<0.05 & data[,7]>0,]
k27down=GRanges(seqnames=paste("chr",as.character(data2[,1]),sep=""),IRanges(start=as.numeric(data2[,2]),end=as.numeric(data2[,3])))
data2=data[data[,9]<0.05 & data[,7]<0,]
k27up=GRanges(seqnames=paste("chr",as.character(data2[,1]),sep=""),IRanges(start=as.numeric(data2[,2]),end=as.numeric(data2[,3])))
#k27up=GRanges(seqnames=paste(as.character(data2[,1]),sep=""),IRanges(start=as.numeric(data2[,2]),end=as.numeric(data2[,3])))
data2=data[data[,9]>0.5,]
k27back=GRanges(seqnames=paste("chr",as.character(data2[,1]),sep=""),IRanges(start=as.numeric(data2[,2]),end=as.numeric(data2[,3])))

data2=data[data[,9]<0.01 & data[,7]>0,]
write.table(data2[,c(1,2,3)],"../results/downstream-analysis/dmr-identification/Diffbind-H3K27ac-Down-Pval-0.01.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
data2=data[data[,9]<0.05 & data[,7]>0,]
write.table(data2[,c(1,2,3)],"../results/downstream-analysis/dmr-identification/Diffbind-H3K27ac-Down-Pval-0.05.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
data2=data[data[,9]<0.1 & data[,7]>0,]
write.table(data2[,c(1,2,3)],"../results/downstream-analysis/dmr-identification/Diffbind-H3K27ac-Down-Pval-0.1.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
data2=data[data[,9]<0.01 & data[,7]<0,]
write.table(data2[,c(1,2,3)],"../results/downstream-analysis/dmr-identification/Diffbind-H3K27ac-Up-Pval-0.01.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
data2=data[data[,9]<0.05 & data[,7]<0,]
write.table(data2[,c(1,2,3)],"../results/downstream-analysis/dmr-identification/Diffbind-H3K27ac-Up-Pval-0.05.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
data2=data[data[,9]<0.1 & data[,7]<0,]
write.table(data2[,c(1,2,3)],"../results/downstream-analysis/dmr-identification/Diffbind-H3K27ac-Up-Pval-0.1.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)



################################################################################################
################################################################################################
################################################################################################
# loading B cell maturation:
l=list.files("../results/wgbs",pattern="_complete_filtered_merged.dupmarked.RData.Beta.RData",recursive="TRUE",full.names=TRUE)
oakeslist=c("GCF","loMBC","NBC","sMGZ","hiMBC","intMBC")
l=l[grep(paste(oakeslist,collapse="|"),l)]
allBetas=list()
for(i in 1:length(l))
{
  print(i)
  load(l[i])
  allBetas[[i]]=data.frame(Beta)
}

lname=basename(l)
lname=sub("dm.combined.","",lname)
lname=sub("_complete_filtered_merged.dupmarked.RData.Beta.RData","",lname)

# loading cpg positions:
source("downstream-analysis/get_wgbs_position.R")
cpgPoslist=get_wgbs_position()
cpg.gr=cpgPoslist[[1]]
cpg.range=cpgPoslist[[2]]
# allBetas2=lapply(allBetas,function(x)data.frame(x))
# betamat=matrix(unlist(rbindlist(allBetas2)),nc=length(lname))
bmatgr=list()
for(i in 1:length(l))
{
  a=cpg.range
  mcols(a)=allBetas[[i]]
  bmatgr[[i]]=a
}

#######################################################################################
#######################################################################################
# loading Chip-Seq peaks:

#chip-seq peaks:
#l=list.files("../../../ChIP-Seq/results_per_pid/",pattern="-W200-G800-FDR1e-3-island.bed$",full.names=TRUE,recursive=TRUE)
chipmark=c("_H3K27ac_macs_peaks.bed","_H3K9ac_macs_peaks.bed","_H3K4me3_macs_peaks.bed","_H3K4me1_macs_peaks.bed","_H3K9me3-W200-G800-FDR1e-3-island.bed","_H3K27me3-W200-G800-FDR1e-3-island.bed","_H3K36me3-W200-G800-FDR1e-3-island.bed")
h3peakl=list()
for(k in 5:7)
{
  print(k)
l=list.files("../../../ChIP-Seq/results_per_pid/",pattern=paste0(chipmark[k],"$"),full.names=TRUE,recursive=TRUE)
nam=basename(l)
nam=gsub("CES_","",nam)
nam=gsub(chipmark[k],"",nam)
samplab=explabel[match(nam,explabel[,1]),3]
l2=l[samplab=="CLL"]
#l2=l
peaksgr=GRanges()
for(i in 1:length(l2))
{
  print(i)
  if(samplab[i]=="CLL")
  {
    peak=read.table(l2[i])
    if(k<=4)
    {
      peaksgr=c(peaksgr,GRanges(seqnames=paste0("chr",peak[,1]),IRanges(start=peak[,2],end=peak[,3])))
    }
    else
    {
      peaksgr=c(peaksgr,GRanges(seqnames=paste0(peak[,1]),IRanges(start=peak[,2],end=peak[,3])))
    }
  }
}
h3peaks=reduce(peaksgr)
h3peaks=h3peaks[countOverlaps(h3peaks,peaksgr)>=10]
h3peaks=h3peaks[as.character(seqnames(h3peaks))%in%paste0("chr",1:22)]
h3peakl[[k]]=h3peaks
}


#save(h3peakl,file="../results/downstream-analysis/dmr-identification/h3peaksl.RData")
load("../results/downstream-analysis/dmr-identification/h3peaksl.RData")
#is reserved for H3K27ac:
h3peaks=h3peakl[[1]]
###########################################################################################
###########################################################################################
#12078#
#########################################################################################
###########################################################################################
load("../results/downstream-analysis/dmr-identification/bsseq-input.RData")
dim(allCov3mat)
dim(allM3mat)
cllm=GRanges(seqnames=a[,1],IRanges(start=a[,2],end=a[,2]))
allCov4mat=allCov3mat
allCov4mat[allCov4mat<5]=NA
mcols(cllm)=allCov4mat

dmrna=getGRNAnumber(cllm,downdmr2)


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
# loading DMVs and UMRs:
l=list.files("../results/wgbs/PMD-DMV-list/defined-dmvs/",pattern=".bed$",full.names=TRUE,recursive=TRUE)
nam=basename(l)
nam=gsub("dm.combined.CES_","",nam)
nam=gsub("dm.combined.","",nam)
nam=gsub("_complete_filtered_merged.dupmarked.DMV.PMD.ARRAY.065-015Parameter..cpgwindow015thresDMV.bed","",nam)
nam=gsub("_complete_filtered_merged.dupmarked.DMV.PMD.ARRAY.065-015Parameter..cpgwindow065thresdmv-healthy-intersection.bed","",nam)
nam=gsub("_complete_filtered_merged.dupmarked.DMV.PMD.ARRAY.065-015Parameter..cpgwindow065thresdmv-patient-intersection.bed","",nam)
samplab=explabel[match(nam,explabel[,1]),3]
#l2=l[samplab=="CLL"]
#l2=l
dmv=GRanges()
for(i in 1:length(l))
{
  print(i)
  peak=read.table(l[i])
  dmv=c(dmv,GRanges(seqnames=paste0(peak[,1]),IRanges(start=peak[,2],end=peak[,3])))
}
dmv2=reduce(dmv)
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
#Integrating DMR with UMRs:
load("../results/downstream-analysis/dmr-identification/dmls.RData")
load("../results/downstream-analysis/dmr-identification/dmrs2.Rdata")

hist(dmrs2[,8])
sum(dmrs2[,8]>0)
df=data.frame(dmrs2[,])
dmr.gr=GRanges(seqnames=dmrs2[,1],IRanges(start=dmrs2[,2],end=dmrs2[,3]))
mcols(dmr.gr)=dmrs2[,c(4,5,6,7,8)]
outdmr.gr=GRanges(seqnames=dmrs2[,1],IRanges(start=dmrs2[,2],end=dmrs2[,3]))
mcols(outdmr.gr)=dmrs2[,c(4,5,6,7,8,9)]
dt=as.matrix(mcols(dmr.gr))
smoothScatter(dt[,3],dt[,5])

# calculating meth. difference between b cell maturation cell types:
bmatmean=lapply(bmatgr,function(x) countBinMean(x, dmr.gr))
# the order of the samples:
# NBC, GCF, loMBC, intMBC, sMGZ, hiMBC
averbmat= cbind(apply(cbind(bmatmean[[9]],bmatmean[[10]]),1,function(x)mean(x,na.rm=TRUE)),apply(cbind(bmatmean[[1]],bmatmean[[2]]),1,function(x)mean(x,na.rm=TRUE)),
      apply(cbind(bmatmean[[7]],bmatmean[[8]]),1,function(x)mean(x,na.rm=TRUE)),apply(cbind(bmatmean[[5]],bmatmean[[6]]),1,function(x)mean(x,na.rm=TRUE)),
      apply(cbind(bmatmean[[11]],bmatmean[[12]]),1,function(x)mean(x,na.rm=TRUE)),apply(cbind(bmatmean[[3]],bmatmean[[4]]),1,function(x)mean(x,na.rm=TRUE)))
bmatrange=apply(averbmat,1,function(x)diff(range(x,na.rm=TRUE))*sign(which.max(x)-which.min(x)))
# missing 211344 238296 270062:
bmatrange[[211344]]=NA
bmatrange[[238296]]=NA
bmatrange[[270062]]=NA
bmatrange=unlist(bmatrange)
bmatrange[!is.finite(bmatrange)]=NA

pdf(paste0("../results/downstream-analysis/dmr-identification/figures/",Sys.Date(),"-MI-diffmethylation-Bmaturation-vs-Cllcomp.pdf"))
smoothScatter(bmatrange,dt[,5],xlab="Diff. meth. during B cell maturation",ylab="Diff. methylation between CLL and B cells",xlim=c(-1,1),ylim=c(-1,1))
dev.off()
odmr.gr=dmr.gr
obmatrange=bmatrange

bmatrange=obmatrange
dmr.gr=odmr.gr
bmatrange=bmatrange[countOverlaps(dmr.gr,pmdgr)==0]
dmr.gr=dmr.gr[countOverlaps(dmr.gr,pmdgr)==0]
dt=as.matrix(mcols(dmr.gr))
dmr.gr=dmr.gr[dt[,1]>200 & (dt[,3]<0.3 |dt[,3]>0.7)]
bmatrange=bmatrange[dt[,1]>200 & (dt[,3]<0.3 |dt[,3]>0.7)]
dt=as.matrix(mcols(dmr.gr))

methvalues=cbind(bmatrange,dt[,5])
bmatrange2=bmatrange
bmatrange2[is.na(bmatrange)]=1
downdmr=dmr.gr[dt[,1]>200 & dt[,3]<0.3 & dt[,5]<(-0.3) & dt[,5]<(bmatrange2-0.2) & dt[,5]<(-0.3)]
selM=which(dt[,1]>200 & dt[,3]<0.3 & dt[,5]<(-0.3) & dt[,5]<(bmatrange2-0.2) & dt[,5]<(-0.3))
bmatrange2[is.na(bmatrange)]=(-1)
updmr=dmr.gr[dt[,1]>200 & dt[,3]>0.7 & dt[,5]>(0.3) & dt[,5]>(bmatrange2+0.2)& dt[,5]>0.3]
selM=c(selM,which(dt[,1]>200 & dt[,3]>0.7 & dt[,5]>(0.3) & dt[,5]>(bmatrange2+0.2)& dt[,5]>0.3))
methvalues2=methvalues[selM,]

pdf(paste0("../results/downstream-analysis/dmr-identification/figures/",Sys.Date(),"-MI-diffmethylation-Bmaturation-vs-Cllcomp-filtered.pdf"))
par(pty="s")
#smoothScatter(bmatrange,dt[,5],nbin=1000,col="deepskyblue",colramp=colorRampPalette(c("white", "deepskyblue")),xlab="Diff. meth. during B cell maturation",ylab="Diff. methylation between CLL and B cells",xlim=c(-1,1),ylim=c(-1,1)); 
#smoothScatter(methvalues2[,1],methvalues2[,2],nbin=1000,col="red",colramp = colorRampPalette(c(rgb(1, 1, 1, 0), rgb(1, 0, 0, 1)), alpha = TRUE),add=T);
plot(bmatrange,dt[,5],col = alpha("blue", 0.4),pch=16,cex=0.3,xlab="Diff. meth. during B cell maturation",ylab="Diff. methylation between CLL and B cells",xlim=c(-1,1),ylim=c(-1,1)); 
points(methvalues2[,1],methvalues2[,2],col = alpha("red", 0.4),pch=16,cex=0.3,xlab="Diff. meth. during B cell maturation",ylab="Diff. methylation between CLL and B cells",xlim=c(-1,1),ylim=c(-1,1)); 
abline(h=0.3)
abline(h=-0.3)
lines(c(-1,1.2),c(-1.2,1))
lines(c(-1.2,1),c(-1,1.2))
dev.off()

#writing the DMR analysis output:
outdmr.gr2=outdmr.gr[countOverlaps(outdmr.gr,downdmr)>0]
outmat=cbind(as.character(seqnames(outdmr.gr2)),start(outdmr.gr2),end(outdmr.gr2),as.matrix(mcols(outdmr.gr2)))
write.table(outmat,file=paste0("../results/downstream-analysis/dmr-identification/figures/",Sys.Date(),"-downDMR-table.txt"),sep="\t",row.names=FALSE,quote=FALSE)

methvalues=cbind(bmatrange,as.matrix(mcols(dmr.gr))[,5],countOverlaps(dmr.gr,h3peaks)>0)
methvalues2=cbind(as.matrix(mcols(dmr.gr))[,5]-bmatrange,countOverlaps(dmr.gr,h3peaks)>0)
dat=c()
num=c()
for(i in 1:30)
{
  s1=c((-15:15)/30)[i]
  s2=c((-15:15)/30)[i+1]
  dat[i]=mean(methvalues2[methvalues2[,1]>s1 & methvalues2[,1]<s2,2],na.rm=TRUE)
  num[i]=length(methvalues2[methvalues2[,1]>s1 & methvalues2[,1]<s2,2])
}
al=cbind(dat,c((-15:15)/30)[-1])
df=data.frame(rat=dat,loc=c((-15:15)/30)[-1])
pdf(paste0("../results/downstream-analysis/dmr-identification/figures/",Sys.Date(),"-MI-Ratio-of-DMRs-as-function-of-methylation-difference.pdf"))
ggplot(dat=df,aes(y=rat,x=loc))+geom_line()+scale_x_continuous(breaks=seq(-0.4,0.4,by=0.1))+
  xlab("Methylation diff. between the comparisons of CLLvsB and Bcell maturation")+ ylab("Ratio of DMRs overlapping with H3K27ac peaks")
dev.off()
al=cbind(dat,c((-15:15)/30)[-1])
df=data.frame(rat=num,loc=c((-15:15)/30)[-1])
pdf(paste0("../results/downstream-analysis/dmr-identification/figures/",Sys.Date(),"-MI-Ratio-of-DMRs-as-function-of-methylation-difference-number-of-DMRs.pdf"))
ggplot(dat=df,aes(y=rat,x=loc))+geom_line()+scale_x_continuous(breaks=seq(-0.4,0.4,by=0.1))+
  xlab("Methylation diff. between the comparisons of CLLvsB and Bcell maturation")+ ylab("Ratio of DMRs overlapping with H3K27ac peaks")
dev.off()

peakratio=
library(ggplot2)
ggplot(df,aes(x=x,y=y))+geom_density2d()



#########################################################################
#########################################################################
dmr3=read.table("../results/downstream-analysis/dmr-identification/downdmr-homer3/downdmr3.bed",sep="\t")
dmr3.gr=GRanges(seqnames=dmr3[,1],IRanges(start=dmr3[,2],end=dmr3[,3]))

#Alternative filtering of DMRs for visualization:
dt=as.matrix(mcols(dmr.gr))
methvalues=cbind(bmatrange,dt[,5])
bmatrange2=bmatrange
bmatrange2[is.na(bmatrange)]=1
#dmr selection:
downdmr=dmr.gr[(dt[,1]>200 & dt[,3]<0.3 & dt[,5]<(-0.3))]
bmatrange2=bmatrange[(dt[,1]>200 & dt[,3]<0.3 & dt[,5]<(-0.3))]

#pmd filtering
downdmr2=downdmr[countOverlaps(downdmr,pmdgr)==0]
bmatrange3=bmatrange2[countOverlaps(downdmr,pmdgr)==0]
# methvalues=cbind(bmatrange,dt[,5])
# selM=which(countOverlaps(dmr.gr,pmdgr)>0)
# methvalues3=methvalues[selM,]

dt2=as.matrix(mcols(downdmr2))
#B cell maturation:
methvalues=cbind(bmatrange3,dt2[,5])
bmatrange4=bmatrange3
bmatrange4[is.na(bmatrange3)]=1
downdmr3=downdmr2[dt2[,5]<(bmatrange4-0.2)]
selM=which(dt2[,5]<(bmatrange4-0.2))
bmatrange4[is.na(bmatrange)]=(-1)
updmr=dmr.gr2[dt2[,5]>(bmatrange4+0.2)& dt2[,5]>0.2]
selM=c(selM,which(dt[,5]>(bmatrange2+0.2)& dt2[,5]>0.2))
methvalues2=methvalues[selM,]

pdf(paste0("../results/downstream-analysis/dmr-identification/figures/",Sys.Date(),"-MI-diffmethylation-Bmaturation-vs-Cllcomp-filtered.pdf"))
smoothScatter(bmatrange,dt[,5],col="blue",colramp=colorRampPalette(c("white", "blue")),xlab="Diff. meth. during B cell maturation",ylab="Diff. methylation between CLL and B cells",xlim=c(-1,1),ylim=c(-1,1)); 
smoothScatter(methvalues2[,1],methvalues2[,2],col="red",colramp = colorRampPalette(c(rgb(1, 1, 1, 0), rgb(1, 0, 0, 1)), alpha = TRUE),add=T);
abline(h=0.2)
abline(h=-0.2)
dev.off()

selM=which(dt[,5]<(bmatrange2-0.2)& dt[,5]<(-0.2))
bmatrange2[is.na(bmatrange)]=(-1)
updmr=dmr.gr[dt[,5]>(bmatrange2+0.2)& dt[,5]>0.2]
selM=c(selM,which(dt[,5]>(bmatrange2+0.2)& dt[,5]>0.2))
methvalues2=methvalues[selM,]


#########################################################################
#########################################################################


####################################
# Filtering with PMDs:
methvalues=cbind(bmatrange,dt[,5])
selM=which(countOverlaps(dmr.gr,pmdgr)>0)
methvalues3=methvalues[selM,]
# pdf(paste0("../results/downstream-analysis/dmr-identification/figures/",Sys.Date(),"-MI-diffmethylation-Bmaturation-vs-Cllcomp-PMDs.pdf"))
# smoothScatter(bmatrange,dt[,5],col="blue",colramp=colorRampPalette(c("white", "blue")),xlab="Diff. meth. during B cell maturation",ylab="Diff. methylation between CLL and B cells",xlim=c(-1,1),ylim=c(-1,1)); 
# smoothScatter(methvalues3[,1],methvalues3[,2],col="red",colramp = colorRampPalette(c(rgb(1, 1, 1, 0), rgb(1, 0, 0, 1)), alpha = TRUE),add=T);
# dev.off()

sum(countOverlaps(dmr.gr,pmdgr)>0)/length(dmr.gr)
sum(countOverlaps(downdmr,pmdgr)>0)/length(downdmr)
sum(countOverlaps(updmr,pmdgr)>0)/length(updmr)
downdmr2=downdmr[countOverlaps(downdmr,pmdgr)==0]
restdowndmr2=downdmr[countOverlaps(downdmr,pmdgr)>0]
updmr2=updmr[countOverlaps(updmr,pmdgr)==0]
restupdmr2=updmr[countOverlaps(updmr,pmdgr)>0]

#H3K27ac ratios:
mean(countOverlaps(downdmr2,h3peaks)>0)
mean(countOverlaps(restdowndmr2,h3peaks)>0)
mean(countOverlaps(updmr2,h3peaks)>0)
mean(countOverlaps(restupdmr2,h3peaks)>0)

mean(countOverlaps(h3peaks,downdmr2)>0)

sum(countOverlaps(downdmr2,k27up)>0)
sum(countOverlaps(downdmr2,k27up)>0)
sum(countOverlaps(dmr.gr,k27up)>0)

write.table(cbind(as.character(seqnames(downdmr2)),start(downdmr2),end(downdmr2),paste0("dmr",1:length(downdmr2)),"A","+"),file="../results/downstream-analysis/dmr-identification/downdmr-homer2/downdmr2.bed",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

dt=as.matrix(mcols(downdmr2))
downdmr3=downdmr2[(dt[,1]>200 & dt[,3]<0.3 & dt[,5]<(-0.3))]
restdowndmr3=downdmr2[!(dt[,1]>100 & dt[,3]<0.3 & dt[,5]<(-0.3))]
#  mean(countOverlaps(downdmr3,h3peaks)>0)
#  mean(countOverlaps(restdowndmr3,h3peaks)>0)
# sum(countOverlaps(downdmr3,k27up)>0)
mean(countOverlaps(downdmr3,k27up)>0)
mean(countOverlaps(restdowndmr3,k27up)>0)

sum(countOverlaps(downdmr3,k27up)>0)

k27up2=k27up[countOverlaps(k27up,pmdgr)==0]
k27down2=k27down[countOverlaps(k27down,pmdgr)==0]
pdf("../results/downstream-analysis/dmr-identification/dmr-vs-k27up.pdf",width=4,height=3)
draw.pairwise.venn(length(downdmr3),length(k27up2),sum(countOverlaps(downdmr3,k27up2)>0),scaled=FALSE)
dev.off()
pdf("../results/downstream-analysis/dmr-identification/dmr-vs-k27down.pdf",width=4,height=3)
draw.pairwise.venn(length(downdmr3),length(k27down2),sum(countOverlaps(downdmr3,k27down2)>0),scaled=FALSE)
dev.off()
h3peaks2=h3peaks[countOverlaps(h3peaks,pmdgr)==0]
pdf("../results/downstream-analysis/dmr-identification/dmr-vs-k27peaks.pdf",width=4,height=3)
draw.pairwise.venn(length(downdmr3),length(h3peaks2),sum(countOverlaps(downdmr3,h3peaks2)>0),scaled=FALSE)
dev.off()
#Filtering based on Coverage and individual CpGs from the test:
dmrna=getGRNAnumber(cllm,downdmr2)
sum(dmrna>2)

dmlsgr=GRanges(seqnames=dmls[,1],IRanges(start=dmls[,2],end=dmls[,2]+1))
mcols(dmlsgr)=dmls[,c(5)]
dmlsgr2=dmlsgr
mcols(dmlsgr2)=dmls[,c(11)]

medianMdiff=getGRCovnumber(dmlsgr,downdmr2)

write.table(cbind(as.character(seqnames(downdmr3)),start(downdmr3),end(downdmr3),paste0("dmr",1:length(downdmr3)),"A","+"),file="../results/downstream-analysis/dmr-identification/downdmr-homer3/downdmr3.bed",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

randdowndmr3 = randomGRanges(downdmr3,chrslength = .chrs())
write.table(cbind(as.character(seqnames(randdowndmr3)),start(randdowndmr3),end(randdowndmr3),paste0("dmr",1:length(randdowndmr3)),"A","+"),file="../results/downstream-analysis/dmr-identification/downdmr-homer3/randdowndmr3.bed",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

write.table(cbind(as.character(seqnames(dmr.gr)),start(dmr.gr),end(dmr.gr),paste0("dmr",1:length(dmr.gr)),"A","+"),file="../results/downstream-analysis/dmr-identification/dmr.gr.bed",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
# finding motifs:
#findMotifsGenome.pl downdmr3.bed hg19 . -bg randdowndmr3.bed -size given -p 48 -S 20

# dt=as.matrix(mcols(downdmr2))
# selM2=countOverlaps(downdmr2,umr
# sum(dt[,1]>150 & dt[,5]<(-0.4) & dt[,3]<0.4)
# downdmr.gr4=dmr.gr2[(dt[,1]>150 & dt[,5]<(-0.4) & dt[,3]<0.4)]
# updmr.gr4=dmr.gr2[(dt[,1]>150 & dt[,5]>(0.3) & dt[,4]<0.3)]
# sum(as.matrix(mcols(dmr.gr2))[,3]<0.4)
# smoothScatter(dt[,3],dt[,5])
# sum(as.matrix(mcols(dmr.gr2))[,4]<(-0.3) & as.matrix(mcols(dmr.gr2))[,4]>0.3)
# #Filtering with B cell maturation
# #Showing enrichment of H3K27ac signal.
# write.table(cbind(as.character(seqnames(downdmr.gr4)),start(downdmr.gr4),end(downdmr.gr4),paste0("dmr",1:length(downdmr.gr4)),"A","+"),file="../results/downstream-analysis/dmr-identification/downdmr-homer/downdmr-gr4.bed",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
# write.table(cbind(as.character(seqnames(updmr.gr4)),start(updmr.gr4),end(updmr.gr4),paste0("dmr",1:length(updmr.gr4)),"A","+"),file="../results/downstream-analysis/dmr-identification/updmr-homer/updmr-gr4.bed",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
# findMotifsGenome.pl downdmr-gr4.bed hg19 . -size given -p 48 -S 10 -mask
# sum(countOverlaps(nonredumr,dmr.gr)>0)
# sigumr=nonredumr[(countOverlaps(nonredumr,dmr.gr)>0)]
# sigumr2=sigumr[(countOverlaps(sigumr,pmdgr)==0)]
# sum(countOverlaps(sigumr,pmdgr)>0)
# dmr.gr2=dmr.gr[countOverlaps(dmr.gr,sigumr2)>0]
# hist(as.matrix(mcols(dmr.gr2))[,1],40)



#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

tfbed=read.table("publicly-available-data/wgEncodeRegTfbsClusteredWithCellsV3.bed",sep="\t")
unTF=unique(tfbed[,4])
encodeTF=list()
size=c()
tfbed2=tfbed[tfbed[,5]>500,]
for(i in 1:length(unTF))
{
  a=tfbed2[tfbed2[,4]==unTF[i],]
  encodeTF[[i]]=GRanges(seqnames=a[,1],IRanges(start=as.numeric(a[,2]),end=as.numeric(a[,3])))
  size[i]=length(encodeTF[[i]])
}


TFenr=list()
TFodds=list()
TFratio=list()
TFratio2=list()
TFcount=list()
for(k in 19:21)
{
  print(k)
   queryGRanges=resize(downdmr3,width=200,fix="center")
   background=resize(dmr.gr[sample(1:length(dmr.gr),length(downdmr2))],width=1000,fix="center")
   #background=randomGRangesFixedWidth(downdmr2,width=2500)
    queryGRanges=downdmr2
    background=dmr.gr
  tfpval=c()
  odds=c()
  count=c()
  tcount=c()
  siteperc=c()
  siteperc2=c()
  for(i in 1:length(encodeTF))
  {
    #     if(i%%100==0)
    #     {
    #       print(c(k,i))
    #     }
    etf=resize(encodeTF[[i]],width=200,fix="center")
    x=sum(countOverlaps(queryGRanges,etf)>0)
    y=length(queryGRanges)
    z=sum(countOverlaps(background,etf)>0)
    t=length(background)
    tfpval[i]=fisher.test(matrix(c(x,y-x,z,t-z),nr=2),alternative="greater")$p.value
    odds[i]=fisher.test(matrix(c(x,y-x,z,t-z),nr=2),alternative="greater")[[3]]
    count[i]=x
    tcount[i]=length(etf)
    siteperc[i]=x/y
    siteperc2[i]=z/t
  }
  #print(c(x,y,z,t))
  adjtfpval=p.adjust(tfpval)
  names(adjtfpval)=unTF
  TFenr[[k]]=adjtfpval
  names(odds)=unTF
  TFodds[[k]]=odds
  TFcount[[k]]=tcount
  TFratio[[k]]=siteperc
  TFratio2[[k]]=siteperc2
}
k=1
write.table((cbind(k, length(queryGRanges),names(sort(TFodds[[k]],decreasing=TRUE)[1:161]), sprintf("%.3f",sort(TFodds[[k]],decreasing=TRUE)[1:161]),sprintf("%.3f",TFratio[[k]][order(TFodds[[k]],decreasing=TRUE)[1:161]]),
            sprintf("%.3f",TFcount[[k]][order(TFodds[[k]],decreasing=TRUE)[1:161]]),TFenr[[k]][order(TFodds[[k]],decreasing=TRUE)[1:161]] )),
            file="../results/downstream-analysis/dmr-identification/TF-enrichment.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
# Separate with Promoters and check the FPKM values for the nearby PMD regions:
#RNA-seq data:
rnaseq=read.table("../results/rna-seq/Sandra-RNAseq/RNA-seq-results-normalized-matrix/coding_results_diff_c1_pAll2_withCoords.csv",sep="\t",header=TRUE)
rownames(rnaseq)=rnaseq[,5]
tss=rnaseq[,2]
tss[rnaseq[,4]=="-"]=rnaseq[rnaseq[,4]=="-",3]
rnagr=GRanges(seqnames=paste("chr",rnaseq[,1],sep=""),IRanges(start=tss,end=tss))
rnagr3=GRanges(seqnames=paste("chr",rnaseq[,1],sep=""),IRanges(start=rnaseq[,2],end=rnaseq[,3]))
mcols(rnagr)=rnaseq[,c(5,6,7,10,12)]
mcols(rnagr3)=rnaseq[,c(5,6,7,10,12)]
#
add=rep(5000,length(tss))
add[rnaseq[,4]=="-"]=add[rnaseq[,4]=="-"]*(-1)

RNAdownst=GRanges(seqnames=paste("chr",rnaseq[,1],sep=""),IRanges(start=apply(cbind(tss,tss+add),1,min),end=apply(cbind(tss,tss+add),1,max)))
RNAupst=GRanges(seqnames=paste("chr",rnaseq[,1],sep=""),IRanges(start=apply(cbind(tss,tss-add),1,min),end=apply(cbind(tss,tss-add),1,max)))
mcols(RNAdownst)=rnaseq[,c(5,6,7,10,12)]
mcols(RNAupst)=rnaseq[,c(5,6,7,10,12)]
rnagb=GRanges(seqnames=paste("chr",rnaseq[,1],sep=""),IRanges(start=rnaseq[,2],end=rnaseq[,3]))
mcols(rnagb)=rnaseq[,c(5,6,7,10,12)]
# Overlap with PMDs and filter if needed:

rnagr2=rnagr[countOverlaps(rnagr,pmdgr)==0]
rnagr4=rnagr3[countOverlaps(rnagr3,pmdgr)==0]
#TAD data:
tad=read.table("publicly-available-data/TAD-IMR90_domains_hg19.bed",sep="\t")
tadgr=GRanges(seqnames=tad[,1],IRanges(start=tad[,2],end=tad[,3]))
#active TADs:
uptad=tadgr[countOverlaps(tadgr,downdmr2)>0]
dmrgenes=as.matrix(mcols(rnagr[countOverlaps(rnagr,uptad)>0]))[,1]

df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr2))[,4]),dmr=as.numeric(countOverlaps(resize(rnagr2,width=5000,fix="center"),downdmr2[1:2500])>0))
pdf("../results/downstream-analysis/dmr-identification/dmr-vs-expression.pdf",width=3,height=4)
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))+theme_bw()
dev.off()
table(df[,2])
df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr2))[,4]),dmr=as.numeric(countOverlaps(resize(rnagr2,width=5000,fix="center"),k27up)>0))
pdf("../results/downstream-analysis/dmr-identification/k27up-vs-expression.pdf",width=3,height=4)
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,5))+theme_bw()
dev.off()
table(df[,2])

df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr2))[,4]),dmr=as.numeric(countOverlaps(resize(rnagr2,width=5000,fix="center"),k27down)>0))
pdf("../results/downstream-analysis/dmr-identification/k27down-vs-expression.pdf",width=3,height=4)
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-5,3))+theme_bw()
dev.off()
table(df[,2])

uptad=tadgr[countOverlaps(tadgr,k27up)>0]
df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr2))[,4]),dmr=as.numeric(countOverlaps(resize(rnagr2,width=5000,fix="center"),uptad)>0))
pdf("../results/downstream-analysis/dmr-identification/k27uptad-vs-expression.pdf",width=3,height=4)
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-5,3))+theme_bw()
dev.off()
table(df[,2])

downdmr3
df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr4))[,4]),dmr=as.numeric(countOverlaps(rnagr4,downdmr3)>0))
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))

df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr4))[,4]),dmr=as.numeric(countOverlaps(rnagr4,k27up)>0))
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))

df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr4))[,4]),dmr=as.numeric(countOverlaps(rnagr4,k27down)>0))
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))


df=data.frame(exp=as.numeric(as.matrix(mcols(RNAdownst))[,4]),dmr=as.numeric(countOverlaps(RNAdownst,downdmr2)>0))
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))

df=data.frame(exp=as.numeric(as.matrix(mcols(RNAupst))[,4]),dmr=as.numeric(countOverlaps(RNAupst,downdmr2)>0))
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))

#active TADs:
uptad=tadgr[countOverlaps(tadgr,k27up)>0]
df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr4))[,4]),dmr=as.numeric(countOverlaps(rnagr4,uptad)>0))
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))

uptad=tadgr[countOverlaps(tadgr,k27up)>0]
df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr2))[,4]),dmr=as.numeric(countOverlaps(resize(rnagr2,5000,fix="center"),k27up)>0))
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))

downtad=tadgr[countOverlaps(tadgr,k27down)>0]
df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr2))[,4]),dmr=as.numeric(countOverlaps(resize(rnagr2,10000,fix="center"),downtad)>0))
ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))

upgen=nearest(downdmr2,rnagr2)
grdistan1=distance(downdmr2,rnagr[upgen])
#df=data.frame(exp=as.numeric(as.matrix(mcols(rnagr2))[,4]),dmr=c(1:length(rnagr2)%in%unique(upgen)))
#ggplot(df,aes(y=as.numeric(exp),x=as.factor(dmr)))+geom_boxplot()+ coord_cartesian(ylim=c(-3,3))
fcexp1=cbind(as.matrix(mcols(downdmr2))[,2],as.matrix(mcols(rnagr2[upgen]))[,3])
df=data.frame(m=round(as.numeric(fcexp1[,1]),digits=1),e=as.numeric(fcexp1[,2]))
ggplot(df,aes(y=as.numeric(e),x=as.factor(m)))+geom_boxplot()+ coord_cartesian(ylim=c(0,5000))
ggplot(df,aes(y=e,x=d))+geom_point(alpha=0.2)+ coord_cartesian(ylim=c(0,2500))

upgen=nearest(downdmr2,rnagr)
grdistan1=distance(downdmr2,rnagr[upgen])
fcexp1=cbind(as.matrix(mcols(downdmr2))[,5],as.matrix(mcols(rnagr[upgen]))[,2],grdistan1)
df=data.frame(m=as.numeric(fcexp1[,1]),e=as.numeric(fcexp1[,2]),d=as.numeric(fcexp1[,3]))
ggplot(df,aes(x=m,y=e))+stat_density2d()+geom_point(alpha=0.2)+scale_y_continuous(limits=c(0,100))
ggplot(df,aes(y=e,x=d))+geom_point(alpha=0.2)+ coord_cartesian(ylim=c(0,2500))

#promoter vs distal:
sum(countOverlaps(resize(rnagr,width=3000,fix="center"),downdmr2)>0)
sum(countOverlaps(downdmr2,nonredumr2)>0)
sum(countOverlaps(nonredumr2,downdmr2)>0)
sum(countOverlaps(downdmr2,dmv2)>0)
sum(countOverlaps(dmv2,downdmr2)>0)
sum(countOverlaps(downdmr2,k27up)>0)
sum(countOverlaps(downdmr3,k27up)>0)

cand=19233
sum(countOverlaps(downdmr2[cand],dmv2)>0)
sum(countOverlaps(downdmr2[cand],nonredumr2)>0)
sum(countOverlaps(downdmr2[cand],resize(rnagr,width=10000,fix="center"))>0)
downdmr2[cand]
              
              
k27upgenes=as.matrix(mcols(rnagr[upgen]))[,1]
grdistan1=distance(downdmr2,rnagr[upgen])

f=as.matrix(mcols(rnagr[upgen]))
sel=which(as.numeric(f[,3])>(2) & as.numeric(f[,4])<0.01)
#k27up=k27up[ sel ]
pdf("../results/downstream-analysis/dmr-tf-enrichment/k27up-gene-FC-hist.pdf")
hist(as.numeric(f[,3]),20)
dev.off()
downgen=nearest(k27down,rnagr)
k27downgenes=as.matrix(mcols(rnagr[downgen]))[,1]
grdistan2=distance(k27down,rnagr[downgen])
exp2=as.matrix(mcols(rnagr[downgen]))[,3]
f=as.matrix(mcols(rnagr[downgen]))
sel=which(as.numeric(f[,3])<(-2) & as.numeric(f[,4])<0.01)
#plot distance vs expression fold change relation for the nearest TSS:
a=rbind(cbind(exp1,grdistan1,"Increased-diffbind"),cbind(exp2,grdistan2,"Decreased-diffbind"))
a[a[,2]==0,2]=10
a=a[is.finite(as.numeric(a[,1])),]
df=data.frame(FoldChange=as.numeric(a[,1]),DistancetoTSS=as.numeric(a[,2]),Direction=as.factor(a[,3]))
ggplot(df,aes(x=DistancetoTSS,y=FoldChange,colour=Direction))+
  stat_smooth( method=lm, formula = y ~ poly(x,1), level=0.95) +
  stat_smooth( method=lm, formula = y ~ poly(x,1), level=0.95) +
  geom_jitter(aes(alpha=0.5),width = 0.1)+scale_x_continuous(trans="log10")+theme_bw()
ggsave(paste("../results/downstream-analysis/dmr-tf-enrichment/",Sys.Date(),"-MI-Diffbind-distancetoTSS-expressionFC.pdf",sep=""),height=12,width=16)



