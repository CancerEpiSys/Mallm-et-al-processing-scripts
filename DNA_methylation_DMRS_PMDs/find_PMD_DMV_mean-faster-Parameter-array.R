## find_PMD_DMV.R
## v0.0.1
## volker hovestadt <git@hovestadt.bio>
## 25 feb 2015
##
## script to detect DMVs and PMDs, as used in hovestadt et al. 2014.
## based on approaches described in xie et al. 2013 and berman et al. 2012, respectively.
##
## DMVs:
## "Average methylation levels within windows of 1 kb were calculated
## (individual CpGs weighted by the distance of both adjacent CpGs), in steps of
## 1 bp. Overlapping 1 kb windows with average methylation levels smaller than 0.15
## were merged and resulting regions larger than 5 kb were termed DMVs."
##
## PMDs:
## "Average methylation levels within windows of 10 kb were calculated 
## (individual CpGs weighted by the distance of both second-next CpGs), in steps of
## 1 bp. Overlapping 10 kb windows with an average methylation level <0.6 were
## merged, and resulting regions larger than 100 kb were termed PMDs."


options(max.print=200)
options(stringsAsFactors=FALSE)


library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
#library(TTR)  # used for running mean calculation, very fast
library(caTools) # used for running mean calculation, faster than TTR, handles NAs

## parameters (hardcoded)

# number of CpGs combined
#n <- 3  # weighting by adjacent CpGs
n <- 5  # weighting by second-next CpGs
# window size
w.dmv <- 1000
w.pmd <- 10000
# thresholds
threshold.dmv <- 0.15
threshold.pmd <- 0.5
threshold.dmv.arr <- seq(15,15,by=1)/100
threshold.pmd.arr <- seq(50,50,by=1)/100


## load data (as produced by bcall2beta.R)

#in.file <- "/icgc/lsdf/mb/analysis/hovestad/ICGC_methylome/final_analysis/pmd/1_dmcombined/Beta_20130118/dm.combined.MB113_20120720_5lanes.RData.Beta.RData"
in.file <- commandArgs(trailingOnly = TRUE)[1]
#in.file= "dm.combined.H038-MCJT_complete_filtered_merged.dupmarked.RData.Beta.RData"
load(in.file)
print(in.file)
#head(Beta, n = 100)

s="label"
in.file.object <- setdiff(ls(), "in.file")

## cpg annotation
load("/icgc/dkfzlsdf/analysis/CancerEpiSys/methylation/results-mem-1000gen/results/wgbs/CES_BK-1/dm.combined.CES_BK-1_complete_filtered_merged.dupmarked.RData")
gr=GRanges(seqnames=dm.chr.combined,IRanges(start=dm.pos.combined, end=dm.pos.combined+1))
# GRanges of CpG positions (width: 2), must be in the same order as Beta!
#hg19.csize <- getChromInfoFromUCSC("hg19") # load once and save
#save(hg19.csize, file="hg19.csize.20130531.RData")
load("/icgc/dkfzlsdf/analysis/CancerEpiSys/methylation/results-mem-1000gen/scripts/publicly-available-data/hg19.csize.20140513.RData")
seqlengths(gr) <- hg19.csize[match(as.vector(seqlevels(gr)), hg19.csize[, 1]), 2]
length(gr) == length(Beta)

#mySession <- browserSession("UCSC")
#genome(mySession) <- "hg19"
#getTable(ucscTableQuery(mySession, track="gap", table="gap"))

#ucsc.gap <- sort(features(makeFeatureDbFromUCSC("hg19", "gap", "gap"))) # load once and save
#save(ucsc.gap.gr, file="ucsc.gap.gr.20130531.RData")
load("/icgc/dkfzlsdf/analysis/CancerEpiSys/methylation/results-mem-1000gen/scripts/publicly-available-data/ucsc.gap.gr.20140513.RData")


## process

# split by chromosome
values(gr)$Beta <- Beta
gr.split <- split(gr, seqnames(gr))
rm(gr)
gc()


grl.window <- GRangesList()

grl.dmv <- GRangesList()
grl.dmv.list <- rep(list(grl.dmv),length(threshold.dmv.arr))
grl.pmd <- GRangesList()
grl.pmd.list <- rep(list(grl.pmd),length(threshold.pmd.arr))
seqlevels(grl.window) <- seqlevels(gr.split)
seqlengths(grl.window) <- seqlengths(gr.split)
seqlevels(grl.dmv) <- seqlevels(gr.split)
seqlengths(grl.dmv) <- seqlengths(gr.split)
seqlevels(grl.pmd) <- seqlevels(gr.split)
seqlengths(grl.pmd) <- seqlengths(gr.split)

# shuffle chromosomes, memory usage depends on chromosome size, ...
chrlist=sample(names(sort(seqlengths(gr.split), decreasing=TRUE)))
# Y chromosome removed from analysis.
#chrlist=chrlist[!chrlist%in%"chrY"]
for(chr in chrlist) {
  gr.chr <- gr.split[[chr]]
  
  message(paste(in.file.object, s, chr, "- combine CpGs"))
  # define regions spanning n CpGs
  gr.window <- GRanges(seqnames = chr, ranges = IRanges(start = start(gr.chr)[1:(length(gr.chr)-n+1)]+1, end = end(gr.chr)[n:length(gr.chr)]-1), Beta.mean = as.numeric(rep(NA, length(gr.chr)-n+1)))
  seqlevels(gr.window) <- seqlevels(gr.split)
  seqlengths(gr.window) <- seqlengths(gr.split)
  
  gr.window <- gr.window[! overlapsAny( gr.window , ucsc.gap.gr)]
  gr.window.overlap <- as.matrix(findOverlaps(gr.window, gr.chr))
  values(gr.window[unique(gr.window.overlap[, 1])])$Beta.mean <- rowMeans(matrix(unlist(split(values(gr.chr)[gr.window.overlap[, 2], "Beta"], gr.window.overlap[, 1]), use.names=FALSE), ncol = 5, byrow = TRUE), na.rm=TRUE)
  grl.window[[chr]] <- gr.window
    
  message(paste(in.file.object, s, chr, "- vector"))
  # put in vector for calculating running mean
  v.gap.sum <- rep(0, seqlengths(gr.chr)[chr])
  v.gap.n <- rep(0, seqlengths(gr.chr)[chr])
  gr.window.start <- start(gr.window)
  gr.window.end <- end(gr.window)
  gr.window.beta <- values(gr.window)$Beta.mean
  
  for(i in which(!is.na(gr.window.beta))) {
    temppos= seq.int(gr.window.start[i], gr.window.end[i]);
    v.gap.sum[temppos] <- v.gap.sum[temppos] + gr.window.beta[i]
    v.gap.n[temppos] <- v.gap.n[temppos] + 1
  }
  v.gap.beta <- v.gap.sum/v.gap.n	
  rm(list=c("gr.window.overlap", "gr.window.start", "gr.window.end", "gr.window.beta", "v.gap.sum", "v.gap.n"))
  gc()
  

  message(paste(in.file.object, s, chr, "- DMV running mean"))
  # WMA doesnt handle NA, use weights 0 and 1
  #v.gap.weight <- rep(1, length(v.gap.beta))
  #v.gap.weight[is.na(v.gap.beta)] <- 0
  #v.gap.beta2=v.gap.beta
  #v.gap.beta2[is.na(v.gap.beta)] <- 1
  v.gap.beta.dmv <- runmean(v.gap.beta,w.dmv,align="right",endrule="NA")[w.dmv:length(v.gap.beta)]
  #v.gap.beta.dmv <- WMA(v.gap.beta2[1000000:1050000], n=w.dmv, wts=v.gap.weight[1000000:1050000])[w.dmv:length(v.gap.beta)]
  # returns -Inf if all wts are 0, replace to NA
  #in case of WMA, this is required
  #v.gap.beta.dmv <- ifelse(is.finite(v.gap.beta.dmv), v.gap.beta.dmv, NA)
  gc()
  
  gr.window.dmv.org <- GRanges(seqnames = chr, ranges = IRanges(start = seq(1, (seqlengths(gr.chr)[chr]-w.dmv+1), 1), end = seq(w.dmv, seqlengths(gr.chr)[chr], 1)))
  message(paste(in.file.object, s, chr, "- DMV window"))
  for(i in 1:length(threshold.dmv.arr))
  {
    threshold.dmv <- threshold.dmv.arr[i]
    gr.window.dmv <- gr.window.dmv.org
    seqlevels(gr.window.dmv) <- seqlevels(gr.split)
    seqlengths(gr.window.dmv) <- seqlengths(gr.split)
  #length(v.gap.beta.dmv) == length(gr.window.dmv)
  
  # only keep windows with mean methylation value below threshold
    gr.dmv <- gr.window.dmv[which(v.gap.beta.dmv < threshold.dmv)]
  # remove windows which contain gaps in the genome assembly
    gr.dmv <- gr.dmv[! overlapsAny(gr.dmv, ucsc.gap.gr)]
  # combine overlapping windows
    gr.dmv.reduce <- reduce(gr.dmv)
  # number of CpGs
    values(gr.dmv.reduce)$score <- countOverlaps(gr.dmv.reduce, gr.chr)
    gr.dmv.length <- gr.dmv.reduce[width(gr.dmv.reduce) >= 5000]
    if(length(gr.dmv.length)) { grl.dmv.list[[i]][[chr]] <- gr.dmv.length }
  }
  message(paste(in.file.object, s, chr, "-", length(gr.dmv.length), "DMVs found\n"))
  rm(list=c("v.gap.beta.dmv", "gr.window.dmv", "gr.dmv", "gr.dmv.reduce")) 
  gc()
  
  message(paste(in.file.object, s, chr, "- PMD running mean"))
  v.gap.beta.pmd <- runmean(v.gap.beta,w.pmd,align="right",endrule="NA")[w.pmd:length(v.gap.beta)]
  #v.gap.beta.pmd <- WMA(v.gap.beta, n=w.pmd, wts=v.gap.weight)[w.pmd:length(v.gap.beta)]
  # returns -Inf if all wts are 0, replace to NA
  #v.gap.beta.pmd <- ifelse(is.finite(v.gap.beta.pmd), v.gap.beta.pmd, NA)
  
  #rm(list=c("v.gap.beta"))
  gc()
  
  gr.window.pmd.org <- GRanges(seqnames = chr, ranges = IRanges(start = seq(1, (seqlengths(gr.chr)[chr]-w.pmd+1), 1), end = seq(w.pmd, seqlengths(gr.chr)[chr], 1)))
  message(paste(in.file.object, s, chr, "- PMD window"))
  for(i in 1:length(threshold.pmd.arr))
  {
    threshold.pmd <- threshold.pmd.arr[i]
    gr.window.pmd <- gr.window.pmd.org
    seqlevels(gr.window.pmd) <- seqlevels(gr.split)
    seqlengths(gr.window.pmd) <- seqlengths(gr.split)
  
    # only keep windows with mean methylation value below threshold
    gr.pmd <- gr.window.pmd[which(v.gap.beta.pmd < threshold.pmd)]
  
    #rm(list=c("v.gap.beta.pmd", "gr.window.pmd"))
    # remove windows which contain gaps in the genome assembly
    gr.pmd <- gr.pmd[! overlapsAny(gr.pmd, ucsc.gap.gr)]
    gc()
    # combine overlapping windows
    gr.pmd.reduce <- reduce(gr.pmd)
    # number of CpGs
    values(gr.pmd.reduce)$score <- countOverlaps(gr.pmd.reduce, gr.chr)
    # apply filter on PMD (berman et al., >= 100kb)
    gr.pmd.length <- gr.pmd.reduce[width(gr.pmd.reduce) >= 100000]
    if(length(gr.pmd.length)) { grl.pmd.list[[i]][[chr]] <- gr.pmd.length }
  }
  message(paste(in.file.object, s, chr, "-", length(gr.pmd.length), "PMDs found\n"))
    rm(list=c("gr.pmd", "gr.pmd.length"))
    gc()
  
  warnings()
}

save(grl.pmd.list,grl.dmv.list,file=sub(".RData.Beta.RData", paste(".DMV.PMD.ARRAY.065-015Parameter.", ".RData", sep=""), in.file))

## save output files

# # as Rdata
# grl.window.unlist <- sort(unlist(grl.window, use.names=FALSE))
# out.file <- sub(".RData.Beta.RData", paste(".cpg", n, ".RData", sep=""), in.file)
# save(grl.window.unlist, file=out.file)

# # as bed
# grl.pmd.unlist <- sort(unlist(grl.pmd, use.names=FALSE))
# out.file <- sub(".RData.Beta.RData", paste(".cpg", n, "window", w.pmd, ".07.PMD.bed", sep=""), in.file)
# if(length(grl.pmd.unlist)>0)
# {
#   export(grl.pmd.unlist, out.file)
# }

# grl.dmv.unlist <- sort(unlist(grl.dmv, use.names=FALSE))
# out.file <- sub(".RData.Beta.RData", paste(".cpg", n, "window", w.dmv, ".015.2.DMV.bed", sep=""), in.file)
# export(grl.dmv.unlist, out.file)

q("no")


## when complete, combine all samples into GRangesList (run manually)

lf.pmd <- list.files(pattern="cpg5window10000.PMD.015.bed")
lf.dmv <- list.files(pattern="cpg5window1000.DMV.07.bed")

grl.pmd <- GRangesList()
grl.dmv <- GRangesList()

for(f in lf.pmd) {
  #message(f)
  f.short <- strsplit(f, "\\.")[[1]][3]
  message(f.short)
  gr <- import(f, asRangedData=FALSE)
  grl.pmd[[f.short]] <- gr
}

for(f in lf.dmv) {
  #message(f)
  f.short <- strsplit(f, "\\.")[[1]][3]
  message(f.short)
  gr <- import(f, asRangedData=FALSE)
  grl.dmv[[f.short]] <- gr
}

grl.dmv.filter <- mendoapply(function(e1, e2) e1[! e1 %in% e2], grl.dmv, grl.pmd[names(grl.dmv)])

sort(sapply(grl.dmv.filter, function(x) sum(width(x))))
sort(sapply(grl.dmv.filter, function(x) length(x[width(x) >= 5000])))

in.file=sub(".RData.Beta.RData","",in.file);
save(list=c("grl.pmd", "grl.dmv", "grl.dmv.filter"), file=paste(in.file,".pmd_dmv.RData",sep=""))

# to decide for the best parameter, we calculate PMD-DMV for a range of parameters. then check the size and length of segments.

options(max.print=200)
options(stringsAsFactors=FALSE)


library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)


threshold.dmv.arr <- seq(10,20,by=1)/100
threshold.pmd.arr <- seq(60,70,by=1)/100
l=list.files(pattern=".RData")

dmvsize=mat.or.vec(nr=21,nc=11)
dmvlength=mat.or.vec(nr=21,nc=11)
pmdsize=mat.or.vec(nr=21,nc=11)
pmdlength=mat.or.vec(nr=21,nc=11)
genomecov=mat.or.vec(nr=21,nc=11)
for(i in 1:length(l))
{
  print(i)
  load(l[i])
  for(k in 1:11)
  {
    b=unlist(grl.pmd.list[[k]])
    #b=b[!overlapsAny(b,ucsc.gap.gr)]
    pmdsize[i,k]=length(b)
    pmdlength[i,k]=mean(width(b))
    genomecov[i,k]=sum(width(b))/sum(as.numeric(hg19.csize[1:22,2]),na.rm=TRUE)
    a=unlist(grl.dmv.list[[k]])
    a=a[!overlapsAny(a,unlist(grl.pmd.list[[1]]))]
    dmvsize[i,k]=length(a)
    dmvlength[i,k]=mean(width(a))
  }
}

dmvplot=cbind(apply(dmvsize,2,median),apply(dmvlength,2,median),c(threshold.pmd.arr))
dmvplot2=data.frame(X=dmvplot[,1],Y=dmvplot[,2],Col=dmvplot[,3])
library(ggplot2)
p <- ggplot(dmvplot2, aes(X, Y))
p + geom_point(aes(colour = Col))
ggsave("dmv-size-length.pdf")


for(i in 1:length(l))
{
load(l[i])
grl.pmd.unlist <- sort(unlist(grl.pmd.list[[6]], use.names=FALSE))
out.file <- sub(".RData", paste(".cpg", "window065thres", "PMD.bed", sep=""), l[i])
out.file <- paste("defined-pmds/",out.file,sep="")
export(grl.pmd.unlist, out.file)
}

for(i in 1:length(l))
{
load(l[i])
grl.dmv.unlist <- sort(unlist(grl.dmv.list[[6]], use.names=FALSE))
out.file <- sub(".RData", paste(".cpg", "window015thres", "DMV.bed", sep=""), l[i])
out.file <- paste("defined-dmvs/",out.file,sep="")
export(grl.dmv.unlist, out.file)
}


gr.window <- GRanges(seqnames = chr, ranges = IRanges(start = start(gr.chr)[1:(length(gr.chr)-n+1)]+1, end = end(gr.chr)[n:length(gr.chr)]-1), Beta.mean = as.numeric(rep(NA, length(gr.chr)-n+1)))

i=5
n=load(l[i])

pmdlist= unlist(grl.pmd.list[[1]])
averMeth=c()
save(pmdlist,file="example-pmdlist.RData")
load("example-pmdlist.RData")

l=list.files("/icgc/dkfzlsdf/analysis/CancerEpiSys/methylation/results-mem-1000gen/beta-values-combined-chr-X-Y/bed-format/sorted",pattern=)
cpgprof=read.table(paste("/icgc/dkfzlsdf/analysis/CancerEpiSys/methylation/results-mem-1000gen/beta-values-combined-chr-X-Y/bed-format/sorted/",l[1],sep=""),sep="\t",quote="")
gr.window <- GRanges(seqnames = cpgprof[,1], ranges = IRanges(start = cpgprof[,2], end = cpgprof[,2]+1), Beta.mean = as.numeric(cpgprof[,3])))
gr.window <- GRanges(seqnames = "chr1", ranges = IRanges(start = 470250, end = 520500), Beta.mean = as.numeric(1)))

load("/icgc/dkfzlsdf/analysis/CancerEpiSys/methylation/results-mem-1000gen/beta-values-combined-chr-X-Y/dm.combined.CES_FG-1_complete_filtered_merged.dupmarked.RData")

l=list.files("/icgc/dkfzlsdf/analysis/CancerEpiSys/methylation/results-mem-1000gen/beta-values-combined-chr-X-Y/",pattern="Beta.RData")
pmdCGBeta=c()
for(i in 1:length(l))
{
  print(i)
  load(paste("/icgc/dkfzlsdf/analysis/CancerEpiSys/methylation/results-mem-1000gen/beta-values-combined-chr-X-Y/",l[i],sep=""))
  gr.window <- GRanges(seqnames = dm.chr.combined, ranges = IRanges(start = dm.pos.combined, end = dm.pos.combined+1), Beta.mean = as.numeric(Beta))
  pmdCGBeta[i]=mean(gr.window[overlapsAny(gr.window,pmdlist)]$Beta.mean,na.rm=TRUE)
}
pdf("PMD-Beta-mean.pdf")
plot(density(pmdCGBeta))
dev.off()

write.table(cbind(l,pmdCGBeta,pmdsize),"PMD-parameter-summary.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)





rnaseq=read.table("RNA-seq-gene-location-unique-list.txt",sep="\t",quote="")
RNAchr=rnaseq[,1]
RNAstart=rnaseq[,2]
RNAend=rnaseq[,3]
names(RNAchr)=rnaseq[,4]
names(RNAstart)=rnaseq[,4]
names(RNAend)=rnaseq[,4]
FC=read.table("RNA-seq-DeSeq.txt",sep="\t",quote="")
RNAseqload=cbind(sapply(RNAchr[as.character(rownames(FC))],function(x)paste("chr",x,sep="")),RNAstart[as.character(rownames(FC))],RNAend[as.character(rownames(FC))],as.numeric(FC[,2]))
RNAseqload=RNAseqload[!is.na(RNAseqload[,1]),]
RNAseqload=RNAseqload[!is.na(RNAseqload[,2]),]
rnaseq.window= GRanges(seqnames = RNAseqload[,1], ranges = IRanges(start = as.numeric(RNAseqload[,2]), end = as.numeric(RNAseqload[,3])), FC = as.numeric(RNAseqload[,4]))

RNApmdin=rnaseq.window[overlapsAny(rnaseq.window,pmdlist)]
RNApmdout=rnaseq.window[!overlapsAny(rnaseq.window,pmdlist)]

rnabox=rbind(cbind(RNApmdin$FC,1),cbind(RNApmdout$FC,2))
rnabox=data.frame(V1=rnabox[,1],V2=rnabox[,2])
pdf("RNAseq-PMD-FC-boxplot.pdf")
boxplot(V1 ~ V2, data = rnabox, 
    outline = FALSE)
  #beeswarm(V1 ~ V2, data = rnabox, 
  #  col = 4, pch = 5, add = TRUE)
dev.off()

#focusing only on significant genes:
RNAseqload=cbind(sapply(RNAchr[as.character(rownames(FC))],function(x)paste("chr",x,sep="")),RNAstart[as.character(rownames(FC))],RNAend[as.character(rownames(FC))],as.numeric(FC[,2]),as.numeric(FC[,6]))
RNAseqload=RNAseqload[!is.na(RNAseqload[,1]),]
RNAseqload=RNAseqload[!is.na(RNAseqload[,2]),]
RNAseqload=RNAseqload[as.numeric(RNAseqload[,5])<1e-5,]
RNAseqloadup=RNAseqload[as.numeric(RNAseqload[,4])>2,]
RNAseqloaddown=RNAseqload[as.numeric(RNAseqload[,4])< (-2),]
rnaseq.window.up= GRanges(seqnames = RNAseqloadup[,1], ranges = IRanges(start = as.numeric(RNAseqloadup[,2]), end = as.numeric(RNAseqloadup[,3])), FC = as.numeric(RNAseqloadup[,4]))
rnaseq.window.down= GRanges(seqnames = RNAseqloaddown[,1], ranges = IRanges(start = as.numeric(RNAseqloaddown[,2]), end = as.numeric(RNAseqloaddown[,3])), FC = as.numeric(RNAseqloaddown[,4]))

rnaseq.window.rest = rnaseq.window[!overlapsAny(rnaseq.window,rnaseq.window.up)]
rnaseq.window.rest = rnaseq.window.rest[!overlapsAny(rnaseq.window.rest,rnaseq.window.down)]

a=c()
b=c()
d=c()
for(i in 1:21)
{
  print(l[i])
  load(l[i])
  pmdlist=grl.pmd.list[[4]]
  if(sum(overlapsAny(rnaseq.window.down,pmdlist))/length(rnaseq.window.down)>0.1)
  {
    a=c(a,sum(overlapsAny(rnaseq.window.up,pmdlist))/length(rnaseq.window.up)) 
    b=c(b,sum(overlapsAny(rnaseq.window.down,pmdlist))/length(rnaseq.window.down))
    d=c(d,sum(overlapsAny(rnaseq.window,pmdlist))/length(rnaseq.window))
  }
  print(sum(overlapsAny(rnaseq.window.up,pmdlist))/length(rnaseq.window.up))
  print(sum(overlapsAny(rnaseq.window.down,pmdlist))/length(rnaseq.window.down))
  print(sum(overlapsAny(rnaseq.window.rest,pmdlist))/length(rnaseq.window.rest))
}
mean(a)
mean(b)
mean(d)




#ChIP-seq analysis:

dat=read.table("chip-seq/multiintersect_patients_H3K27ac.txt",sep="\t")
chipseq.window.dis.H3K27ac= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_healthy_H3K27ac.txt",sep="\t")
chipseq.window.heal.H3K27ac= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_patients_H3K9ac.txt",sep="\t")
chipseq.window.dis.H3K9ac= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_healthy_H3K9ac.txt",sep="\t")
chipseq.window.heal.H3K9ac= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_patients_H3K9me3.txt",sep="\t")
chipseq.window.dis.H3K9me3= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_healthy_H3K9me3.txt",sep="\t")
chipseq.window.heal.H3K9me3= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_patients_H3K4me1.txt",sep="\t")
chipseq.window.dis.H3K4me1= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_healthy_H3K4me1.txt",sep="\t")
chipseq.window.heal.H3K4me1= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_patients_H3K4me3.txt",sep="\t")
chipseq.window.dis.H3K4me3= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_healthy_H3K4me3.txt",sep="\t")
chipseq.window.heal.H3K4me3= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_patients_H3K27me3.txt",sep="\t")
chipseq.window.dis.H3K27me3= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_healthy_H3K27me3.txt",sep="\t")
chipseq.window.heal.H3K27me3= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_patients_H3K36me3.txt",sep="\t")
chipseq.window.dis.H3K36me3= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))
dat=read.table("chip-seq/multiintersect_healthy_H3K36me3.txt",sep="\t")
chipseq.window.heal.H3K36me3= GRanges(seqnames = dat[,1], ranges = IRanges(start = as.numeric(dat[,2]), end = as.numeric(dat[,3])))


PMD.chip.dis=c(sum(overlapsAny(chipseq.window.dis.H3K27ac,pmdlist)),sum(overlapsAny(chipseq.window.dis.H3K9ac,pmdlist)),sum(overlapsAny(chipseq.window.dis.H3K9me3,pmdlist)),sum(overlapsAny(chipseq.window.dis.H3K4me1,pmdlist)),sum(overlapsAny(chipseq.window.dis.H3K4me3,pmdlist)),sum(overlapsAny(chipseq.window.dis.H3K27me3,pmdlist)),sum(overlapsAny(chipseq.window.dis.H3K36me3,pmdlist)))
PMD.chip.heal=c(sum(overlapsAny(chipseq.window.heal.H3K27ac,pmdlist)),sum(overlapsAny(chipseq.window.heal.H3K9ac,pmdlist)),sum(overlapsAny(chipseq.window.heal.H3K9me3,pmdlist)),sum(overlapsAny(chipseq.window.heal.H3K4me1,pmdlist)),sum(overlapsAny(chipseq.window.heal.H3K4me3,pmdlist)),sum(overlapsAny(chipseq.window.heal.H3K27me3,pmdlist)),sum(overlapsAny(chipseq.window.heal.H3K36me3,pmdlist)))

unPMD.chip.dis=c(sum(!overlapsAny(chipseq.window.dis.H3K27ac,pmdlist)),sum(!overlapsAny(chipseq.window.dis.H3K9ac,pmdlist)),sum(!overlapsAny(chipseq.window.dis.H3K9me3,pmdlist)),sum(!overlapsAny(chipseq.window.dis.H3K4me1,pmdlist)),sum(!overlapsAny(chipseq.window.dis.H3K4me3,pmdlist)),sum(!overlapsAny(chipseq.window.dis.H3K27me3,pmdlist)),sum(!overlapsAny(chipseq.window.dis.H3K36me3,pmdlist)))
unPMD.chip.heal=c(sum(!overlapsAny(chipseq.window.heal.H3K27ac,pmdlist)),sum(!overlapsAny(chipseq.window.heal.H3K9ac,pmdlist)),sum(!overlapsAny(chipseq.window.heal.H3K9me3,pmdlist)),sum(!overlapsAny(chipseq.window.heal.H3K4me1,pmdlist)),sum(!overlapsAny(chipseq.window.heal.H3K4me3,pmdlist)),sum(!overlapsAny(chipseq.window.heal.H3K27me3,pmdlist)),sum(!overlapsAny(chipseq.window.heal.H3K36me3,pmdlist)))

PMD.chip.dis2=PMD.chip.dis/sum(PMD.chip.dis)
PMD.chip.heal2=PMD.chip.heal/sum(PMD.chip.heal)
unPMD.chip.dis2=unPMD.chip.dis/sum(unPMD.chip.dis)
unPMD.chip.heal2=unPMD.chip.heal/sum(unPMD.chip.heal)


plotdat= rbind(cbind(cbind(PMD.chip.dis2,c("H3K27ac","H3K9ac","H3K9me3","H3K4me1","H3K4me3","H3K27me3","H3K36me3")),"CLL-PMD"),cbind(cbind(PMD.chip.heal2,c("H3K27ac","H3K9ac","H3K9me3","H3K4me1","H3K4me3","H3K27me3","H3K36me3")),"Healthy-PMD"),cbind(cbind(unPMD.chip.dis2,c("H3K27ac","H3K9ac","H3K9me3","H3K4me1","H3K4me3","H3K27me3","H3K36me3")),"CLL-notPMD"),cbind(cbind(unPMD.chip.heal2,c("H3K27ac","H3K9ac","H3K9me3","H3K4me1","H3K4me3","H3K27me3","H3K36me3")),"Healthy-notPMD"))
plotdat=data.frame(V1=as.numeric(plotdat[,1]),V2=as.factor(plotdat[,2]),V3=as.factor(plotdat[,3]))

ggplot(plotdat, aes(x=V3,y=(V1))) + geom_bar(stat="identity", aes(fill=V2),width=0.7)  +ylab("")+xlab("") +theme_bw()+theme(legend.justification=c(1,1), legend.position=c(1,1))
ggsave("PMD-overlap-with-chip-seq-peaks.pdf",width=8,height=10)

#total number
plotdat= rbind(cbind(cbind(PMD.chip.dis,c("H3K27ac","H3K9ac","H3K9me3","H3K4me1","H3K4me3","H3K27me3","H3K36me3")),"CLL-PMD"),cbind(cbind(PMD.chip.heal,c("H3K27ac","H3K9ac","H3K9me3","H3K4me1","H3K4me3","H3K27me3","H3K36me3")),"Healthy-PMD"),cbind(cbind(unPMD.chip.dis,c("H3K27ac","H3K9ac","H3K9me3","H3K4me1","H3K4me3","H3K27me3","H3K36me3")),"CLL-notPMD"),cbind(cbind(unPMD.chip.heal,c("H3K27ac","H3K9ac","H3K9me3","H3K4me1","H3K4me3","H3K27me3","H3K36me3")),"Healthy-notPMD"))
plotdat=data.frame(V1=as.numeric(plotdat[,1]),V2=as.factor(plotdat[,2]),V3=as.factor(plotdat[,3]))

ggplot(plotdat, aes(x=V3,y=(V1))) + geom_bar(stat="identity", aes(fill=V2),width=0.7)  +ylab("")+xlab("") +theme_bw()+theme(legend.justification=c(1,1), legend.position=c(1,1))
ggsave("PMD-overlap-with-chip-seq-peaks-total-counts.pdf",width=8,height=10)
