###########################################################################
###          Rsript: Annotation of differential ATAC-seq peaks          ###
###              (author: Lara Klett    date: 2017-10-12)               ###
###########################################################################

########
writeLines("\n***\nCreated on 2017-10-12\n\nauthor: Lara Klett\n***\n"); writeLines("Hello User!"); writeLines("\n Assigns peaks to certain functional states (takes output produced by annotatebed in script Annotation_DiffBind_ATAC_regions_ChromHMM.sh) ...\n")
########


folder="your_path"
setwd(folder)



#################
#### all
#################
set="all"
anno=read.table(file="DBA_DiffBind_ATAC_Regions_fdr0.0038-fold1_-2.3_fold2_1.13_all_forHomer_annotated-ChromHMM.bed")
colnames(anno)=c("chromosom", "start", "end","peaknumber","NA","strand", "state01", "state02","state03","state04","state05","state06","state07","state08","state09","state10","state11","state12")

tab1 = table(colnames(anno[,c(7:18)])[apply(anno[,c(7:18)],1,which.max)]) #if overlap percentage is equal for several states it is counted in following order state1 before 2, 2 before3,...
barplot(cbind(tab1,c(rep(0,times=12))), cex.names=0.6, col=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow","snow2"), main="Annotation of differential accessible ATAC-seq regions", beside=F, horiz=F, names=c(set,""))
legend("topright", fill=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow","snow2"), legend=c("Enhancer (genic)","Transcription","Repressed (ZNF)","Repressed (HET)","Repressed","Repressed (polycomb)","Bivalent","Enhancer (active2)","Enhancer (active1)", "TSS (active)","Enhancer (weak)","Quiescent"))

#################
#### gain
#################

set="gain"
annoGain=read.table(file=paste0("DBA_DiffBind_ATAC_Regions_fdr0.0038-fold1_-2.3_fold2_1.13_",set,"_forHomer_annotated-ChromHMM.bed"))
colnames(annoGain)=c("chromosom", "start", "end","peaknumber","NA","strand", "state01", "state02","state03","state04","state05","state06","state07","state08","state09","state10","state11","state12")

tab2 = table(colnames(annoGain[,c(7:18)])[apply(annoGain[,c(7:18)],1,which.max)]) #if overlap percentage is equal for several states it is counted in following order state1 before 2, 2 before3,...
barplot(cbind(tab2,c(rep(0,times=12))), cex.names=0.6, col=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow","snow2"), main="Annotation of differential accessible ATAC-seq regions", beside=F, horiz=F, names=c(set,""))
legend("topright", fill=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow","snow2"), legend=c("Enhancer (genic)","Transcription","Repressed (ZNF)","Repressed (HET)","Repressed","Repressed (polycomb)","Bivalent","Enhancer (active2)","Enhancer (active1)", "TSS (active)","Enhancer (weak)","Quiescent"))

#################
#### loss
#################

set="loss"
annoLoss=read.table(file=paste0("DBA_DiffBind_ATAC_Regions_fdr0.0038-fold1_-2.3_fold2_1.13_",set,"_forHomer_annotated-ChromHMM.bed"))
colnames(annoLoss)=c("chromosom", "start", "end","peaknumber","NA","strand", "state01", "state02","state03","state04","state05","state06","state07","state08","state09","state10","state11","state12")

tab3 = table(colnames(annoLoss[,c(7:18)])[apply(annoLoss[,c(7:18)],1,which.max)]) #if overlap percentage is equal for several states it is counted in following order state1 before 2, 2 before3,...
barplot(cbind(tab3,c(rep(0,times=12))), cex.names=0.6, 
        col=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow","snow2"), 
        main="Annotation of differential accessible ATAC-seq regions", beside=F, horiz=F, names=c(set,""))
legend("topright", fill=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow","snow2"), legend=c("Enhancer (genic)","Transcription","Repressed (ZNF)","Repressed (HET)","Repressed","Repressed (polycomb)","Bivalent","Enhancer (active2)","Enhancer (active1)", "TSS (active)","Enhancer (weak)","Quiescent"))



#################
#### combined (no quiescent state)
#################
tabcomb = cbind(tab1,tab2,tab3)
colnames(tabcomb)=c("all","gained","lost")
tabcomb_red=tabcomb[1:11,]

pdf("Annotation_DiffBind_ATAC_regions_ChromHMM_barplot_noQuiescent.pdf", width=5, height=8)
barplot(tabcomb_red, cex.names=2, ylim=c(0,25000), 
        col=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow"), 
        main="Annotation of differential accessible ATAC-seq regions", beside=F, horiz=F, names=colnames(tabcomb_red))
legend("topright", cex= 0.8,fill=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow"), legend=c("Enhancer (genic)","Transcription","Repressed (ZNF)","Repressed (HET)","Repressed","Repressed (polycomb)","Bivalent","Enhancer (active2)","Enhancer (active1)", "TSS (active)","Enhancer (weak)"))

dev.off()

tabcomb_red2=tabcomb_red
tabcomb_red2[,1]=tabcomb_red2[,2]+tabcomb_red2[,3]

pdf("Annotation_DiffBind_ATAC_regions_ChromHMM_barplot_added_noQuiescent.pdf", width=5, height=8)

barplot(tabcomb_red2, cex.names=2, ylim=c(0,25000), 
        col=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow"), 
        main="Annotation of \ndifferential accessible ATAC-seq regions", beside=F, horiz=F, names=colnames(tabcomb_red))
legend("topright", cex= 0.8,fill=c("green3","darkgreen","lightblue","blue","lightgrey","darkgrey","brown","gold2","darkorange","red","yellow"), legend=c("Enhancer (genic)","Transcription","Repressed (ZNF)","Repressed (HET)","Repressed","Repressed (polycomb)","Bivalent","Enhancer (active2)","Enhancer (active1)", "TSS (active)","Enhancer (weak)"))

dev.off()
