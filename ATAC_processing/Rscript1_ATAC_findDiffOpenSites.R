###########################################################################
###           Rsript: Find fold change threshold for DiffBind           ###
###              (author: Lara Klett    date: 2017-04-04)               ###
###########################################################################

########
writeLines("\n***\nCreated on 2017-04-04\n\nauthor: Lara Klett\n***\n"); writeLines("Hello User!"); writeLines("\nFind fold change threshold for DiffBind ...\n")
########


folder="your_path"

fdrvalue=1

fname=paste0("DBA_DiffBind_all_fdr",fdrvalue)

#############################################
### Read in Diff bound sites list
setwd(folder)
diffs=read.csv(file=paste0(fname,".csv"), header=T, sep=",")
write.table(diffs, file=paste0(fname,".bed"), row.names=FALSE, quote=FALSE,col.names=FALSE, sep = "\t")

############################
### Prepare info matrix  ###
############################
info=read.csv("sample-names.txt", sep="\t", header=F)
head(info)

infomat = sapply(strsplit(colnames(diffs)[10:dim(diffs)[2]],"_"),unlist)
rownames (infomat) = c("individual","replicate")
colnames (infomat) = colnames (diffs[10:dim(diffs)[2]])


######################
### Prepare Matrix ###
######################

# Remove first 9 columns
diffs1=diffs[,10:dim(diffs)[2]]

# Remove Rep1 and Rep2 from column name
type = gsub("_Rep[1-3]","",colnames (diffs1))
colnames (diffs1)= type

# Combine replicates (addition)
diffs2 = sapply(unique(type),function(i){rowSums(diffs1[type %in% i])})

# Display how often a colname appears
table=table(colnames (diffs1))

### Order colname frequency table according to matrix
table3=as.vector(table[unique(type)[]])

### Divide all columns with corresponding number of combined replicates to ensure correct normalization
diffs3=t(t(diffs2)/table3)

### Calculate fold change correctly with replicates
diffs3H = diffs3[,colnames(diffs3) %in% info[info[,3]=="H",1]]
diffs3CLL = diffs3[,colnames(diffs3) %in% info[info[,3]=="CLL",1]]

Conc_control = log2(rowSums(diffs3H)/ncol(diffs3H))
Conc_CLL = log2(rowSums(diffs3CLL)/ncol(diffs3CLL))

Conc = log2(rowSums(diffs3)/ncol(diffs3))

Fold=log2(2 ** Conc_CLL/2 ** Conc_control)

######

diffs4 = cbind(diffs[,1:3],Conc,Conc_CLL,Conc_control,Fold,diffs3)



##########################################
### Fit three gaussian distributions to fold change distribution
#

ta=as.data.frame(table(round(x=(diffs4$Fold), digits=2)))
ta[,1]=as.numeric(levels(ta$Var1))
tacum=cbind(ta$Var1,cumsum(ta$Freq))

write.csv(ta,paste0("DBA_DiffBind_all_fdr",fdrvalue,"_number-diffRegs-vs-absolute-fold-change-threshold.csv"))
write.csv(tacum,paste0("DBA_DiffBind_all_fdr",fdrvalue,"_number-diffRegs-with-this-fold-threshold_vs_absolute-fold-change-threshold.csv"))


name="three gaussian distributions"

x <- ta[,"Var1"]
y <- ta[,"Freq"] ##range(ta$Freq)[2]=1162


### define fit function
fit <- function(x)
  function(a1,width1,pos1,a2,width2,pos2,a3,width3,pos3)
    a1*exp(-((x-pos1)/width1)^2) + a2*exp(-((x-pos2)/width2)^2) + a3*exp(-((x-pos3)/width3)^2)

### perform fit
ctrl = list(maxiter = 500000, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = TRUE)

init = list(a1=200, width1=1, pos1=0, a2=200, width2=2, pos2=2, a3=150, width3=2, pos3=-3)
lw = list(a1=0, width1=0, pos1=-1, a2=0, width2=0, pos2=0,a3=0, width3=0, pos3=-5)
up = list(a1=300, width1=2, pos1=1, a2=500, width2=30, pos2=12, a3=150, width3=30, pos3=0)

fitres <- nls(y ~ fit(x)(a1,width1,pos1,a2,width2,pos2,a3,width3,pos3), alg="port", start=init, control=ctrl, lower=lw, upper=up)
### plot data and fit functions
plot(x,y, type="p", pch="x", main=paste("Fit result for",name), xlab="Fold change", ylab="Frequency")
lines(x, predict(fitres), col="blue", lwd=4)

### output fit result
print(summary(fitres))


a1=summary(fitres)$coef[1]
width1=summary(fitres)$coef[2]
pos1=summary(fitres)$coef[3]
a2=summary(fitres)$coef[4]
width2=summary(fitres)$coef[5]
pos2=summary(fitres)$coef[6]
a3=summary(fitres)$coef[7]
width3=summary(fitres)$coef[8]
pos3=summary(fitres)$coef[9]

plot(ta[,"Var1"],ta[,"Freq"], main=paste("Fit result for",name), xlab="Fold change", ylab="Frequency")

points(x,a1*exp(-((x-pos1)/width1)^2),col="grey")
points(x,a2*exp(-((x-pos2)/width2)^2),col="green")
points(x,a3*exp(-((x-pos3)/width3)^2),col="red")
lines(x, predict(fitres), col="blue", lwd=4)

### Display formula of fitting functions
print("Fitted function  for regions with no change");paste0(a1,"*exp(-((x-",pos1,")/",width1,")^2)")
print("Fitted function for regions with positive fold change");paste0(a2,"*exp(-((x-",pos2,")/",width2,")^2)")
print("Fitted functionfor regions with negative fold change");paste0(a3,"*exp(-((x-",pos3,")/",width3,")^2)")


# determine intersection of the distributions manually (e.g. with Wolfram alpha)
thresh1=-2.30328
thresh2=1.13122

foldchange1=thresh2
foldchange2=abs(thresh1)

colv=c("#225ea8", "#41b6c4", "#7fcdbb")

abline(v=thresh1, lwd=2)
text(x=thresh1-3, y=350, labels=paste0("fold change \nthreshold 1=",round(thresh1, digits=5)))
abline(v=thresh2, lwd=2)
text(x=thresh2+4, y=350, labels=paste0("fold change \nthreshold 2=",round(thresh2, digits=5)))
text(x= 0.0000001,y= c(130000,150000), labels=c(paste0("pos fold change threshold = ",round(foldchange1,digits=2)),paste0("neg fold change threshold = ",round(foldchange2,digits=2))),col=colv, pch=1, cex=2)

#################      Plot in pdf      #################


threshfold1=round(thresh1,digits=2)
threshfold2=round(thresh2,digits=2)


colv=c("#225ea8", "#41b6c4", "#7fcdbb")

pdf(paste0(folder, "/DBA_DiffBind_ATAC_Regions_fdr",fdrvalue,"-fold1_",threshfold1,"_fold2_",threshfold2,"_fold-change-threshold_gaussian-fit_LaraRecalcFold_lines.pdf"))
plot(ta[,"Var1"],ta[,"Freq"], type="n", pch="x", main=paste("Fit result for",name), xlab="Fold change", ylab="Frequency", cex.lab=1.5, cex.axis=1.5)

lines(x,a1*exp(-((x-pos1)/width1)^2),col="grey", lwd=4)
lines(x,a2*exp(-((x-pos2)/width2)^2),col="green", lwd=4)
lines(x,a3*exp(-((x-pos3)/width3)^2),col="red", lwd=4)

lines(x, predict(fitres), col="black", lwd=6)


abline(v=thresh1, lwd=2)
text(x=thresh1-3, y=350, labels=paste0("fold change \nthreshold 1=",round(thresh1, digits=5)), cex=1.2)
abline(v=thresh2, lwd=2)
text(x=thresh2+4, y=350, labels=paste0("fold change \nthreshold 2=",round(thresh2, digits=5)), cex=1.2)
text(x= 0.0000001,y= c(130000,150000), labels=c(paste0("pos fold change threshold = ",round(foldchange1,digits=2)),paste0("neg fold change threshold = ",round(foldchange2,digits=2))),col=colv, pch=1, cex=2)

dev.off()

#
##########################






#############################################
### Calculate inflexion point of fdr scan with two different fold change thresholds for negative and positive fold changes
#
#



library(DiffBind)
ta2 <- dba.load('dba.count')
ta3 <- dba.load('dba.contrast')
ta4 <- dba.load('dba.analyze')



foldchange1=thresh2
foldchange2=abs(thresh1)

### Create matrix containing calculated overall and group specific "Conc", and diffbind calculated p.value and fdr
diffs5 = cbind(diffs[,1:3],Conc,Conc_CLL,Conc_control,Fold,diffs$p.value,diffs$FDR,diffs3)
colnames(diffs5)[8]="p.value"
colnames(diffs5)[9]="FDR"

diffs5_reps = cbind(diffs[,1:3],Conc,Conc_CLL,Conc_control,Fold,diffs$p.value,diffs$FDR,diffs[10:ncol(diffs)])
colnames(diffs5_reps)[8]="p.value"
colnames(diffs5_reps)[9]="FDR"

### Check how many differential sites are gained and lost with chosen fold thresholds
dim(diffs5_reps[diffs5_reps$Fold>=thresh2,]) [1] # result: 24369
dim(diffs5_reps[diffs5_reps$Fold<=thresh1,]) [1] # result: 18848

### Save subsets as 
gained = diffs5[diffs5$Fold>=thresh2,]
lost = diffs5[diffs5$Fold<=thresh1,]

gained_reps = diffs5_reps[diffs5_reps$Fold>=thresh2,]
lost_reps = diffs5_reps[diffs5_reps$Fold<=thresh1,]



diffs6 = rbind(gained,lost)
diffs6_reps = rbind(gained_reps,lost_reps)





fold = matrix(nrow=0, ncol=6)
colnames(fold) = c("fdr","DiffBoundPeaks","DiffBoundPeaks_pos", "fold_threshold1","DiffBoundPeaks_neg", "fold_threshold2")


fdr = 0.00000000001

repeat{
  
  ### Report number of differentially gained sites with chosen fdr threshold
  ndiffpeaks2=nrow(gained[gained$FDR<=fdr,])
  
  ### Report number of differentially gained sites with chosen fdr threshold
  ndiffpeaks1=nrow(lost[lost$FDR<=fdr,])
  
    
  ########
  
  ndiffpeaks=ndiffpeaks1+ndiffpeaks2
  ndiffpeaks
    
  fold=rbind(fold,c(fdr,ndiffpeaks,ndiffpeaks1,thresh2,ndiffpeaks2,thresh1))
  
  fdr = fdr * 2
    
  if(fdr>1){
    break
  }
}
  

foldnew=cbind(fold,c(diff(fold[,2]),1))
colnames(foldnew)[dim(foldnew)[2]]="slope"
numeric_fdr_thresh=((foldnew[grep(max(diff(fold[,2])), x=foldnew[,"slope"]),1])+(foldnew[grep(max(diff(fold[,2])), x=foldnew[,"slope"])+1,1]))/2
numeric_fdr_thresh2=(log2(foldnew[grep(max(diff(fold[,2])), x=foldnew[,"slope"]),1])+log2(foldnew[grep(max(diff(fold[,2])), x=foldnew[,"slope"])+1,1]))/2
numeric_fdr_thresh3=2 ** numeric_fdr_thresh2
numeric_fdr_thresh3




### Data.frame with all differentially open sites
diffs7=diffs6[diffs6$FDR<=numeric_fdr_thresh3,]
diffs7_reps=diffs6_reps[diffs6_reps$FDR<=numeric_fdr_thresh3,]

### Data.frames with all gained or lost differentially open sites
gained=diffs7[diffs7$Fold>0,]
gained_reps=diffs7_reps[diffs7_reps$Fold>0,]

lost=diffs7[diffs7$Fold<0,]
lost_reps=diffs7_reps[diffs7_reps$Fold<0,]


#########
writeLines("\n#############################################\n"); writeLines("Your numeric fdr threshold for DiffBind is:"); writeLines(paste0("\n",numeric_fdr_thresh3,"\n")); writeLines("\n#############################################\n")
#########

######### number of all differential sites
writeLines("\n#############################################\n"); writeLines(paste0("\n",nrow(diffs7)," sites are differentially open\n"));writeLines(paste0("\nof this\n",nrow(gained)," sites are gained")); writeLines(paste0("\n",nrow(lost)," sites are lost\n")); writeLines("\n#############################################\n")
#########


### Plot fdr scan
thresh=round(numeric_fdr_thresh3,digits=5)

pdf(paste0("Result_diffbind_fdr-scan-th",thresh,"-fold1_", thresh1, "_fold2_",thresh2,"_LaraFoldCalc.pdf"))
par(mar=c(5.1,6.1,4.1,2.1))

   plot((fold[,1]), fold[,2], xlab="FDR", ylab="Number of differentially accessible regions", col=colv[1], log="x",ylim=c(0,ndiffpeaks+1/5*ndiffpeaks), cex.axis=2, cex.lab=2, cex.main=2,main="FDR scan")
   abline(v=numeric_fdr_thresh3, lty=2)
   legend("topleft", legend=paste0("neg. fold change threshold1 = ",round(thresh1,digits=2),"\npos. fold change threshold2 = ",round(thresh2,digits=2)),col=colv, pch=1, cex=1.5)
   text(x=numeric_fdr_thresh3+0.01, y=ndiffpeaks-40000, labels=paste0("fdr threshold=\n",round(numeric_fdr_thresh3,digits=5)))


dev.off()



pdf(paste0("Result_diffbind_fdr-scan-th",thresh,"-fold1_", thresh1, "_fold2_",thresh2,"_LaraFoldCalc_derivative.pdf"))
par(mar=c(5.1,6.1,4.1,2.1))

   # Shift slope points to middle between points they were calculated from as difference in y-value
   mat=cbind(log2(foldnew[-37,1]),log2(foldnew[-1,1]))
   mat2=2 ** (rowSums(mat)/2)

   plot(mat2,foldnew[-37,"slope"],log="x", main="Slope of fdr scan\n(Determination of inflexion point\nfrom maximum turning point of derivative)", xlab="fdr threshold", ylab="slope")
   abline(v=numeric_fdr_thresh3)
   text(x=numeric_fdr_thresh3, y=500, labels=paste0("fdr threshold=",round(numeric_fdr_thresh3,digits=5)))

dev.off()

#############################################
### Save "fold" in a file

nameFDR = paste0(fname,"_number-diffRegs-vs-fdr-threshold-fold1_", thresh1, "_fold2_",thresh2)
write.csv(foldnew, file=paste0(nameFDR,"_inflexionPoint.csv"), row.names=TRUE, quote=FALSE)

#
#
#######################


##############################
### PCA (only diffbound sites)
##############################
info=read.csv("sample-names.txt", sep="\t", header=F)
colnames(info)=c("individual","ID","health_status", "mut_status")

pcl=c("red","orange","grey40")

col=pcl[factor(info[,4])]
info[,5]=col
info3=info[match(colnames(diffs7[,-c(1:9)]), info[,1]),]


library(ggbiplot)
library("ggplot2")
## install.packages("ggrepel") 
library("ggrepel")

###Calculate Principle components
pca.diffs7 <- prcomp(t(diffs7[,-c(1:9)]), scale. = T)
summary(pca.diffs7)
# Plot how much variance is explained by which PC
plot(pca.diffs7, type = "l")

p <- ggbiplot::ggbiplot(pca.diffs7,
                        #choices=c(3,4), # plot different PCs
                        var.axes = F, var.scale = -1, varname.size = 2, varname.adjust = 4, varname.abbrev = F, 
                        circle = F, ellipse =F, labels.size=10) + 
  theme_bw()+
  geom_point(aes(size=2,colour = factor(info3[match(colnames(diffs7[,-c(1:9)]), info3[,1]),"mut_status"]))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #geom_text_repel(aes(label=colnames(diffs7[,-c(1:9)])), size = 5, box.padding=unit(0.47, "lines"), point.padding = unit(0.5, "lines")) +
  scale_color_manual(breaks = c("1", "2", "3"),labels = c("CLL-mut", "CLL-unmut", "H"), values = pcl, c(0.1,1))


pdf(paste0("Result_diffbind_ATAC_PCA_CLLmu-CLLunmut-healthy_onlyDiffSites.pdf"), height=5)
print(p + xlim(-2.7, 2.7) + ylim(-2.7, 2.7))
dev.off()







###################################
### Save differentially bound peaks with different fold thresholds for up and down regulated in csv file (as well as separated by gain)

thresh=round(numeric_fdr_thresh3,digits=5)

threshfold1=round(thresh1,digits=2)
threshfold2=round(thresh2,digits=2)


write.csv(lost_reps, file=paste0("DBA_DiffBind_ATAC_Regions_fdr",thresh,"-fold1_", threshfold1, "_fold2_",threshfold2,"_loss.csv"), row.names=FALSE, quote=FALSE)
write.table(lost_reps, file=paste0("DBA_DiffBind_ATAC_Regions_fdr",thresh,"-fold1_", threshfold1, "_fold2_",threshfold2,"_loss.bed"), row.names=FALSE, quote=FALSE,col.names=FALSE, sep = "\t")

write.csv(gained_reps, file=paste0("DBA_DiffBind_ATAC_Regions_fdr",thresh,"-fold1_", threshfold1, "_fold2_",threshfold2,"_gain.csv"), row.names=FALSE, quote=FALSE)
write.table(gained_reps, file=paste0("DBA_DiffBind_ATAC_Regions_fdr",thresh,"-fold1_", threshfold1, "_fold2_",threshfold2,"_gain.bed"), row.names=FALSE, quote=FALSE,col.names=FALSE, sep = "\t")

######

# Combine both data sets and save as all differential peaks
all_reps=rbind(gained_reps,lost_reps)
all=rbind(gained,lost)
write.table(all_reps, file=paste0("DBA_DiffBind_ATAC_Regions_fdr",thresh,"-fold1_", threshfold1, "_fold2_",threshfold2,"_all.bed"), row.names=FALSE, quote=FALSE,col.names=FALSE, sep = "\t")



#
#
#####################################


##############################################
#### Plot correlation (manually calculated fold change)
##
library(gplots) #install.packages("gplots", dependencies = TRUE)

### Clustered (patients) plotting of data:
library("amap", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")

library("corrplot", lib.loc="/Users/klett/Library/R/3.3/library")
library(RColorBrewer)   #install.packages("RColorBrewer", dependencies = TRUE)


all2_reps=all_reps[,10:ncol(all_reps)]
all2=all[,10:ncol(all)]




### Colors:
col1 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7", # reds
                           "#FFFFFF",                                            # white
                           "#FFFFE5","#D9F0A3", "#ADDD8E", "#41AB5D", "#238443", "#004529")) # greens( similar to: brewer.pal(9, "YlGn") )




#############################################
### Correlation with all regions

# Rename columns of "diffs3
diffs3=t(t(diffs2)/table3)
colnames(diffs3)=gsub("P","CLL",info3[match(info3$individual,colnames(diffs3)),2])

diffs_cor=diffs[10:ncol(diffs)]

names=gsub("*\\_Rep[0-9]","",colnames(diffs_cor))
names2<-info3[match(names, info3$individual),"ID"]
vec=gsub("P","CLL",names2)
colnames(diffs_cor)=paste0(vec,"_", gsub(".*_","",colnames(diffs_cor))) # replace everything behind slash



###### Cluster according to correlation

#Replicates combined
corVals = round(cor(diffs3, method="spearman"),2)
pdf(paste0(folder, "/DiffBind1_CorrelationHeatmapOccupancy_LaraRecalcFold.pdf"))
corrplot(corVals, order = "hclust", hclust.method = "average", addrect = 2, tl.col = "black",col=col1(200))
mtext("spearman correlation - all ATAC-seq regions",line=2, cex=2) 
dev.off()

# Replicates separate
corVals = round(cor(diffs_cor, method="spearman"),2)
pdf(paste0(folder, "/DiffBind1_CorrelationHeatmapOccupancy_Reps_LaraRecalcFold.pdf"))
corrplot(corVals, order = "hclust", hclust.method = "average", addrect = 2, tl.col = "black",col=col1(200),tl.cex=0.6)
mtext("spearman correlation - all ATAC-seq regions",line=2, cex=2) 
dev.off()

# Revert renaming of columns of "diffs3"
diffs3=t(t(diffs2)/table3)


#############################################
### Correlation with all differential regions

### Rename samples

# Save original sample names
originalcolnames_all2_reps=colnames(all2_reps)

# Change sample names according to info3 ID list
names=gsub("*\\_Rep[0-9]","",colnames(all2_reps))
names2<-info3[match(names, info3$individual),"ID"]
vec=gsub("P","CLL",names2)
colnames(all2_reps)=paste0(vec,"_", gsub(".*_","",colnames(all2_reps))) # replace everything behind slash


### Correlation with Replicates separated

corVals = round(cor(all2_reps, method="spearman"),2)

pdf(paste0(folder, "/DiffBind1_CorrelationHeatmapOccupancy_fdr",thresh,"-fold1_",thresh1,"_fold2_",thresh2,"_Reps_LaraRecalcFold.pdf"))
corrplot(corVals, order = "hclust", hclust.method = "average", addrect = 2, tl.col = "black", tl.cex=0.6,col=col1(200))
mtext("spearman correlation - differentially open regions",line=2, cex=1.5) 
dev.off()

# Restore original sample names
colnames(all2_reps)=originalcolnames_all2_reps

##############################################
### Correlation with Replicates combined 

### Rename samples

# Save original sample names
originalcolnames_all2=colnames(all2)

# Change sample names according to info3 ID list
names=gsub("*\\_Rep[0-9]","",colnames(all2))
names2<-info3[match(names, info3$individual),"ID"]
vec=gsub("P","CLL",names2)
colnames(all2)=vec # replace everything behind slash


corVals = round(cor(all2, method="spearman"),2)

pdf(paste0(folder, "/DiffBind1_CorrelationHeatmapOccupancy_fdr",thresh,"-fold1_",thresh1,"_fold2_",thresh2,"_LaraRecalcFold.pdf"))
corrplot(corVals, order = "hclust", hclust.method = "average", addrect = 2, tl.col = "black", tl.cex=0.8,col=col1(200))
mtext("spearman correlation - differentially open regions",line=2, cex=1.5) 
dev.off()

# Restore original sample names
colnames(all2)=originalcolnames_all2


#
#
#####################################
