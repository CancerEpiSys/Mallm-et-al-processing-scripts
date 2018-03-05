###########################################################################
###                    Rsript: Cluster ATAC-seq data                    ###
###              (author: Lara Klett    date: 2017-04-04)               ###
###########################################################################

########
writeLines("\n***\nCreated on 2017-07-24\n\nauthor: Lara Klett\n***\n"); writeLines("Hello User!"); writeLines("\nCluster ATAC-seq data (produced by diffbind) ...\n")
########


folder="my_path"
fdrvalue=1
fname=paste0("DBA_DiffBind_all_fdr",fdrvalue)

#############################################
### Read in Diff ATAC sites list
#############################################
setwd(folder)
diffs=read.csv(file=paste0(fname,".csv"), header=T, sep=",")

info=read.csv("sample-names.txt", sep="\t", header=F)
colnames(info)=c("individual","ID","health_status", "mut_status")
############################
### Prepare info matrix  ###
############################

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
sites=nrow(diffs3)

#
#
#########################

##############################
### Load needed libraries  ###
##############################

library(gplots) #install.packages("gplots", dependencies = TRUE)
library(RColorBrewer)   #install.packages("RColorBrewer", dependencies = TRUE)

library("amap", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")


##############################
### Define cluster method  ###
##############################


clust_meth="spearman" #"spearman" or "pearson"
clust_link="average"  # "average" or "complete" or "single" or ...

##############################
###       Clustering       ###
##############################


hComCor <-hcluster(t(log2(diffs3[1:sites,])), method=clust_meth, lin=clust_link)

### plot dendrogram
plot(as.dendrogram(hComCor))

### print clusters
n=2
cutree(hComCor,n)





###########################################
###########################################
###      Produce Dendrogram Plots       ###
###########################################
###########################################


#install.packages('dendextend')
library(dendextend)


##########################
###  Add colored bar  ###
##########################
pcl=c("red","orange","grey40")

col=pcl[factor(info[,4])]
info[,5]=col
info3=info[match(colnames(diffs3), info[,1]),]

#### Plot dendrogram
op <- par(mfrow = c(1, 1))
par(mar=c(15,5,3,3))
dend=as.dendrogram(hComCor, dLeaf=10)
dend %>% set("labels") %>% plot
colored_bars(colors = info3[hComCor$order,5], rowLabels="ATAC-seq")


### Save plotted dendrogram in pdf
pdf(paste0("Result_diffbind_ATAC_clusteringTree_",clust_link,"Link_",clust_meth,"_CLL-healthy.pdf"), height=5)

par(mar=c(15,5,3,3))
dend=as.dendrogram(hComCor, dLeaf=10)
dend %>% set("labels",paste0(info3[hComCor$order,2],"-B-",info3[hComCor$order,1])) %>% plot
colored_bars(colors = info3[hComCor$order,5], rowLabels="ATAC-seq", y_scale=0.1, y_shift=-0.3)

dev.off()


