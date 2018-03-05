
########################################
#                Dotplot               #
#              2017-09-15              #
########################################


########
writeLines("\n***\nCreated on 2017-07-31\n\nauthor: Lara Klett\n***\n"); writeLines("Hello User!"); writeLines("\nPlotting Homer results ...\n")
########

library("ggplot2")
## install.packages("ggrepel") 
library("ggrepel")


##############################################################
### PLOT 1 (loss at enhancers)
##############################################################


folder="your_path"

name="knownResults_ATAC_loss_enhancers"

ttl="Lost ATAC signal at enhancers"

##########################

setwd(folder)

### Read in file
ta=read.csv("knownResults.txt", sep="\t")

Gene=gsub("\\(.*","",ta$Motif.Name)
Perc_of_Targets=as.numeric(gsub("%","",ta$X..of.Target.Sequences.with.Motif))/100
Perc_of_Background=as.numeric(gsub("%","",ta$X..of.Background.Sequences.with.Motif))/100
Enrichment=Perc_of_Targets/Perc_of_Background

ta1=cbind(Gene,ta,Perc_of_Targets,Perc_of_Background,Enrichment,log10(Enrichment))
##########################
### Renaming to achieve consistent naming (e.g. E2A is identical to Irf3)
ta1$Gene=as.vector(ta1$Gene)

ta1$Gene=gsub("EWS:ERG-fusion","ERG",ta1$Gene)
ta1$Gene=gsub("RUNX-AML","RUNX",ta1$Gene)
ta1[grep("NFkB",ta1$Gene),]$Gene="NFkB"
ta1$Gene=gsub("Mef2c","MEF2C",ta1$Gene)
ta1$Gene=gsub("Mef2d","MEF2D",ta1$Gene)
ta1$Gene=gsub("Oct","OCT",ta1$Gene)
ta1$Gene=gsub("Etv","ETV",ta1$Gene)
ta1$Gene=gsub("MafK","MAFK",ta1$Gene)
ta1$Gene=gsub("Brn1","BRN1",ta1$Gene)
ta1$Gene=gsub("SpiB","SPIB",ta1$Gene)
ta1$Gene=gsub("Nrf2","NRF2",ta1$Gene)

ta1$Gene=gsub("Fra","FRA",ta1$Gene)
ta1$Gene=gsub("Fosl2","FOSL2",ta1$Gene)

ta1$Gene=gsub("Atf3","ATF3",ta1$Gene)
ta1$Gene=gsub("Ets1","ETS1",ta1$Gene)
ta1$Gene=gsub("Jun","JUN",ta1$Gene)
##########################



############################################################################
### Color according to similar motifs
############################################################################

### Combine similar motifs as same motif type and save in vector "Motif.type"
ta2=ta1[1:35,]
Motif.type=ta2$Gene

Motif.type[grep("EBF",Motif.type)]="EBF"
Motif.type[grep("BORIS",Motif.type)]="CTCF"
Motif.type[grep("RUNX",Motif.type)]="RUNX"
Motif.type[Motif.type %in% c("JUNB","JUN-AP1","FOSL2","FRA2","FRA1","ATF3","AP-1","Bach2","Bach1","NRF2","NF-E2","MAFK","BATF")]="JUN"
Motif.type[Motif.type %in% c("ERG","PU.1","ELF3","ETS1-distal","ETV2","SPIB")]="ETS"
Motif.type[Motif.type %in% c("OCT2","OCT4","BRN1")]="OCT"
Motif.type[Motif.type %in% c("MEF2C","MEF2D")]="MEF"

ta2=cbind(ta1[1:35,],Motif.type)

### Plot

pdf(paste0(name,"_dotplot5_names_top35_motif-type_coloring.pdf"))

#my_bluegreen = rbind(brewer.pal(n = 9, "Blues")[5:9],brewer.pal(n = 11, "PiYG")[8:11])
#my_blue = brewer.pal(n = 9, "YlGnBu")[9:2]
#my_blue = brewer.pal(n = 9, "Blues")[9:2]

col1 <- colorRampPalette(c("#67001F", "#D6604D",  # reds
                           "#67001F","#ACD0E8",
                           "#ADDD8E", "#41AB5D", "#004529")) # greens( similar to: brewer.pal(9, "YlGn") )

my_blue = col1(length(unique(Motif.type)))

p3=ggplot(ta2,
          aes(x = ta2$"log10(Enrichment)", y = -1 * log10(ta2$P.value), color = (ta2$Motif.type), size=ta2$Perc_of_Targets*100)) + 
  #scale_colour_gradient(limits=c(0,1)) + 
  scale_colour_manual(values=my_blue) +
  #scale_colour_brewer(direction=-1) + 
  coord_cartesian(xlim = c(-0.2,0.4), ylim=c(0,90)) + 
  scale_size_continuous(limits=c(0,100)) + 
  #scale_size_area() + 
  #scale_size(limits = c(10, 100)) + 
  geom_point(aes(color = (ta2$Motif.type)))


p3 +
  geom_point() + 
  geom_text_repel(data = ta2,aes(label=ta2$Gene), size = 5, box.padding=unit(0.2, "lines"), point.padding = unit(0.2, "lines")) + 
  labs(x="log10(enrichment)") + 
  labs(y="-log10(p-value)") +
  labs(title=ttl) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(colour = "Motif type", size="target sequence\nwith motif [%]" ) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17)
        #, panel.background = element_rect(fill = 'grey')
        )

dev.off()



##############################################################
### PLOT 2 (loss at promoters)
##############################################################


folder="your_path2"

name="knownResults_ATAC_loss_promoters"

ttl="Lost ATAC signal at promoters"



##########################

setwd(folder)

### Read in file
ta=read.csv("knownResults.txt", sep="\t")

Gene=gsub("\\(.*","",ta$Motif.Name)
Perc_of_Targets=as.numeric(gsub("%","",ta$X..of.Target.Sequences.with.Motif))/100
Perc_of_Background=as.numeric(gsub("%","",ta$X..of.Background.Sequences.with.Motif))/100
Enrichment=Perc_of_Targets/Perc_of_Background

ta1=cbind(Gene,ta,Perc_of_Targets,Perc_of_Background,Enrichment,log10(Enrichment))

##########################
### Renaming to achieve consistent naming (e.g. E2A is identical to Irf3)
ta1$Gene=as.vector(ta1$Gene)

ta1$Gene=gsub("Sp","SP",ta1$Gene)
ta1$Gene=gsub("Klf","KLF",ta1$Gene)
ta1$Gene=gsub("Fli","FLI",ta1$Gene)
ta1$Gene=gsub("Elk","ELK",ta1$Gene)
ta1$Gene=gsub("Etv","ETV",ta1$Gene)
ta1$Gene=gsub("Atf","ATF",ta1$Gene)
ta1$Gene=gsub("Rfx","RFX",ta1$Gene)
ta1$Gene=gsub("Maz","MAZ",ta1$Gene)


##########################





############################################################################
### Color according to similar motifs
############################################################################

### Combine similar motifs as same motif type and save in vector "Motif.type"

ta2=ta1[1:35,]
Motif.type=ta2$Gene

Motif.type[grep("KLF",Motif.type)]="KLF-SP"
Motif.type[Motif.type %in% c("SP1","SP5","MAZ")]="KLF-SP"

Motif.type[Motif.type %in% c("ELK4","FLI1","ETV1","GABPA","ETV2","ELK1","EWS:FLI1-fusion","ETS","ERG","ELF1","EHF")]="ETS"
Motif.type[Motif.type %in% c("Ronin","GFY-Staf","GFY")]="GFY"
Motif.type[grep("ATF",Motif.type)]="ATF"
Motif.type[Motif.type %in% c("CRE")]="ATF"
Motif.type[grep("RFX",Motif.type)]="RFX"
Motif.type[Motif.type %in% c("X-box")]="RFX"
Motif.type[grep("NRF",Motif.type)]="NRF"

Motif.type[Motif.type %in% c("YY1","ZNF143|STAF")]="others"

ta2=cbind(ta1[1:35,],Motif.type)

### Plot

pdf(paste0(name,"_dotplot5_names_top35_motif-type_coloring.pdf"))

#my_bluegreen = rbind(brewer.pal(n = 9, "Blues")[5:9],brewer.pal(n = 11, "PiYG")[8:11])
#my_blue = brewer.pal(n = 9, "YlGnBu")[9:2]
#my_blue = brewer.pal(n = 9, "Blues")[9:2]

col1 <- colorRampPalette(c("#67001F", "#D6604D",  # reds
                           "#67001F","#ACD0E8",
                           "#ADDD8E", "#41AB5D", "#004529")) # greens( similar to: brewer.pal(9, "YlGn") )

my_blue = col1(length(unique(Motif.type)))


p3=ggplot(ta2, 
          aes(x = ta2$"log10(Enrichment)", y = -1 * log10(ta2$P.value), color = (ta2$Motif.type), size=ta2$Perc_of_Targets*100)) + 
  #scale_colour_gradient(limits=c(0,1)) + 
  scale_colour_manual(values=my_blue) +
  coord_cartesian(xlim = c(-0.1,0.4), ylim=c(0,60)) + 
  scale_size_continuous(limits=c(0,100)) + 
  #scale_size_area() + 
  #scale_size(limits = c(10, 100)) + 
  geom_point(aes(color = (ta2$Motif.type)))


p3 +
  geom_point() + 
  geom_text_repel(data = ta2,aes(label=ta2$Gene), size = 5, box.padding=unit(0.2, "lines"), point.padding = unit(0.5, "lines")) + 
  labs(x="log10(enrichment)") + 
  labs(y="-log10(p-value)") +
  labs(title=ttl) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(colour = "Motif type", size="target sequence\nwith motif [%]" ) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17)
        #, panel.background = element_rect(fill = 'grey')
  )


dev.off()






##############################################################
### PLOT 3 (gain at enhancers)
##############################################################


folder="your_path3"

name="knownResults_ATAC_gain_enhancers"

ttl="Gained ATAC signal at enhancers"

##########################

setwd(folder)

### Read in file
ta=read.csv("knownResults.txt", sep="\t")

Gene=gsub("\\(.*","",ta$Motif.Name)
Perc_of_Targets=as.numeric(gsub("%","",ta$X..of.Target.Sequences.with.Motif))/100
Perc_of_Background=as.numeric(gsub("%","",ta$X..of.Background.Sequences.with.Motif))/100
Enrichment=Perc_of_Targets/Perc_of_Background

ta1=cbind(Gene,ta,Perc_of_Targets,Perc_of_Background,Enrichment,log10(Enrichment))

### Renaming to achieve consistent naming (e.g. E2A is identical to Irf3)
ta1$Gene=as.vector(ta1$Gene)

ta1$Gene=gsub("E2A","TCF3/4",ta1$Gene)
ta1[grep("IRF",ta1$Gene),]$Gene="IRF"
ta1[grep("Tcf4",ta1$Gene),]$Gene="TCF3/4"
ta1[grep("Tcf3",ta1$Gene),]$Gene="TCF3/4"
ta1$Gene=gsub("NFAT","NFATC1",ta1$Gene)
ta1$Gene=gsub("Tcf","TCF",ta1$Gene)

ta1$Gene=gsub("Ascl1","ASCL1",ta1$Gene)
ta1$Gene=gsub("Ap4","AP4",ta1$Gene)
ta1$Gene=gsub("Myf5","MYF5",ta1$Gene)
ta1$Gene=gsub("MyoD","MYOD",ta1$Gene)
ta1$Gene=gsub("Lhx2","LHX2",ta1$Gene)
ta1$Gene=gsub("Egr1","EGR1",ta1$Gene)
ta1$Gene=gsub("Sox6","SOX6",ta1$Gene)
ta1$Gene=gsub("Atoh1","ATOH1",ta1$Gene)
ta1$Gene=gsub("Olig2","OLIG2",ta1$Gene)
ta1$Gene=gsub("Foxo3","FOXO3",ta1$Gene)
ta1$Gene=gsub("Ptf1a","PTF1a",ta1$Gene)
ta1$Gene=gsub("Slug","SLUG",ta1$Gene)
##########################




############################################################################
### Color according to similar motifs
############################################################################

### Combine similar motifs as same motif type and save in vector "Motif.type"

ta2=ta1[1:35,]
Motif.type=ta2$Gene

#Motif.type[grep("TCF",Motif.type)]="TCF"
Motif.type[Motif.type %in% c("TCF3/4","HEB","MyoG","PTF1a","TCF21","MYF5","AP4","ASCL1","ATOH1","MYOD","TCF12","NeuroG2","ZEB1","SLUG","ZBTB18","OLIG2","NeuroD1","SCL","TCFL2","BMAL1")]="TCF"
Motif.type[Motif.type %in% c("NFATC1")]="NFAT"
Motif.type[Motif.type %in% c("EGR1")]="EGR"
Motif.type[Motif.type %in% c("FOXK2","FOXK1","FOXO3")]="FOX"


ta2=cbind(ta1[1:35,],Motif.type)

### Plot

pdf(paste0(name,"_dotplot5_names_top35_motif-type_coloring.pdf"))

#my_bluegreen = rbind(brewer.pal(n = 9, "Blues")[5:9],brewer.pal(n = 11, "PiYG")[8:11])
#my_blue = brewer.pal(n = 9, "YlGnBu")[9:2]
#my_blue = brewer.pal(n = 9, "Blues")[9:2]

col1 <- colorRampPalette(c("#67001F", "#D6604D",  # reds
                           "#67001F","#ACD0E8",
                           "#ADDD8E", "#41AB5D", "#004529")) # greens( similar to: brewer.pal(9, "YlGn") )

my_blue = col1(length(unique(Motif.type)))


p3=ggplot(ta2, 
          aes(x = ta2$"log10(Enrichment)", y = -1 * log10(ta2$P.value), color = (ta2$Motif.type), size=ta2$Perc_of_Targets*100)) + 
  #scale_colour_gradient(limits=c(0,1)) + 
  scale_colour_manual(values=my_blue) +
  coord_cartesian(xlim = c(0,0.3), ylim=c(0,65)) + 
  scale_size_continuous(limits=c(0,100)) + 
  #scale_size_area() + 
  #scale_size(limits = c(10, 100)) + 
  geom_point(aes(color = (ta2$Motif.type)))


p3 +
  geom_point() + 
  geom_text_repel(data = ta2,aes(label=ta2$Gene), size = 5, box.padding=unit(0.3, "lines"), point.padding = unit(0.3, "lines")) + 
  labs(x="log10(enrichment)") + 
  labs(y="-log10(p-value)") +
  labs(title=ttl) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(colour = "Motif type", size="target sequence\nwith motif [%]" ) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17)
        #, panel.background = element_rect(fill = 'grey')
  )


dev.off()



