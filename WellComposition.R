###########################################################################
###### Script for summarising MySeq reads into per well composition  #####
###########################################################################


### Clear the workspace
rm(list=ls())

# .libPaths(c("C:\\rlib", .libPaths()))

### load the libraries you need
j <- c("reshape2","ggplot2","akima","plyr","dplyr","grid","RColorBrewer","scales","lme4")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

for (i in j){
  library(j,character.only=TRUE)
}
 
 
library(reshape2)
library(ggplot2)
library(akima)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(scales)

########################################################################################
################### stacked bar charts of well composition   ###########################
########################################################################################

### read in the data - this is the read count output from metaBEAT, I have manually removed the .biom header line and the taxonomy column as this make live a lot easier in R
my.reads<-read.csv(file=paste("data/metaBEAT-processed.txt",sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)

### read in the sample by plate data - this is the querymap file used in metaBEAT with additional columns for PCR plate
my.plates<-read.csv(file="data/sample_metadata.txt",
                    sep="\t", stringsAsFactors=FALSE, header=TRUE)

### trim the plate data to the necessary columns
my.plates<-my.plates[,1:3]


### set a minimum occurance for an assignment to be trusted - i.e. at this point you can decide that taxa in only one sample might be considered unreliable, in which case set occurance=2
occurance=1

### subset the data frame to drop all assignments occuring fewer than the frequency specified above
my.reads.subs<-subset(my.reads, subset = rowSums(my.reads[2:ncol(my.reads)] > 0) >= occurance)

### transpose the read data
my.reads.trans<-recast(my.reads.subs, variable~OTU_ID)
colnames(my.reads.trans)[1]<-"sample"
my.reads.trans$sample<-as.character(my.reads.trans$sample)

### use match to add the plate data to the read data
my.reads.trans$plate<-my.plates$plate[match(my.reads.trans$sample,my.plates$sample)]
my.reads.trans$plate.numeric<-my.plates$plate.numeric[match(my.reads.trans$sample,my.plates$sample)]


### use grep to add nesting data

my.reads.trans$type<-ifelse(grepl("POS1",my.reads.trans$sample),"DNApositive",
                            ifelse(grepl("POS2",my.reads.trans$sample),"DNApositive",
                                   ifelse(grepl("POS3",my.reads.trans$sample),"PCRpositive",
                                          ifelse(grepl("POS4",my.reads.trans$sample),"PCRpositive",
                                                 ifelse(grepl("NEG",my.reads.trans$sample),"Negative",
                                                        ifelse(grepl("PROC",my.reads.trans$sample),"Procedural","Sample"))))))
                                   

### make sample and plate factors for faceting and ordering
my.reads.trans$sample<-as.factor(my.reads.trans$sample)
my.reads.trans$plate<-as.factor(my.reads.trans$plate)



### total all the reads
my.reads.trans$total<-rowSums(my.reads.trans[c(2:81)])

### calculate the percentage of reads in each well that are Buddleia
#my.reads.trans$Percent.Budd<-(my.reads.trans$Buddleja_davidii/my.reads.trans$total)*100

### calculate the percentage of reads in each well that are Fraxinus (=main environmental contamination)
my.reads.trans$Percent.Frax<-(my.reads.trans$Fraxinus_excelsior/my.reads.trans$total)*100


### set a minimum read coverage for accepting an assignment
cutoff<-20

### replace all the species that are less than the cutoff of the total reads with zero
my.reads.trans[,2:81][my.reads.trans[,2:81] < cutoff] <- 0

### subset the data frame to drop all coloumns containing only zeros
my.reads.trans.nozero<-cbind(my.reads.trans[,c(1, 82:86)], subset(my.reads.trans[,c(2:81)], select = colSums(my.reads.trans[,c(2:81)]) !=0))

### remove the total variable
my.reads.trans.drop<-my.reads.trans.nozero[,c(1:4,6:68)]

### to remove the positive and negative samples for OTU counting in a barplot, I need to iteratively subset them out using the nest column.
### I've not found a more elegant way of doing this yet
my.reads.trans.samps.only<-my.reads.trans.drop[my.reads.trans$type != "Negative", ]
my.reads.trans.samps.only<-my.reads.trans.samps.only[my.reads.trans.samps.only$type != "DNApositive", ]
my.reads.trans.samps.only<-my.reads.trans.samps.only[my.reads.trans.samps.only$type != "PCRpositive", ]
my.reads.trans.samps.only<-my.reads.trans.samps.only[my.reads.trans.samps.only$type != "Procedural", ]

### order the type for ordering the samples in the plot
my.reads.trans.drop$type <- factor(my.reads.trans.drop$type,
                                   levels=c("Sample","DNApositive","PCRpositive","Procedural","Negative"))




### reorder my.reads.melt by both percentage Fraxinus, and plate
my.reads.trans.drop<-my.reads.trans.drop[order(my.reads.trans$Percent.Frax,my.reads.trans.drop$type),]


### make a panel factor after setting the decreasing Fraxinus order
my.reads.trans.drop$panel<-as.factor(c(rep(1,times=336/3),
                                       rep(2,times=336/3),
                                       rep(3,times=(336/3)-1)))

### melt the data into long format
my.reads.melt<-melt(my.reads.trans.drop, id.vars=c("sample","plate", "plate.numeric","Percent.Frax","type","panel"))
colnames(my.reads.melt)<-c("Sample","Plate","Plate.numeric","Percent.Frax","Type","Panel","Species","Reads")


### Label hypothesized contamination
my.reads.melt$contam<-ifelse(grepl("Holcus",my.reads.melt$Species),"Environmental",
                             ifelse(grepl("Hordeum",my.reads.melt$Species),"Environmental",
                                    ifelse(grepl("Zea",my.reads.melt$Species),"Environmental",
                                           ifelse(grepl("Fraxinus",my.reads.melt$Species),"Environmental",
                                                  ifelse(grepl("Plantago",my.reads.melt$Species),"Environmental",
                                                         ifelse(grepl("Avena",my.reads.melt$Species),"Environmental",
                                                                ifelse(grepl("Festuca",my.reads.melt$Species),"Environmental",
                                                                       ifelse(grepl("Poa",my.reads.melt$Species),"Environmental",
                                                                              ifelse(grepl("Triticum",my.reads.melt$Species),"Environmental",
                                                                                     ifelse(grepl("Rumex",my.reads.melt$Species),"Environmental",
                                                                                            ifelse(grepl("Populus",my.reads.melt$Species),"Environmental",
                                                                                                   ifelse(grepl("Salix",my.reads.melt$Species),"Environmental",
                                                                                                          ifelse(grepl("Urtica",my.reads.melt$Species),"Environmental",
                                                                                                                 ifelse(grepl("Poaceae",my.reads.melt$Species),"Environmental",
                                                                                                                        ifelse(grepl("Salicaceae",my.reads.melt$Species),"Environmental",
                                                                                                                               ifelse(grepl("Saxifragales",my.reads.melt$Species),"Laboratory",
                                                                                                                                      ifelse(grepl("Crassulaceae",my.reads.melt$Species),"Laboratory",
                                                                                                                                             ifelse(grepl("Asparagaceae",my.reads.melt$Species),"Laboratory",
                                                                                                                                                    ifelse(grepl("Kalanchoe",my.reads.melt$Species),"Laboratory",
                                                                                                                                                           ifelse(grepl("Rhodiola",my.reads.melt$Species),"Laboratory",
                                                                                                                                                                  ifelse(grepl("Dracaena",my.reads.melt$Species),"Laboratory",
                                                                                                                                                                         ifelse(grepl("Liliopsida",my.reads.melt$Species),"unassigned",
                                                                                                                                                                                ifelse(grepl("Streptophyta",my.reads.melt$Species),"unassigned",
                                                                                                                                                                                       ifelse(grepl("unassigned",my.reads.melt$Species),"unassigned","Sample"))))))))))))))))))))))))

          

### order by Type
my.reads.melt <- my.reads.melt[order(my.reads.melt$Type, my.reads.melt$Sample, my.reads.melt$contam),]

### create an ID variable for the order
my.reads.melt$level.order<-seq(1, length(my.reads.melt$Type),1)

### see what species remain after filtering
levels(my.reads.melt$Species)
levels(my.reads.melt$contam)

### order the species for plotting the legend
#my.reads.melt$Sample <- reorder(my.reads.melt$Sample, my.reads.melt$level.order)

###order contamination
my.reads.melt$contam <- factor(my.reads.melt$contam,
                                levels=c("Sample",
                                         "Environmental",
                                         "Laboratory",
                                         "unassigned"))

### order the type
my.reads.melt$Type <- factor(my.reads.melt$Type,
                               levels=c("Sample",
                                        "DNApositive",
                                        "PCRpositive",
                                        "Procedural",
                                        "Negative"))


### order the colours for plotting the legend
my.colours <- c("mediumorchid4",
                "cornflowerblue",
                "gold1",
                "gray80")


### order the samples putting controls last and in increasing plate number
my.reads.melt <- my.reads.melt[order(my.reads.melt$level.order),]


### create a nice colour scale using colourRampPalette and RColorBrewer
vivid.colours2<-brewer.pal(length(unique(my.reads.melt$contam)),"Set1")

### count the columns greater than zero and write to a new data frame - this is used for the barplot below the composition diagram
hit.hist<-data.frame(OTUs = rowSums(my.reads.trans.samps.only[c(7:12)] != 0), Type=my.reads.trans.samps.only$type)

########################################################################################
################### make a plot of % composition by PCR plate  ######################################
########################################################################################

### aggregate
my.reads.agg <- aggregate(x = my.reads.melt$Reads, by = list(my.reads.melt$Type,my.reads.melt$contam), FUN = sum)

colnames(my.reads.agg)[1] <- "Type"
colnames(my.reads.agg)[2] <- "contam"
colnames(my.reads.agg)[3] <- "Reads"

### set up the ggplot - well type
well.type<-ggplot(data=my.reads.agg, aes(x=Type, y=Reads, fill=contam)) +
  ### make it a stacked barplot and set the bars to be the same height
  geom_bar(position="fill", stat="identity") +
  ### give it a percentage scale
  scale_y_continuous(labels = percent_format()) +
  ### set the colours
  # scale_fill_manual(name="Species",
  #                    values = jet.colors3(length(unique(my.reads.melt$Species)))) +
  scale_fill_manual(name="contamination", labels=c("Sample",
                                                   "Environmental",
                                                   "Laboratory",
                                                   "unassigned"),
                    values = my.colours) +
  ### add a sensible y axis label
  labs(y = "% of reads per well", x="PCR wells") +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_text(size = rel(1.1), colour="black"),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = rel(1.1), colour="black"),
        axis.title.y = element_text(size = rel(1), vjust=2),
        axis.title.x = element_text(size = rel(1), vjust=-1.3),
        legend.text = element_text(size = rel(1), face="italic"),
        legend.title = element_text(size = rel(1)),
        strip.text.x = element_blank(),
        strip.background=element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.1),
        plot.margin=unit(c(0.1, 0.1, 1, 1), "lines"))

well.type


### set up the ggplot - wells
well.composition<-ggplot(data=my.reads.melt, aes(x=Type, y=Reads, fill=contam)) +
  ### make it a stacked barplot and set the bars to be the same height
  geom_bar(position="fill", stat="identity") +
  ### wrap the plot by plate
  facet_wrap(~Panel, scales="free_x", nrow=3, ncol=1) +
  ### give it a percentage scale
  scale_y_continuous(labels = percent_format()) +
  ### set the colours
  # scale_fill_manual(name="Species",
  #                    values = jet.colors3(length(unique(my.reads.melt$Species)))) +
  scale_fill_manual(name="contamination", labels=c("Sample",
                                             "Environmental",
                                             "Laboratory",
                                             "unassigned"),
                    values = my.colours) +
  ### add a sensible y axis label
  labs(y = "% of reads per well", x="PCR wells") +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = rel(1.1), colour="black"),
        axis.title.y = element_text(size = rel(1), vjust=2),
        axis.title.x = element_text(size = rel(1), vjust=-1.3),
        legend.text = element_text(size = rel(1), face="italic"),
        legend.title = element_text(size = rel(1)),
        strip.text.x = element_blank(),
        strip.background=element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.1),
        plot.margin=unit(c(0.1, 0.1, 1, 1), "lines"))

well.composition

### Make a ggplot object of our OTU counts
c<- ggplot(hit.hist, aes(factor(reorder(OTUs, -OTUs)))) +
  geom_bar(alpha=0.5, fill="#3399ff") + labs(y = "Frequency", x="OTUs per well") +
  coord_flip() +
  ### rotate the x-axis labels and resize the text for the svg
  theme(axis.text = element_text(size = rel(1.1), colour="black"),
        axis.title.y = element_text(size = rel(1), vjust=2),
        axis.title.x = element_text(size = rel(1), vjust=-1),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        plot.margin=unit(c(0.1, 0.1, 1, 1), "lines"))

grid.arrange(well.composition, c, heights=c(3/4, 1/4), ncol=1)

dev.off()

########################################################################################
# work out percentage parasitism with carcelia #
########################################################################################

### divide the number of well containing Carcelia by all well that are not +ve or -ve to get percentage
percent.cacelia<-(colSums(my.reads.trans[3] > 0)/919)*100