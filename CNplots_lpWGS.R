library(ggplot2)

dna_ng <- read.table("dna_ng.txt", header=F)

names(dna_ng) <- c("Type", "DNA_ng", "Depth")

dna_ng$Type <- as.factor(dna_ng$Type)

theme_set(theme_gray(base_size = 18))

ggplot(dna_ng, aes(x=Depth, y=DNA_ng, shape=Type)) +
  geom_point(aes(colour = Type), size = 2.7)+
  ylab("DNA ng")

load("dna_ng.Rdata")


#stacked barplot of the number and types of aberrations for samples in each group
#first get the aberration type - use the 50Kb segmentation



require(QDNAseq)
require(CGHcall)
require(reshape)

load("/Users/salpienowinski/Desktop/lpWGS/plots/for_thesis/calls_50_5.Rdata")


calls$location <- rownames(calls)
calls <- calls[,c(31,1:30)]


setwd("/Users/salpienowinski/Desktop/lpWGS/plots/for_thesis")
load("calls_50_5.Rdata")

names(calls) <- gsub(".bwa.sam.bam.merged.filtered.fixed_mate.marked_duplicates_picard", "", names(calls))

total_homozygous_loss <- data.frame()
total_heterozygous_loss <- data.frame()
total_gain <- data.frame()
total_amplification <- data.frame()

#all aberrations
for (i in 1:30) {
	homozygous_loss <- length(subset(calls[,i], calls[,i] == -2))
	heterozygous_loss <- length(subset(calls[,i], calls[,i] == -1))
	gain <- length(subset(calls[,i], calls[,i] == 1))
	amplification <- length(subset(calls[,i], calls[,i] == 2))
	total_homozygous_loss <- rbind(total_homozygous_loss, homozygous_loss)
	total_heterozygous_loss <- rbind(total_heterozygous_loss, heterozygous_loss)
	total_gain <- rbind(total_gain, gain)
	total_amplification <- rbind(total_amplification, amplification)
}

total <- cbind(total_homozygous_loss, total_heterozygous_loss, total_gain, total_amplification)
names(total) <- c("homozygous_loss", "heterozygous_loss", "gain", "amplification")
type <- read.table("type.txt", header=F)
total$type <- type$V1
total$sampleName <- 1:30

data <- melt(total, id=c("sampleName", "type"))

#the separate into gains and losses
names(data)[3] <- "Type_of_aberration"

data <- data[c(91:120, 61:90, 31:60, 1:30),]

ggplot(data=data, aes(x=sampleName, y=value, fill=Type_of_aberration, order=Type_of_aberration)) +
    geom_bar(stat="identity")+
    facet_grid(. ~ type) +
      ylab("Number aberrant bins")+
      xlab("Sample ID")+
      scale_fill_manual(values=c("darkblue", "lightblue", "red", "pink"))

number <- read.table("number.txt", header=T)

theme_set(theme_gray(base_size = 18))

levels(number$Type_of_aberration) <- levels(number$Type_of_aberration)[c(4,3,2,1)]


number[order(number$Type_of_aberration)

ggplot(data=number, aes(x=type, y=Number_aberrant_bins, fill=Type_of_aberration)) +
    geom_bar(stat="identity")+
      xlab("Sample Type")+
      scale_fill_manual(values=c("pink", "red", "lightblue", "darkblue"))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))












### similarity score and bin stuff

setwd("/Users/salpienowinski/Desktop/lpWGS/plots/for_thesis")


load("/Users/salpienowinski/Desktop/lpWGS/plots/for_thesis/calls_50_5.Rdata")


compareTo2 <- function(copynumber1, copynumber2) {
  as.numeric(((copynumber1 == -1 & copynumber2 == -1) | (copynumber1 == -2 & copynumber2 == -2) | (copynumber1 == 1 & copynumber2 == 1) | (copynumber1 == 2 & copynumber2 == 2)))
}

subset_location <- function(pair) {
	subset(pair, pair[,1] == -1 & pair[,2] == -1 | pair[,1] == -2 & pair[,2] == -2 | pair[,1] == 1 & pair[,2] == 1 | pair[,1] == 2 & pair[,2] == 2)
}

#load("500_all_samples.Rdata")
unstable_pairs <- read.table("/Users/salpienowinski/Desktop/lpWGS/plots/tumour/unstable_pairs.txt")
stable_pairs <- read.table("/Users/salpienowinski/Desktop/lpWGS/plots/tumour/stable_pairs.txt")


unstable_pairs_info <- read.table("/Users/salpienowinski/Desktop/lpWGS/plots/tumour/unstable_pairs_info.txt")
stable_pairs_info <- read.table("/Users/salpienowinski/Desktop/lpWGS/plots/tumour/stable_pairs_info.txt")


calls$location <- rownames(calls)
calls <- na.omit(calls)
ca <- calls

total_proportion <- data.frame()
total_proportion_time <- data.frame()
total_location <- data.frame()

time_to_biopsy <- 

for (i in 1:nrow(unstable_pairs)) {
	pair <- ca[,c(unstable_pairs$V1[i], unstable_pairs$V2[i])]
	pair <- na.omit(pair)
	a <- length(subset(pair[,1], pair[,1]!=0))
	b <- length(subset(pair[,2], pair[,2]!=0))
	#proportion <- (sum(compareTo2(pair[,1], pair[,2]))/(a+b))*nrow(calls)
	proportion <- (sum(compareTo2(pair[,1], pair[,2]))/(a+b))
	total_proportion <- rbind(total_proportion, proportion)
	total_proportion_time <- rbind(total_proportion_time, (proportion/unstable_pairs_info$V1[i]))
	location <- subset_location(pair)
	if (nrow(location)>1) {
		location$names <- paste(c(names(location)[1:2]), collapse="-")
		location$location <- rownames(location)
		names(location) <- c("sample1", "sample2", "names", "location")
		rownames(location) <- 1:nrow(location)
		total_location <- rbind(total_location, location)
	}
}
total_proportion_unstable <- total_proportion
total_location_unstable <- total_location
total_proportion_time_unstable <- total_proportion_time

total_proportion <- data.frame()
total_proportion_time <- data.frame()
total_location <- data.frame()


for (i in 1:nrow(stable_pairs)) {
	pair <- ca[,c(stable_pairs$V1[i], stable_pairs$V2[i])]
	pair <- na.omit(pair)
	a <- length(subset(pair[,1], pair[,1]!=0))
	b <- length(subset(pair[,2], pair[,2]!=0))
	#proportion <- (sum(compareTo2(pair[,1], pair[,2]))/(a+b))*nrow(calls)
	proportion <- (sum(compareTo2(pair[,1], pair[,2]))/(a+b))
	total_proportion <- rbind(total_proportion, proportion)
	total_proportion_time <- rbind(total_proportion_time, (proportion/stable_pairs_info$V1[i]))
	location <- subset_location(pair)
	if (nrow(location)>1) {
		location$names <- paste(c(names(location)[1:2]), collapse="-")
		location$location <- rownames(location)
		names(location) <- c("sample1", "sample2", "names", "location")
		rownames(location) <- 1:nrow(location)
		total_location <- rbind(total_location, location)
	}
}


total_proportion_stable <- total_proportion
total_location_stable <- total_location
total_proportion_time_stable <- total_proportion_time

##boxplots
total_proportion_stable$status <- "stable"
total_proportion_unstable$status <- "upgrading"
names(total_proportion_stable)[1] <- "proportion"
names(total_proportion_unstable)[1] <- "proportion"

tot <- rbind(total_proportion_stable, total_proportion_unstable)

tot$patient_ID <- c("Patient_2", "Patient_10", "Patient_9", "Patient_15", "Patient_16", "Patient_17", "Patient_5", "Patient_8", "Patient_3", "Patient_1", "Patient_11", "Patient_12", "Patient_13")


ggplot(tot, aes(x=status, y=proportion, fill=status)) + 
geom_boxplot(fill=c("darkgrey", "red")) + 
      geom_jitter(width = 0.2)+
      ylab("proportion of similar aberrant bins")
      theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))

#barplot - so can add in copy number and location


ggplot(data=tot, aes(x=reorder(patient_ID, -proportion), y=proportion, fill=status)) +
    geom_bar(stat="identity")+
     theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual("legend", values = c("stable" = "darkgrey", "upgrading" = "red"))+
      ylab("Proportion of similar aberrant bins")+
      xlab("Patient ID")





#recurrent regions of scnas in same location in more than one patient

#recurrent

tot_stab <- total_location_stable[duplicated(total_location_stable$location),]
tot_unstab <- total_location_unstable[duplicated(total_location_unstable$location),]


##number of aberrant regions in initial stable vs initial unstable

total_bins <- data.frame()
for (i in 1:nrow(unstable_pairs)) {
	sample <- calls[,c(unstable_pairs$V1[i])]
	aberrant_bins <- length(subset(sample, sample==-1|sample==-2|sample==1|sample==2))
	total_bins <- rbind(total_bins, aberrant_bins)
}
total_bins_unstable <- total_bins

total_bins <- data.frame()
for (i in 1:nrow(stable_pairs)) {
	sample <- calls[,c(stable_pairs$V1[i])]
	aberrant_bins <- length(subset(sample, sample==-1|sample==-2|sample==1|sample==2))
	total_bins <- rbind(total_bins, aberrant_bins)
}
total_bins_stable <- total_bins



###add in other stables
stab <- stable_pairs$V1
stab <- c(stab, 3, 8, 22)

total_bins <- data.frame()
for (i in 1:length(stab)) {
	sample <- calls[,c(stab[i])]
	aberrant_bins <- length(subset(sample, sample==-1|sample==-2|sample==1|sample==2))
	total_bins <- rbind(total_bins, aberrant_bins)
}
total_bins_stable <- total_bins


total_bins_loss <- data.frame()
for (i in 1:length(stab)) {
	sample <- calls[,c(stab[i])]
	aberrant_bins <- length(subset(sample, sample==-1|sample==-2))
	total_bins_loss <- rbind(total_bins_loss, aberrant_bins)
}
total_bins_stable_loss <- total_bins_loss


total_bins_loss <- data.frame()
for (i in 1:length(unstable_pairs)) {
	sample <- calls[,c(unstable_pairs$V1[i])]
	aberrant_bins <- length(subset(sample, sample==-1|sample==-2))
	total_bins_loss <- rbind(total_bins_loss, aberrant_bins)
}
total_bins_unstable_loss <- total_bins_loss

###independant segments of loss....
