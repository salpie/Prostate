group_one <- (Sys.glob("/Users/salpienowinski/Documents/prostate_cancer/TOOLS/apt-1.17.0-x86_64-apple-lion/bin/finals/*output.txt"))
group_two <- Sys.glob("/Users/salpienowinski/Documents/prostate_cancer/TCGA/SNP6/gleason_7/*output.txt") # path to *segmentCN.txt files produced from TAPS 2.0


findGain_Loss <- function(Cn, mCn) {
  as.numeric(Cn > 2 | Cn < 2 | Cn == 2 & mCn == 1)
}


findGain_Loss <- function(Cn, mCn) {
  as.numeric(Cn < 2)
}


 number_group_one <- list()
 group_one_file <- read.table(group_one[1], header=T)
 names <- colnames(group_one_file)
 group_one_file$Cn <- group_one_file$nA + group_one_file$nB

 calc_mCn_forcnloh <- function(nA, nB) {
 	as.numeric(nA == 0 & nB > 1 | nB == 0 & nA > 1)
 }

for (i in 1:length(group_one)) {
	group_one_file <- read.table(group_one[i], header=T)
	colnames(group_one_file) <- names
	group_one_file$Cn <- group_one_file$nA + group_one_file$nB
	group_one_file$mCn <- calc_mCn_forcnloh(group_one_file$nA, group_one_file$nB)
	group_one_file$scna <- findGain_Loss(group_one_file$Cn, group_one_file$mCn)
	zero <- length(group_one_file$scna)
	one <- length(group_one_file$scna[group_one_file$scna == 1])
	number_group_one <- rbind(number_group_one, one/zero)

}


 number_group_two <- list()
 group_two_file <- read.table(group_two[1], header=T)
 names <- colnames(group_two_file)
 group_two_file$Cn <- group_two_file$nA + group_two_file$nB

 calc_mCn_forcnloh <- function(nA, nB) {
 	as.numeric(nA == 0 & nB > 1 | nB == 0 & nA > 1)
 }

for (i in 1:length(group_two)) {
	group_two_file <- read.table(group_two[i], header=T)
	colnames(group_two_file) <- names
	group_two_file$Cn <- group_two_file$nA + group_two_file$nB
	group_two_file$mCn <- calc_mCn_forcnloh(group_two_file$nA, group_two_file$nB)
	group_two_file$scna <- findGain_Loss(group_two_file$Cn, group_two_file$mCn)
	zero <- length(group_two_file$scna)
	two <- length(group_two_file$scna[group_two_file$scna == 1])
	number_group_two <- rbind(number_group_two, two/zero)

}
boxplot(as.numeric(number_group_one), names="PRAD_GLEASON_6", main="number of SCNA break points in PRAD GLEASON 6", col=c("gold"), ylim=c(0,1))
boxplot(as.numeric(number_group_one)*100, as.numeric(number_group_two)*100, names=c("PRAD_GLEASON_6", "PRAD_GLEASON_7"), beside=T, main="Percentage of SCNA in PRAD GLEASON 6 and 7", col=c("gold", "forest green"), ylim=c(0,100), ylab="% of SCNAs")
