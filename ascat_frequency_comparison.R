#####################################################
############To produce Frequency Table in R##########
#################Salpie Nowinski#####################
###########Date: 25/03/15############################



library(dplyr)
library(data.table)
library(biomaRt)




####### Group ONE and TWO location
group_one <- (Sys.glob("/Users/salpienowinski/Documents/prostate_cancer/TOOLS/apt-1.17.0-x86_64-apple-lion/bin/finals/*output.txt")) # path to *segmentCN.txt files produced from TAPS 2.0
group_two <- Sys.glob("/Users/salpienowinski/Documents/prostate_cancer/TCGA/SNP6/gleason_7/*output.txt") # path to *segmentCN.txt files produced from TAPS 2.0



library(dplyr)
library(data.table)
library(biomaRt)

findChromothripsis <- function(Cn) {
  if (length(Cn[Cn >= 5]) > 2) { #if more than 2 segements in chromosome have copy number 5
    y1 <- as.numeric(Cn >= 7 & lead(Cn) %in% 1:2) #if there are at least 2 oscillation of from copy number 1/2 to 7 or more then chromothripsis then score 1
    as.numeric(Cn > 7 & lead(Cn) %in% 1:2) #then score 1
  }
  else {
    as.numeric(!(Cn))
  }
}

findAmplicon_single <- function(Cn) {
  if (length(Cn[Cn >= 5]) > 2) {
    y1 <- as.numeric(Cn > 7 & lead(Cn) %in% 1:2)
    as.numeric(Cn >= 7 & (lead(Cn) <= 4 | lag(Cn) <= 4) & !y1)     #if chromothripsis still look for amplicons
  }
  else {
    as.numeric(Cn >= 7 & (lead(Cn) <= 4 | lag(Cn) <= 4))
  }  
}

findAmplicon_stepwise <- function(Cn) {
    as.integer((Cn >= 5 & (lead(Cn) > 7 | lag(Cn) > 7) )| (Cn > 7 & lead(Cn) >= 5 & lag(Cn) >= 5)) #if amplification is 'stepwise' then score 1
}


findGain <- function(Cn) {
  as.numeric(Cn > 2)
}

findGain_not_amp <- function(Cn) {
  as.numeric(Cn > 2 & (!y1 | !y2 | !y3))
}

findLoss <- function(Cn) {  
  as.numeric(Cn < 2)
}

findCnLoh <- function(nA, nB) {
  as.numeric(nA == 2 & nB == 0 | nB == 2 & nA == 0)
}

findAmp <- function(Cn) { 
  as.numeric(Cn >= 5)
}

rm(y)
rm(x)
rm(d)
rm(r)
rm(i)


  chromothrip <- as.numeric()
  amplicon_single <- as.numeric()
  amplicon_step <- as.numeric()


group_one_first <- read.table(group_one[1], header=TRUE)
names <- names(group_one_first)
group_one_first$SampleID <- NULL

group_one_first$Cn <- group_one_first$nA + group_one_first$nB

group_one_first$Cn[group_one_first$Cn == 0 ] <- 1
ds_group_one <- group_one_first

group_one_first$number_samples <- 1

group<- ds_group_one[which(is.na(ds_group_one)), ]

for (m in 1:23) {
     group_one_chr <- subset(group_one_first, group_one_first$Chr == m)
     if (nrow(group_one_chr) > 1) {
      for (y in 1:nrow(group_one_chr)-1) {
        a <- group_one_chr$End[y] == group_one_chr$Start[y+1]
        print(as.numeric(a))
        a <- a+1
        }
       to_iterate <- 1:((nrow(group_one_chr)*2)-2-a)
       for (j in to_iterate) {
         print(j)
         if (group_one_chr$End[j] != group_one_chr$Start[j+1]) {
          number_rows <- nrow(group_one_chr)
           rownames(group_one_chr) <- 1:number_rows
          print(group_one_chr)
           newrow <- cbind(m, group_one_chr$End[j], group_one_chr$Start[j+1], 0, 0, 0, 0, 0)
           print(nrow(group_one_chr))
           number_rows <- nrow(group_one_chr)
           rownames(group_one_chr) <- 1:number_rows
           group_one_chr[seq(j+1,nrow(group_one_chr)+1),] <- group_one_chr[seq(j,nrow(group_one_chr)),]
           group_one_chr[j+1,] <- newrow
           j <- j
         }
       }
      } 
   group <- rbind(group, group_one_chr)
   }


group_one_first$nA[group_one_first$nA==-1] <- 0
group_one_first$nA[group_one_first$nB==-1] <- 0



group_one_first <- group
group_one_first$loss <- findLoss(group_one_first$Cn)
group_one_first$loss[group_one_first$nProbes == 0] <- 0
group_one_first$gain <- findGain(group_one_first$Cn)
group_one_first$amp <- findAmp(group_one_first$Cn)
group_one_first$cnloh <- findCnLoh(group_one_first$nA, group_one_first$nB)


ds_group_one <- group_one_first


for (d in 2:length(group_one)) {

  print(group_one[d])



  group_one_file <- read.table(group_one[d], header=TRUE)
 group_one_file$SampleID <- NULL
 
group_one_file$Cn <- group_one_file$nA + group_one_file$nB



group_one_file$number_samples <- 1

group<- group_one_file[which(is.na(group_one_file)), ]

for (m in 1:23) {
    a <- 0
     group_one_chr <- subset(group_one_file, group_one_file$Chr == m)
     if (nrow(group_one_chr) > 1) {
      for (y in 1:nrow(group_one_chr)-1) {
        a <- group_one_chr$End[y] == group_one_chr$Start[y+1]
        print(as.numeric(a))
        a <- a+1
        }
      to_iterate <- 1:((nrow(group_one_chr)*2)-2 - a)
      for (j in to_iterate) {
         print(j)
         if (group_one_chr$End[j] != group_one_chr$Start[j+1]) {
          number_rows <- nrow(group_one_chr)
           rownames(group_one_chr) <- 1:number_rows
          print(group_one_chr)
           newrow <- cbind(m, group_one_chr$End[j], group_one_chr$Start[j+1], 0, 0, 0, 0, 0)
           print(nrow(group_one_chr))
           number_rows <- nrow(group_one_chr)
           rownames(group_one_chr) <- 1:number_rows
           group_one_chr[seq(j+1,nrow(group_one_chr)+1),] <- group_one_chr[seq(j,nrow(group_one_chr)),]
           group_one_chr[j+1,] <- newrow
           j <- j
         }
       }
      } 
   group <- rbind(group, group_one_chr)
   }


for (y in 1:nrow(group_one_chr)) {
  a <- group_one_chr$End[y] == group_one_chr$Start[y+1]
  print(as.numeric(a))
}




group_one_file <- group


group_one_file$nA[group_one_file$nA==-1] <- 0
group_one_file$nA[group_one_file$nB==-1] <- 0


group_one_file$loss <- findLoss(group_one_file$Cn)
group_one_file$loss[group_one_file$nProbes == 0] <- 0
  group_one_file$gain <- findGain(group_one_file$Cn)
  group_one_file$amp <- findAmp(group_one_file$Cn)
  group_one_file$cnloh <- findCnLoh(group_one_file$nA, group_one_file$nB)


  lst2 <- split(ds_group_one, ds_group_one$Chr)
  lst1 <- split(group_one_file, group_one_file$Chr)


require(survival)


 # merge
 df <- do.call(rbind, mapply(FUN = function(x, y) {
 
   x$event <- y$event <- 0
   idc.spl <- survSplit(x, cut=y$End, start='Start', end='End', event='event')
   dcis.spl <- survSplit(y, cut=x$End, start='Start', end='End', event='event')
   mrg <- merge(idc.spl, dcis.spl, 
                by=c('Chr', 'Start', 'End'), 
                #all=TRUE, 
                suffixes = c("_first","_second"))
   mrg[c('Chr', 'Start', 'End', 'loss_first', 'loss_second', 'gain_first', 'gain_second', 'cnloh_first', 'cnloh_second', 'amp_first', 'amp_second', 'number_samples_first', 'number_samples_second')]
 },
 lst1, lst2, SIMPLIFY=FALSE))

df$loss <- df$loss_second + df$loss_first
 df$gain <- df$gain_second + df$gain_first
 df$cnloh <- df$cnloh_second + df$cnloh_first
 df$amp <- df$amp_second + df$amp_first
 df$number_samples <- df$number_samples_first + df$number_samples_second
#df$amplicon_step <- df$amplicon_step_first + df$amplicon_step_second
#df$amplicon_single <- df$amplicon_single_first + df$amplicon_single_second
#df$chromothrip <- df$chromothrip_first + df$chromothrip_second

ds_group_one <- df[,c(1:3,14:18)]


}

rm(y)
rm(x)
rm(d)
rm(r)
rm(i)


  chromothrip <- as.numeric()
  amplicon_single <- as.numeric()
  amplicon_step <- as.numeric()


group_two_first <- read.table(group_two[1], header=TRUE)
names <- names(group_two_first)
group_two_first$SampleID <- NULL

group_two_first$Cn <- group_two_first$nA + group_two_first$nB

group_two_first$Cn[group_two_first$Cn == 0 ] <- 1
ds_group_two <- group_two_first

group_two_first$number_samples <- 1

group<- ds_group_two[which(is.na(ds_group_two)), ]

for (m in 1:23) {
     group_two_chr <- subset(group_two_first, group_two_first$Chr == m)
     if (nrow(group_two_chr) > 1) {
      for (y in 1:nrow(group_two_chr)-1) {
        a <- group_two_chr$End[y] == group_two_chr$Start[y+1]
        print(as.numeric(a))
        a <- a+1
        }
       to_iterate <- 1:((nrow(group_two_chr)*2)-2-a)
       for (j in to_iterate) {
         print(j)
         if (group_two_chr$End[j] != group_two_chr$Start[j+1]) {
          number_rows <- nrow(group_two_chr)
           rownames(group_two_chr) <- 1:number_rows
          print(group_two_chr)
           newrow <- cbind(m, group_two_chr$End[j], group_two_chr$Start[j+1], 0, 0, 0, 0, 0)
           print(nrow(group_two_chr))
           number_rows <- nrow(group_two_chr)
           rownames(group_two_chr) <- 1:number_rows
           group_two_chr[seq(j+1,nrow(group_two_chr)+1),] <- group_two_chr[seq(j,nrow(group_two_chr)),]
           group_two_chr[j+1,] <- newrow
           j <- j
         }
       }
      } 
   group <- rbind(group, group_two_chr)
   }






group_two_first <- group
group_two_first$loss <- findLoss(group_two_first$Cn)
group_two_first$loss[group_two_first$nProbes == 0] <- 0
group_two_first$gain <- findGain(group_two_first$Cn)
group_two_first$amp <- findAmp(group_two_first$Cn)
group_two_first$cnloh <- findCnLoh(group_two_first$nA, group_two_first$nB)


ds_group_two <- group_two_first


for (d in 2:length(group_two)) {

  print(group_two[d])



  group_two_file <- read.table(group_two[d], header=TRUE)
 group_two_file$SampleID <- NULL
 
group_two_file$Cn <- group_two_file$nA + group_two_file$nB



group_two_file$number_samples <- 1

group<- group_two_file[which(is.na(group_two_file)), ]

for (m in 1:23) {
    a <- 0
     group_two_chr <- subset(group_two_file, group_two_file$Chr == m)
     if (nrow(group_two_chr) > 1) {
      for (y in 1:nrow(group_two_chr)-1) {
        a <- group_two_chr$End[y] == group_two_chr$Start[y+1]
        print(as.numeric(a))
        a <- a+1
        }
      to_iterate <- 1:((nrow(group_two_chr)*2)-2 - a)
      for (j in to_iterate) {
         print(j)
         if (group_two_chr$End[j] != group_two_chr$Start[j+1]) {
          number_rows <- nrow(group_two_chr)
           rownames(group_two_chr) <- 1:number_rows
          print(group_two_chr)
           newrow <- cbind(m, group_two_chr$End[j], group_two_chr$Start[j+1], 0, 0, 0, 0, 0)
           print(nrow(group_two_chr))
           number_rows <- nrow(group_two_chr)
           rownames(group_two_chr) <- 1:number_rows
           group_two_chr[seq(j+1,nrow(group_two_chr)+1),] <- group_two_chr[seq(j,nrow(group_two_chr)),]
           group_two_chr[j+1,] <- newrow
           j <- j
         }
       }
      } 
   group <- rbind(group, group_two_chr)
   }


for (y in 1:nrow(group_two_chr)) {
  a <- group_two_chr$End[y] == group_two_chr$Start[y+1]
  print(as.numeric(a))
}




group_two_file <- group

group_two_file$loss <- findLoss(group_two_file$Cn)
group_two_file$loss[group_two_file$nProbes == 0] <- 0
  group_two_file$gain <- findGain(group_two_file$Cn)
  group_two_file$amp <- findAmp(group_two_file$Cn)
  group_two_file$cnloh <- findCnLoh(group_two_file$nA, group_two_file$nB)


  lst2 <- split(ds_group_two, ds_group_two$Chr)
  lst1 <- split(group_two_file, group_two_file$Chr)


require(survival)


 # merge
 df <- do.call(rbind, mapply(FUN = function(x, y) {
 
   x$event <- y$event <- 0
   idc.spl <- survSplit(x, cut=y$End, start='Start', end='End', event='event')
   dcis.spl <- survSplit(y, cut=x$End, start='Start', end='End', event='event')
   mrg <- merge(idc.spl, dcis.spl, 
                by=c('Chr', 'Start', 'End'), 
                #all=TRUE, 
                suffixes = c("_first","_second"))
   mrg[c('Chr', 'Start', 'End', 'loss_first', 'loss_second', 'gain_first', 'gain_second', 'cnloh_first', 'cnloh_second', 'amp_first', 'amp_second', 'number_samples_first', 'number_samples_second')]
 },
 lst1, lst2, SIMPLIFY=FALSE))

df$loss <- df$loss_second + df$loss_first
 df$gain <- df$gain_second + df$gain_first
 df$cnloh <- df$cnloh_second + df$cnloh_first
 df$amp <- df$amp_second + df$amp_first
 df$number_samples <- df$number_samples_first + df$number_samples_second
#df$amplicon_step <- df$amplicon_step_first + df$amplicon_step_second
#df$amplicon_single <- df$amplicon_single_first + df$amplicon_single_second
#df$chromothrip <- df$chromothrip_first + df$chromothrip_second

ds_group_two <- df[,c(1:3,14:18)]


}



lst1 <- split(ds_group_one, ds_group_one$Chr)
lst2 <- split(ds_group_two, ds_group_two$Chr)


df <- do.call(rbind, mapply(FUN = function(x, y) {
 
   x$event <- y$event <- 0
   idc.spl <- survSplit(x, cut=y$End, start='Start', end='End', event='event')
   dcis.spl <- survSplit(y, cut=x$End, start='Start', end='End', event='event')
   mrg <- merge(idc.spl, dcis.spl, 
                by=c('Chr', 'Start', 'End'), 
                #all=TRUE, 
                suffixes = c("_group_one","_group_two"))
   mrg[c('Chr', 'Start', 'End', 'loss_group_one', 'loss_group_two', 'gain_group_one', 'gain_group_two', 'cnloh_group_one', 'cnloh_group_two', 'amp_group_one', 'amp_group_two', 'number_samples_group_one', 'number_samples_group_two')]
 },
 lst1, lst2, SIMPLIFY=FALSE))



df$group_one_loss <- df$number_samples_group_one - df$loss_group_one
df$group_two_loss <- df$number_samples_group_two - df$loss_group_two
df$group_one_gain <- df$number_samples_group_one - df$gain_group_one
df$group_two_gain <- df$number_samples_group_two - df$gain_group_two
df$group_one_cnloh <- df$number_samples_group_one - df$cnloh_group_one
df$group_two_cnloh <- df$number_samples_group_two - df$cnloh_group_two


df$pvalue_loss <- apply(as.matrix(df[,c(4,14,5,15)]), 1, function(x) 
         fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)


df$pvalue_gain <- apply(as.matrix(df[,c(6,16,7,17)]), 1, function(x) 
         fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)

df$pvalue_cnloh <- apply(as.matrix(df[,c(8,18,9,19)]), 1, function(x) 
         fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)


newdata <- df[order(df$Chr, df$Start),]


we <- setDT(newdata)[, .ind:= cumsum(c(TRUE,Start[-1]!=End[-.N])),
        list(loss_group_one, loss_group_two, gain_group_one, gain_group_two, cnloh_group_one, cnloh_group_two, amp_group_one, amp_group_two, chromothrip_group_one, chromothrip_group_two, amplicon_single_group_one, amplicon_single_group_two, amplicon_step_group_one, amplicon_step_group_two,pvalue_gain, pvalue_loss, pvalue_cnloh)][,
       list(chr=Chromosome[1], start=Start[1], stop=End[.N]),
       list(loss_group_one, loss_group_two, gain_group_one, gain_group_two, cnloh_group_one, cnloh_group_two, amp_group_one, amp_group_two, chromothrip_group_one, chromothrip_group_two, amplicon_single_group_one, amplicon_single_group_two, amplicon_step_group_one, amplicon_step_group_two,pvalue_gain, pvalue_loss, pvalue_cnloh, .ind)][,.ind:=NULL][]

write.table(newdata, file="gleason_6_gleason_7.txt", sep="\t", quote=F)


for (i in 1:23) {
  grp1 <- subset(ds_group_one, ds_group_one$Chromosome == i)
  grp2 <- subset(ds_group_two, ds_group_two$Chromosome == i)

  print(min(grp1$Start))
  print(min(grp2$Start))
  print(max(grp1$Start))
  print(max(grp2$Start))

}
