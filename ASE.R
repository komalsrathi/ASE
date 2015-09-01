# load libraries
library(data.table)
library(plyr)
library(doMC)
library(qdapDictionaries)
library(qdap)
registerDoMC(32) # register 32 cores for the analysis

# read genotype file
dt = fread('genotype.csv')

# create function for doing allele specific expression
myfunction <- function(arg1)
{
  # create output bed file name
  p0 <- paste("~/allele_specific_expression/",paste(unique(arg1$sample_name),".bed",sep=""),sep="")

  write.table(arg1[,c(9,11,11)],paste(p0),row.names=F,quote=F,sep="\t",col.names=F)
  
  # print sample name
  cat(noquote(paste("\nAllele Specific Expression",unique(arg1$sample_name),sep="...")))
  
  # print Processing Counts
  cat(noquote(paste("\nProcessing Counts",unique(arg1$sample_name),sep="...")))
  
  # call bam-readcount program
  # filter reads that have base quality less than 20
  # p2 to p6 will generate a system command
  p2 <- paste("bam-readcount -b 20 -f hg19.fa -l",p0,sep=" ")
  p3 <- paste("~/allele_specific_expression/",paste(unique(arg1$sample_name),"_filter_sort.bam",sep=""),sep="")
  p4 <- paste(p2,p3,sep=" ")
  p5 <- paste("~/allele_specific_expression/",paste(unique(arg1$sample_name),"_output.txt",sep=""),sep="")
  p6 <- paste(p4,p5,sep=" > ")
  cat("\n")
  system(p6) # call the command generated above
  
  # write output
  cat(noquote(paste("\nWriting Output...",unique(arg1$sample_name),sep="")))  
  
  # different columns per line bc of indels
  no_col <- max(count.fields(paste(p5), sep = "\t")) 
  
  # read in bamreadcount output
  bamreadcounts <- read.table(paste(p5),sep="\t",fill=TRUE,col.names=1:no_col,check.names=F) 
  
  # separate contents of columns 6-9 using based on ':'
  x <- colsplit2df(bamreadcounts[,c(1:4,6:9)],c(5:8),sep=":") 
  
  # extract relevant cols
  x <- x[,c(1:4,10:11,24:25,38:39,52:53)] 
  
  # change column names
  colnames(x) <- c("chrom","chrom_end","ref","depth","A_pos","A_neg","C_pos","C_neg","G_pos","G_neg","T_pos","T_neg") 
  
  # merge the two files
  gt.merge<-unique(merge(arg1,x,by=c("chrom","chrom_end")))
  
  # process the counts
  gt.merge <- data.frame(gt.merge,as.data.frame(sapply(gt.merge[c("allelea", "alleleb")],
                                                       function(x) as.integer(gt.merge[cbind(seq_len(nrow(gt.merge)),
                                                                                             match(paste(x, "pos", sep = "_"), names(gt.merge)))]))))
  
  gt.merge <- data.frame(gt.merge,as.data.frame(sapply(gt.merge[c("allelea", "alleleb")],
                                                       function(x) as.integer(gt.merge[cbind(seq_len(nrow(gt.merge)),
                                                                                             match(paste(x, "neg", sep = "_"), names(gt.merge)))]))))
  
  # summation
  alleleA_freq <- gt.merge[,24]+gt.merge[,26]
  alleleB_freq <- gt.merge[,25]+gt.merge[,27]
  gt.merge <- cbind(gt.merge,alleleA_freq,alleleB_freq)
  
  # remove rows where both values are 0s
  gt.merge <- gt.merge[!apply(gt.merge[,c(28,29)],1,function(x) x[1]==0 & x[2]==0),]
  
  # calculate chisquare/p-val
  gt.merge <- cbind(gt.merge,t(apply(gt.merge[,c(28,29)],1,function(x) with(chisq.test(x[1:2]),c(statistic,p.value=p.value))))) 
  
  # filter out depth below 20
  gt.merge <- gt.merge[which(gt.merge$depth>=20),] 
  
  # calculate fdr adjusted p-value
  gt.merge <- cbind(gt.merge,p.adjust(gt.merge$p.value,method="BH"))
  colnames(gt.merge)[32] <- "p_adjust"
  
  # write output
  write.csv(gt.merge[,c(1:10,12:15,28:32)],paste("~/allele_specific_expression/",paste(unique(arg1$sample_name),"_results.csv",sep=""),sep=""),row.names=F,quote=F) #write output to file
  
  cat(noquote(paste("\nDone...",unique(arg1$sample_name),sep="")))
  cat("\n")
}

# use ddply to call the function in a parallel manner
# dt is the genotype file
# split dt by sample_name (in my case 64 samples)
ddply(dt,.(sample_name),myfunction,.parallel=TRUE)

