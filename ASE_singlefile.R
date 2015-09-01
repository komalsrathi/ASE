#######################################################################
cat(noquote("\nAllele Specific Expression:\n"))
cat(noquote("\nUsage:Rscript ASE_singlefile.R <infile-name.bam> <outfile-name.txt>"))

library(qdap)
library(data.table)    

# use arguments instead of hard coded file names
args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]

#read in genotype file
gt<-read.csv("genotype.csv")

#output bed file for input
write.table(gt[,c(9,11,11)],"input_test.bed",row.names=F,quote=F,sep="\t",col.names=F) 

#run bam-readcount program, filtered by base quality > 20
p <- paste("bam-readcount -b 20 -l input_test.bed -f hg19.fa",infile,">",outfile,sep=" ")
print(p)
system(p) 

#different columns per line bc of indels
no_col <- max(count.fields("bamreadcount_output_test.txt", sep = "\t")) 

#read in bamreadcount output
bamreadcounts <- read.table("bamreadcount_output_test.txt",sep="\t",fill=TRUE,col.names=1:no_col,check.names=F) 

#separate contents of columns 6-9 using based on ':'
x <- colsplit2df(bamreadcounts[,c(1:4,6:9)],c(5:8),sep=":") 

#extract relevant cols
x <- x[,c(1:4,10:11,24:25,38:39,52:53)] 

#change colnames
colnames(x) <- c("chrom","chrom_end","ref","depth","A_pos","A_neg","C_pos","C_neg","G_pos","G_neg","T_pos","T_neg") 

#merge the two files
gt.merge<-unique(merge(gt,x,by=c("chrom","chrom_end")))

#process the counts
gt.merge <- data.frame(gt.merge,as.data.frame(sapply(gt.merge[c("allelea", "alleleb")],
            function(x) as.integer(gt.merge[cbind(seq_len(nrow(gt.merge)),
            match(paste(x, "pos", sep = "_"), names(gt.merge)))]))))

gt.merge <- data.frame(gt.merge,as.data.frame(sapply(gt.merge[c("allelea", "alleleb")],
            function(x) as.integer(gt.merge[cbind(seq_len(nrow(gt.merge)),
            match(paste(x, "neg", sep = "_"), names(gt.merge)))]))))

#summation
alleleA_freq <- gt.merge[,24]+gt.merge[,26]
alleleB_freq <- gt.merge[,25]+gt.merge[,27]
gt.merge <- cbind(gt.merge,alleleA_freq,alleleB_freq)

#remove rows where both values are 0s
gt.merge <- gt.merge[!apply(gt.merge[,c(28,29)],1,function(x) x[1]==0 & x[2]==0),]

#calculate chisquare/p-val
gt.merge <- cbind(gt.merge,t(apply(gt.merge[,c(28,29)],1,function(x) with(chisq.test(x[1:2]),c(statistic,p.value=p.value))))) 

#filter out depth below 20
gt.merge <- gt.merge[which(gt.merge$depth>=20),] 

#calculate fdr adjusted p-value
gt.merge <- cbind(gt.merge,p.adjust(gt.merge$p.value,method="BH"))
colnames(gt.merge)[32] <- "p_adjust"

#write output
write.csv(gt.merge[,c(1:10,12:15,28:32)],"output_results_test.csv",row.names=F,quote=F)


