# after merging output from all samples into one file called all_results_merged.csv
# use the below script to summate chi-square value and get a pvalue

library(plyr)
ase <- read.csv("all_results_merged.csv") # read in merged output files obtained from bamreadcount program
ase <- ase[order(ase$snpprobeset_id),] # order by snp id
ase <- ase[which(ase$alleleA_freq>=10 & ase$alleleB_freq>=10),] # both alleles should have atleast 10 reads

# summate the chi-square value and get a p-value
chisq.sum <- function(arg1)
{
  X.squared.sum <- sum(arg1$X.squared)
  n  <- length(unique(arg1$sample_name))
  if(n>1)
  {
    p.value <- pchisq(q=X.squared.sum,df=n-1,lower.tail=FALSE)
    dt <- rbind(data.frame(n,X.squared.sum,p.value))
    return(dt)
  }
}

# apply function by snp id
ase.res <- ddply(.data=ase,.variables="snpprobeset_id",.fun=chisq.sum)

# order by decreasing chi-square sum
ase.res <- ase.res[order(ase.res$X.squared.sum,decreasing=T),]

# write out
write.csv(ase.res,"ase_pchisq.csv",quote=F,row.names=F)

