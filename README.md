# ASE 
ASE is a R script for allele specific expression in RNASeq data using bam-readcount.
After getting allele specific counts using bam-readcount, it performs a chisquare test to measure allelic imbalance.
It requires a BAM file and a genotype file as input.
Output is a text file with allele specific counts & chisq p-value determining if there is an allelic imbalance or not.
