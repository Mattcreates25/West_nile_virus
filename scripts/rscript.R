#go into the correct working directory
setwd("~/WNV/results")

#load the library seqinr which has the function kaks
library(seqinr)


all_alignment= read.alignment(file = "all_alignment.fasta", format = "fasta", whole.header= T)

kaks(all_alignment, verbose = F, debug = F, forceUpperCase = T, rmgap = T)


