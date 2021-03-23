#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
aln_file <- args[1]
seq_name <- args[2]
out_file <- args[3]

library("seqinr")
library("phangorn")

aln <- read.alignment(file = aln_file, format = "fasta")
dat <- as.phyDat(aln, type="AA")
dm <- dist.hamming(dat)
tree <- NJ(dm)
score<-sankoff(tree, dat)

write.table( t(c(seq_name, score)),  
             file=out_file, 
             append = T, 
             sep=';',
             quote = F,
             row.names=F, 
             col.names=F )
