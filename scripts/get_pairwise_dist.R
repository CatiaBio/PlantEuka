#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

library(ape, quietly=TRUE)
library(Biostrings, quietly=TRUE)


dna <- readDNAStringSet(sprintf("/home/projects/animalia_mito/data/%s/%s_prank.best.fas", args[1], args[1]))

dist_mat <- dist.dna(as.DNAbin(dna), model = "raw", as.matrix = T)



mu <- mean(dist_mat)
print(mu)
