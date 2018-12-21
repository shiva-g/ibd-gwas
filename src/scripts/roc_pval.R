# roc pval
require(verification)

args <- commandArgs(trailingOnly = TRUE)
dat_file = args[1]
out = args[2]

results = read.delim(assoc_file, sep='\t', header=T)
p = roc.area(df$y, df$PRS)$p.value