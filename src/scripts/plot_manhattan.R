# manhattan plot
require(qqman)

args <- commandArgs(trailingOnly = TRUE)
assoc_file = args[1]
png_out = args[2]

results = read.delim(assoc_file, sep='\t', header=T)
#results = na.omit(results)
#head(results)
png(png_out, width=2000, height=500)
manhattan(results)
dev.off()
