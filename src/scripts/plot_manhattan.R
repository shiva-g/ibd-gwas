# manhattan plot
require(qqman)

args <- commandArgs(trailingOnly = TRUE)
assoc_file = args[1]
man_png_out = args[2]
qq_png_out = args[3]

results = read.delim(assoc_file, sep='\t', header=T)
png(man_png_out, width=2000, height=500)
manhattan(results)
dev.off()

png(qq_png_out, width=2000, height=500)
qq(results$P)
dev.off()