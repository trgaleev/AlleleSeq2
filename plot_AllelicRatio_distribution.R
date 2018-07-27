library(ggplot2)

args = (commandArgs(TRUE))

df = read.delim(args[1], header=TRUE, sep="\t", check.names=FALSE)

if (length(args) > 3){
	colmn = args[4]
	label = gsub('_', ' ', args[4]) 
}else{ 
	colmn = "ref_allele_ratio"
	label = "reference allele ratio"
}

p0 = ggplot(df, aes(get(colmn))) + 
	geom_histogram(bins=41) + 
	xlab(label) + 
#	ggtitle(paste0(gsub('_', ' ', args[2]), gsub('_', ' ', args[3]), sep=" ")) +
	ggtitle(paste0(args[2], "\n", gsub('\\.', ', ', gsub('_', ' ',args[3])))) + 
	theme_minimal()
	


#ggsave(paste0(args[2], "_", colmn, "s.", args[3], ".pdf"), p0, width=5, height=5)
pdf(file=paste0(args[2], "_", colmn, "s.", args[3], ".pdf"), width=5, height=5)
print(p0)
#dev.off()
