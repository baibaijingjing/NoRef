#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript count_fpkm_cor.r in.txt out_prefix")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) pair.txt: the pairs of sample id for cor analysis")
	print("3) out_prefix: the prefix of output")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 ) {
	print(args)
	usage()
	stop("the length of args != 2")
}


# load library
require(ggplot2)


# read data
print("read  data ...")
df <- read.table(args[1], header=TRUE,sep="\t")
head(df)
colnames(df) <- c("Sample", "x", "y")
head(df)
print("read data is over")


# plot
#p <- ggplot(df, aes(x=x, y=y, group=Sample))
#p <- p + geom_line(aes(colour = Sample)) + labs(colour="Sample") + xlab("Relative Position in Genes(5'-3')")+ylab("Percent of Reads")
p <- qplot(x, y, data=df, geom="line", group=Sample, colour=Sample)
p <- p + labs(colour="Sample") + xlab("Relative Position in Genes(5'-3')")+ylab("Percent of Reads")
size <- 18
p <- p + theme(axis.title.x = element_text(face="bold", size=size),
	axis.text.x  = element_text(face="bold", size=size),
	axis.title.y = element_text(face="bold", size=size),
	axis.text.y  = element_text(face="bold", size=size) )
p <- p + theme(legend.title = element_text(face="bold", size=size),
	legend.text = element_text(size=size) )
# legend col
levels_num <- length( levels(df$Sample) )
legend_col <- as.integer( levels_num / 18 ) + 1
if( levels_num %% 18 == 0 ) {
	legend_col <- legend_col - 1
}
p <- p + guides(col = guide_legend(ncol = legend_col))


# output plot
png(filename=args[2], height = 3000, width = 6000 + (legend_col-1)*400, res = 500, units = "px")
print(p)
dev.off()



