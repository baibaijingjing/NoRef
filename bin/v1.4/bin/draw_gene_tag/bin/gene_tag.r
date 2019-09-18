#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript gene_tag.r Total.gene_tag.list Total.gene_tag.png")
	print("1) Total.gene_tag.list: Basic_Analysis/geneExpression/Total.gene_tag.list")
	print("2) Total.gene_tag.png: output file")
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
#p <- p + geom_line(aes(colour = Sample)) + labs(colour="Sample") + xlab("Total Tag Number (M)")+ylab("Total Gene Number (K)")+ggtitle("Saturation Curve(FPKM>=0.1)")
p <- qplot(x, y, data=df, geom="line", group=Sample, colour=Sample)
p <- p + labs(colour="Sample") + xlab("Total Reads Number (M)")+ylab("Total Genes Number (K)")+ggtitle("Saturation Curve")
size <- 18
p <- p + theme(axis.title.x = element_text(face="bold", size=size),
	axis.text.x  = element_text(face="bold", size=size),
	axis.title.y = element_text(face="bold", size=size),
	axis.text.y  = element_text(face="bold", size=size) )
p <- p + theme(legend.title = element_text(face="bold", size=size),
	legend.text = element_text(size=size) )

# remove background, horizontal & vertical grid lines 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# legend col
levels_num <- length( levels(df$Sample) )
legend_col <- as.integer( levels_num / 18 ) + 1
if( levels_num %% 18 == 0 ) {
	legend_col <- legend_col - 1
}
p <- p + guides(col = guide_legend(ncol = legend_col))


# output plot
png(filename=args[2], height = 3000, width = 5000 + (legend_col-1)*400, res = 500, units = "px")
print(p)
dev.off()


