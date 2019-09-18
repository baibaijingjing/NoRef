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
df <- read.table(args[1], header=TRUE)
head(df)
print("read data is over")


# plot
p <- ggplot(df, aes(x=Site, y=Percent, group=Type))
p <- p + geom_line(aes(colour = Type)) + labs(colour="Type")
size <- 18
p <- p + theme(axis.title.x = element_text(face="bold", size=size),
	axis.text.x  = element_text(face="bold", size=size),
	axis.title.y = element_text(face="bold", size=size),
	axis.text.y  = element_text(face="bold", size=size) )
p <- p + theme(legend.title = element_text(face="bold", size=size),
	legend.text = element_text(size=size) )
p <- p + theme(legend.position=c(.7, .85))
# remove background, horizontal & vertical grid lines 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# output plot
png(filename=args[2], height = 3000, width = 6000, res = 500, units = "px")
print(p)
dev.off()



