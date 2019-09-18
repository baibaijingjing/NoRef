#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript fpkm_density_plot_func.r read_count.txt outdir")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) outdir: the dir for output")
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
require(edgeR)
require(ggplot2)


# read count data
print("read count data ...")
count_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names = FALSE)
head(count_data)
print("read count data is over")


# check geneLength key
if( !("geneLength" %in% colnames(count_data)) ) {
	stop("geneLength is not in read count file")
}


# calculate FPKM
sam_name <- colnames(count_data)[ !(colnames(count_data)%in%c("geneLength")) ]
fpkm <- count_data[ ,sam_name ]
log10_fpkm <- fpkm
for ( i in 1:(dim(fpkm)[2]) ){
	fpkm[,i] <- 10^9 * fpkm[,i] / sum(fpkm[,i]) / count_data[,"geneLength"]
	log10_fpkm[,i] <- log10(fpkm[,i])
}
print("calculate FPKM is over")


# init 
all <- NULL
all_sam <- NULL

# iter plot fpkm density
for( i in 1:length(sam_name) ){
	# fpkm
	c <- count_data[,sam_name[i]]
	keep <- c > 0
	r <- 10^9 * c[keep] / sum(c[keep]) / count_data[,"geneLength"][keep]
	log10fpkm <- data.frame(log10fpkm=log10(r))
	
	# update all
	all <- c(all, log10(r))
	all_sam <- c(all_sam, rep(sam_name[i], length(r)))
}


# create data.frame
log10fpkm <- data.frame(log10fpkm=all, sample=all_sam)
Sample <- factor(all_sam)
# plot fpkm box for all
# plot
m <- ggplot(log10fpkm, aes(factor(sample), log10fpkm))
p <- m + geom_boxplot(aes(fill=Sample)) + xlab("Sample") + ylab("log10(FPKM)")
p <- p + theme(axis.title.x = element_text(face="bold", size=14),
		axis.text.x  = element_text(face="bold", size=12),
		axis.title.y = element_text(face="bold", size=14),
		axis.text.y  = element_text(face="bold", size=12) )
#p <- p + theme(legend.title = element_text(face="bold", size=14),
#	legend.text = element_text(size=12),axis.text.x=element_text(angle = 30,vjust=1,hjust=1))
p=p+theme(legend.position='none',axis.text.x=element_text(angle = 30,vjust=1,hjust=1))
# remove background, horizontal & vertical grid lines 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# output
filepng <- paste(args[2], "/all", ".fpkm_box.png", sep="")
filepdf <- paste(args[2],"/all",".fpkm_box.pdf",sep="")
png(filename=filepng, height = 3000, width = 3400, res = 500, units = "px")
print(p)
dev.off()
pdf(file=filepdf,height=6,width=6)
print(p)
dev.off()







