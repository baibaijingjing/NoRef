#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript

### include DGE plot functions
source("dge_utils_funcs.r")


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript ebseq_analysis.r read_count.txt S1,S2,S3 out.de")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) S1,S2,S3: the sample id for condition")
        print("3) 0.01: FDR cutoff")
        print("4) 0: the cutoff of fold change(FC)")
        print("5) out.de: the output of DE")
        print("6) fpkm.txt: the fpkm file for RNA-SEQ")
	print("-------------------------------------------------------------------------------")
}


# only process factor without replicates
get_sample_name <- function(str=NULL) {
	# check 
	if( is.null(str) ) stop("str is NULL")
 
	# abstract factor
	s <- unlist( strsplit(str, ",") )

	# check length
	if( length(s) < 2 ) stop(paste("length(s) < 2: ", str, sep=""))

	# return
	return(s)
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 6 ) {
	print(args)
	usage()
	stop("the length of args != 6")
}



# load library
require(edgeR)
require(EBSeq)
require(ggplot2)


# abstract sample name
sample_name <- get_sample_name(args[2])


# read count data
print("read count data ...")
ori_data <- read.delim(args[1], row.names = 1, header=TRUE,check.names =F,comment.char="")
ori_fpkm <- read.delim(args[6], row.names = 1, header=TRUE,check.names =F,comment.char="")

fpkm_data<-ori_fpkm[,sample_name]
head(ori_data)
count_data <- ori_data[,sample_name]
head(count_data)
print("read count data is over")
head(ori_fpkm)
head(fpkm_data)

# create DEG_list object
de <- new_DGE_list(counts=count_data, fpkm=fpkm_data,condition=as.factor(sample_name), geneLength=ori_data[,"geneLength"],filter="cpm")
print("create DEG_list object is over")
cat("the number of origin gene is ", dim(de$counts)[1], "\n")
cat("the number of gene (after filtering) is ", dim(de$filterCounts)[1], "\n")

# calculate the FC cutoff
fc_cutoff <- ifelse( as.numeric(args[4])<=1, 0, log2(as.numeric(args[4])) )
print(fc_cutoff)

# Library size factor
Sizes=MedianNorm(de$filterCounts)

######## two level
if( length(sample_name) == 2 ) {
# gene expression estimates
# do estimate
#EBOut <- EBTest(Data=de$filterCounts, Conditions=as.factor(sample_name), sizeFactors=Sizes, maxround=5)
EBOut <- EBTest(Data=de$filterCounts, Conditions=as.factor(sample_name), sizeFactors=Sizes, maxround=5, Qtrm=0.5, QtrmCut=0)
# get pp matrix
PP <- GetPPMat(EBOut)
# get post FC
GeneFC <- PostFC(EBOut)
# print
head(GeneFC$PostFC)

# select DEG
isDGE <- (PP[,"PPDE"]>=1-as.numeric(args[3])) & (abs(log2(GeneFC$PostFC)) > fc_cutoff)
cat("de: ", length(isDGE), "\n")
cat("PPDE = ", length(PP[,"PPDE"]), "\n")
cat("GeneFC$PostFC = ", length(GeneFC$PostFC), "\n")


# check compare direction
C1_name <- strsplit(GeneFC$Direction, " ")[[1]][1]
log2_PostFC <- NULL
if( C1_name == sample_name[1] ){
	log2_PostFC <- -log2(GeneFC$PostFC)
} else {
	log2_PostFC <- log2(GeneFC$PostFC)
}

# update DGE_list object
de <- update_DGE_list(de, FDR=1-PP[,"PPDE"], log2FC=log2_PostFC,isDGE=isDGE)


# output DGE
out <- output_DGE_list(de, de_file=paste0(args[5],".DEG_final.xls"), 
        all_file=paste0(args[5],".final.xls"), cluster_file=paste0(args[5],".DEG_final.xls.forheatmap"))
} 



