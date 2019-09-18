#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript


# load library
library('getopt');
library(pheatmapm)
library(RColorBrewer)
library(NbClust)
library(cluster)
#require(ggplot2)


#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'outfile' , 'o', 1, "character",
	'file.type' , 'f', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'show_rownames' , 'a', 0, "logical",
	'show_colnames' , 'b', 0, "logical",
	'cluster_rows' , 'c', 0, "logical",
	'cluster_cols' , 'd', 0, "logical",
	'treeheight_row' , 'e', 1, "integer",
	'treeheight_col' , 'g', 1, "integer",
	'legend' , 'j', 0, "logical",
	'fontsize' , 'k', 1, "integer",
	'fontsize_row' , 'm', 1, "integer",
	'fontsize_col' , 'n', 1, "integer",
	'cellwidth' , 'p', 1, "integer",
	'cellheight' , 'q', 1, "integer",
	'is.log' , 'J', 0, "logical",
	'color.type' , 'r', 1, "integer",
	'scale' , 's', 1, "character",
	'div.clust' , 'v', 0, "logical"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
1) Rscript pheatmap.r --infile in_pheatmap.data --outfile out_pheatmap.png \\
	--file.type 1 --show_rownames --show_colnames --cluster_rows --cluster_cols \\
	--legend
2) Rscript pheatmap.r --infile in_pheatmap.data --outfile out_pheatmap.png \\
	--file.type 1 --show_colnames --cluster_rows --cluster_cols \\
	--legend
3) Rscript pheatmap.r --infile in_pheatmap.data --outfile out_pheatmap.png \\
	--file.type 1 --show_rownames --show_colnames --cluster_cols \\
	--legend --treeheight_row 50 --treeheight_col 20
4) Rscript pheatmap.r --infile in_pheatmap.data --outfile out_pheatmap.png \\
	--file.type 1 --show_rownames --show_colnames --cluster_rows --cluster_cols \\
	--legend --fontsize_row 3 --fontsize_col 9
5) Rscript pheatmap.r --infile in_pheatmap.data --outfile out_pheatmap.png \\
	--file.type 1 --show_rownames --show_colnames --cluster_rows --cluster_cols \\
	--legend --cellwidth 16 --cellheight 5
6) Rscript pheatmap.r --infile in_pheatmap.data --outfile out_pheatmap.png \\
	--file.type 1 --show_rownames --show_colnames --cluster_rows --cluster_cols \\
	--legend --height 5000 --width 3000

Options: 
--help		-h 	NULL 		get this help
--infile 		character 	the input file [forced]
--outfile 		character 	the filename for output graph [forced]
--file.type 		integer 	the type of file, must be 1 or 2 [optional, default: 1]
--height 		integer 	the height of graph [optional, default: 5000]
--width 		integer 	the width of graph [optional, default: 3000]
--show_rownames		logical 	show rownames [optional, default: FALSE]
--show_colnames		logical 	show colnames [optional, default: FALSE]
--cluster_rows		logical 	show row cluster tree [optional, default: FALSE]
--cluster_cols		logical 	show col cluster tree [optional, default: FALSE]
--treeheight_row 	integer 	the height of row cluster tree [optional, default: 50]
--treeheight_col 	integer 	the height of col cluster tree [optional, default: 50]
--legend		logical 	show legend [optional, default: FALSE]
--fontsize 		integer 	the size of font [optional, default: NULL]
--fontsize_row 		integer 	the size of font for rowname [optional, default: NULL]
--fontsize_col 		integer 	the size of font for colname [optional, default: NULL]
--cellwidth 		integer 	the width of cell [optional, default: NA]
--cellheight 		integer 	the height of cell [optional, default: NA]
--color.type 		integer 	the type of color [optional, default: 1]
--is.log		logical 	log(data) [optional, default: FALSE]
--scale			character	row or column or none [optional, default: none]
--div.clust		logical		output div result of all cluster [optional, default: FALSE]
\n")
	q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$file.type ) )		{ opt$file.type = 1 }
if ( is.null(opt$height ) )		{ opt$height = 5000 }
if ( is.null(opt$width ) )		{ opt$width = 3000 }
if ( is.null(opt$show_rownames ) )	{ opt$show_rownames = FALSE }
if ( is.null(opt$show_colnames ) )	{ opt$show_colnames = FALSE }
if ( is.null(opt$cluster_rows ) )	{ opt$cluster_rows = FALSE }
if ( is.null(opt$cluster_cols ) )	{ opt$cluster_cols = FALSE }
if ( is.null(opt$treeheight_row ) )	{ opt$treeheight_row = ifelse(opt$cluster_rows, 50, 0) }
if ( is.null(opt$treeheight_col ) )	{ opt$treeheight_col = ifelse(opt$cluster_cols, 50, 0) }
if ( is.null(opt$legend) )		{ opt$legend = FALSE }
if ( is.null(opt$fontsize ) )		{ opt$fontsize = 8 }
if ( is.null(opt$fontsize_row ) )	{ opt$fontsize_row = NULL }
if ( is.null(opt$fontsize_col ) )	{ opt$fontsize_col = NULL }
if ( is.null(opt$cellwidth ) )		{ opt$cellwidth = NA }
if ( is.null(opt$cellheight ) )		{ opt$cellheight = NA }
if ( is.null(opt$color.type ) )		{ opt$color.type = 1 }
if ( is.null(opt$is.log) )		{ opt$is.log = FALSE }
if ( is.null(opt$scale) )		{ opt$scale = "none" }
if ( is.null(opt$div.clust))		{ opt$div.clust = FALSE }

# check file type
if( (opt$file.type != 1) && (opt$file.type != 2) ){
	cat("Final Error: the type of file must be 1 or 2\n")
	print_usage(spec)
}
# check scale
if( (opt$scale != "column") && (opt$scale != "row") && (opt$scale != "none") ){
	cat("Final Error: the scale must be column or row or none\n")
	print_usage(spec)
}


#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
if( opt$file.type == 1 ){
	data <- read.delim(opt$infile, row.names = 1, header=TRUE ,check.names = F)
} else {
	data <- read.table(opt$infile, row.names = 1, header=TRUE ,check.names = F)
}
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
# log
if( opt$is.log ){
	# check value < 0
	if( sum(data<0) > 0 ) {
		cat("Final Error: there exist some values in data < 0, option is.log ERROR\n")
		print_usage(spec)
	}
	# log
	data <- log10(data+0.000001)
}
if( opt$div.clust){
	out_file = paste(opt$outfile,".clustered.data",sep="")
} else {
	out_file = NA
}


#-----------------------------------------------------------------
# plot
# output plot
#-----------------------------------------------------------------
# graph size
height <- opt$height*2/1000
width <- opt$width*2/1000
# color
# check
if ( (opt$color.type < 1) || (opt$color.type > 6) ) {
	cat("Final Error: color.type must be 1-6\n")
	print_usage(spec)
}
# set
color.1 <- colorRampPalette(rev(c("#ff0000", "#000000", "#00ff00")))(100)
color.2 <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", 
	"#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4", rep("#4575B4",2))))(100)
color.3 <- colorRampPalette(brewer.pal(3,"Set1"))(100)
color.4 <- colorRampPalette(brewer.pal(3,"Set2"))(100)
color.5 <- colorRampPalette(brewer.pal(3,"Set3"))(100)
color.6 <- colorRampPalette(rev(c("#ff0000", "#ffffff",  "#0000ff")))(100)
color.set <- list(color.1, color.2, color.3, color.4, color.5, color.6)

# output pdf
pheatmapm(data, scale=opt$scale, height=height, width=width, color=color.set[[opt$color.type]], 
	show_rownames=opt$show_rownames, show_colnames=opt$show_colnames,
	cluster_rows=opt$cluster_rows, cluster_cols=opt$cluster_cols,
	treeheight_row=opt$treeheight_row, treeheight_col=opt$treeheight_col,
	legend=opt$legend, fontsize=opt$fontsize, 
	fontsize_row=opt$fontsize_row, fontsize_col=opt$fontsize_col,
	cellwidth=opt$cellwidth, cellheight=opt$cellheight,
	filename=paste(opt$outfile,".pdf",sep=""))
# output png
png(filename=opt$outfile, height=opt$height, width=opt$width, res=500, units="px")
pheatmapm(data, out_file=out_file, scale=opt$scale, height=height, width=width, color=color.set[[opt$color.type]],
	show_rownames=opt$show_rownames, show_colnames=opt$show_colnames,
	cluster_rows=opt$cluster_rows, cluster_cols=opt$cluster_cols,
	treeheight_row=opt$treeheight_row, treeheight_col=opt$treeheight_col,
	legend=opt$legend, fontsize=opt$fontsize, 
	fontsize_row=opt$fontsize_row, fontsize_col=opt$fontsize_col,
	cellwidth=opt$cellwidth, cellheight=opt$cellheight)
#file_path=paste(opt$outfile,sep="",collapse = "/")
#file_path
dir_file=as.vector(strsplit(opt$outfile,split = "/"))
dir_file_f=dir_file[[1]]
dir_file_f2=dir_file_f[-length(dir_file_f)]
dir_file_path=paste(dir_file_f2,sep="",collapse = "/")
dir_file_path
cmd=paste("cp ", dirname(get_Rscript_filename()), "/readme.txt ", dir_file_path, sep="")
#cmd=paste("cp /share/nas1/licq/testing/toolkits/20_tookits/heatmap_draw/v1.0/readme.txt",dir_file_path)
dev.off()
system(cmd)






