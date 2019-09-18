#!/share/nas2/genome/biosoft/R/test/3.1.1/bin/Rscript

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
# load library
library('getopt');

spec = matrix(c(
                'help',   'h', 0, 'logical',  'print help document, and exit',
                'defile', 'd', 1, 'character','file, different expression analysis result, required',
                'outfile','o', 2, 'character','file, volcano plot output file, [./volcano_plot.png]',
                'fc',     'f', 2, 'integer',  'integer, fold change threshold, [2]',
                'fdr',    'q', 2, 'double',   'double, false discovery rate threshold, [0.05]',
                'height', 't', 2, 'integer',  'integer, height of plot, [3000]',
                'width',  'w', 2, 'integer',  'integer, width of plot, [3000]',
                'color',  'c', 2, 'character','colors of point, split by comma, [green,blue,red]'
                ), byrow=TRUE, ncol=5);
# getopt usage
# Col3: 0=no argument, 1=required argument, 2=optional argument.
# Col4: logical, integer, double, complex, character.
opt = getopt(spec);

# define usage function
print_usage <- function(spec=NULL){
    script <- get_Rscript_filename()
#   script <- unlist(strsplit(get_Rscript_filename(), "/"))[length(unlist(strsplit(get_Rscript_filename(), "/")))]
    cat(getopt(spec, usage=TRUE));
    cat("Example: \nRscript", script, "-d T1_T2_vs_T3_T4.final.xls \n")
    cat("Rscript", script, "-d T1_T2_vs_T3_T4.final.xls -o T1_T2_vs_T3_T4.FC_FDR.png -c 4 -q 0.01 -t 800 -w 800 -c green,black,red \n")
    q(status=1);
}

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }

# check required args 
if ( is.null(opt$defile) )	{ print_usage(spec) }

# check & init optional args 
if ( is.null(opt$outfile) ) { opt$outfile=paste(getwd(), "/volcano_plot.png", sep="")}
if (is.null(opt$fc)) { opt$fc=2 }
if (is.null(opt$fdr)) { opt$fdr=0.05 }
if (opt$fc < 1 || opt$fdr >= 1) { cat("ERROR: Illegal option --fc or --fdr !\n"); print_usage(spec) }

if (is.null(opt$height)) { opt$height = 3000 }
if (is.null(opt$width)) { opt$width = 3000 }
if (opt$height <= 0 || opt$width <= 1) { cat("ERROR: Illegal option --height or --width !\n"); print_usage(spec) }

if (is.null(opt$color)) { opt$color = 'green,blue,red' }
opt$color <- unlist(strsplit(opt$color, ",")) # split opt$color by comma
if (length(opt$color) < 3) { cat("ERROR: Illegal option --color !\n"); print_usage(spec) }

cat("[", format(Sys.time()), "]", script<-unlist(strsplit(get_Rscript_filename(), "/"))[length(unlist(strsplit(get_Rscript_filename(), "/")))], "start...\n")
#-----------------------------------------------------------------
# main 
#-----------------------------------------------------------------
# read data
der <- read.delim(opt$defile, row.names = 1, header=TRUE,check.names = FALSE)
der[,ncol(der)]<-as.character(der[,ncol(der)])

# check input
if ( (colnames(der)[ncol(der)] != "regulated") || (colnames(der)[ncol(der)-1] != "log2FC") || (colnames(der)[ncol(der)-2] != "FDR") ) {
    cat("ERROR: Illegal input file, please check opt$defile !\n");
    print_usage(spec)
}

# switch column "regulated" 
for (i in 1:nrow(der) ) {
    if ( as.double(der[i,ncol(der)-1]) >= opt$fc && as.double(der[i,ncol(der)-2]) <= opt$fdr ) {
        der[i,ncol(der)] <- "up"
    } else if ( as.double(der[i,ncol(der)-1]) <= 0-opt$fc && as.double(der[i,ncol(der)-2]) <= opt$fdr ) {
        der[i,ncol(der)] <- "down"
    } else {
        der[i,ncol(der)] <- "normal"
    }
}

# load library 
library(ggplot2)

# draw volcano plot 
png(filename=opt$outfile, height = opt$height, width = opt$width, res = opt$height/6, units = "px")

p <- ggplot(der, aes(x=log2FC, y=-log10(FDR), colour=regulated) ) + geom_point(size=1)
p <- p + scale_colour_manual(values=opt$color)

# add dash threshold line
p <- p + geom_vline(xintercept=c(-opt$fc,opt$fc), linetype="longdash", size=0.2)
p <- p + geom_hline(yintercept=c(-log10(opt$fdr)), linetype="longdash", size=0.2) 

# add title
p <- p + labs(list(title="Volcano Plot", x="log2(FC)"))

print(p)
dev.off()

# 
cat("[", format(Sys.time()), "]", script<-unlist(strsplit(get_Rscript_filename(), "/"))[length(unlist(strsplit(get_Rscript_filename(), "/")))], "done.\n")

