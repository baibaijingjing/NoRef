#!/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript

# load library

options(digits=3, width=95)
library('getopt')

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'outdir' , 'o', 2, "character",
  'flag' , 'f', 1, "character",
  'rep' , 'X', 2, "character",
  'q' , 'q', 2, "double",
  'norm' , 'n', 2, "character",
  'width' , 'W', 2, "integer",
  'hight' , 'H', 2, "integer",
  'bio','p',2,"logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      Rscript NOISeqAnalysis.R --infile bayesian.value.xls --outfile manhattanplot.png 
      
      Options: 
      --help		-h 	NULL 		get this help
      --infile 	-i 	character 	the input file [forced]
      --flag 	  -f 	character 	sample to analyze use the format sample11_sample12_vs_sample21_sample22 [forced]
      --outdir 	-o 	character 	the directory for output graph [optional, default:current dir]
      --height 	-H 	integer 	the height of graph [optional, default: 600]
      --width 	-W 	integer 	the width of graph [optional, default:800]
      --rep 	  -r 	character 	In this argument, the type of replicates to be used is defined: 'technical', 'biological' or 'no' replicates. [optional,default: 'technical'].
      --q      	-q	curoff to filter [optional, default: 0.9]
      --bio     -b  use NOISeqBIO to analyze data [optional,default:True when data condition more than two, rep is biological]
      --norm    -n  Normalization method. It can be one of 'rpkm', 'uqua' (upper quartile), 'tmm' (trimmed mean of M) or 'n' (no normalization). [optional,default: 'n']
      \n")
  q(status=1);
}
opt = getopt(spec);

#setwd("E:\\R_workplace\\20151204NOIseq")
#opt<- data.frame(infile="All_gene_fpkm.list",outfile="test",sample="S.vs.NS",stringsAsFactors = F,flag="EK-1_EK-2_vs_WOX5-1_WOX5-2",rep="n")

# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$flag) )	{ print_usage(spec)}



#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$outdir) )	{ opt$outdir="NOISeq"}
#opt$outdir=trimws(opt$outdir)
if(!file.exists(opt$outdir)){
  dir.create(opt$outdir,recursive = TRUE)
}
if ( is.null(opt$height ) )		{ opt$height =600 }
if ( is.null(opt$width ) )		{ opt$width = 800}
if ( is.null(opt$norm ) )		{ opt$norm ='n'}

labsplit<-function(x,split="_"){
  tmp<-unlist(strsplit(x,split = split))
  tmp=tmp[tmp!=""]
  tmp=sapply(tmp, function(x){x=sub("^_*",replacement = "",x = x);x=sub("_*$",replacement = "",x = x)},simplify = "vector")
  return(tmp)
}

group=labsplit(opt$flag,split = "vs")
names(group)<-group
group_lab<- sapply(group,labsplit,simplify = F)
sample<-as.vector(unlist(group_lab))
valCol=sample
geneCol=1
library(methods)
rt1 <- read.table(opt$infile, header = T, sep = "\t", row.names = NULL,comment.char = "",check.names = F)
exp_values <- as(rt1[valCol], "matrix")
exp_values[is.na(exp_values)] <- 0
Gene_names <- as(rt1[geneCol], "matrix")
dimnames(exp_values) <- list(as.character(Gene_names), dimnames(exp_values)[[2]])
if (!is.numeric(exp_values)) {
  cat("Some invalid values in the matrix!\n")
  cat("Please check the input file!\n")
}

mylength=""
if(is.element("geneLength",colnames(rt1))){
  mylength="geneLength"
  mylength<- as(rt1[mylength], "matrix")
  mylength[is.na(mylength)] <- 0
  Gene_names <- as(rt1[geneCol], "matrix")
  dimnames(mylength) <- list(as.character(Gene_names), dimnames(mylength)[[2]])
  if (!is.numeric(exp_values)) {
    cat("Some invalid values in the matrix!\n")
    cat("Please check the input file!\n")
  }
}

if(opt$norm!='n'){
  k=NULL
}else{
  k=0.5
}
###################################################
### output degenes to file
###################################################
output<-function(x,file){
  if(ncol(x)>0){
    output=data.frame(ID=rownames(x),x)
  }else{
    output=c("ID",colnames(x))
  }
  write.table(output,file = file,col.names = TRUE,row.names =FALSE,quote = FALSE,sep="\t")
}
library(NOISeq)
###################################################
### code chunk number 3: NOISeq.Rnw:88-89
###################################################
group<-names(group_lab)
if(length(group_lab)<2)stop("[NOISeq Error]:Group number less than 2!")
for(i in length(group_lab):2){
  sample1=as.vector(group_lab[[i]])
  lab1=rep(group[i],length(sample1))
  for(j in (i-1):1){
    sample2<-as.vector(group_lab[[j]])
    lab2=rep(group[j],length(sample2))
    lab=c(lab1,lab2)
    mycounts=exp_values[,c(sample1,sample2)]
    myfactors=data.frame(Group =lab,
                         GroupRun = c(sample1,sample2))
   
    ###################################################
    ### code chunk number 6: readData
    ###################################################
    mydata <- readData(data=mycounts, factors=myfactors)
    if(length(mylength)>0||mylength!=""){
      mydata <- addData(mydata, length=mylength)
    }
    if(length(sample1)>1&&length(sample2)>1){
      if(is.null(opt$rep))opt$rep="technical"
      mynoiseq = noiseq(mydata, k = k, norm = opt$norm, factor="Group",
                        nss = 0, lc = 1, replicates = opt$rep)
      if(is.null(opt$q))opt$q=0.9
    }else{
      mynoiseq = noiseq(mydata, k = k, norm = opt$norm, factor="Group", pnr = 0.2, 
                        nss = 5, v = 0.02, lc = 1, replicates = "no")
      if(is.null(opt$q))opt$q=0.9
    }
    outdir=paste(c(opt$outdir,"/",group[i],"_vs_",group[j]),collapse = "")
    if (!file.exists(outdir)) {
      dir.create(outdir)
    }
    mynoiseq.deg = degenes(mynoiseq, q = opt$q, M = NULL)
    mynoiseq.up.deg = degenes(mynoiseq, q = opt$q, M = "up")
    mynoiseq.down.deg = degenes(mynoiseq, q = opt$q, M = "down")
    output(x=mynoiseq.deg,file = paste(outdir,"deg.list",sep="/"))
    output(x=mynoiseq.up.deg,file = paste(outdir,"deg.up.list",sep="/"))
    output(x=mynoiseq.down.deg,file = paste(outdir,"deg.down.list",sep="/"))
    
    ###################################################
    ### code chunk number 36: fig_summ_expr
    ###################################################
    png(paste(outdir,"deg.expr.png",sep="/"),width = opt$width,height=opt$height)
    DE.plot(mynoiseq, q = opt$q, graphic = "expr", log.scale = FALSE)
    dev.off()
    
    ###################################################
    ### code chunk number 37: fig_summ_MD
    ###################################################
    png(paste(outdir,"deg.MD.png",sep="/"),width = opt$width,height=opt$height)
    DE.plot(mynoiseq, q =  opt$q, graphic = "MD")
    dev.off()
    if(is.null(opt$bio)){
      if(opt$rep=="biological"){
        mynoiseq =  noiseqbio(mydata, k = NULL, norm = "n", factor="Group", lc = 1,
                              r = 50, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345,
                              filter = 1)
        outdir=paste(c(opt$outdir,"/NOISeqBIO/",group[i],"_vs_",group[j]),collapse = "")
        if (!file.exists(outdir)) {
          dir.create(outdir,recursive = TRUE)
        }
        mynoiseq.deg = degenes(mynoiseq, q = opt$q, M = NULL)
        mynoiseq.up.deg = degenes(mynoiseq, q = opt$q, M = "up")
        mynoiseq.down.deg = degenes(mynoiseq, q = opt$q, M = "down")
        output(x=mynoiseq.deg,file = paste(outdir,"deg.list",sep="/"))
        output(x=mynoiseq.up.deg,file = paste(outdir,"deg.up.list",sep="/"))
        output(x=mynoiseq.down.deg,file = paste(outdir,"deg.down.list",sep="/"))
        
      }
    }
    
  }
}







