### load library
require(edgeR)

#################################################################
### create DGE_list object
### here it will filter low count gene
#################################################################
new_DGE_list <- function(counts=NULL,fpkm=NULL, condition=NULL, geneLength=NULL,filter=filter){
	# check NULL
	if( is.null(counts) ) stop("counts is NULL")
	if( is.null(fpkm) ) stop("fpkm is NULL")
	if( is.null(condition) ) stop("condition is NULL")
	if( is.null(geneLength) ) stop("geneLength is NULL")

	# check length
	if( dim(counts)[2] != length(condition) )
		stop(paste("dim(counts)[2] != length(condition): ", 
			dim(counts)[2], " != ", length(condition), sep=""))
	if( dim(counts)[1] != length(geneLength) )
		stop(paste("dim(counts)[1] != length(geneLength): ", 
			dim(counts)[1], " != ", length(geneLength), sep=""))
	if( dim(counts)[2] < 2 )
		stop(paste("dim(counts)[2] < 2: ", dim(counts)[2], sep=""))

	# create DGE_list object
	de <- list(counts=as.matrix(counts),FPKM=as.matrix(fpkm), condition=condition, geneLength=geneLength,
		filterCounts=NULL, filterLength=NULL, fpkm=NULL, FDR=NULL, log2FC=NULL, isDGE=NULL,
		regulated=NULL)
    # filter low expression gene
    if(filter=="cpm"){
        # filter low expression gene
        keep <- rowSums(cpm(de$counts)>1) >= min( table(condition) )
        }else if(filter=="count"){
        keep <- rowSums(de$counts)>=2
        }else{print("You must choose one filter method!")}

	
	
	de$filterCounts <- de$counts[keep,]
	de$filterLength <- de$geneLength[keep]

	# calculate the fpkm
	print(head(de$FPKM))
	de$fpkm <- de$FPKM[keep,]
	print(head(de$fpkm))
	# return
	return(de)
}


#################################################################
### update DGE_list object
### here it update the dge information according to DGE analysis result
#################################################################
update_DGE_list <- function(de=NULL, isDGE=NULL, FDR=NULL, log2FC=NULL){
	# check null
	# note: log2FC doesnot need to check, which is only for two condition
	if( is.null(de) ) stop("de is NULL")
	if( is.null(isDGE) ) stop("isDGE is NULL")
	if( is.null(FDR) ) stop("FDR is NULL")
	# update
	de$isDGE = isDGE;
	de$FDR = FDR;
	de$log2FC = log2FC;

	# update the regulated
	if( is.null(log2FC) ){
		de$regulated <- rep("up", length(de$isDGE))
                names(de$regulated) <- ( rownames(de$fpkm) )
                }else{
		de$regulated <- ifelse(de$log2FC < 0, "down", "up")
                names(de$regulated) <- ( rownames(de$fpkm) )
                de$regulated[which(!isDGE)] <- "normal"
        }
	# set name
	# NOTE: there may be only one DEG, so must first: ( rownames(de$fpkm) ), then [de$isDGE]
	#names(de$regulated) <- rownames(de$fpkm[de$isDGE,]) 	# NOTE: error for one DEG
        #names(de$regulated) <- ( rownames(de$fpkm) )
	

	# return
	return(de)
}


#################################################################
### output DGE_list object
### here it will output the result of DGE
#################################################################
output_DGE_list <- function(de=NULL, de_file=NULL, all_file=NULL, cluster_file=NULL){
	# check null
	if( is.null(de) ) stop("de is NULL")
	if( is.null(de_file) ) stop("de_file is NULL")
	if( is.null(all_file) ) stop("all_file is NULL")
	if( is.null(cluster_file) ) stop("cluster_file is NULL")

	# check the number of DEG
	if( sum(de$isDGE) == 0 ) {
		print("The number of DEG is 0")
		df <- data.frame(zero=c("The number of DEG is 0"))
		write.table(df, file=paste(de_file, ".zero", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
		return(0)
	}

	# output DEG fpkm for cluster
	df <- NULL
	if( length(de$isDGE) == 1 ) { 	# only get one DEG
		print("##################")
		print(names(de$regulated[de$isDGE]))
		print("##################")
		df <- data.frame(Gene_lib=names(de$regulated[de$isDGE]),matrix(de$fpkm[de$isDGE,],1,))
	} else {
		df <- data.frame(Gene_lib=names(de$regulated[de$isDGE]),de$fpkm[de$isDGE,])
	}
        count_name<-paste(colnames(de$fpkm),"Count",sep="_")
	fpkm_name<-paste(colnames(de$fpkm),"FPKM",sep="_")
	colnames(df) <- c("#ID",colnames(de$fpkm))
	write.table(df, file=cluster_file, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

	# for de
	if( is.null(de$log2FC) ) {

		if( length(de$isDGE) == 1 ) { 	# only get one DEG
			df <- data.frame(id=names(de$regulated[de$isDGE]),matrix(de$filterCounts[de$isDGE,],1,),matrix(de$fpkm[de$isDGE,],1,),FDR=matrix(de$FDR[de$isDGE],1,),regulated=de$regulated[de$isDGE])
			colnames(df) <- c("#ID",count_name,fpkm_name,"FDR", "regulated")
		} else {
			df <- data.frame(id=names(de$regulated[de$isDGE]),de$filterCounts[de$isDGE,],de$fpkm[de$isDGE,],FDR=de$FDR[de$isDGE], regulated=de$regulated[de$isDGE])
			colnames(df) <- c("#ID",count_name,fpkm_name,"FDR", "regulated")
		}

	} else {
		if( length(de$isDGE) == 1 ) { 	# only get one DEG
			df <- data.frame(id=names(de$regulated[de$isDGE]),matrix(de$filterCounts[de$isDGE,],1,),matrix(de$fpkm[de$isDGE,],1,),FDR=matrix(de$FDR[de$isDGE],1,),log2FC=matrix(de$log2FC[de$isDGE],1,),regulated=de$regulated[de$isDGE])
			colnames(df) <- c("#ID", count_name,fpkm_name,"FDR", "log2FC", "regulated")
		} else {
			df <- data.frame(id=names(de$regulated[de$isDGE]),de$filterCounts[de$isDGE,],de$fpkm[de$isDGE,],FDR=de$FDR[de$isDGE], log2FC=de$log2FC[de$isDGE],regulated=de$regulated[de$isDGE])
			colnames(df) <- c("#ID", count_name,fpkm_name,"FDR", "log2FC", "regulated")
		}
	}
	write.table(df, file=de_file, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

	# for all
	if( is.null(de$log2FC) ){
		df <- data.frame(id=rownames(de$fpkm),de$fpkm, FDR=de$FDR,regulated=de$regulated)
		colnames(df) <- c("#ID",colnames(de$fpkm),"FDR","regulated")
	}else{
		df <- data.frame(id=rownames(de$fpkm),de$fpkm,FDR=de$FDR, log2FC=de$log2FC,regulated=de$regulated)
		colnames(df) <- c("#ID", colnames(de$fpkm),"FDR", "log2FC","regulated")
	}
	write.table(df, file=all_file, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}




