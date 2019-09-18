#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($All_GO,$DEG_list,$key,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"All_GO:s"=>\$All_GO,
				"DEG_list:s"=>\$DEG_list,
				"key:s"=>\$key,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($All_GO and $DEG_list and $key and $od);

mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);

my %H;
open (IN,"$DEG_list") or die $!;
while (<IN>) {
	next if /^\#/;
	my @A=split/\s+/,$_;
	my $gene=$A[0];
	my $val=(@A>=4)?$A[-3]:0.01;
	$H{$gene}=$val;
}
close IN;

open (IN,"$All_GO") or die $!;
open (OUT,">$od/topGO.map") or die $!;
open (OUT2,">$od/topGO.list") or die $!;
print OUT2 "#ID\tValue\n";
while (<IN>) {
	chomp;
	my ($name,$info)=(split/\t/,$_,2)[0,1];
	next unless defined $info;
	$info=~s/\t/,/g;
	print OUT "$name\t$info\n";
	print OUT2 "$name\t$H{$name}\n" if exists $H{$name};
	print OUT2 "$name\t1\n" unless exists $H{$name};
}
close IN;
close OUT;

open (OUT,">$od/topGO.R") or die $!;
print OUT '#!usr/bin/env Rscript',"\n";
print OUT 'library(topGO)'."\n";
print OUT 'library(ALL)'."\n";
print OUT 'data(ALL)'."\n";
print OUT 'geneID2GO<-readMappings(file = "'."$od/topGO.map".'")'."\n";
print OUT 'data<-read.table("'."$od/topGO.list".'", row.names = 1, header=TRUE)'."\n";
print OUT 'geneList<-data[,1]'."\n";
print OUT 'names(geneList) <- rownames(data)'."\n";
print OUT 'topDiffGenes<-function(allScore){return(allScore<0.05)}'."\n";
print OUT 'sampleGOdata <- new("topGOdata",nodeSize = 1,ontology="BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,geneSel=topDiffGenes)'."\n";
print OUT 'resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")'."\n";
print OUT 'allRes <- GenTable(sampleGOdata,KS = resultKS.elim,ranksOf = "classic", topNodes = attributes(resultKS.elim)$geneData[4])'."\n";
print OUT 'ug <- usedGO(sampleGOdata)'."\n";
print OUT 'ann.genes<- genesInTerm(sampleGOdata,ug)'."\n";
#print OUT 'write.table(ug, file="'."$od/$key".'.topGO_BP_terms.xls", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)'."\n";
#print OUT 'length(ug)';
#print OUT 'usedGO';
print OUT 'for(i in 1:length(ann.genes)){'."\n";
print OUT 'write.table(matrix(c(names(ann.genes[i]),ann.genes[[i]]),nrow=1),file="'."$od/$key".'.topGO_BP_gene.xls",col.names =F,append=T,sep="\t",row.names=F,quote=F)'."\n";
print OUT '}'."\n";
#print OUT 'test<-lapply(ann.genes,matrix,nrow=1)'."\n";
##print OUT 'lapply(test,write.table,append=T,"test1.txt",col.names=F)'."\n";
#print OUT 'printGenes(sampleGOdata, whichTerms=ug, file="test.xls",)';
print OUT 'write.table(allRes, file="'."$od/$key".'.topGO_BP.xls", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)'."\n";
print OUT 'pdf("'."$od/$key".'.topGO_BP.pdf")'."\n";
print OUT 'showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 10, useInfo = "all")'."\n";
print OUT 'dev.off()'."\n";
print OUT 'png("'."$od/$key".'.topGO_BP.png")'."\n";
print OUT 'showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 10, useInfo = "all")'."\n";
print OUT 'dev.off()'."\n";
print OUT 'sampleGOdata <- new("topGOdata",nodeSize = 1,ontology="MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,geneSel=topDiffGenes)'."\n";
print OUT 'resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")'."\n";
print OUT 'allRes <- GenTable(sampleGOdata,KS = resultKS.elim,ranksOf = "classic", topNodes = attributes(resultKS.elim)$geneData[4])'."\n";
print OUT 'ug <- usedGO(sampleGOdata)'."\n";
print OUT 'ann.genes<- genesInTerm(sampleGOdata,ug)'."\n";
print OUT 'for(i in 1:length(ann.genes)){'."\n";
print OUT 'write.table(matrix(c(names(ann.genes[i]),ann.genes[[i]]),nrow=1),file="'."$od/$key".'.topGO_MF_gene.xls",col.names =F,append=T,sep="\t",row.names=F,quote=F)'."\n";
print OUT '}'."\n";
#print OUT 'test<-lapply(ann.genes,matrix,nrow=1)'."\n";
#print OUT 'lapply(test,write.table,append=T,"test2.txt",col.names=F)'."\n";
print OUT 'write.table(allRes, file="'."$od/$key".'.topGO_MF.xls", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)'."\n";
print OUT 'pdf("'."$od/$key".'.topGO_MF.pdf")'."\n";
print OUT 'showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 10, useInfo = "all")'."\n";
print OUT 'dev.off()'."\n";
print OUT 'png("'."$od/$key".'.topGO_MF.png")'."\n";
print OUT 'showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 10, useInfo = "all")'."\n";
print OUT 'dev.off()'."\n";
print OUT 'sampleGOdata <- new("topGOdata",nodeSize = 1,ontology="CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,geneSel=topDiffGenes)'."\n";
print OUT 'resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")'."\n";
print OUT 'allRes <- GenTable(sampleGOdata,KS = resultKS.elim,ranksOf = "classic", topNodes = attributes(resultKS.elim)$geneData[4])'."\n";
print OUT 'ug <- usedGO(sampleGOdata)'."\n";
print OUT 'ann.genes<- genesInTerm(sampleGOdata,ug)'."\n";
print OUT 'for(i in 1:length(ann.genes)){'."\n";
print OUT 'write.table(matrix(c(names(ann.genes[i]),ann.genes[[i]]),nrow=1),file="'."$od/$key".'.topGO_CC_gene.xls",col.names =F,append=T,sep="\t",row.names=F,quote=F)'."\n";
print OUT '}'."\n";
#print OUT 'test<-lapply(ann.genes,matrix,nrow=1)'."\n";
#print OUT 'lapply(test,write.table,append=T,"test3.txt",col.names=F)'."\n";
print OUT 'write.table(allRes, file="'."$od/$key".'.topGO_CC.xls", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)'."\n";
print OUT 'pdf("'."$od/$key".'.topGO_CC.pdf")'."\n";
print OUT 'showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 10, useInfo = "all")'."\n";
print OUT 'dev.off()'."\n";
print OUT 'png("'."$od/$key".'.topGO_CC.png")'."\n";
print OUT 'showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 10, useInfo = "all")'."\n";
print OUT 'dev.off()'."\n";
close OUT;

#`ssh compute-0-14 /share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $od/topGO.R`;
my $Rscript= ${&readconf("$Bin/./../../../../../config/sys.cfg")}{Rscript};
`$Rscript $od/topGO.R`;

my $svg2xxx= ${&readconf("$Bin/./../../../../../config/sys.cfg")}{svg2xxx};
chdir $od;
system("$svg2xxx $key.topGO_BP.svg");
system("$svg2xxx $key.topGO_MF.svg");
system("$svg2xxx $key.topGO_CC.svg");
#`rm $od/topGO.R`;

#######################################################################################
&timeLog("$Script Done. Total elapsed time : ".time()-$BEGIN_TIME."s");
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub max{#&max(lists or arry);
	#���б��е����ֵ
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#���б��е���Сֵ
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#��ȡ�ַ������еķ��򻥲����У����ַ�����ʽ���ء�ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Program Date:   2013.8.21
Usage:
  Options:
  -All_GO    <file>  input file,forced 
  
  -DEG_list  <file>  input file,forced 
  
  -key       <str>   index of output files,forced 
  
  -od        <dir>   output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
