#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Encode;
use Spreadsheet::WriteExcel;
use Spreadsheet::ParseExcel;  
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$deg,$key,$od,$map);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"deg:s"=>\$deg,
				"k:s"=>\$key,
				"od:s"=>\$od,
				"map:s"=>\$map,
				) or &USAGE;
&USAGE unless ($fIn and $deg and $key and $od);

mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);
$fIn=&ABSOLUTE_DIR($fIn);

#&runOrDie("ln -sf $fIn/*Kegg.ko $od/") ;
#my $mid_file_ko=(glob("$od/*Kegg.ko"))[0];

#&runOrDie("ln -sf $fIn/*Kegg.pathway $od/") ;
#my $mid_file_pathway=(glob("$od/*Kegg.pathway"))[0];
#die "Err:file *.Kegg.ko or *.Kegg.pathway does not exist\n" if(!-e $mid_file_ko or !-e $mid_file_pathway);
open CFG,"$Bin/../../../../../config/sys.cfg" or die "$Bin/../../../../../config/sys.cfg";
my $kegg_png;
while(<CFG>){
	chomp;
	if(/^kegg_png_file/){
		$kegg_png=(split/\s+/)[1];
		last;
	}
}
close CFG;
`perl $Bin/KeggGo_enrich_map_web.pl -d $deg -k $key -i $fIn -o $od -func kegg -map $kegg_png`;
#`rm $mid_file_ko`;
#`rm $mid_file_pathway`;

&timeLog("perl $Bin/kegg_enrichment_plot.pl -enrich_file $od/pathway/kegg_enrichment/$key.KEGG.stat -od $od/Graph -key $key");
system "perl $Bin/kegg_enrichment_plot.pl -enrich_file $od/pathway/kegg_enrichment/$key.KEGG.stat -od $od/Graph -key $key";

&timeLog("perl $Bin/draw_KEGG_histogram.pl --ipf $od/pathway/kegg_enrichment/$key.KEGG.xls --opd $od/pathway/kegg_enrichment/ --prf $key.KEGG ");
system "perl $Bin/draw_KEGG_histogram.pl --ipf $od/pathway/kegg_enrichment/$key.KEGG.xls --opd $od/pathway/kegg_enrichment/ --prf $key.KEGG";

#######################################################################################
&timeLog("$Script Done");
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
Program Date:   2013.10.17
Usage:
  Options:
  -i     <file>  input file dirname where All_Database_annotation.xls is ,forced 
  
  -deg   <file>  deg file,forced 
  
  -k     <str>   keywords of output file,forced 
  
  -od    <file>  output dir,forced 
  -map   <path>  the kegg pathway png path
  
  -h         Help

USAGE
	print $usage;
	exit;
}
