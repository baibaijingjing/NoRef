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
my ($fIn,$deg,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"deg:s"=>\$deg,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($fIn and $deg and $o);

$fIn=&ABSOLUTE_DIR($fIn);
my $head_line;
my $limit=0;
my %DEG;
open (IN,"$deg") or die $!;
while (<IN>) {
	$_=~s/\s+$//;
	$head_line=$_ if /^\#/;
	next if /^\#/;
	my ($name,$val)=split/\t/,$_,2;
	$DEG{$name}=$val;
}
close IN;
open (OUT,">$o") or die $!;

open I,"$fIn/Integrated_Function.annotation.xls";
while(<I>){
	chomp;
	my @tmp=split/\t/;
	if( $.==1){
		
		print OUT "$head_line\t".join("\t",@tmp[1..$#tmp])."\n";
		next;	
	}
	next if (/^\s*$/); 
	if(exists $DEG{$tmp[0]}){
		print OUT "$tmp[0]\t$DEG{$tmp[0]}\t".join("\t",@tmp[1..$#tmp])."\n";	
	}
}

close(OUT);


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
	#求列表中的最大值
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
	#求列表中的最小值
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
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
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
  -i    <file>  input file dirname where All_Database_annotation.xls is,forced 
  
  -deg  <file>  deg file,forced

  -o    <file>  output file,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
