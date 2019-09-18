#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fK);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"od:s"=>\$fK,
				) or &USAGE;
&USAGE unless ($fIn and $fK);

mkdir $fK unless (-d $fK) ;
$fK=&ABSOLUTE_DIR($fK);
$fK=~/\/([^\/]+)$/;
my $Key=$1;

chdir $fK;

print "Get insertSize :$fIn\n";
my %hash;
open (IS,"<","$fIn") or die $!;
while (<IS>) {
	my ($Mode,$num,$D1,$D2,$D3,$D4,$D5,$D6)=(split "\t",$_)[1,3,4..9];
	if ($Mode eq "PE") {
			$hash{$D4-$D1}++;
			print "Warning:This Record's InsertSize is less than 0, The progarm will ignore it!\n$_" if ($D4-$D1<0) ;
	}elsif($Mode eq "SO" or $Mode eq "ST"){
		$hash{$D4-$D1-($D6-$D5)}++;
		print "Warning:This Record's InsertSize is less than 0, The progarm will ignore it!\n$_" if ($D4-$D1-($D6-$D5)<0) ;
	}
}
close (IS) ;

open (OUT,">","$Key.insertSize") or die $!;
foreach my $insertSize (sort {$a <=> $b} keys %hash) {
	print OUT $insertSize,"\t",$hash{$insertSize},"\n";
}
close (OUT) ;

##########modify by xugl 2015-8-11#########
open OUT,">$fK/$Key.insertSize.list";
foreach my $key (sort {$a <=> $b} keys %hash) {
    next if $key<=145;
    print OUT $key,"\t",$hash{$key},"\n";
    last if $key>=800;
}
close (OUT);
my $Rscript = "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
system "$Rscript $Bin/../../share_script/plot_insertsize.R infile=$fK/$Key.insertSize.list outfile=$fK/$Key.insertSize.png bg=F";


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Program Date:   2013.11.21
Usage:
  Options:
  -i    <file>  input file,rmap format,forced 
  
  -od   <dir>   output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
