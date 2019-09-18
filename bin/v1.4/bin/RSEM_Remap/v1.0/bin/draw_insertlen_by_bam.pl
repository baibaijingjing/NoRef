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
my ($fIn,$index,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"index:s"=>\$index,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn and $index and $od);

mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);

my %hash;
open (IN,"samtools view $fIn|") or die $!;
while (<IN>) {
	my ($name,$match,$type,$insert)=(split/\t/,$_)[0,5,6,8];
	next unless $type eq '=';
#	next unless ($match eq '101M' or $match eq '100M');
	next if $insert<=0;
#	$hash{$name}=$insert;
	if (exists $hash{$name}) {
		$hash{$name}=($hash{$name}>$insert)?$hash{$name}:$insert;
	}
	else{
		$hash{$name}=$insert;
	}
}
close (IN) ;

my %Num;
foreach my $key (keys %hash) {
	$Num{$hash{$key}}++;
}

open (OUT,">","$od/$index.insertSize") or die $!;
foreach my $insertSize (sort {$a <=> $b} keys %Num) {
    print OUT $insertSize,"\t",$Num{$insertSize},"\n";
}
close (OUT) ;


##########modify by xugl 2015-8-11#########

open OUT,">$od/$index.insertSize.list";
foreach my $key (sort {$a <=> $b} keys %Num) {
    next if $key<=145;
    print OUT $key,"\t",$Num{$key},"\n";
    last if $key>=800;
}
close (OUT);
my $Rscript = ${&readconf("$Bin/../../../../../../../config/sys.cfg")}{Rscript};
system "$Rscript $Bin/../../../share_script/plot_insertsize.R infile=$od/$index.insertSize.list outfile=$od/$index.insertSize.png bg=F";

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
Usage:
  Options:
  -i        <file>         input file,bam format,forced 
  
  -index    <str>          index of outfile,forced 
  
  -od       <dir>          output dir,forced 
  
  -h                       Help

USAGE
	print $usage;
	exit;
}
