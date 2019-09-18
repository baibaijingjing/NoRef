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
my ($id,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($id and $od);

mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);

my @Files=glob "$id/*/*.randcheck.list";

my %H;
foreach my $file (@Files) {
	$file=~/([^\/]+)\.randcheck\.list/;
	my $name=$1;
	open (IN,$file) or die $!;
	while (<IN>) {
		next unless /^\d/;
		my ($site,$num)=(split/:/,(split/\s+/,$_)[0])[0,1];
		$H{$name}{$site}=$num;
	}
	close IN;
}

open (OUT,">$od/Total.randcheck.list") or die $!;
print OUT "Sample\tRelative Position in Genes(5'-3')\tPercent of Reads\n";
foreach my $sam (keys %H) {
	foreach my $s (sort {$a<=>$b} keys %{$H{$sam}}) {
		print OUT "$sam\t$s\t$H{$sam}{$s}\n";
	}
}
close OUT;

#`ssh compute-0-14 /share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/bin/random.r $od/Total.randcheck.list $od/Total.randcheck.png`;
#`/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $Bin/bin/random.r $od/Total.randcheck.list $od/Total.randcheck.png`;
my $Rscript=${&readconf("$Bin/../../../../../config/sys.cfg")}{Rscript};
&timeLog("$Rscript $Bin/../share_script/plot_random.R infile=$od/Total.randcheck.list od=$od/");
system "$Rscript $Bin/../share_script/plot_random.R infile=$od/Total.randcheck.list od=$od/";

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
Usage:
  Options:
  -id   <dir>  include */*.randcheck.list,forced 
  
  -od   <dir>  output dir,forced 

  -h         Help

USAGE
	print $usage;
	exit;
}
