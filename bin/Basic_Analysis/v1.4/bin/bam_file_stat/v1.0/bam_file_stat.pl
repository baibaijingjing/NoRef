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
my ($total_num,$fIn,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"t:s"=>\$total_num,
				"i:s"=>\$fIn,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($total_num and $fIn and $o);

	my %B;
	open (IN,"samtools view $fIn|") or die $!;
	while (<IN>) {
		next if (split/\t+/,$_)[2] eq '*';
		my $name=(split/\s+/,$_)[0];
		$name=~s/\/\d$// if $name=~/\/\d$/;
		$B{$name}++;
	}
	close IN;
	my ($Map_num,$Uniq_num,$Multi_num);
	foreach my $key (keys %B) {
		$Map_num++;
		$Uniq_num++ if $B{$key}==2;
		$Multi_num++ if $B{$key}>2;
	}
	my $Map_per=sprintf "%.2f",100*$Map_num/$total_num;
	my $Uniq_per=sprintf "%.2f",100*$Uniq_num/$Map_num;
	my $Multi_per=sprintf "%.2f",100*$Multi_num/$Map_num;
	open (OUT,">$o") or die $!;
	print OUT "Total Reads\t$total_num\t100%\n";
	print OUT "Mapped Reads\t$Map_num\t$Map_per%\n";
	print OUT "Uniq mapped Reads\t$Uniq_num\t$Uniq_per%\n";
	print OUT "Multi mapped Reads\t$Multi_num\t$Multi_per%\n";
	close OUT;

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
Program Date:   2013.10.12
Usage:
  Options:
  -t <num>   Total number of Reads,forced 
  
  -i <file>  input file,bam format,forced 
  
  -o <file>  output file,forced 
 
  -h         Help

USAGE
	print $usage;
	exit;
}
