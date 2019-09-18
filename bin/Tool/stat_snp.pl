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
my ($id);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
				) or &USAGE;
&USAGE unless ($id);

if (defined $id) {
	my %SNP;
	my @SNP_dir=glob "$id/*";
	print "\n#Type\tS1.homo.S2.homo\tS1.hete.S2.homo\tS1.homo.S2.hete\tS1.hete.S2.hete\tTotal\n";
	foreach my $dir (@SNP_dir) {
		if (-d $dir) {
			$dir=~m/.*\/(\S+)/;
			my $nam=$1;
			my $file;
			foreach my $f (glob "$dir/*.homo.snp.xls") {
				unless ($f=~/\.hete_/) {
					$file=$f;
					last;
				}
			}
			my $mm=`wc -l $file`;
			$SNP{$nam}{mm}=(split/\s+/,$mm)[0]-1;
			my $tm=`wc -l $dir/*.hete_*.homo.snp.xls`;
			$SNP{$nam}{tm}=(split/\s+/,$tm)[0]-1;
			my $mt=`wc -l $dir/*.homo_*.hete.snp.xls`;
			$SNP{$nam}{mt}=(split/\s+/,$mt)[0]-1;
			foreach my $f (glob "$dir/*.hete.snp.xls") {
				unless ($f=~/\.homo_/) {
					$file=$f;
					last;
				}
			}
			my $tt=`wc -l $file`;
			$SNP{$nam}{tt}=(split/\s+/,$tt)[0]-1;
			$SNP{$nam}{total}=$SNP{$nam}{mm}+$SNP{$nam}{tm}+$SNP{$nam}{mt}+$SNP{$nam}{tt};
		}
	}
	foreach my $name (keys %SNP) {
		print "$name\t$SNP{$name}{mm}\t$SNP{$name}{tm}\t$SNP{$name}{mt}\t$SNP{$name}{tt}\t$SNP{$name}{total}\n";
	}
	print "\n";
}

#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
  -id <dir>  input dir,forced 
  
 
  -h         Help

USAGE
	print $usage;
	exit;
}
