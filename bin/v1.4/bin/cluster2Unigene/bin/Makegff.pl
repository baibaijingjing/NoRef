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
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);


my $line="";
my $id;
my $length;
open(IN,$fIn);
open(OUT,">$fOut");
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		if (length($line)==0) {
		}else{
			$length=length($line);
			print  OUT $id,"\t","blat","\t","mRNA","\t","1","\t",$length,"\t","1.0000","\t","+","\t",".","\t","ID=$id","\n";
			print OUT $id,"\t","blat","\t","CDS","\t","1","\t",$length,"\t","100","\t","+","\t",".","\t","Parent=$id","\n";
		}
		$id=$1;
		$line=0;
	}else {
		$line.=$_;
	}
}
close IN;
$length=length($line);
print  OUT $id,"\t","blat","\t","mRNA","\t","1","\t",$length,"\t","1.0000","\t","+","\t",".","\t","ID=$id","\n";
print OUT $id,"\t","blat","\t","CDS","\t","1","\t",$length,"\t","100","\t","+","\t",".","\t","Parent=$id","\n";
$line=0;
close OUT;




#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact: wu shuang <wus\@biomarker.com.cn> 

Usage:
  Options:
  -i <file>   Input EST file, forced
  -o <file>   Output EST.gff file, forced
  -h         Help

USAGE
	print $usage;
	exit;
}
