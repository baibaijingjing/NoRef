#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################

# ==============================================================
# Get Options
# ==============================================================
my ($fIn1,$fIn2,$cfg,$od);

GetOptions(
				"help|?" =>\&USAGE,
				"snp:s"=>\$fIn1,
#				"chi:s"=>\$fIn2,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($fIn1);
#===============================================================
# Default optional value 
#===============================================================
$od||="./stat";
&MKDIR($od);
&begining;

my %type = ('AG'=>'R', 'GA'=>'R', 'CT'=>'Y', 'TC'=>'Y', 'GT'=>'K', 'TG'=>'K', 'AC'=>'M', 'CA'=>'M', 
	'CG'=>'S', 'GC'=>'S', 'AT'=>'W', 'TA'=>'W', 'A'=>'A', 'C'=>'C', 'G'=>'G', 'T'=>'T', 
	'CGT'=>'B', 'AGT'=>'D', 'ACT'=>'H', 'ACG'=>'V', 'ACGT'=>'N',
);
my %Retype;
foreach my $key (keys(%type)) {
	$Retype{$type{$key}}=$key;
}
#===============================================================
# Process
#===============================================================
my @header;
my %base;
my %char;
my %count;
open IN,"$fIn1" or die $!;
while (<IN>) {
	chomp;
	my @line=split /\t/,$_;
	if ($_=~/\#/){
		@header=@line;
		next;
	}
	my $str="$line[0]:$line[1]";
	$base{$str}="$line[0]\t$line[1]\t$line[2]\t$line[3]";
	for (my $i=4;$i<@line ;$i=$i+3) {
		my $x=$line[$i];
		next if($header[$i]=~/effect/i);
		next if($x eq "N");
		$count{$header[$i]}{"homo"}=0 if(!exists $count{$header[$i]}{"homo"}); 
		$count{$header[$i]}{"hete"}=0 if(!exists $count{$header[$i]}{"hete"});
		$count{$header[$i]}{"homo"}++ if($x eq 'A'  or $x eq 'T'  or $x eq 'G'  or $x eq 'C');
		$count{$header[$i]}{"hete"}++ if($x ne 'A' and $x ne 'T' and $x ne 'G' and $x ne 'C');
		$char{$header[$i]}{$str}=$line[$i]."\t".$line[$i+1]."\t".$line[$i+2];
	}
}
close IN;

open STAT,">$od/AllSample.snp.stat";
print STAT "#Samples\tHomoSNP\tHeteSNP\tAllSNP\n";
foreach my $sam (sort keys %count) {
	my $sum=$count{$sam}{homo}+$count{$sam}{hete};
	print STAT "$sam\t$count{$sam}{homo}\t$count{$sam}{hete}\t$sum\n";
	open OUT,">$od/$sam.snp.list" or die $!;
	print OUT "#GeneID\tPos\tRef\tAlt\t$sam\tDepth\tAlleDepth\n";
	foreach my $str (sort keys %{$char{$sam}}) {
		my ($chr,$pos)=split /:/,$str;
		print OUT "$base{$str}\t$char{$sam}{$str}\n";
	}
	close OUT;
}
close STAT;

&end;
# ==============================================================
# sub function
# ==============================================================
sub MKDIR{
	my ($dirname)=@_;
	if (-e $dirname) {
		print "$dirname exists!";
	}else{
		mkdir $dirname;
	}
}
sub begining{
	my $cmd=$0." ";
	$cmd.=join(" ",@Original_ARGV);
	my $time=&GetTime;
	print STDOUT "$cmd\n START:$time\n";
}
sub end{
	print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
}
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub ReadConfig{
	my($cfg)=@_;
	open (IN,"$cfg") || die "$!";
	my %para;
	my %sample;
	while (<IN>) {
		chomp;
		s/\s+//;
		s/\s+$//;
		next if (/\#/);
		next if (/$/);
		my @tmp=split /\s+/,$_;
		if ($tmp[0]=~m/Sample/) {
			my $fq1=<IN>;chomp $fq1;
			my $fq2=<IN>;chomp $fq2;
			my @fq_1=split /\s+/,$fq1;
			$sample{$tmp[1]}{FQ1}=$fq_1[1];
			my @fq_2=split /\s+/,$fq2;
			$sample{$tmp[1]}{FQ2}=$fq_2[1];
		}
		$para{$tmp[0]}=$tmp[1];
	}
	close IN;

	return (\%para,\%sample);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Shi Tongwei <shitw\@biomarker.com.cn> 
=====================================================================================================
Discription:
=====================================================================================================

Usage:
  Options:
  -snp	<file>	required	parent snplist file
  -od	<str>	optional	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}

