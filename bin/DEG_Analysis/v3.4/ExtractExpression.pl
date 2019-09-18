#!/usr/local/bin/perl -w
# Copyright (c) BMK 2013
# Writer:         lium <lium@biomarker.com.cn>
# Program Date:   2013.
# Modifier:       baij <baij@biomarker.com.cn>
# Last Modified:  2017-01-12
my $ver="1.0.0";

use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $programe_dir=basename($0);
my $path=dirname($0);


# get option
our %opts;
GetOptions(\%opts,"indir=s","out=s","h" );
# check option
if(!defined($opts{indir}) || !defined($opts{out}) || defined($opts{h})){
	&help();
	exit;
}
# get indir
my $indir=&ABSOLUTE_DIR($opts{indir});
my $out=&ABSOLUTE_DIR($opts{out});

###############Time
&timeLog("$programe_dir Start.");


##########################################################
##
##           extract read count and fpkm from dir
##
##########################################################

############### Load geneExpression.xls files
my @geneExpress = glob "$indir/*.geneExpression.xls" ;
my %count;
my %GeneLength;
my %fpkm;
my @sample;

foreach my $file (@geneExpress) {
	my $name=basename($file);
	$name=~/(.*)\.geneExpression\.xls/;
	my $key=$1;
	push @sample,$key;
	open (IN,"$file") or die $!;
	while (<IN>) {
		s/^\s+//;s/\s+$//;s/\r+$//;
		next if (/^$/ || /^\#/ || /^Gene/);
		my @tmp=split/\s+/,$_;
		$count{$tmp[0]}{$key}=int($tmp[6]);
                $fpkm{$tmp[0]}{$key}=$tmp[4];
		$GeneLength{$tmp[0]} = $tmp[1];
	}
	close (IN) ;
}


############## combinate read count and output them into one file ############
open (COUNT, ">$out/All_gene_counts.list") or die $!;
print COUNT "#ID\t";
print COUNT join("\t",(sort @sample)),"\tgeneLength\n";
foreach my $id (sort keys %count) {
	print COUNT "$id";
	my $str;
	foreach my $key (sort @sample) {
		if (!defined $count{$id}{$key}) {
			$str.="\t"."0";
			next;
		}
		$str.="\t"."$count{$id}{$key}";
	}
	print COUNT "$str\t$GeneLength{$id}\n";
}
close (COUNT) ;

##########################################################
##
##           count to expression
##
##########################################################
open (FPKM, ">$out/All_gene_fpkm.list") or die $!;
print FPKM "#ID\t";
print FPKM join("\t",(sort @sample)),"\n";
foreach my $id (sort keys %fpkm) {
        print FPKM "$id";
        my $str;
        foreach my $key (sort @sample) {
                if (!defined $fpkm{$id}{$key}) {
                        $str.="\t"."0";
                        next;
                }
                $str.="\t"."$fpkm{$id}{$key}";
        }
        print FPKM "$str\n";
}
close (FPKM) ;


##########################################################
##
##           produce all expression list
##
##########################################################
open EXP,">$out/All_gene_expression.list" or die $!;
my @sam;
foreach my $s (@sample) {
	push @sam,$s."_count";
	push @sam,$s."_FPKM";
}
my $exp_line=join "\t",@sam;
print EXP "\#ID\t$exp_line\n";

foreach my $id (sort keys %count) {
        print EXP "$id";
        my $str;
        foreach my $key (sort @sample) {
                if (!defined $fpkm{$id}{$key}) {
                        $fpkm{$id}{$key}=0;
                        next;
                }
                $str.="\t"."$count{$id}{$key}"."\t"."$fpkm{$id}{$key}";
        }
        print EXP "$str\n";
}


###############Time
&timeLog("$programe_dir end.");
&Runtime($BEGIN);



sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";

	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}

	chdir $cur_dir;
	return $return;
}

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


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub help {
print <<"Usage End.";
Description: Extract  read count and fpkm  from the directory;
Version: $ver
Writer by lium <lium\@biomarker.com.cn>
Usage:
-indir		geneExpress files dir [ forced ]
-out		the dir name of output [ forced ]
-h		for help
Usage End.
exit;
}
