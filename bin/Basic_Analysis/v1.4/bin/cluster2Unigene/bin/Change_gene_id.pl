#!/usr/bin/perl -w
# 
#Copyright (c) 
#Writer         baij <baij@biomarker.com.cn>
#Program Date   2016101018
my $ver="1.0.0";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;


my %opts;
GetOptions(\%opts,"i=s","od=s","h");
if (!defined($opts{i}) || !defined($opts{od}) || defined($opts{h})) {
	&help;
	exit;
}

my $in=&ABSOLUTE_DIR($opts{i});
mkdir $opts{od} unless -d $opts{od};
my $odir=&ABSOLUTE_DIR($opts{od});
my %seq;



open (IN,"$in") || die;
$/='>';
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/) ;
	my ($head,$seq_1)=split /\n+/,$_,2;
	my $id=(split/\s+/,$head)[0];
	$seq_1=~s/\s+//g;
	$seq{$id}=$seq_1;
}
close (IN);


open OUT1,">$odir/trans_TGICL_results_ChangeId.fa" || die ;
open OUT2,">$odir/TGICL_BMK_ID.list" || die ;

my $i=1;
my $num = "%0".length(scalar keys %seq)."d";
foreach my $comp (sort keys %seq) {
        my $idnum =sprintf $num, $i;
        my $gene_id ="BMK_Unigene_".$idnum;	
        print OUT1 ">$gene_id\n$seq{$comp}\n";
	print OUT2 ">$gene_id\t$comp\n";
        $i++;
}

close OUT1;
close OUT2;
############# subs

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
		die "Warning just for file and dir $in\n";
	}
	
	chdir $cur_dir;
	return $return;
}

sub help
{
	print <<"	Usage End.";
	Description:
		Function : change ID from trans_TGICL_results.fa to BMK ID;
		Version  : $ver
		Writer   : baij <baij\@biomarker.com.cn>
		Usage    :
		-i
		    fasta file of trans_TGICL_results;
		-od
		    fasta file with BMK ID dir;
		-h
		    Help document
		
	Usage End.

	exit;
}

