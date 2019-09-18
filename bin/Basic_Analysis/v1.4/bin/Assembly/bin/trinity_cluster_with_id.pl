#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $ver="1.0.0";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;

##############################����д����֮ǰ��һ��д��ʱ�䡢������;������˵����ÿ���޸ĳ���ʱ��Ҳ������ע�͹�����
my %opts;
GetOptions(\%opts,"i=s","od=s","key=s","h");
if (!defined($opts{i}) || !defined($opts{od}) || !defined($opts{key}) || defined($opts{h})) {
	&help;
	exit;
}

######################
my $in=&ABSOLUTE_DIR($opts{i});
&MKDIR($opts{od});
my $odir=&ABSOLUTE_DIR($opts{od});
my $key=$opts{key};

my %seq;
my %ID;

open (IN,"$in") || die;
$/='>';
while (<IN>) {
	chomp;
	next if (/^$/ || /^\#/) ;
	my ($head,$seq_1)=split /\n+/,$_,2;
	my $id=(split/\s+/,$head)[0];
	$seq_1=~s/\s+//g;
	$seq{$id}=$seq_1;
	$id=~/(.*_i)\d+/;
	my $path=$1;
	push @{$ID{$path}},$id;
}

close (IN);

my %AS_Stat;
my $AS_sum;
open OUT1,">$odir/$key.Unigenes.fa" || die ;
open OUT2,">$odir/$key.Unigenes.Cluster.xls" || die ;
open OUT3,">$odir/$key.Unigene_Trans_ID.list" || die ;

my $i=1;
foreach my $comp (sort keys %ID) {
	my $gene_id="$key"."_Unigene_BMK.$i";
	print OUT2 ">$gene_id\n";
	my $gene_seq;
	my $max_length=0;
	my $max_id;
	my $AS=0;
	foreach my $Trans_id (sort @{$ID{$comp}}) {
		my $length=length $seq{$Trans_id};
		print OUT2 "$Trans_id\t$length\t$seq{$Trans_id}\n";
		if ($length>$max_length) {
			$gene_seq=$seq{$Trans_id};
			$max_length=$length;
			$max_id=$Trans_id;
		}
		$AS++;
	}
	my $AS_num=$AS-1;
	print OUT1 ">$gene_id\n$gene_seq\n";
	print OUT3 "$gene_id\t$max_id\n";
	if ($AS_num<=5) {
		$AS_Stat{$AS_num}++;
		$AS_sum++;
	}
	else {
		$AS_Stat{"5"}++;
		$AS_sum++;
	}
	$i++;
}

close OUT1;
close OUT2;
close OUT3;

open OUT4,">$odir/$key.Unigene.Cluster.stat.xls" || die ;

print OUT4 "AS_number\tGene_number\tPercent\n";

foreach my $AS_num (sort{$a<=>$b} keys %AS_Stat) {
	my $percent=$AS_Stat{$AS_num}/$AS_sum;
	print OUT4 "$AS_num\t$AS_Stat{$AS_num}\t$percent\n";
}
print OUT4 "All\t$AS_sum\t1\n";
close OUT4;


############# subs
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
		die "Warning just for file and dir $in\n";
	}
	
	chdir $cur_dir;
	return $return;
}

sub help
{
	print <<"	Usage End.";
	Description:
		Function : use Trinity OUT fasta ID to Cluster Alter Splice Transcriptome to Uingene;
		Version  : $ver
		Writer   : mengf <mengf\@biomarker.com.cn>
		Usage    :
		-i
		    Denovo Transcripts file;
		-od
		    Unigene Cluster with Trinity ID dir;
		-key
		    prefix of Unigene ID ��work shell (The Sample Name is Recommend eg. Tea_XXY);
		-h
		    Help document
		
	Usage End.

	exit;
}
