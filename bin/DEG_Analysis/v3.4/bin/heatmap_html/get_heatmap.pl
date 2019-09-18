#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);
use Term::ANSIColor qw(:constants);
my $version="1.0";
my ($infile,$id_list,$odir,$index,$column);

my @anno;

GetOptions(
		"i:s"=>\$infile,
		"id:s"=>\$id_list,
		"col:s"=>\$column,
		"p:s"=>\$index,
		"od:s"=>\$odir,
		"help"
	) or &USAGE;

&USAGE if (!defined $infile || !defined $odir || !defined $index|| !defined $id_list) ;

system "mkdir -p  $odir" unless (-d $odir);

my @b;
my @col ;
my $head =`head -n 1 $infile`; 
$head =~s/#ID\t//;
my @col_list =split "\t\+",$head;
if ($column){
$column=~s/(_vs_)|_/,/g;
$column=~s/\s+//g;
open COL,">$odir/col.txt" or die $!;
print COL "$column";
`sed -ri  's/,/\\\n/g ' \"$odir\/col\.txt\"`;
close COL;
open IN,"$odir/col.txt" or die $!;
	while (<IN>){
		chomp;
		for (my $j=0;$j<$#col_list;++$j){
#			print "das\t$j\t$#col_list\n";
			next if ($col_list[$j]=~/_/);
			next if ($col_list[$j]=~/ /);
			 if ("$col_list[$j]" eq "$_"){
			 push (@b,$j+2);
				}
			}
	}
close IN;
}else{
	for (my $j=0;$j<$#col_list;$j++){
		next if ($col_list[$j]=~/_/);
		next if ($col_list[$j]=~/ /);
		push (@b,$j+2);
		}
	}
my $a=join ",",@b;
print "$a\n";
`cut -f  1,$a  $infile >$odir/$index.txt `;

my $cmd_sh="cd $odir && perl  /share/nas2/genome/cloud_soft/app_pipeline/plug-in/pheatmap_modified_by_liux/v1.0/heatmap_draw.pl  -infile  $odir/$index.txt  --outfile $odir/$index.png  -cluster  row   -color_type  1 -is_log -id $id_list  -scale  none -anno 2>$odir/$index.error";
$cmd_sh.=" -scale  row" if (@b>=3);
$cmd_sh.=" -scale  none" if (@b<3);
$cmd_sh.="&& rm $odir/$index.txt ";
print "$cmd_sh\n";
system "$cmd_sh";

sub USAGE {
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Usage:
  Options:
	-i		input file, must be separated as a matrix by space or tab
	-od		output png file, must end with .png
	-id		file contains id to draw, separated by tab, firse column
	-i		scale by column or row: column or row or none
	-col		sample name. eq:T01,T02_vs_T03 
	-p		index name
	-h		Help

USAGE
	print $usage;
	exit;
}



