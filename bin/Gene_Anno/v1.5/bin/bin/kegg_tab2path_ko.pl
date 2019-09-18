#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use autodie;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0";
my @Times = localtime();
#######################################################################################
&timeLog("$Script start.");
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($ftab,$foutdir,$fkey,$gene2pathway,$gene2KO,$anno);
GetOptions(
				"help|h|?" =>\&USAGE,
				"tab:s"=>\$ftab,
				"od:s"=>\$foutdir,
				"key:s"=>\$fkey,
				"g2p:s"=>\$gene2pathway,
				"g2k:s"=>\$gene2KO,
				"anno:s"=>\$anno,
				) or &USAGE;
&USAGE unless ($ftab and $foutdir and $fkey);
$ftab = Cwd::realpath($ftab);
$foutdir = Cwd::realpath($foutdir);
mkdir $foutdir if (!-d $foutdir);

my %globemap = (
		"ko01100" => 1,
		"ko01110" => 1,
		"ko01120" => 1,
		);

########database file############
my $keggsmap;
my %sys_cfg=%{&selectconf("$Bin/../../../../../config")};

$keggsmap||=$sys_cfg{kegg_png_file};
$gene2KO||=$sys_cfg{gene2KO_file};
$gene2pathway||=$sys_cfg{gene2pathway_file};
$anno||=$sys_cfg{kegg_anno_file};

my %result_path;
my %globe_path;
my %for_pie;
my %gene2path;
my %path_title;
my %genes_ko;
open I,"$gene2pathway";
while(<I>){
	chomp;
	my @tmp=split/\t/;
	$gene2path{$tmp[0]}{$tmp[1]}=1;
	$path_title{$tmp[1]}=$tmp[2];
}
close I;
open GK,"$gene2KO";
my %gene2K;
while(<GK>){
	chomp;
	my @tmp=split/\s+/;
	$gene2K{$tmp[0]}=$tmp[1];
}
close GK;
my %annotation;
open I,"$anno";
while(<I>){
	chomp;
	next if /^\s*$/ or /^\#/;
	my ($id,$des)=split/\s+/,$_,2;
	$id=~s/^>//;
	$annotation{$id}=$des;
	$annotation{$id}="" unless $des ;
}
close I;
open TAB,"$ftab" or die "cannot open $ftab \n";
open OUTKO,">$foutdir/$fkey.ko" or die "cannot open $foutdir/$fkey.ko ,$!\n";
print OUTKO "#Gene_id\tKO|e_value|Database_Genes|Anno\n";
#Pumpkin_171-T1_Unigene_BMK.35822        243     30      242     rcu:RCOM_1032360        441     1       69      0.51    0.66    0.01    72    68.9    4e-13   --      nucleic acid binding protein, putative
while (<TAB>) {
	chomp;
	next if (/^\s*$/);
	my @temp = split /\t/, $_;
	#$temp[-1] =~ /(.*?);?.*/;
	if (defined $gene2K{$temp[-1]}) {
        print OUTKO "$temp[0]\t$gene2K{$temp[-1]}|$temp[-3]|$temp[-1]|$annotation{$temp[-1]}\n";
	}	
	if (defined $gene2path{$temp[-1]}) {
		foreach my $pa (sort keys %{$gene2path{$temp[-1]}}) {
			next unless (defined $gene2K{$temp[-1]});
			if (defined $globemap{$pa}) {
				$globe_path{$pa}{$temp[0]} = $gene2K{$temp[-1]};
			}else{
				$result_path{$pa}{$temp[0]} = $gene2K{$temp[-1]};
			}
		}
	}
}
close TAB;
close OUTKO;

open OUTPATH,">$foutdir/$fkey.pathway" or die "cannot open $foutdir/$fkey.pathway ,$!\n";
print OUTPATH "#pathway\tpathway_id\tGene_number\tGene_id\tKOs\n";
foreach my $key (sort keys %result_path) {
	print OUTPATH "$path_title{$key}\t$key\t";
	my ($num, $realgenes,$kos);
	$num = 0;
	foreach my $ge (sort keys %{$result_path{$key}}) {
		$num++;
		$realgenes .= "$ge;";
		$kos .= "$result_path{$key}{$ge}\+";
	}
	chop $kos;
	$for_pie{$key} = $num;
	print OUTPATH "$num\t$realgenes\t$kos\n";
}
close OUTPATH;

open OUTG,">$foutdir/$fkey.globe" or die "cannot open $foutdir/$fkey.globe ,$!\n";
print OUTG "#pathway\tpathway_id\tGene_number\tGene_id\tKOs\n";
foreach my $key (sort keys %globe_path) {
	print OUTG "$path_title{$key}\t$key\t";
	my ($num, $realgenes,$kos);
	$num = 0;
	foreach my $ge (sort keys %{$globe_path{$key}}) {
		$num++;
		$realgenes .= "$ge;";
		$kos .= "$globe_path{$key}{$ge}\+";
	}
	chop $kos;
	print OUTG "$num\t$realgenes\t$kos\n";
}
close OUTG;

#open PIE,">$foutdir/$fkey.pie.stat" or die;
#print PIE "#ko\tnum\n";
#my $pathway_num = scalar (keys %for_pie);
#my $n = ($pathway_num >=20) ? 20 : 10;
#
#my $k = 1;
#my $other = 0;


################pick up kegg map##############    add this part     :modified by sunhy 2013/4/15
my @map_picture=glob "$keggsmap/*";
my %pathway_map;

foreach my $map_dir (@map_picture) {
	if ($map_dir=~m/.*(ko\d+)\.png/) {
		$pathway_map{$1}=$map_dir;
	}
}
mkdir("$foutdir/Kegg_map") if (!-d "$foutdir/Kegg_map");
foreach my $path (keys %result_path) {
	if (defined $pathway_map{$path}) {
		system("cp $pathway_map{$path} $foutdir/Kegg_map/");
	}
}

#################draw pie graph##############
#foreach my $key (sort {$for_pie{$b}<=>$for_pie{$a}} keys %for_pie) {
#	if ($k < $n ) {
#		print PIE "$key\t$for_pie{$key}\n";
#		$k++;
#	}else{
#		$other += $for_pie{$key};
#	}
#}
#print PIE "Other\t$other\n";
#close PIE;
#if ($n == 20) {
#	system("perl $Bin/Just_Pie.pl -i $foutdir/$fkey.pie.stat -o $foutdir/$fkey.pie.svg -w 600 -h 450 -css $Bin/pie20.css -note 'KEGG pathway distribution'");
#}else{
#	system("perl $Bin/Just_Pie.pl -i $foutdir/$fkey.pie.stat -o $foutdir/$fkey.pie.svg -css $Bin/pie12.css -note 'KEGG pathway distribution'");
#}
#######################################################################################
#&timeLog("$Script Done.\t\tTotal elapsed time : ".time()-$BEGIN_TIME."s");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE {
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[2013/3/26]
	Contact:Sun Huaiyu <sunhy\@biomarker.com.cn>
	Description:	Extract KEGG Pathway and KOs file From Blast tab file;

	Usage:
		-tab         The blast_tab Out of seq with KEGG           must be given

		-od          Output dir                                   must be given

		-key         Prefix of OUT files (eg. key.path key.ko)    must be given
		
		-g2p         the file contains the gene with pathway
		-g2k         the file  contains the gene with KEGG Orthology
		-anno        The annotaion file for database.fa           must be given
		-h          Help document
USAGE
	print $usage;
	exit;
}
