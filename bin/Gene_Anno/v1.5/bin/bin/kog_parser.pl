#!/usr/bin/perl

=head1 Name

  kog_parser.pl  --  extract the kog number.

=head1 Description

  This program is designed for extracting the kog number after blast.

=head1 Version

  Author: Simon Young, simonyoung8824@gmail.com
  Version: 1.0,  Date: 2014-09-28

=head1 Usage
	
  perl kog_parser.pl <blast_tab>
  --nohead      do not show the first instruction line.
  --verbose     output verbose information to screen  
  --help        output help information to screen  

=head1 Exmple

  perl kog_parser.pl PLASMID.fasta.ori.glimmer3.pep.blast.tab

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
my ($Nohead,$Verbose,$Help);
GetOptions(
	"nohead"=>\$Nohead,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my ($kog,$fun);
open CFG,"$Bin/../../../../../config/sys.cfg" or die "Error:Can not open $Bin/../../../../../config/sys.cfg";
while(<CFG>){
	chomp;
	if(/^kog_fun_file/){
		$fun=(split/\s+/)[1];
	}
	elsif(/^kog_file/){
		$kog=(split/\s+/)[1];
	}
}
close CFG;
my $blast_tab = shift;

my $seq_name = $1 if ($blast_tab=~/([^\/]+)\.Kog\.blast/);
my $outdir = dirname($blast_tab);
#print "$seq_name\n";

my (%KOG,%fun);
my $category;

##read fun.txt
open FUN,$fun || die "fail $fun";
while (<FUN>) {
	if (/(\[\w\])\s(.+?)\n$/) {
		my ($class,$function) = ($1,$2);
		$fun{$class}{defination} = $function;
		$fun{$class}{category} = $category;
	}else {
		chomp $_;
		$category = $_;
	}
}
close FUN;

##read kog
open KOG,$kog || die "fail $kog";
$/="\n[";
while (<KOG>) {
    chomp;
    s/^\[// if ($.==1);
	my ($class,$kog_num,$kog_anno) = ($1,$2,$3) if (/(\S+\])\s+(KOG\d+)\s+(.+?)\n/s);
    $class = '['.$class;
	$kog_anno =~ s/\s+/ /g;
	$_ =~ s/.+?\[.+?\n//s;
	my @line = split (/\n/,$_);
	shift @line;
	my $org;
	for (my $i=0;$line[$i];$i++) {
		$org = $1 if ($line[$i]=~/(\S+)\:/);
		$line[$i] =~ s/.+?\://;
		$line[$i] =~ s/^\s+//;
		my @protein = split (/ /,$line[$i]);
		for (my $j=0;$protein[$j];$j++) {
			$fun{$class}{category} = "--" unless (exists $fun{$class}{category});
			$fun{$class}{defination} = "--" unless (exists $fun{$class}{defination});
			$KOG{$protein[$j]} = "$org\t$kog_num\t$kog_anno\t$class\t$fun{$class}{category}\t$fun{$class}{defination}";
		}
	}
}
$/="\n";
close KOG;

##read the tab file and create a file including kog info
open IN,$blast_tab || die "fail $blast_tab";
open OUT,">$outdir/$seq_name.Kog.class" || die "fail $seq_name.Kog.class";
print OUT "#Gene name\tPortein_name_in_KOG\tE_value\tIdentity\tScore\tOrganism\tKOG_id\tKOG_class_defination\tFunction_code\tFunctional_categories\tFunction_class_defination\n" unless (defined $Nohead);
while (<IN>) {
	next if /^\#/;
	my @t = split /\t/;
	print OUT "$t[0]\t$t[4]\t$t[13]\t$t[8]\t$t[12]\t$KOG{$t[4]}\n" if (exists $KOG{$t[4]});
}
close OUT;
close IN;

