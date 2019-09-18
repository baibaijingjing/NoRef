#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.1.0";
my @Original_ARGV=@ARGV;

# ==============================================================
# Get Options
# ==============================================================
my ($cfg,$dOut,$qphred);

GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"=>\$cfg,
				"od:s"=>\$dOut,
                                "qphred:s"=>\$qphred,
				) or &USAGE;
&USAGE unless ($cfg);
#===============================================================
# Default optional value 
#===============================================================
###### software dir ###########
my $_mapping_	=	"$Bin/STAR_align/v1.1/STAR_align.pl";
my $_SNP_		=	"$Bin/bin/GATK_calling_Trans_noRef.pl";
$dOut||="./SNP_Trans";
###### global value ###########
$qphred||=33;
my %para;
my %sample;
my $cmd;
$dOut=abs_path($dOut);
#===============================================================
# pipeline
#===============================================================
log_current_time("$Script start...");
&read_config($cfg);
&MKDIR("$dOut");
&MKDIR("$dOut/STAR");
$cmd = "perl $_mapping_ -cfg $cfg -od $dOut/STAR >/dev/null ";
log_current_time("STAR alignment start...");
&runOrDie($cmd);
log_current_time("STAR alignment done.");

&MKDIR("$dOut/aligndir");
my @bam=glob("$dOut/STAR/Alignment/*/*.bam");
foreach my $bam (@bam) {
	next if ($bam=~/work_sh/);
	$bam=~/Alignment\/([^\/]*)\//;
	system "ln -s $bam $dOut/aligndir/$1.bam";
}
$cmd = "perl $_SNP_ -ref $para{'genome'} -aligndir $dOut/aligndir -ploidy $para{'ploidy'} -win $para{'window'} -clu $para{'cluster'} -QD $para{'QD'} -FS $para{'FS'} -od $dOut/SNP -doRecal $para{'Recal'} -doIndel $para{'ReAlignIndel'} -vf $para{vf} -queue $para{queue} ";
$cmd.= " --qphred " if ($qphred==64);
$cmd.= " >/dev/null ";

log_current_time("GATK calling start...");
&runOrDie($cmd);
log_current_time("GATK calling done.");

##SNP sites stat
log_current_time("pairwised SNP abstrct and SNP density plot start...");
system "perl $Bin/util/SNP_stat.pl -snp $dOut/SNP/final.snp.list -od $dOut/SNP/stat/ ";
$cmd= "perl $Bin/util/pairwised_SNP_abstrct_density_plot.pl --ref $para{'genome'} --snp $dOut/SNP/final.snp.list --od $dOut/ >/dev/null ";
&runOrDie($cmd);
log_current_time("pairwised SNP abstrct and SNP density plot done.");

#######################################################################################
my $elapsed_time = (time()-$BEGIN_TIME).'s';
log_current_time("$Script done. elapsed time: $elapsed_time.");
####################################################################################################
sub read_config {#
	my($cfg)=@_;
	open (IN,"$cfg") || die "$!";
	while (<IN>) {
		chomp;
		s/\r$//;
		s/\s+$//;
		next if (/\#/ || /^$/);
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
}

####################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

####################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	mkdir($dir) if(!-d $dir);
}
sub USAGE {#
	my $usage=<<"USAGE";
    Program: $0
    Version: $version
    Contact: Shi Tongwei <shitw\@biomarker.com.cn> 
Discription:
      Usage:
        Options:
        -cfg    <file>  required, config file, (template: $Bin/conf_main.cfg)
        -od     <path>  optional, directory where output file produced, default [./SNP_Trans]
        -help           help

USAGE
	print $usage;
	exit;
}
