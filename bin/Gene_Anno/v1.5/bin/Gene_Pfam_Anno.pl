#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($idir, $pfam_db, $cpu_threads_num, $odir,$queue);

GetOptions(
    "idir:s"=>\$idir,
    "pfam:s"=>\$pfam_db,
    "cpu:i" =>\$cpu_threads_num,
    "odir:s"=>\$odir,
    "help|?" =>\&USAGE,
    "queue:s"=>\$queue,
    ) or &USAGE;
&USAGE unless ($idir);
#######################################################################################
&timeLog("$Script start ...");
#######################################################################################
# ------------------------------------------------------------------
# preparation
# ------------------------------------------------------------------
chomp(my $wd=`pwd`);
$odir            ||= $wd;
$cpu_threads_num ||= 4;

system "mkdir -p $odir" unless (-d $odir);
$odir = &ABSOLUTE_DIR($odir);
$idir = &ABSOLUTE_DIR($idir);
$pfam_db = &ABSOLUTE_DIR($pfam_db);

my $notename = `hostname`; chomp $notename;
my $Q_name = (glob "$idir/*.fa")[0];
$Q_name = basename($Q_name);

if ($Q_name=~/(.+)\.\d+\.fa$/) {
    $Q_name=$1;
} else {
    die "Your file name is Wrong!\n";
}

my $Result_dir="$odir/Result";
my $Tab_dir="$odir/02.gene-annotation";
my $sh_dir="$odir/work_sh/Pfam_sh";
my $cds_predict_dir="$odir/Unigene_CDS_Predict";

&MKDIR("$odir/Pfam_Dir");
&MKDIR("$odir/work_sh");
&MKDIR("$Result_dir");
&MKDIR("$Tab_dir");
&MKDIR("$sh_dir");
&MKDIR("$cds_predict_dir");

my $TransDecoder = "$Bin/TransDecoder/TransDecoder.pl";
# ------------------------------------------------------------------
# creat and run shell file
# ------------------------------------------------------------------
my $align_shell_file = "$sh_dir/Pfam.align.sh";
my @subfiles = glob("$idir/*.fa");

open (OUT,">$align_shell_file") or die "fail $align_shell_file";

foreach my $subfile (@subfiles) {
	my $name = &basename($subfile);
    print OUT "perl $TransDecoder --trans $subfile --pref $name --odir $odir/Pfam_Dir --pfam $pfam_db --cup $cpu_threads_num --hmmn 200 && \n";
}

close OUT;

$cpu_threads_num = ($cpu_threads_num < 50)? 50 : $cpu_threads_num ;
&qsubOrDie("$align_shell_file",$queue,$cpu_threads_num,"6G");

# ------------------------------------------------------------------
# cds and pep prediction combine
# ------------------------------------------------------------------
my $Q_basename = $Q_name;
$Q_basename =~s/.fa$//;

my @cds_file = glob "$odir/Pfam_Dir/$Q_name*.best_candidates.final.cds.fa";
my @pep_file = glob "$odir/Pfam_Dir/$Q_name*.best_candidates.final.pep.fa";
open (ACDS, ">$cds_predict_dir/$Q_basename.cds.fa") or die;
$/ = ">";

for my $cds_file (@cds_file) {
    open (CDS, $cds_file) or die;
    <CDS>;

    while (<CDS>) {
        chomp;
        my ($id,$seq) = (split /\n/,$_,2);
        if ($id =~ /^(\S+) \S+\s+(.+)$/ ) {
            my ($mid, $inf) = ($1, $2);
            $mid =~ s/\|m\.\w+$//;
            $id = "$mid $inf";
        }
        print ACDS ">$id\n$seq";
    }

    close CDS;
}

close ACDS;
open (APEP, ">$cds_predict_dir/$Q_basename.pep.fa") or die;

for my $pep_file (@pep_file) {
    open (PEP, $pep_file) or die;
    <PEP>;

    while (<PEP>) {
        chomp;
        my ($id,$seq) = (split /\n/,$_,2);
        if ($id =~ /^(\S+) \S+\s+(.+)$/ ) {
            my ($mid, $inf) = ($1, $2);
            $mid =~ s/\|m\.\w+$//;
            $id = "$mid $inf";
        }
        print APEP ">$id\n$seq";
    }
    close PEP;
}

$/ = "\n";
close APEP;

# len stat & plot
system "perl $Bin/util/len_stat_plot.pl -i $Q_basename.cds -fa $cds_predict_dir/$Q_basename.cds.fa -od $cds_predict_dir/ ";

# ------------------------------------------------------------------
# pfam annotation combine
# ------------------------------------------------------------------

my @pfam_anno_file = glob "$odir/Pfam_Dir/$Q_name*.best_candidates.final.pfam.anno";
my @pfam_tab_file = glob "$odir/Pfam_Dir/$Q_name*.best_candidates.final.pfam.tab";

open (AANNO, ">$Tab_dir/$Q_name.Pfam.align.anno") or die;
print AANNO "#Gene_ID\tPfam_IDs\tPfam_Description\n";

for my $pfam_anno_file (@pfam_anno_file) {
    open (ANNO, $pfam_anno_file) or die;

    while (<ANNO>) {
        next if (/^\s+$/ or /^#/);
        print AANNO $_;
    }

    close ANNO;
}

close AANNO;

open (ATAB, ">$Tab_dir/$Q_name.Pfam.align.tab") or die;
my $pfam_tab_head = <<"_HEAD_";
#                                                                   --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name               accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ----------     -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
_HEAD_
print ATAB $pfam_tab_head;

for my $pfam_tab_file (@pfam_tab_file) {
    open (TAB, $pfam_tab_file) or die;

    while (<TAB>) {
        next if (/^\s+$/ or /^#/);
        $_ =~ s/(\|orf\d+)\|m\.\d+/$1/;
        print ATAB $_;
    }
    close TAB;
}

$/ = "\n";
close ATAB;

system "cp $Tab_dir/$Q_name.Pfam.align.tab $Result_dir/$Q_name.Pfam.anno.details";
system "cp $Tab_dir/$Q_name.Pfam.align.anno $Result_dir/$Q_name.Pfam.anno.txt";

#######################################################################################
&timeLog("$Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir=`pwd`; chomp($cur_dir);
    my ($in)=@_;
    my $return="";
    if(-f $in){
        my $dir=dirname($in);
        my $file=basename($in);
        chdir $dir; $dir=`pwd`; chomp $dir;
        $return="$dir/$file";
    }elsif(-d $in){
        chdir $in; $return=`pwd`; chomp $return;
    }else{
        warn "Warning just for file and dir\n";
        exit;
    }
    chdir $cur_dir;
    return $return;
}

################################################################################################################
sub MKDIR { # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if (-d $dir);
	mkdir($dir) if (!-d $dir);
}

################################################################################################################
sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	chomp(my $line = `wc -l $shell|cut -f 1 -d" "`);
	if ($line<=1000) {
		#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $cpu --resource vf=$vf --independent --reqsub $shell ";
	}
	if ($line>1000) {
		my @div=glob "$shell.div*sh";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		my $div_index=1;
		my $line_num=1;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index.sh" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=0;
				close OUT;
			}
		}
		if ($line_num!=0) {
			close OUT;
		}
		@div=glob "$shell.div*sh";
		foreach my $div_file (@div) {
			#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $cpu --resource vf=$vf --reqsub  $div_file";
		}
	}
}

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
#------------------------------------------------------------------------------
Program: $Script
Version: $version
Contact: Simon Young <simonyoung8824\@gmail.com>
   Data: 2014-09-28
Fuction: the script is used to search Pfam domain.
  Usage:
        --idir  <STR>   input transcript dir.
        --pfam  <STR>   path to Pfam database, *.hmm.                    [/share/nas19/yangxh/Database/pfam/27.0/Pfam-AB.hmm]
        --cpu   <INT>   number of threads to use to searche Pfam domain. [4]
        --odir  <STR>   output directory.                                [./]

Example:
    perl $Script --idir mid/ --odir Uni_Anno/

#------------------------------------------------------------------------------
USAGE
    print $usage;
    exit;
}