#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME = time();
my $version = "1.0.0";
##########################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($analysis_issue_id, $indir, $odir, $backup_fin, $backup_med, $check);

GetOptions(
    "help|?" =>\&USAGE,
    "ai:s" => \$analysis_issue_id,
    "id:s" => \$indir,
    "od:s" => \$odir,
    "bf:s" => \$backup_fin,
    "bm"   => \$backup_med,
    "ck"   => \$check,
) or &USAGE;
&USAGE unless ($analysis_issue_id and $indir and $odir);
unless ($analysis_issue_id =~/^BMK\d{6}-\w{3}(-\d{2})?(-[a-zA-Z])?$/) {
     &writeLog("ERROR: Analysis Issue is illegal.");
    &USAGE;
}

system "mkdir -p $odir" unless (-d $odir);
$indir = ABSOLUTE_DIR($indir);
$odir = ABSOLUTE_DIR($odir);
$backup_fin ||= 'both';

#print STDOUT "\n[".&GetTime($BEGIN_TIME)."] $Script start ...\n";
&timeLog("$Script start ...");
# --------------------------------------------------------------------
# check analysis directory
# --------------------------------------------------------------------

if ($check) {
    my $error_files_print='';
    open (IN, "ls `find $indir -name *.error`|") or die;

    while (<IN>) {
        chomp;
        my @stat=stat($_);
        next if ($stat[7] == 0);
        $error_files_print .= "$_\n";
    }

    close IN;

    &writeLog( "The following files are not empty, related task may be unfinished. please check!\n$error_files_print\n" ) unless ($error_files_print eq '');

    &timeLog("analysis directory check done.");
}

opendir (DIR, "$indir/Web_Report") or die $!;
my @file = grep { !/^\./ } readdir DIR;
die "$indir/Web_Report is a empty directionary.\n" if (scalar @file == 0);

# --------------------------------------------------------------------
# backup all final analysis result
# --------------------------------------------------------------------
unless ($backup_fin eq 'part') {
    &mkdirOrDie ("$odir/${analysis_issue_id}_Transcriptome_final") unless (-d "$odir/${analysis_issue_id}_Transcriptome_final");
    system "cp -r $indir/Web_Report/* $odir/${analysis_issue_id}_Transcriptome_final/";
	&timeLog("all final analysis result backup done.");
    #print "[".&GetTime()."] all final analysis result backup done.\n";
}

# --------------------------------------------------------------------
# backup one third of final analysis result
# --------------------------------------------------------------------
unless ($backup_fin eq 'all') {
    mkdir "$odir/${analysis_issue_id}_Transcriptome_partial" unless (-d "$odir/${analysis_issue_id}_Transcriptome_partial");
    system "cp -r $indir/Web_Report/cleandata $odir/${analysis_issue_id}_Transcriptome_partial/";
    system "cp -r $indir/Web_Report/geneExpression $odir/${analysis_issue_id}_Transcriptome_partial/";
    print "[".&GetTime()."] one third of final analysis result backup done.\n";
}

# --------------------------------------------------------------------
# backup essential intermediate result.
# --------------------------------------------------------------------
if ($backup_med) {
    mkdir "$odir/Analysis" unless (-d "$odir/Analysis");

    if (-d "$indir/work_sh") {
        system "cp -r $indir/work_sh $odir/Analysis/";
    } else {
        &writeLog ("WARNNING: $indir/work_sh is lost.");
    }

    if (-d "$indir/Config") {
        system "cp -r $indir/Config $odir/Analysis/";
    } else {
        &writeLog ("WARNNING: $indir/Config is lost.");
    }

    &mkdirOrDie ("$odir/Analysis/Unigene_Annotation") unless (-d "$odir/Analysis/Unigene_Annotation");
    if (-d "$indir/Uni_Anno/02.gene-annotation" and -d "$indir/Uni_Anno/Result") {
        system "cp -r $indir/Uni_Anno/02.gene-annotation/ $odir/Analysis/Unigene_Annotation/";
        system "cp -r $indir/Uni_Anno/Result/ $odir/Analysis/Unigene_Annotation/";
    } else {
        &writeLog ("WARNNING: Unigene Annotation result is lost.");
    }

    my @files = (glob "$indir/QC_Report/*.QC.stat");
    if (@files!=0) {
        system "cp $indir/QC_Report/*.QC.stat $odir/Analysis/";
    } else {
        &writeLog ("WARNNING: Analysis Quality Control result is lost.");
    }

    @files = (glob "$indir/*.rtf");
    if (@files!=0) {
        system "cp $indir/*.rtf $odir/Analysis/";
    } else {
        &writeLog ("WARNNING: Analysis Report is lost.");
    }

    &timeLog("essential intermediate result backup done.");
}

# --------------------------------------------------------------------
# 
# --------------------------------------------------------------------



##########################################################################################
&timeLog("$Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s");
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
####################################################################################################
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir=`pwd`;chomp($cur_dir);
    my ($in)=@_;
    my $return="";
    if(-f $in){
        my $dir=dirname($in);
        my $file=basename($in);
        chdir $dir;$dir=`pwd`;chomp $dir;
        $return="$dir/$file";
    }elsif(-d $in){
        chdir $in;$return=`pwd`;chomp $return;
    }else{
        warn "Warning just for file and dir\n";
        exit;
    }
    chdir $cur_dir;
    return $return;
}

####################################################################################################
sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

####################################################################################################
sub USAGE {
    my $usage =<<"USAGE";
#-------------------------------------------------
Program: $Script
Version: $version
Contact: Simon Young <simonyoung8824\@gmail.com>
   Data: 2014-09-10
Fuction: the script is used to backup no ref trans analysis result.
  Usage: 
    --ai  <STR>    no ref trans analysis issue id.
    --id  <STR>    no ref trans final analysis directory.
    --od  <STR>    output directory.

    --bf           way to backup final result.                   ['both']
          part     one third of final analysis result
          all      all final analysis result
          both     [default]
    --bm           backup essential intermediate result.
    --ck           check result anaysis.

Example:
    perl $Script --ai BMK140606-C67 --id Data_Analysis/BMK140606-C67_Mink/ --od Backup/BMK140606-C67_Mink/ --md --ck
    perl $Script --ai BMK140606-C67 --id BMK140606-C67_Mink/ --od ../Backup/BMK140606-C67_Mink/ --br all

#-------------------------------------------------
USAGE
    print $usage;
    exit;
}