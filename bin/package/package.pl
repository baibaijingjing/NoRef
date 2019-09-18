#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my $indir;

#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut);

GetOptions(
                "help|?" =>\&USAGE,
                "in:s"=>\$indir,
                ) or &USAGE;
&USAGE unless ($indir);
# ------------------------------------------------------------------
#
# ---indir is program folder--------------------------------------------------------------
$indir = ABSOLUTE_DIR("$indir");

system "mkdir -p $indir/../Needed_Data" unless -d "$indir/../Needed_Data";

my $Results_zip = -e  "$indir/../Needed_Data/biomarker_Web_Report.zip";
my $html_zip = -e  "$indir/../Needed_Data/biomarker_htmlReport.zip";

if (! $html_zip) {
   system "cd $indir && zip  -r  -q $indir/../Needed_Data/biomarker_htmlReport.zip index.html ";
   system "cd $indir && zip  -r  -q $indir/../Needed_Data/biomarker_htmlReport.zip src/  HTML/";
}

if (! $Results_zip) {
    system "cd $indir && zip  -r  -q $indir/../Needed_Data/biomarker_Web_Report.zip $indir";
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

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

################################################################################################################

sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   Yang nan <yangn\@biomarker.com.cn> 
Program Date:   2016.03.08
      Modify:   
 Description:   This program is used to package the zip file to Need_Data......
       Usage:
        Options:
        -in <dir>   input BMK_Results directory,forced

        -h      help

USAGE
    print $usage;
    exit;
}
