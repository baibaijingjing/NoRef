#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME = time();
my $version = "1.0.0";
##########################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($ref_seq, $prefix);

GetOptions(
    "help|?" =>\&USAGE,
    "i:s" => \$ref_seq,
    "p:s" => \$prefix,
) or &USAGE;
&USAGE unless ($ref_seq and $prefix);

print STDOUT "\n[".&GetTime($BEGIN_TIME)."] $Script start ...\n";
# --------------------------------------------------------------------
# get refseq length
# --------------------------------------------------------------------
my $ref_len = 0;

open (REF, $ref_seq) or die $!;
$/="\n>";

while (<REF>) {
    chomp;
    s/^>//g if ($.==1);
    next if (/^\s+$/);
    my ($id,$seq) = (split /\n/,$_,2);
    $seq =~ s/\n//g;
    $ref_len += length $seq;
}

close REF;

$/ = "\n";
# --------------------------------------------------------------------
# SSR density stat
# --------------------------------------------------------------------
my %SSR_type_num;

open (IN, "$prefix.misa") or die;

while (<IN>) {
    next if (/^\s+$/ or /^#/);
    my $SSR_type = (split /\t/)[2];

    if (exists $SSR_type_num{$SSR_type}) {
        $SSR_type_num{$SSR_type}++;
    } else {
        $SSR_type_num{$SSR_type} = 1;
    }
}

close IN;

open (STAT,">$prefix.ssr.density.xls") or die;
print STAT "#type\tdensity\n";

for my $type (sort {$a cmp $b} keys %SSR_type_num) {
    my $density = $SSR_type_num{$type} * 1000000 / $ref_len; #SSR number per Mb
    print STAT "$type\t$density\n";
}

close STAT;

my $Rscript = "/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
my $x_lab = "Type of SSRs";
my $y_lab = "Number of SSRs per Mb";
my $title_lab = "SSR Density";

# plot
system "$Rscript $Bin/simpleBar.r --infile $prefix.ssr.density.xls --outfile $prefix.ssr.density.png --x.col 1 --y.col 2 --x.lab \"$x_lab\" --y.lab \"$y_lab\" --title.lab \"$title_lab\" --no.grid";

##########################################################################################
print STDOUT "\n[".&GetTime(time())."] $Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n";
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
   Data: 2014-10-15
Fuction: draw SSR density plot.
  Usage: 
    -i  <STR>    input file, FASTA format.
    -p  <STR>    input file prefix.

Example:
    perl $Script -i prefix.1000.fa -p prefix.1000.fa

#-------------------------------------------------
USAGE
    print $usage;
    exit;
}