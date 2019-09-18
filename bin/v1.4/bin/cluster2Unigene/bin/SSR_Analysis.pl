#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $version="1.1.0";
my $BEGIN_TIME=time();
##############################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作；
my %sys_cfg=%{&readconf("$Bin/../../../../../../config/sys.cfg")};
my ($fa, $od, $len);

GetOptions(
    "help|?" => \&USAGE,
    "fa:s" => \$fa,
    "od:s" => \$od,
    "len=i" => \$len,
) or &USAGE;
&USAGE unless ($fa and $od);

`mkdir $od ` if(!-d $od);
$od=&ABSOLUTE_DIR($od);
$fa=&ABSOLUTE_DIR($fa);
$len ||= 1000;

print STDOUT "\n[".&date_time_format($BEGIN_TIME)."] $Script start ...\n";
#######################
# --------------------------------------------------------------------
# SSR Analysis
# --------------------------------------------------------------------
open (IN,"$fa") || die "$!";
my $name = basename($fa);
$name =~/(.*)\.fa/;
my $prefix=$1;

open (OUT,">$od/$prefix.$len.fa") || die "$!";
$/='>';
<IN>;
while (<IN>) {
	chomp;
	my ($id,$seq)=split /\n+/,$_,2;
	my $seq_id=(split /\s+/,$id)[0];
	$seq=~s/\s+//g;
	my $length=length $seq;
	if ($length>=$len) {
		print OUT ">$seq_id\n$seq\n";
	}
}
close IN;
close OUT;

print STDOUT "\n[".&date_time_format(time())."] filting of short sequence that less than $len nt is done! ...ready for MISA analysis ...\n";

chdir $od;
`cp $Bin/MISA/misa.ini ./ `;
`perl $Bin/MISA/misa.pl $prefix.$len.fa `;
`perl $Bin/MISA/SSR.pl -key $prefix.$len.fa `;
`perl $Bin/MISA/SSR_density_plot.pl -i $prefix.$len.fa -p $prefix.$len.fa`;

print STDOUT "\n[".&date_time_format(time())."] MISA analysis done. start primer design...\n";
# --------------------------------------------------------------------
# SSR primer design by using primer3
# --------------------------------------------------------------------
#my $primer3_core = "/share/nas2/genome/biosoft/primer3/primer3_core"; # 2014-12-10 ~ 
my $primer3_core = $sys_cfg{"primer3_core"};
`perl $Bin/MISA/p3_in.pl $prefix.$len.fa.misa`;
`$primer3_core <$prefix.$len.fa.ssr.p3in >$prefix.$len.fa.ssr.p3out`;
`perl $Bin/MISA/p3_out.pl $prefix.$len.fa.ssr.p3out $prefix.$len.fa.misa`;
`rm $prefix.$len.fa.ssr.p3in`;

system "perl $Bin/combine_SSR_result.pl -p $prefix.$len.fa.ssr.primer.xls -s $prefix.$len.fa.misa -o $prefix.$len.fa.SSR.result.xls";

###############Time
print STDOUT "\n[".&date_time_format(time())."] $Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n";
###########subs
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
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR { # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub USAGE {#
    my $usage =<<"USAGE";
#-------------------------------------------------
 Program: $Script
 Version: $version
 Contact: Meng Fei <mengf\@biomarker.com.cn>
    Data: 2012-00-00
Modifier: Simon Young <simonyoung8824\@gmail.com>
    Data: 2014-10-15
 Fuction: SSR analysis & SSR primer design by using primer3
   Usage:
        -fa     <STR>    unigene seqence file, FASTA format.
        -od     <STR>    output directory.
        -len    <INT>    length cutoff for SSR Analysis         [1000]

#-------------------------------------------------
USAGE
    print $usage;
    exit;
}
