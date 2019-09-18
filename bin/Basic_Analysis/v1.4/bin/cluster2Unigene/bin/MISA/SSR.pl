#!/usr/bin/perl -w
#
# Copyright (c) BMK 2011
# Writer:         He hua <heh@biomarker.com.cn>
# Program Date:   2011.
# Modifier:       He hua <heh@biomarker.com.cn>
# Last Modified:  2011.
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"key=s","h");
if (!defined($opts{key})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
###############
my $cur_pwd=`pwd`;
chomp($cur_pwd);
my $key=&ABSOLUTE_DIR($opts{key},$cur_pwd);

my $ssrfile=$key.".misa";
my $fa=$key;

my %stat;
my %types;
my %pos;
my %end;
my $total;
open SSR,"$ssrfile"||die "$!";
while (<SSR>) {
	chomp;
	next if(!/\d+/);
	my ($id,$type,$ssr,$start,$end)=(split/\s+/,$_)[0,2,3,5,-1];
	my $num;
	if ($type=~m/p/) {
#		$type=~s/p/perfect_SSR_/;
		$num=1;
	}elsif($type=~m/c/){
#		$type=~s/c/compound_SSR/;
#		$num=()=$ssr=~/\(/g;
        $num=1;
	}
	$types{$id}{$start}=$type;
	$pos{$id}{$start}{$ssr}=$end;
	$stat{$type}+=$num;
	$total+=$num;
}

my $out_key=basename($key);
open STAT,">$cur_pwd/$out_key.stat.xls"||die "$!";
print STAT "#type\tnumber\n";

foreach my $new_type (sort keys %stat) {
	print STAT "$new_type\t$stat{$new_type}\n";
}
print STAT "Total\t$total\n";
undef %stat;
close STAT;

open OUT,">$cur_pwd/$out_key.detail.xls"||die "$!";
print OUT "#Gene_ID\tLength\tSSR_type\tSSR\tStart\tEnd\tSequence\n";
$/=">";
open FA,$fa||die "$!";
while (<FA>) {
	chomp;
	next if(/^$/);
	my ($id,$seq)=split/\s+/,$_,2;
	$seq=~tr/AGCTN/agctn/;
	$seq=~s/\n//;
	if (exists $types{$id}) {
		foreach my $start (sort {$a<=>$b} keys %{$pos{$id}}) {
			foreach my $ssr (keys %{$pos{$id}{$start}}) {
				print OUT "$id\t",length($seq),"\t$types{$id}{$start}\t$ssr\t$start\t$pos{$id}{$start}{$ssr}\t";
				my $pos1=$start>300 ? $start-300 : 1;
				my $pos2=$pos{$id}{$start}{$ssr};
				my $pos2end=$pos2+300>=length($seq) ? length($seq) : $pos2+300;
				my $long=$start-$pos1;
				print OUT substr($seq,$pos1-1,$long);
##-----------------------------------------------------------------------------
#				$ssr=~s/\(//;
#				my @SSR=split/\(|\)/,$ssr;
#				if ($ssr=~m/\*/){
#					if (length($SSR[-4])>1) {
#						$SSR[-1]=~s/\*//;
#						$SSR[-1]=$SSR[-1]-1;
#					}else{
#						$SSR[-3]=$SSR[-3]-1;
#					}
#				}
#				for (my $i=0;$i<@SSR-1;$i+=2) {
#					$SSR[$i+1]=~m/(\d+)/;
#					print OUT $SSR[$i]x$1;
#					$SSR[$i+1]=~s/\d+|\*//g;
#					print OUT $SSR[$i+1];
#				}
##=================== modified by Simon Young on 2014-10-31 ===================
                $ssr=~s/\*//;
                while ($ssr =~/\(|\)/) {
                    my ($abbr,$unit,$unit_n) = $ssr =~ /(\(([ATCGN]+)\)(\d+))/;
                    my $unit_x_n = $unit x $unit_n;
                    $abbr =~s/\(/\\\(/; $abbr =~s/\)/\\\)/;
                    $ssr =~s/$abbr/$unit_x_n/;
                }
                print OUT $ssr;
##-----------------------------------------------------------------------------
				print OUT substr($seq,$pos2,$pos2end-$pos2),"\n";
			}
		}
	}
}
close OUT;

################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";

###############Subs
sub ABSOLUTE_DIR
{
        my ($in,$cur_dir)=@_;
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
                warn "Warning just for file and dir\n";
                exit;
        }
        chdir $cur_dir;
        return $return;
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime
{ # &Runtime($BEGIN);
        my ($t1)=@_;
        my $t=time()-$t1;
        print "\nTotal elapsed time: ${t}s\n";
}
sub help{
	print << "	Usage End.";
	Description:
		version:$ver
	Usage:

		-key      infile(key key.misa) and outfile(key.detail key.stat at current dir)          must be given;

	Usage End.
		exit;
}
