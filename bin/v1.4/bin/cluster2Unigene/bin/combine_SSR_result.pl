#!/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use File::Basename;

my ($primerFile,$SSRFile,$outFile,$help,$outFile);
GetOptions(
	'p=s'=>\$primerFile,
	's=s'=>\$SSRFile,
	'o=s'=>\$outFile,
	'h|help'=>\$help,
);
usage() and exit unless ($primerFile and $SSRFile and $outFile);
$primerFile=ABSOLUTE_DIR($primerFile);
$SSRFile   =ABSOLUTE_DIR($SSRFile);

open P,"$primerFile";
my (%primer,$headTitle);
while(<P>){
	chomp;
	$headTitle=$_ and next if ($.==1);
	my ($geneID)=(split/\s+/)[0];
	$_=~s/$geneID\s+//;
	$primer{$geneID}{$_}=1;
}
close P;
$headTitle=~s/\#Gene_ID//;
$headTitle=~s/Start/PStart/g;
$headTitle=~s/End/PEnd/g;
$headTitle=~s/^\s+//;
$headTitle=~s/\s+$//;

open S,"$SSRFile";
open O,">$outFile";
while(<S>){
	chomp;
	if ($.==1){
		$_=~s/Start/SSR_Start/g;
		$_=~s/End/SSR_End/g;
		print O "$_\t$headTitle\n";
		next;
	}
	my ($geneID,$start,$end)=(split/\s+/)[0,5,6];
	if(exists $primer{$geneID}){
		my $hasPrimer=0;
		my $primer="";
		for my $info(sort {(split/\s+/,$a)[8] <=> (split/\s+/,$b)[8]} keys %{$primer{$geneID}}) {
			my ($s,$e)=(split/\s+/,$info)[7,8];
			if(&is_in("$s,,$e","$start,,$end") and !$hasPrimer){
				$primer.="\t$info";
				$hasPrimer=1;
			}
		}
		if($hasPrimer){
			print O $_.$primer."\n";
		}
		else{
			print O "$_\t--\n";
		}
	}
	else{
		print O "$_\t--\n";
	}
}
close S;close O;

##########Sub Function#########
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

sub is_in(){
	my ($mark1,$mark2)=@_;
	my ($s1,$e1)=split/,,/,$mark1;
	my ($s2,$e2)=split/,,/,$mark2;
	if($s2>=$s1 and $e2<=$e1 ){
		return 1;
	}
	return 0;
	
}

sub usage(){
	print STDOUT <<USAGE
Usage:$0 <options> 
	-p <string>			primer file out of MISA primer out file 
	-s <string>			SSR result of MISA
	-o <string>			out file
	-h -help			print usage information
USAGE
}
