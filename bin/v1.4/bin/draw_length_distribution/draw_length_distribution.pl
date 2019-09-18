#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Find;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
use newPerlBase;
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($in,$n,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$in,
				"n:s"=>\$n,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($in and $od);

$n||=0;
mkdir $od unless -d $od;
$od=&ABSOLUTE_DIR($od);


my @sample;
if ($n == 1) {
	@sample= glob "$in/Cluster/Unigene/*.stat.info.xls";
	push @sample,"$in/Assembly/Trinity_assembly/All_Combination/All_Combination_Result/Transcripts/All_Combination.Transcripts.stat.info.xls";
}else {
	my $statname;
	find(\&wanted,$in);
	sub wanted {
	         if ($_ =~ /\.Unigenes?\.stat\.info\.xls$/)
		{
				  $statname.="$File::Find::name";
				  $statname=$statname."YYY";
			 }
}
 @sample=split /YYY/,$statname;
 push @sample,(glob "$in/Cluster/Unigene/*.stat.info.xls")[0];
 }
my %H;
foreach my $file (@sample) {
	$file=~/([^\/]+)\.stat\.info\.xls$/;
	my $name=$1;
	open (IN,$file) or die $!;
	<IN>;
	while (<IN>) {
		my @A=split/\s+/,$_;
		if ($A[0]=~/>/){
			$A[0]=~s/>//;
			$A[0]+=50;
			$H{$name}{$A[0]}=$A[2];
		}
		else{
			my @B=split/\D/,$A[0];
			my $site=($B[1]+$B[0])/2;
			$H{$name}{$site}=$A[2];
		}
	}
	close IN;
}
open (OUT,">$od/Total.list") or die $!;
print OUT "Type\tSite\tPercent\n";
foreach my $sam (keys %H) {
	foreach my $s (sort {$a<=>$b} keys %{$H{$sam}}) {
		print OUT "$sam\t$s\t$H{$sam}{$s}\n";
	}
}
close OUT;

my $Rscript=${&readconf("$Bin/../../../../../config/sys.cfg")}{Rscript};
`$Rscript $Bin/bin/draw_length_distribution.r $od/Total.list $od/AssemblyLengthDis.png`;



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

sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Usage:
  Options:
		-i					Basic_Dir		must be given
		-od					output dir		must be given
		-n					only one group(1) or more than one group(0) ,default 0
		-h					Help document


USAGE
	print $usage;
	exit;
}
