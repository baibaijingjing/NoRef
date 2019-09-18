#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$limit_max);

GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"m:s"=>\$limit_max,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
$limit_max||=10;
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
&nr_pie_stat($fIn,$fOut);
my $pre=$fOut;
$pre=~s/.stat//;
`perl $Bin/Just_Pie.pl -i $fOut -o $pre.svg  -w 800 -css $Bin/pie12.css -note "Nr Homologous Species Distribution" `;

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
sub nr_pie_stat{#&nr_pie_stat($in,$out)
	my ($in,$out)=@_;
	my %H;
	my $total;
	open (IN,$in) or die $!;
	<IN>;
	while (<IN>) {
		chomp;
		my @A=split/\t/,$_;;
		my $name=$A[-1];
		$name=~s/\s*$//;
		$total++;
                next unless $name=~/\[([^\]]+)\]$/;
		$H{$1}=0  unless exists $H{$1};
		$H{$1}++ if exists $H{$1};
	}
	# foreach my $key (keys %H) {
		# my $ration=sprintf("%.2f",$H{$key}/$total*100);
		# $H{$key}="$H{$key}";
	# }
	close (IN) ;
	open (OUT,">$out") or die $!;
    print OUT "#Species_Name\tHomologous_Number\tFull_Name\n";
	my @name=sort {$H{$b} <=> $H{$a}} keys %H;
	my $ref=&str_name(@name);
	my $limit=keys %H;
	if ($limit<=$limit_max+1) {
		if ($limit == 1) {
			my ($key,$value) = each %H;
			my $str = &cut_str($key);
			print OUT "$str\t$value\t$key\n";
		}else{
			foreach my $key (sort {$H{$b}<=>$H{$a}} keys %H) {
				#my $str = &cut_str($key);
				my $str=$$ref{$key};
				print OUT "$str\t$H{$key}\t$key\n";
			}
		}
	}
	else {
		my $n=0;
		my $other=0;
		my $ration_total=0;
		foreach my $key (sort {$H{$b} <=> $H{$a}} keys %H) {
			$n++;
			if($n<$limit_max+1){
				#my $str = &cut_str($key);
				my $str=$$ref{$key};
				print OUT "$str\t$H{$key}\t$key\n";
			}
		$other+=$H{$key} unless $n<$limit_max+1;
		}
	print OUT "Other\t$other\tOther species\n";
	}
	close OUT;
}
#####################
sub str_name{
	my @name=@_;
	my %nameHash;
	my %lastName;
	my $len;
	if (@name>$limit_max+1){
		$len=$limit_max;
	}
	else{
		$len=@name;
	}
	my $str_num=2;
	for (my $i=0;$i<$len;$i++){
		my @str = split /\s+/,$name[$i];
		if(@str>$str_num){
			my $nameStr=join(" ",@str[0..$str_num-1]);
			if(exists $nameHash{$nameStr}){
				$i=-1;
				$str_num++;
				undef(%nameHash);
				next;
			}
			$nameHash{$nameStr}=1;
			$lastName{$name[$i]}=$nameStr;
		}
		else{
			$nameHash{$name[$i]}=1;
			$lastName{$name[$i]}=$name[$i];
		}
	}
	return (\%lastName);
}
################################################################################################################
sub cut_str {
	my $string = shift;
	my @str = split /\s+/,$string;
	if (@str > 2) {
		return "$str[0] $str[1]"
	}else{
		return $string;
	}
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
     Version:	$version
     Contact:	Simon Young <yangxh\@biomarker.com.cn> 
Program Date:	2012.07.02
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		-i <file>	Annotation/Result/*.nr.anno.txt

		-o <file>	Annotation/Result/*.nr.lib.stat

		-m <int>	maximum species number to show [10]

		-h		help

USAGE
	print $usage;
	exit;
}
