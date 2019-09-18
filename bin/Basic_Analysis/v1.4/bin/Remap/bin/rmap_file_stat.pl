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

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fF1,$file1Format,$fK);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"file1:s"=>\$fF1,
				"file1Format:s"=>\$file1Format,
				"od:s"=>\$fK,
				) or &USAGE;
&USAGE unless ($fIn and $file1Format and $fK);

mkdir $fK unless (-d $fK) ;
$fK=&ABSOLUTE_DIR($fK);
$fF1=&ABSOLUTE_DIR($fF1);
$fK=~/\/([^\/]+)$/;
my $Key=$1;

chdir $fK;

print "Statistics maped info :$fIn\n";
	
my $Total_Reads;
my $Total_line;
if ($file1Format eq "fastq") {
	$Total_line=`less -S $fF1| wc -l`;chomp $Total_line;
	$Total_Reads=$Total_line/4;
}elsif($file1Format eq "fasta"){
	$Total_Reads=`less -S $fF1|grep >|wc -l`;chomp $Total_Reads;
}else{
	die "Error in Total_Reads\n";
}
	
my %STAT;

open (ST,"<","$fIn") or die $!;
my $Line=<ST>;chomp $Line;
while ($Line) {
	chomp $Line;
	my ($id,$misMatch,$inDel,$hits)=(split "\t",$Line)[0,11,12,13];
	$STAT{"Perfect"}++ if ($misMatch==0 and $inDel==0 and $hits==1) ;
	$STAT{"MisMatch"}{$misMatch}++;
	$STAT{"InDel"}{$inDel}++;
	$hits==1?$STAT{"Uniq"}++:$STAT{"Multi"}++;
	while ($Line=<ST>) {
		chomp $Line;
		my ($Newid)=split "\t",$Line;
		last if ($Newid ne $id) ;
	}
}
close (ST) ;

my $TotalMappedReads=$STAT{"Uniq"}+$STAT{"Multi"};
open (STOUT,">","$Key.Mapped.stat.xls") or die $!;
print STOUT "Total Reads","\t",$Total_Reads,"\t","100%\n";
print STOUT "Mapped Reads","\t",$TotalMappedReads,"\t",int(10000*$TotalMappedReads/$Total_Reads)/100,"%\n";
#print STOUT "Perfect Mapped Reads","\t",$STAT{"Perfect"},"\t",int(10000*$STAT{"Perfect"}/$TotalMappedReads)/100,"%\n";
#print STOUT "\nMisMatch\n";
#foreach my $elem (sort {$a <=> $b} keys %{$STAT{"MisMatch"}}) {
	#print STOUT $elem,"\t",$STAT{"MisMatch"}{$elem},"\t",int(10000*$STAT{"MisMatch"}{$elem}/$TotalMappedReads)/100,"%\n";
#}

#print STOUT "\nInDel\n";
#foreach my $elem (sort {$a <=> $b} keys %{$STAT{"InDel"}}) {
#	print STOUT $elem,"\t",$STAT{"InDel"}{$elem},"\t",int(10000*$STAT{"InDel"}{$elem}/$TotalMappedReads)/100,"%\n";
#}

#print STOUT "\nUniq\/Multi\n";
$STAT{"Uniq"}||=0;
$STAT{"Multi"}||=0;
print STOUT "Uniq mapped Reads\t",$STAT{"Uniq"},"\t",int(10000*$STAT{"Uniq"}/$TotalMappedReads)/100,"%\n";
print STOUT "Multi mapped Reads\t",$STAT{"Multi"},"\t",int(10000*$STAT{"Multi"}/$TotalMappedReads)/100,"%\n";

close (STOUT) ;


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
	#���б��е����ֵ
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
	#���б��е���Сֵ
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
	#��ȡ�ַ������еķ��򻥲����У����ַ�����ʽ���ء�ATTCCC->GGGAAT
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
Program Date:   2013.11.21
Usage:
  Options:
  -i              <file>  input file,rmap format,forced 
  
  -file1          <file>  cleandata file,forced
  
  -file1Format    <str>   fastq or fasta,forced
  
  -od             <dir>   output dir,forced 
  
  -h         Help

USAGE
	print $usage;
	exit;
}
