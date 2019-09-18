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
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOrf,$exp,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"orf:s"=>\$fOrf,
				"exp:s"=>\$exp,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($fIn and $exp and $o);

my %U;
my %M;
open (IN,"samtools view $fIn|") or die $!;
while (<IN>) {
	my ($tag,$trans,$site)=(split/\t/,$_)[0,2,3];
#-----------------------------------------------------------
#	my ($tag,$trans,$start,$seq)=(split/\t/,$_)[0,2,3,9];
#   my $site = $start + length($seq)/2;

	next if $trans=~/\*/;
	next if exists $M{$tag};
	$U{$tag}{$trans}=$site unless exists $U{$tag}{$trans};
	$U{$tag}{$trans}=($U{$tag}{$trans}+$site)/2 if exists $U{$tag}{$trans};
	my $limit=keys %{$U{$tag}};
	if ($limit>1){
		$M{$tag}=1;
		delete $U{$tag};
		next;
	}
}
close IN;

my %ST;
if (defined $fOrf) {
	$/=">";
	open (IN,"$fOrf") or die $!;
	<IN>;
	while (<IN>) {
		chomp;
		my ($genename,$strand)=(split/\s+/,$_)[0,2];
		if ($strand=~/-/) {
			$ST{$genename}=1;
		}
	}
	close IN;
	$/="\n";
}

my %G;
my %T_G;
open (IN,"$exp") or die $!;
while (<IN>) {
	next if /^\#/;
	my ($gene,$len,$trans)=(split/\s+/,$_)[0,2,5];
	foreach my $t (split/,/,$trans) {
		$T_G{$t}=$gene;
	}
	$G{$gene}{len}=$len;
}
close IN;


my %num;
foreach my $tag (keys %U) {
	foreach my $trans (keys %{$U{$tag}}) {
		my $gene=$T_G{$trans};
		next if $G{$gene}{len}<500;
		my $site=(exists $ST{$gene})?(100-100*$U{$tag}{$trans}/$G{$gene}{len}):(100*$U{$tag}{$trans}/$G{$gene}{len});
		$site=int $site;
		$num{$gene}{total}++;
		$num{$gene}{site}{$site}++;
	}
}

my $out=$o.".randcheck";

open (OUT,">$out.list")||die "$!";

my %per;
foreach my $key1 (keys %num) {
	foreach my $key2 (keys %{$num{$key1}{site}}) {
		$per{$key2}+=$num{$key1}{site}{$key2}/$num{$key1}{total};
	}
}
my $gene_num=keys %num;
foreach my $key (keys %per) {
	$per{$key}=sprintf "%.2f",$per{$key}/$gene_num*100;
}
my $max_y=0;
foreach (values %per) {
	if ($_>$max_y) {
		$max_y=$_;
	}
}
print "$max_y\n";
my ($y,$ystep)=&define_axis ($max_y);
print OUT << "Usage End.";
Type:Line
Width:600
Height:400
WholeScale:0.9
FontSize:25
X:Relative Position in Genes(5'-3')
Y:Percent of Reads
XStep:10
YStep:$ystep
XStart:0
YStart:0
XEnd:100
YEnd:$y

Color:red
Usage End.
foreach  (sort{$a<=>$b}keys %per) {
	next if $_>100;
	next if $_<0;
	print OUT "$_:$per{$_}\n";
}
close OUT;
my $distributing_svg = ${&readconf("$Bin/../../../../../../../config/sys.cfg")}{distributing_svg};
system("perl $distributing_svg $out.list $out.svg");
my $svg_file="$out.svg";
my $svg_dir=dirname $svg_file;
my $svg_name=basename $svg_file;
chdir $svg_dir;
my $svg2xxx = ${&readconf("$Bin/../../../../../../../config/sys.cfg")}{svg2xxx};
system("$svg2xxx $svg_name");


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub define_axis () {
        my $i=0;
        my ($max)=@_;
        my $time=1;
        my @ret=();
        my @limit=(1,2,3,4,5,6,8,10,12,15,16,20,24,25,30,40,50,60,80,100,120);
        my @unlim=(0,1,2,3,4,5,6,8 ,10,11,14,15,18.5,21,23,29,37,47,56,76 ,92 ,110);

        while ($max >$unlim[21]) {
                 $max=$max/10;
                 $time=$time*10;
        }
        for ($i=0;$i<=20 ;$i++) {
                 if ($max>$unlim[$i] && $max<=$unlim[$i+1]) {
                         $ret[0]=$limit[$i]*$time;

                         if ($i==2 || $i==5 || $i==9 || $i==14) {
                                 $ret[1]=$ret[0]/3;
                         }
                         elsif ($i==4 || $i==7 || $i==13 || $i==16){
                                 $ret[1]=$ret[0]/5;
                         }
                         else {
                                 $ret[1]=$ret[0]/4;
                         }

                 }
        }
        return @ret;
}

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
	#���б��е�����ֵ
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
Usage:
  Options:
  -i    <file>  input file,bam format,forced

  -exp  <file>  expression file,.geneExpression.xls,forced

  -orf  <file>  orf result,*.Unigene.cds_pep.stat.xls

  -o    <file>  output file,forced

  -h         Help

USAGE
	print $usage;
	exit;
}
