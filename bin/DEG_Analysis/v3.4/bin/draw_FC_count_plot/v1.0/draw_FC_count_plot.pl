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
my ($t,$c,$all,$de,$o);
GetOptions(
				"help|?" =>\&USAGE,
				"treated:s"=>\$t,
				"control:s"=>\$c,
				"all:s"=>\$all,
				"de:s"=>\$de,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($t and $c and $all and $de and $o);
&USAGE unless ($o=~/\.png$/);

$all=&ABSOLUTE_DIR($all);
$de=&ABSOLUTE_DIR($de);

###### set env variables
my $temp = `echo \$PATH`;
print "PATH=$temp\n";
$ENV{"PATH"} = ${&readconf("$Bin/../../../../../../config/sys.cfg")}{PATH}.":" . $ENV{"PATH"};
$temp = `echo \$PATH`;
print "PATH=$temp\n";

my $str=$t.'_'.$c;
$str=~s/_/","/g;
my %DE;
open (IN,$de) or die $!;
my %num;my $total;
while (<IN>) {
	next if /^\#/ or /^\s*$/;
	$_=~s/\s*$//;
	my ($name,$regulated)=(split/\s+/,$_)[0,-1];
	$DE{$name}=$regulated;
	$num{$regulated}++;
}
close IN;
my @TF;
chomp($total=`less -S $all|wc -l`);
$num{up}=0 if (!exists $num{up});
$num{down}=0 if (!exists $num{down});
$num{unchange}=$total-1-$num{up}-$num{down};
open (IN,$all) or die $!;
while (<IN>) {
	next if /^\#/;
	my $name1=(split/\s+/,$_)[0];
	my $type=(exists $DE{$name1})?"$DE{$name1}":"unchange";
	push @TF,"$type:$num{$type}" if $type ne "unchange";
	push @TF,"$type:$num{unchange}" if $type eq "unchange";
}
close IN;
my $TF=join '","',@TF;


open (OUT,">$o.r") or die $!;
print OUT "#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript\n";
print OUT <<END;
library(ggplot2)
ALL<- read.delim("$all", row.names = 1, header=TRUE,check.names=F)
Significant<-c("$TF")
v<-c("$str")
png(filename="$o", height = 3000, width = 3000, res = 500, units = "px")
df=data.frame(x=log2(rowMeans(ALL[,v])), y=ALL[,c("log2FC")],lab=factor(Significant,levels=unique(Significant)))
png(filename="$o", height = 3000, width = 3000, res = 500, units = "px")
p=ggplot(data=df,mapping=aes(x=log2(rowMeans(ALL[,v])), y=ALL[,c("log2FC")],color=lab))
p=p+geom_point(size=1)
p=p+xlab("log2(FPKM)")+ylab("log2(FC)")+ggtitle("MA plot")+theme(plot.title=element_text(face="bold",size=14))
p=p+theme_classic()
p=p+scale_color_manual(name="Significant",values=c("red","green","black"),limits=c("up:$num{up}","down:$num{down}","unchange:$num{unchange}"))
# add panel border 
p=p+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
print(p)
#qplot(log10(rowMeans(ALL[,v])), ALL[,c("log2FC")], xlab="log2(FPKM)", ylab="log2(FC)", main="MA plot", size=I(1.0), colour=Significant)
dev.off()
END

my $cmd;
#$cmd = "ssh compute-0-14 /share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript $o.r";

$cmd =${&readconf("$Bin/../../../../../../config/sys.cfg")}{Rscript}." $o.r";
print "$cmd\n";
$temp = `$cmd`;





#######################################################################################
&timeLog("$Script Done. Total elapsed time : ".time()-$BEGIN_TIME."s");
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
Program Date:   2013.10.14
Usage:
  Options:
  -treated    <str>    names of treated samples,"T1_T2",forced 
  
  -control    <str>    names of control samples,"T3_T4",forced 
  
  -all        <file>   All genes expression list of Group,forced 
  
  -de         <file>   different genes expression list of Group,forced 
  
  -o          <file>   output file,*.png,forced 
 
  -h           Help

USAGE
	print $usage;
	exit;
}
