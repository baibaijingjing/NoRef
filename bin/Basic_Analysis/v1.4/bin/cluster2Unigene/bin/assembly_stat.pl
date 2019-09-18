#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $version="1.2.0";
my $BEGIN=time();

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $programe_dir=basename($0);
my $path=dirname($0);
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($Index,$fa,$region,$od,$upper_len,$syscfg);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$Index,
				"s=s"=>\$region,
				"fa=s"=>\$fa,
				"len=s"=>\$upper_len,
				"od=s"=>\$od,
				"syscfg=s"=>\$syscfg,
				) or &USAGE;
&USAGE unless ($Index and $fa and $od) ;

$syscfg||="$Bin/../../../../../../config/";
my $svg2xxx=${&selectconf("$syscfg")}{svg2xxx};
$fa=&ABSOLUTE_DIR($fa);
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);

$region = $region || "300,500,1000,2000";

my @range=split/,/,$region;
$upper_len = $upper_len || 30000;

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";

################################ Program
my %length;
my %gap;
my %length_region_dis;
my %upper_length;
my @sub_dis;

my $total_length;
my $count=0;

open IN,"$fa" || die $!;
open OUT,">$od/${Index}.Unigene.fa" || die $!;
$/='>';
while (<IN>) {
	chomp;
	s/^\s+//;s/\s+$//;s/\r+$//;
	next if (/^$/ || /^\#/);
	my ($head,$seq)=split/\n+/,$_,2;
	$count++;
	my $id=(split/\s+/,$head)[0];
	$seq=~s/\s+//g;
	my $len=length($seq);
	if ($len >= $upper_len) {
		$upper_length{$id}=$seq;
		next;
	}
	
	print OUT ">$id\n$seq\n";
	$length{$len}++;
	my $region_site=int($len/100);
	if ($region_site>=30) {
		$length_region_dis{30}++;
	}
	else {
		$length_region_dis{$region_site}++;
	}
	$total_length+=$len;
	if ($len<$range[0]) {
		$sub_dis[0]++;
	}
	elsif ($len>=$range[-1]) {
		$sub_dis[$#range+1]++;
	}
	my $i=1;
	while($i<scalar(@range)){
		if($range[$i-1]<=$len && $range[$i]>$len){
			$sub_dis[$i]++;
			$i++;
			last;
		}
		else{
			$i++;
			next;
		}
	}
}

close IN;
close OUT;

#####################
#Stat output
################################
open OUT,">$od/$Index.Unigene.exceed_upper_length.fa" || die $!;
foreach my $id (sort keys %upper_length) {
	my $id_length =length $upper_length{$id};
	print OUT ">$id\t$id_length\n$upper_length{$id}\n";
}
close OUT;

system "rm $od/$Index.Unigene.exceed_upper_length.fa" if (-e "$od/$Index.Unigene.exceed_upper_length.fa" and -z "$od/$Index.Unigene.exceed_upper_length.fa"); # delete empty file

my $n40=$total_length*0.4;
my $n50=$total_length*0.5;
my $n60=$total_length*0.6;
my $n70=$total_length*0.7;
my $n80=$total_length*0.8;
my $n90=$total_length*0.9;

my $min_len=(int((sort {$a <=> $b} keys %length)[0]/100))*100;
open (OUT0,">$od/$Index.Unigene.stat.xls");
print OUT0 "$Index Unigene Length\tTotal Number\tPercentage\n";
print OUT0 "$min_len"."-$range[0]\t$sub_dis[0]\t";
my $per=$sub_dis[0]/$count*100;
print OUT0 "$per%\n";
for(my $i=1;$i<scalar(@range);$i++) {
	print OUT0 "$range[$i-1]-$range[$i]\t$sub_dis[$i]\t";
	$per=$sub_dis[$i]/$count*100;
	print OUT0 "$per%\n";
}
print OUT0 "$range[-1]+\t$sub_dis[-1]\t";
$per=$sub_dis[-1]/$count*100;
print OUT0 "$per%\n";
print OUT0 "Total Number\t$count\n";
print OUT0 "Total Length\t$total_length\n";

my ($n40_len,$n50_len,$n60_len,$n70_len,$n80_len,$n90_len);
my $n40_index=0;
my $n50_index=0;
my $n60_index=0;
my $n70_index=0;
my $n80_index=0;
my $n90_index=0;
open (OUT1,">$od/$Index.Unigene.n50");
my $sum_length=0;
foreach  my $num(sort {$b <=> $a} keys %length) {
	$sum_length += $num*$length{$num};
	if($sum_length>=$n40 && $n40_index==0){
		$n40_len=$num;
		$n40_index=1;
	}
	if($sum_length>=$n50 && $n50_index==0){
		$n50_len=$num;
		$n50_index=1;
	}
	if($sum_length>=$n60 && $n60_index==0){
		$n60_len=$num;
		$n60_index=1;
	}
	if($sum_length>=$n70 && $n70_index==0){
		$n70_len=$num;
		$n70_index=1;
	}
	if($sum_length>=$n80 && $n80_index==0){
		$n80_len=$num;
		$n80_index=1;
	}
	if($sum_length>=$n90 && $n90_index==0){
		$n90_len=$num;
		$n90_index=1;
	}
}
print OUT1 "N40\tN50\tN60\tN70\tN80\tN90\n";
print OUT1 "$n40_len\t$n50_len\t$n60_len\t$n70_len\t$n80_len\t$n90_len\n";
close OUT1;
print OUT0 "N50 Length\t$n50_len\n";
$per=$total_length/$count;
print OUT0 "Mean Length\t$per\n";
close OUT0;

open OUT2,">$od/$Index.Unigene.stat.info.xls" || die $!;
print OUT2 "Length_span\tNumbers\tPercent\n";
foreach (0..30) {
	my $site=$_;
	my $min=$site*100;
	my $max=($site+1)*100;
	my $len_span="$min"."~"."$max";
	my $percent;
	if ($site==30 && defined $length_region_dis{30}) {
		$len_span=">3000";
		$percent=($length_region_dis{30}/$count)*100;
		printf OUT2 "%s\t%s\t\t%.2f\n",$len_span,$length_region_dis{30},$percent;
	}
	elsif (!defined $length_region_dis{$site}) {
		print OUT2 "$len_span\t0\t0\n";
	}
	else {
		$percent=($length_region_dis{$site}/$count)*100;
		printf OUT2 "%s\t%s\t%.2f\n",$len_span,$length_region_dis{$site},$percent;
	}
}
close OUT2;


#####################
#Length Distribution Graph
###############################

open OUT,">$od/$Index.Unigene.distribution.svg" || die $!;
my $svg=&Length_Distribution_svg(\%length_region_dis,$Index);
print OUT "$svg";
close OUT;
chdir $od;
&runOrDie("perl $svg2xxx $Index.Unigene.distribution.svg");  #svg-->png

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);


####################subs
sub Length_Distribution_svg {#
	my ($len_dis,$lable) = @_;
	my $svg;
	my %length_dis=%$len_dis;
	my $max_value=(sort {$a <=> $b} values %length_dis)[-1];
	my $y_max=int(&log10($max_value))+1;
	my $width=720;
	my $height=480;
	my $left_pos=100;my $right_pos=$width-35;
	my $low_pos=$height-100;my $high_pos=50;
	my $yx_coordinate=$left_pos-5;my $yy_coordinate=$low_pos;
	my $xx_coordinate=$left_pos;my $xy_coordinate=$low_pos+5;
	my $x_iner=3;my $x_span=9;
	my $rect_color = "rgb(0,165,210)";

	###############�������� ������̶� ������
	$svg.=&svg_paper($width,$height);
	$svg.=&svg_line($left_pos,$low_pos,$left_pos,$high_pos,"#000000");
	$svg.=&svg_txt_middle($left_pos/2,$height/2,16,"#000000","$lable Unigene Number",3);
	$svg.=&svg_line($left_pos,$low_pos,$right_pos,$low_pos,"#000000");
	$svg.=&svg_txt_middle($width/2,$height-40,16,"#000000","Length (nt)");
	$svg.=&svg_txt_middle($width/2,$high_pos/2+10,18,"#000000","$lable Unigene Length Distribution");

	my $y_step=($low_pos-$high_pos)/$y_max;
	my $x_step=($right_pos-$left_pos)/31;
	for (0..$y_max) {
		if ($_==0) {
			$svg.=&svg_line($yx_coordinate,$low_pos-($y_step*$_),$xx_coordinate,$low_pos-($y_step*$_),"#000000");
			$svg.=&svg_txt_end($xx_coordinate-10,$low_pos-($y_step*$_)+5,10,"#000000",0);
		}
        elsif ($_<=5) {
			$svg.=&svg_line($yx_coordinate,$low_pos-($y_step*$_),$xx_coordinate,$low_pos-($y_step*$_),"#000000");
			$svg.=&svg_txt_end($xx_coordinate-10,$low_pos-($y_step*$_)+5,10,"#000000",10**$_);
		}
		else {
			$svg.=&svg_line($yx_coordinate,$low_pos-($y_step*$_),$xx_coordinate,$low_pos-($y_step*$_),"#000000");
			$svg.=&svg_txt_end($xx_coordinate-10,$low_pos-($y_step*$_)+5,10,"#000000","1e+$_");
		}
	}
	for (1..30) {
		$svg.=&svg_line($left_pos+($x_step*$_),$yy_coordinate,$left_pos+($x_step*$_),$xy_coordinate,"#000000");
        my $lab = 100*($_-1)."-".100*$_;
		$svg.=&svg_txt_end($left_pos+($x_step*$_)+5-0.5*$x_step,$xy_coordinate+13-0.5*$x_step,8,"#000000",$lab,4);
#		$svg.=&svg_txt_middle($left_pos+($x_step*$_)+5,$xy_coordinate+13,10,"#000000",$lab,4);
	}
	$svg.=&svg_txt_end($right_pos-0.4*$x_step,$xy_coordinate+15-0.5*$x_step,8,"#000000",">3000",4);

	##############�����ȷֲ�����
	for my $site (0..30) {
		next if (!defined $length_dis{$site}) ;
		my $x=$left_pos+($site)*$x_step+$x_iner;
		my $log=&log10($length_dis{$site});
		my $y=($log/$y_max)*($low_pos-$high_pos);
		$svg.=&svg_rect_event($x,$low_pos-$y,$x_span,$y,$rect_color,$rect_color,1);
	}
	$svg.=&svg_end();
	return ($svg);
}

sub svg_txt (){
	#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);						# д�ַ��ĵ���x���꣬y���꣬�ַ���С����ɫ���ַ��ĵ�������
	#$anchor ��xyȷ���ĵ�Ϊ׼��0��text�Ե���ʼ��1���Ե�Ϊ���ģ�2���Ե������
	#$svg_x[6] ��ת��0:0��;1:90��;2:180��;3:270��;
	my @svg_x=@_;
	if (!defined $svg_x[6]) {
		$svg_x[6]=0;
	}
	my $svg_matrix='';
	if ($svg_x[6]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[6]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[6]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[6]==3) {
		$svg_matrix="0 -1 1 0";
	}
	if (!defined $svg_x[5] || $svg_x[5] == 0) {
		my $line="<text fill=\"$svg_x[3]\"  transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"Arial\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
		return $line;
	}else{
		my $anchor="";
		if ($svg_x[5]==1) {
			$anchor="middle";
		}
		if ($svg_x[5]==2) {
			$anchor="end";
		}
		my $line="<text fill=\"$svg_x[3]\" text-anchor=\"$anchor\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"Arial\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
		return $line;
	}
}

sub svg_txt_middle (){#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\" text-anchor=\"middle\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_txt_end (){#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
    if ($svg_x[5]==4) {
		$svg_matrix="1 -0.75 0.75 1";
    }
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\" text-anchor=\"end\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_line (){#&svg_line(x1,y1,x2,y2,color,[width])
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[5]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_paper (){#&svg_paper(width,height,[color])
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\" viewBox=\"0 0 $svg_x[0] $svg_x[1]\">\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub svg_rect_event () {#&svg_rect(x,y,width,height,fill_color,stroke_color,strock_width);
	my @svg_x=@_;
	my $line ;
	$line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6];\" />\n";
	if (defined $svg_x[7]){
		$line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6];fill-opacity:$svg_x[7]\" />\n";
	}
	if (defined $svg_x[8]){
		$line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6];fill-opacity:$svg_x[7]\" onclick=\"alert('$svg_x[8]')\" onmousemove=\"window.status='$svg_x[8]'\" />\n";
	}
	if (defined $svg_x[9]){
		$line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" style=\"fill:$svg_x[4];stroke:$svg_x[5];stroke-width:$svg_x[6];fill-opacity:$svg_x[7]\" onmouseover=\"changeText(evt,'$svg_x[8]','$svg_x[9]')\" onmouseout=\"changeTextNotOver(evt,'$svg_x[8]')\" onclick=\"changeClick(evt)\" onmousemove=\"window.stsatus='$svg_x[8]'\"  attributeName=\"fill-opacity\" from=\"0\" to=\"0.4\" begin=\"mouseover\" end=\"mouseout\"/> \n";
	}
	return $line;
}

sub svg_end (){#
	return "</svg>\n";
}

sub log10 {#
	my ($n) = @_;
	return log($n)/log(10);
}

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
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
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}
sub USAGE {#
    my $usage=<<"USAGE";
     Program: $0
     Version: $version
     Contact: Meng Fei <mengf\@biomarker.com.cn>
Program Date: 2012-07-02
      Modify: Simon Young <simonyoung8824\@gmail.com> 
 Modify Date: 2014-08-08
 Description: 
        Usage:
          -i       <str>    Index of inFiles and outFiles                             must be given;
          -fa      <str>    Infile                                                    must be given;
          -len     <int>    max length of transcripts                                 optional, default 30,000;
          -s       <str>    Length Region for Fasta file Stat(300,500,1000,2000)      optional;
          -od      <str>    OUT file DIR                                              must be given;
		  -syscfg  <file>   system config for the third party softwares                         optional

USAGE
	print $usage;
	exit;
}