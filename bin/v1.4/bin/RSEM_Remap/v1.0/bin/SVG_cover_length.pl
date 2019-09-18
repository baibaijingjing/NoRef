#!/usr/local/bin/perl -w
# Copyright (c) BMK 2011
# Writer:         Guoxd <guoxd@biomarker.com.cn>
# Program Date:   2012/01/04
# Modifier:       Guoxd <guoxd@biomarker.com.cn>
# Last Modified:  2012/01/04
my $ver="1.0";

use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $programe_dir=basename($0);
my $path=dirname($0);

#############################################参数格式

my %opts;
GetOptions(\%opts,"i=s","o=s","h" );

if(! defined($opts{i}) ||! defined($opts{o}) || defined($opts{h})){
	&USAGE;
}

###############Time
&timeLog("$programe_dir Start.");
################

################perl geneExp_diff_svg.pl -in Cabbage_Leaf.geneExpression.xls -o myout
my $in=$opts{i};
my $out=$opts{o};
open (IN,"$in")||die "$!";
open (OUT,">".$out.".Reads2Gene.svg") || die "$!";

my $max_x=0;
my $max_y=0;
my $min_x=10000;
my $min_y=10000;
while(<IN>){
	chomp;
	s/^\s+//;s/\s+$//;s/\r+$//;
	next if (/^\#/ || /^$/) ;
	next if $.==1 ;
	my ($len,$reads);
	$len=(split/\s+/,$_)[3];
	$reads=(split/\s+/,$_)[4];
	if($reads > $max_x){$max_x=$reads;}
	if($len > $max_y){$max_y=$len;}
	next if ($reads==0);
	if($reads < $min_x){$min_x=$reads;}
	if($len < $min_y){$min_y=$len;}
}
close IN;
$max_x=&log_base(10,$max_x);
$max_y=&log_base(10,$max_y);
$min_x=&log_base(10,$min_x);
$min_y=&log_base(10,$min_y);
my $countx=int($max_x)+1;
my $county=int($max_y)+1;
print "max_x\t$countx\tmax_y\t$county\n";
print "min_x\t$min_x\tmin_y\t$min_y\n";
my $width=100*$countx+100;
my $height=100*$county+100;
my $left_pos=50;my $low_pos=$height-50;my $high_pos=50;my $right_pos=$width-50;
my $yx_coordinate=$left_pos-5;my $yy_coordinate=$low_pos;
my $xx_coordinate=$left_pos;my $xy_coordinate=$low_pos+13;
###############画坐标轴 坐标轴刻度 网格线
print OUT &svg_paper($width,$height);
print OUT &svg_line($left_pos,$low_pos,$left_pos,$high_pos,"#000000");
print OUT &svg_txt($left_pos/2,$height/2+100,18,"#000000","Unigene length (log10)",3);
print OUT &svg_line($left_pos,$low_pos,$right_pos,$low_pos,"#000000");
print OUT &svg_txt($width/2-100,$low_pos+30,18,"#000000","Reads Number (log10)");
my $j=0;
my $i;
for ($i=2;$i<=$county+1;$i+=0.5) {
	my $y_plus=$j*100-5;
	my $y=$low_pos-$y_plus;
	print OUT &svg_txt($yx_coordinate-10,$y,10,"#000000","$i");
	$j++;
}
for ($i=1;$i<=$county;$i++) {
	print OUT &svg_line($left_pos,$low_pos-100*$i,$left_pos+5,$low_pos-100*$i,"#000000");
}

$j=1;
for ($i=1;$i<=$countx+1;$i++) {
	my $x_plus=$j*100-5;
	my $x=$left_pos+$x_plus;
	print OUT &svg_txt($x,$xy_coordinate,10,"#000000","$i");
	$j++;
}
for ($i=1;$i<=$countx;$i++) {
	print OUT &svg_line($left_pos+100*$i,$low_pos,$left_pos+100*$i,$low_pos-5,"#000000");
}

################画点
open (IN,"$in")||die "$!";
while(<IN>){
	next if $.==1;
	chomp;
	s/^\s+//;s/\s+$//;s/\r+$//;
	next if (/^\#/ || /^$/) ;
	my ($len,$reads);
	$len=(split/\s+/,$_)[3];
	$reads=(split/\s+/,$_)[4];
	next if ($reads eq "0.00");
	my $pos_x=&log_base(10,$reads);
	next if ($pos_x eq "0.00");
    my $pos_y=&log_base(10,$len)-2;
	print OUT &svg_circle2($left_pos+$pos_x*100,$low_pos-$pos_y*200,3,"none","#0000FF",0.5);
}
close IN;

###################
print OUT &svg_end();
close(OUT);

my $svg2xxx=${&readconf("$Bin/../../../../../../../config/sys.cfg")}{svg2xxx};
system "$svg2xxx $out.Reads2Gene.svg";

#############subs
sub svg_txt (){#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
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
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
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

sub svg_circle2 () {#&svg_circle(x,y,r,fill-color,stroke-color,stroke-width)
         my @svg_x=@_;
         my $line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" fill=\"$svg_x[3]\" />\n";
         return $line;
}


sub svg_circle () {#&svg_circle(x,y,r,color,[info])
	my @svg_x=@_;
	my $line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" fill=\"$svg_x[3]\" />\n";
	if (defined $svg_x[4]) {
		$line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" fill=\"$svg_x[3]\" onclick=\"alert('$svg_x[4]')\" onmousemove=\"window.status='$svg_x[4]'\" />\n";
	}
	return $line;
}

sub svg_end (){#
	return "</svg>\n";
}

sub log_base {
    my ($base, $value) = @_;
    return log($value)/log($base);
}

###############Time
&timeLog("$programe_dir End.");
&Runtime($BEGIN);

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total $programe_dir elapsed time: ${t}s\n\n";
}
sub USAGE{
	print <<"	Usage End.";
	Description:
		Version: $ver
		Writer by Guoxd <guoxd\@biomarker.com.cn>
	Usage:
		-i          infile (geneExpression)         must be given

		-o          outfile prefix                  must be given

		-h          Help document
	Usage End.
	exit;
}
