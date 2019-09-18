#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2012
# Writer:         Zhangyh <zhangyh@biomarker.com.cn>
# Program Date:   2012/8/20 10:00
# Modifier:       Zhangyh <zhangyh@biomarker.com.cn>
# Last Modified:  2012/9/13 18:19
# Modifier:       Xugl <xugl@biomarker.com.cn>
# Last Modified:  2015/05/25 16:00
# Modifier:       baij <baij@biomarker.com.cn>
# Last Modified:  2016/12/14 11:00

my $ver="1.3.0";

use GD;
use GD::Polyline;
use Cwd;
use strict;
use Getopt::Long;
use FindBin qw/$Bin $Script/;
use Data::Dumper;
use File::Basename qw(basename dirname);
use newPerlBase;

my %opts;
GetOptions(\%opts,
	"i=s",
	"o=s",
	"d=s",
	"v",
	"h",
	"map=s",#system dir
);
&help() if ( defined($opts{h}) || !defined($opts{i}) ||
			!defined($opts{o}) || !defined($opts{d}) );

my $BEGIN=time();
my $Time_Start = sub_format_datetime(localtime(time()));
&timeLog("$Script Start.");

###### program begin#############
if (!-d $opts{o}) {
	mkdir $opts{o} || die "Can't creat $opts{o},$!\n";
}
#/share/nas1/baij/TestData/Result/Maize_v2.2.2/Uni_Anno/Result

my $kegg_map = "$opts{o}/Kegg_map";#Kegg_map
$opts{i} = ABSOLUTE_DIR ($opts{i});
$opts{d} = ABSOLUTE_DIR ($opts{d});
mkdir $kegg_map unless -d $kegg_map;



my $path=(glob "$opts{i}/*.Kegg.pathway")[0];
my $ko = (glob "$opts{i}/*.ko")[0];


if ( !defined($path) ) {
		die "\n\tPlease check input files!\n\n";
	}

################# build reference ###############
my %uni;
$\=">";
open (UNI,"$opts{d}") || die "Can't open $opts{d},$!\n";
while (<UNI>) {
	chomp;
	next if /^\s*$/ or /^\#/;
	if(/^>/){
	my $gene=(split/\n/,$_)[0];
        $gene =~s/^>//;
	$uni{$gene}=1;
	}
}


#################################################
#################### kegg map  ##################
#################################################


#��������������  relate-to-gene  ����������������
	my %png_new;
	open (LI,"$path") || die "Can't open $path,$!\n";
	while (<LI>) {
		next if $.<2;
		chomp;
		my @tmp_1 = split /\t/;
		my @tmp_2 = split /;/,$tmp_1[3];
		my @tmp_3 = split /\+/,$tmp_1[4];
		for my $i (0..$#tmp_2) {
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{gname}.="$tmp_2[$i];";
		}
		for my $i (0..$#tmp_2) {
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{gname} =~ s/;$//;
		}
	}
	close (LI);

	#����������������  modify png  ������������������
	#����������������& creat html  ������������������
	my @lost;
	for my $koxx (sort keys %png_new) {
		if (-e "$opts{map}/$koxx.png") {
			my @K_ID = keys %{$png_new{$koxx}};
#			print "@K_ID\n",'=' x 10,"\n"; #tag
			my %html;
			open (PNG,"$opts{map}/$koxx.png") || die "Can' open $opts{map}/$koxx.png,$!\n";
			open (CONF,"$opts{map}/$koxx.conf") || die "Can't open $opts{map}/$koxx.conf,$!\n";
			open (RESULT,">$kegg_map/$koxx.png") || die "Can't creat $kegg_map/$koxx.png,$!\n";
			open (HTML,">$kegg_map/$koxx.html") || die "Can't creat $kegg_map/$koxx.html,$!\n";

			my $im = newFromPng GD::Image(\*PNG,1) ||die;		#truecolorģʽ
			my $red = $im->colorAllocate(255,0,0);
			my $green = $im->colorAllocate(0,255,0);
			my $blue = $im->colorAllocate(0,0,255);
			my $black = $im->colorAllocate(0,0,0);
                        my $purple = $im->colorAllocate(171,130,225);
			while (my $rect=<CONF>) {
				my ($type,$p1,$p2,$p3,$p4) = $rect =~ /^(\w+)\t\((\d+),(\d+),(\d+),(\d+)\)/;
				foreach my $K (@K_ID) {
					if ($rect =~ /$K/){
						if($type=~/rect/){
                                                my $poly=new GD::Polygon;
                                                $poly->addPt($p1,$p2);$poly->addPt($p1,$p4);$poly->addPt($p3,$p4);$poly->addPt($p3,$p2);
                                                $im->openPolygon($poly,$black);
                                                $im->fillToBorder($p1+1,$p2+1,$black,$purple);
                                                for my $x($p1+1..$p3-1){                    ##look for the close area and fill color
                                                        for my $y ($p2+1..$p4-1){
                                                                my $colorIndex=$im->getPixel($x,$y);
                                                                my ($r,$g,$b)=$im->rgb($colorIndex);
                                                                if($r==255 and $g ==255 and $b==255){
                                                                        $im->fillToBorder($x,$y,$black,$purple);
                                                                }
                                                        }
                                                }
                                                        $html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{color}="AB82FF";
							$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{info}.="<ul><li>$K $png_new{$koxx}{$K}{gname} </li></ul>";
							$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{KO}.="+$K";
							$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{html_coord}="$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2";
						}
						elsif($type=~/line/){
							my @temp=split/\t/,$rect;
							for(my $i=1;$i<=$#temp-1;$i++){
								my ($p1,$p2,$p3,$p4) = $temp[$i] =~ /\((\d+),(\d+),(\d+),(\d+)\)/;
								my ($x1,$x2,$x3,$x4,$y1,$y2,$y3,$y4);
								if($p2==$p4){
									$y1=$p2-1;$y2=$p2+1;$y3=$p4-1;$y4=$p4+1;
									$html{$koxx}{"$p1,$p2,$p3,$p4"}{png_coord}="$p1,$p2,$p3,$p4,$p3,$y1,$p1,$y1,$p1,$p2";  ##png_coord for bold the png line 
									$html{$koxx}{"$p1,$p2,$p3,$p4"}{html_coord}="$p1,$y1,$p1,$y2,$p3,$y4,$p3,$y3,$p1,$y1";   ###html_coord for view the href of html
								}
								else{
									$x1=$p1-1;$x2=$p1+1;$x3=$p3-1;$x4=$p3+1;
									$html{$koxx}{"$p1,$p2,$p3,$p4"}{png_coord}="$p1,$p2,$p3,$p4,$x3,$p4,$x1,$p2,$p1,$p2";
									$html{$koxx}{"$p1,$p2,$p3,$p4"}{html_coord}="$x1,$p2,$x2,$p2,$x4,$p4,$x3,$p4,$x1,$p2";

								}
								$html{$koxx}{"$p1,$p2,$p3,$p4"}{info}.="<ul><li>$K $png_new{$koxx}{$K}{gname} </li></ul>";
								$html{$koxx}{"$p1,$p2,$p3,$p4"}{KO}.="+$K";
                                                         my $poly=new GD::Polygon;
                                                        my @new_coord=split/,/,$html{$koxx}{"$p1,$p2,$p3,$p4"}{png_coord};
                                                        for (my $i=0;$i<=$#new_coord ;$i+=2) {
                                                          $poly->addPt($new_coord[$i],$new_coord[$i+1]);
                                                        }
                                                        $im->openPolygon($poly,$purple);
                                                        $html{$koxx}{"$p1,$p2,$p3,$p4"}{color}="AB82FF";
							}
						}
					}
				}
				
             }
			##########################
			######	HTML BEGIN	######
			##########################
			print HTML <<____________HTML;
<html>
<meta http-equiv="content-type" content="text/html; charset=utf-8">
<head>
<title>
$koxx
</title>
<style type="text/css">
<!--

area {cursor: pointer;}

-->
</style>
<link rel="stylesheet" href="/css/kegg.css" type="text/css" />
<script language="JavaScript" src="/js/dhtml.js"></script>
<script type="text/javascript">
<!---

function showInfo(info) {
	obj = document.getElementById("result");
	obj.innerHTML = "<div style='cursor: pointer; position: absolute; right: 5px; color: #000;' onclick='javascript: document.getElementById(\\"result\\").style.display = \\"none\\";' title='close'>X</div>" + info;
	obj.style.top = document.body.scrollTop;
	obj.style.left = document.body.scrollLeft;
	obj.style.display = "";
}

//--->
</script>

</head>
<body>
<img src="$koxx.png" usemap="#mapdata" border="0" />
<map name="mapdata">
____________HTML
			foreach my $coords ( keys %{$html{$koxx}} ) {
				my $color = $html{$koxx}{$coords}{color};
				my $info = $html{$koxx}{$coords}{info};
				print "$koxx  $coords\n" and exit if (!defined $info);
				print HTML <<____________HTML;
<area shape='poly' coords='$html{$koxx}{$coords}{html_coord}' href='http://www.kegg.jp/dbget-bin/www_bget?$html{$koxx}{$coords}{KO}' target="_blank" onmouseover='javascript: showInfo("<ul><li style=\\"color: #$color;\\">$info</li></ul>");' />
____________HTML
				
			}

			print HTML <<____________HTML;
</map>
<div id='result' style='position: absolute; width: 50%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;' onmouseover="javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;" onmouseout="javascript: this.style.filter = 'alpha(opacity=95)'; this.style.opacity = 0.95;"></div>
</body>
</html>
____________HTML
			##########################
			######	HTML END	######
			##########################
			close (HTML);
			close (PNG);
			close (CONF);
			binmode RESULT;
			print RESULT $im->png;
			close (RESULT);
			print $koxx,"\n" if (defined $opts{v});
			print "job done!\n" if (defined $opts{v});
			print '=' x 40,"\n" if (defined $opts{v});
		} else {
			print "$koxx map not found in kegg datebase!\n" if (defined $opts{v});
			print '=' x 40,"\n" if (defined $opts{v});
			push @lost,"$koxx\n";
		}
		
	}
	print "[notice]:These maps are not in datebase\n",@lost;

#��������������  program  end  ������������������

my $Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
&timeLog("$Script End.");

#������������������������������������������������
#                  Sub Routines                ��
#������������������������������������������������


sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "\nTotal elapsed time: ${t}s\n";
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
sub help {#usage
	print <<"	Usage End.";
	Description:
		KEGG GO enrichment & maps & web pages
		Version: $ver

	Usage:
		
		-d    <str>      DEG file
		-i    <str>      input dir
		-o    <str>      output dir
		-v               [option]   view process on screen
	
	Example��

		perl kegg_map_web.pl -d Maize.Unigene.fa  -i input -o output 

	[notice]input dir should at least contain:
			
		xxx.fa.Kegg.path                        
		xxx.fa.Kegg.ko

	Usage End.

	exit;
}
