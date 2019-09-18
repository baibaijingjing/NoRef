#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0";

my($dir);
GetOptions(
	"help|?"	=>\&USAGE,
	"dir:s"		=>\$dir,
) or &USAGE;
&USAGE unless (-d $dir);


################################################################################################################
#main
################################################################################################################

#读入聚类数高度信息，每个cell的长宽。treeheight_row，treeheight_col，cellwidth，cellheight
my %coords;
my @coords = glob "$dir/*.coordinate";
if (@coords == 1) {
	open COORDS,"$coords[0]";
	while (<COORDS>) {
		chomp;
		next if /^\s*$/;
		my @line = split "\t";
		$line[1] = $line[1] * 6.94571399553 * 0.14; #图片点位的单位不同，需要把pt转换为px
		$coords{$line[0]} = $line[1];
	}
	close COORDS;
}
else {
	die "The FILE does not EXISTS: $dir/*.coordinate";
}

#从输入目录中寻找必要文件。
my $heatmap_file;
my %heatmap;
my %go;
my ($x1,$y1,$x2,$y2) = ($coords{'treeheight_row'},0,0,$coords{'treeheight_col'});
foreach my $data_png (glob "$dir/*.clustered.data_subcluster_*.data.png") {
	my ($num) = $data_png =~ /.*\.clustered\.data_subcluster_(\d+)\.data\.png$/;
	
	my $basename = basename($data_png);
	$heatmap{$num}{'data_png'} = $basename;
	
	$heatmap_file = $basename;
	$heatmap_file =~ s/\.clustered\.data_subcluster_$num\.data\.png$//;
	$heatmap_file .= "\.png";
	die "The FILE does not EXISTS: $dir/$heatmap_file" unless (-f "$dir/$heatmap_file");
}
#my ($pic_width) = (`file $dir/$heatmap_file` =~ /\, (\d+) x \d+\,/);
my ($pic_width) = (`file $dir/$heatmap_file|cut -f 2 -d","` =~ /(\d+) x \d+/);
print $pic_width;

foreach my $num (sort {$a<=>$b} keys %heatmap) {
	my $basename = $heatmap{$num}{'data_png'};
	$basename =~ s/\.data\.png$//;
	
	
	my $cluster_num = (split " ",`wc -l $dir/$basename.data`)[0] - 1;
	my $sample_num = (`head -n 1 $dir/$basename.data` =~ tr/\t/\t/) + 1;
#	$x2 = $x1 + $sample_num * $coords{'cellwidth'};
	$x2 = $pic_width;
	$y1 = $y2;
	$y2 = $y1 + $cluster_num * $coords{'cellheight'};
	$heatmap{$num}{'coords'} = "0,$y1,$x2,$y2";
	
	if (-f "$dir/$basename.png") {
		$heatmap{$num}{'png'} = "$basename.png";
	}
	else {
		die "The FILE does not EXISTS: $dir/$basename.png";
	}
	
	my ($go_list) = glob "$dir/$basename.*.go_enrichment.list";
	if (defined $go_list && -f $go_list) {
		open IN,"$go_list";
		while (<IN>) {
			chomp;
			next if (/^\s*$/ or $. == 1);
			my ($term,$pvalue) = split "\t",$_;
			$pvalue = 0 if ($pvalue == 0);
			$go{$num} .=<<"GO_COORDS";
                    <tr>
                        <td>$term</td>
                        <td>$pvalue</td>
                    </tr>
GO_COORDS
		}
		close IN;
		$go{$num} ||= <<"NO_GO_RESULT";
                    <tr>
                        <td>There in no result of GO enrichment.</td>
                        <td>There in no result of GO enrichment.</td>
                    </tr>
NO_GO_RESULT
	}
	else {
		die "The FILE does not EXISTS: $dir/$basename.*.go_enrichment.list";
	}
}

#html代码
my $html = <<"MYHTML";
<!DOCTYPE html>
<html>
<head lang="zh-cn">
    <meta http-equiv="content-type" content="text/html; charset=UTF-8">
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="Keywords" content="百迈客云、云科技、大数据、百迈客云生物科技、生物云、生物信息云平台、生信云、生信云平台">
    <title>Heatmap数据挖掘</title>

    <link href="css/bootstrap.min.css" rel="stylesheet">
    <link href="css/font-awesome/css/font-awesome.css" rel="stylesheet">

    <link href="css/animate.css" rel="stylesheet">
    <link href="css/style.css" rel="stylesheet">

    <style>
        body{
            min-width: 900px;
            /*overflow: auto !important;*/
            /*margin: 0 auto;*/
        }
        td {
            word-wrap: break-word;
            word-break: break-all;
            /*padding:2px 3px !important;*/
        }

        .map-content-right tr td:first-child{
            max-width: 285px;
        }
        .map-content-right td{
            font-weight:700;

        }
        .map-head{
            margin-bottom: 30px;
            margin-right: 20px;
        }
        .map-head>span{
            margin-right: 20px;
            font-weight: 700;
            font-size: 20px;
        }
        .map-head img{
            margin-top: -10px;
        }
        .map-content{
            padding: 20px 40px 0;
            /*min-width: 1250px;*/
            overflow: hidden;
        }
        .map-content-right{
            display: none;
            margin-left: 30px;
            width: 400px;
            max-height:1500px;
        }
        .map-content-right img{
            display: block;
            width: 100%;
        }
        .map-content-right div{
            width:100%;
        }
        .map-content-left{
            /*margin-left:200px;*/
        }
        #area_img_mask{
            display:none;
            position:absolute;
            background-color:#000;
            opacity: .2;
        }
        table{
            height:100%;
            max-height: 600px;;
        }
        table tr{
            max-width: 450px;;
        }
    </style>
</head>
<Script>
    function reurl(){

        var url = location.href; //把当前页面的地址赋给变量 url

        var subStr=url.substr(url.length-3,3);
        if(subStr != "??1"){ //如果?后的值不等于1表示没有刷新
            url += "??1"; //把变量 url 的值加入 ?1

            self.location.replace(url); //刷新页面
        }
    }
</script>
<body class="white-bg">

<div class="map-content">
    <div class="row text-center map-head">
        <span>Heatmap数据挖掘</span>
        <img src="logo_zh.png" alt="">
    </div>
    
    <div class="map-content-div">
    
        <div class="map-content-left pull-left">
            <div>
                <img src="$heatmap_file" usemap="#mapMap" id="sourceImg">
            </div>
            <map  name="mapMap" id="mapMap">
MYHTML

foreach my $num (sort {$a<=>$b} keys %heatmap) {
	$html .=<<"HTMLCOORDS";
                <area class="problem_area"
                      coords="$heatmap{$num}{'coords'}"
                      href="#$num" >
HTMLCOORDS
}



$html .= "            </map>\n            <div id=\"area_img_mask\"></div>\n        </div>\n";


foreach my $num (sort {$a<=>$b} keys %heatmap) {
	$html .= <<"HTML_FILE_AND_TABLE";
        <div class="map-content-right text-center pull-left" id="${num}map-result">
            <div>
                <img src="$heatmap{$num}{'png'}" alt="" class="img-container">
            </div>
            <div>
                <img src="$heatmap{$num}{'data_png'}" alt="" class="img-container">
            </div>
            <div class="text-left">
                <table class="table table-bordered">
                    <tr>
                        <td>GO Term</td>
                        <td>P-value</td>
                    </tr>
$go{$num}
                </table>
            </div>
        </div>
HTML_FILE_AND_TABLE
}


$html .=<<"MAINLY";
    </div>
</div>

<!-- Mainly scripts -->
<script src="js/jquery-2.1.1.js"></script>
<script src="js/bootstrap.min.js"></script>

<script>
    \$(document).ready(function(){
        \$("body").width(\$(window).innerWidth()*0.8+'px');
        \$(".map-content").width(\$(window).innerWidth());
        var imgHeight=\$(".map-content-left>div>img").height();
        var imgWidth=\$(".map-content-left>div>img").width();
        \$(".map-content-left").css("width",imgWidth);
        \$(".map-content-div").css("margin-left",(\$(window).innerWidth()-imgWidth-490)*0.5+'px');
        var url=\$("#sourceImg").attr("src");
	function loadImage(url, callback) {
	var img = new Image(); //创建一个Image对象，实现图片的预下载
	img.src = url;

	if(img.complete) { // 如果图片已经存在于浏览器缓存，直接调用回调函数
		callback.call(img);
		return; // 直接返回，不用再处理onload事件
	}
	img.onload = function () { //图片下载完毕时异步调用callback函数。
		callback.call(img);//将回调函数的this替换为Image对象
	};
	};
	loadImage(url,heatMap);
        function heatMap(){
            showBorder(\$("map").find("area:first"));
            \$("map").find("area").each(function(){
                \$(this).click(function(){
                    showBorder(\$(this));
                })
            });

            \$(".map-content-right").children("div").not(":last-child").css("max-height",imgHeight*0.35+'px');
            \$(".map-content-right>div>img").css("max-height",imgHeight*0.35+'px');
            \$(".map-content-right").find("div").find('table').css("max-height",imgHeight*0.26+'px');
            \$("#1map-result").show();

            \$(".map-content").css("min-width","1000px");
            \$(".map-content-right").css("min-height",imgHeight+50+'px');
        }
        \$(window).resize(function(){
            \$(".map-content-right").each(function(){
                \$(".map-content-div").css("margin-left",(\$(window).innerWidth()-imgWidth-490)*0.5+'px');
                if(\$(window).innerWidth()<(460+imgWidth)){
                    \$(".map-content-div").css("margin-left",(460+imgWidth-\$(window).innerWidth())*0.5+'px');
                }
                if(\$(this).is(":visible")){
                    var num=parseInt(\$(this).attr("id"));
                    \$("area").each(function(){
                        if(\$(this).attr("href").slice(1)==num){
                            showBorder(\$(this));
                        }
                    })
                }
            })
        });
        \$("map").find("area").each(function(){
            \$(this).click(function () {
                var href=\$(this).attr("href");
                console.log(href);
                \$(".map-content-right").hide();
                \$(href+"map-result").css("display","block");
                showBorder(\$(this));
                console.log(\$(this).attr("coords"));
            })
        });
        function getOffset(obj){
            var x = obj.offsetLeft, y = obj.offsetTop;
            while(obj.offsetParent){
                obj = obj.offsetParent;
                x += obj.offsetLeft;
                y += obj.offsetTop;
            }
            return {x : x, y : y};
        }
        function showBorder(obj){
            var img = document.getElementById("sourceImg");
            var div = document.getElementById("area_img_mask");
            var offset = getOffset(img);
            var coords = obj.attr("coords").split(",");

            div.style.display = "block";
            div.style.left = offset.x + parseInt(coords[0]) + "px";
            div.style.top = offset.y + parseInt(coords[1]) + "px";
            div.style.width = parseInt(coords[2]) - parseInt(coords[0]) + "px";
            div.style.height = parseInt(coords[3]) - parseInt(coords[1]) + "px";
        }
        \$(window).resize(function(){
            var wwww = \$(window).width();
            if(wwww<=1200){
                \$("body").css({
                    "overflow":"auto"
                });
            }
        })
    })



</script>

</body>
</html>
MAINLY


`cp -r $Bin/css $Bin/js $Bin/logo_zh.png $dir`;
open OUT,">$dir/heatmap.html";
print OUT "$html";
close OUT;



################################################################################################################
#sub functions
################################################################################################################

sub GetTime {
        my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
        return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

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
                warn "Warning just for file and dir $in\n";
                exit;
        }
        chdir $cur_dir;
        return $return;
}

################################################################################################################

sub cmd_call {
        print "@_\n";
        system(@_) == 0 or die "system @_ failed: $?";
}

################################################################################################################

sub LOG_FORMAT {
        my $info = shift;
        my $star = '*'x80;
        my $time = &GetTime;
        my $error = $star."\n$time\n$info\n".$star."\n";
        return $error;
}

################################################################################################################

sub ERROR_AND_DIE {
        my $info = shift;
        my $error = &LOG_FORMAT($info);
        die "$error";
}

################################################################################################################

sub MAKE_DIR {
        my $directory = shift;
        if (-f $directory) {
                &ERROR_AND_DIE("$directory is a file!");
        }
        elsif (-d $directory) {
#               &ERROR_AND_DIE("Directory $directory exists!");
        }
        else {
                &cmd_call("mkdir -p $directory");
        }
        $directory = &ABSOLUTE_DIR($directory);
        return $directory;
}

################################################################################################################

sub check_array_no_same {#&check_array_no_same("array_name",@array_name); 检查数组中是否有相同的元素，有则die。
	my $array_name = shift;
	my @array = @_;
	foreach my $i (0..$#array-1) {
		foreach my $j ($i+1..$#array) {
			if ($array[$i] eq $array[$j]) {
				print "ERROR: \"$array[$i]\" appears twice in \@$array_name at least. Please check your input.\n";
				die "Illegal parameter input: -$array_name $array[$i]\n";
#				print "$i,$j:\tSame\n";
			}
		}
	}
}

################################################################################################################

sub USAGE {
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Usage:
  Options:
	-dir		input dir

	-h		Help

USAGE
	print $usage;
	exit;
}
