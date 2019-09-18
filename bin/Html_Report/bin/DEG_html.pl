#!/bin/env perl
use strict;
use autodie;
use File::Basename qw(basename);
use Getopt::Long;
use FindBin qw($Bin);
my ($inputdir,$help);
GetOptions(
	"id=s"=>\$inputdir,
	"h|help"=>\$help,
);
mkdir "$inputdir/HTML" unless(-d "$inputdir/HTML");
my $outdir=$inputdir."/HTML";
foreach my $g (glob "$inputdir/DEG_Analysis/*_vs_*") {
my $group=basename($g);
####DEG####
open DEG,"$g/$group.DEG_final.xls";
my @deg;

my $i=0;
while(<DEG>){
	chomp;
	next if (m/^\s*$/);
	my @tmp=split/\s+/;
	map{push @{$deg[$i]},$_} @tmp;
	$i++;
}
close DEG;

open HTML,">$outdir/$group.DEG.html";
print HTML <<XGL;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=gb2312" />
<meta NAME="Author" CONTENT="xugl\@biomarker.com.cn" /> 
<meta NAME="Version" CONTENT="2015821v1.0" /> 
<title>DEG</title>
<link href="src/base.css" type="text/css" rel="stylesheet">
</head>
<body>
<div>
<p><a name="home"><img class="normal" src="src/logo.png" /></a></p>
</div><br />
<h1>Different Expression Unigenes Analysis</h1>
<table class="gy">
<tr>
XGL
for (my $i=0;$i<=$#{$deg[0]} ;$i++) {
	print HTML '<th>'.$deg[0][$i].'</th>'."\n";
}
print HTML '</tr>'."\n";
for (my $i=1;$i<@deg ;$i++) {
	print HTML "\t".'<tr>'."\n\t";
	for (my $j=0;$j<@{$deg[$i]} ;$j++) {
		print HTML '<td>'.$deg[$i][$j].'</td>';
	}
	print HTML "\t".'</tr>'."\n";
}
print HTML <<XGL;
</table>
</body>
</html>
XGL
close HTML;
#######DEG_KEGG#######
my @deg_kegg;
open KEGG,"$g/pathway/kegg_enrichment/$group.KEGG.stat";
my (%KEGG);
while(<KEGG>){
	chomp;
	next if(m/^\#/ or m/^\s*$/);
	my @tmp=split/\t+/;
	$KEGG{$tmp[0]}{ko}=$tmp[1];
	$KEGG{$tmp[0]}{"P-value"}=$tmp[4];
	$KEGG{$tmp[0]}{"Corrected_P-value"}=$tmp[5];
}
close KEGG;
open KEGG,"$g/pathway/kegg_enrichment/$group.KEGG.xls";
while(<KEGG>){
	chomp;
	next if (m/^\#/ or m/^\s*$/);
	my @tmp=split/\t+/;
	$tmp[6]=~s/;/ ;/g;
	$tmp[7]=~s/\+/ ;/g;
	$KEGG{$tmp[0]}{deg}=$tmp[2];
	$KEGG{$tmp[0]}{gene}=$tmp[3];
	$KEGG{$tmp[0]}{geneName}=$tmp[6];
	$KEGG{$tmp[0]}{KO}=$tmp[7];
}
close KEGG;
open HTML,">$outdir/$group.DEG_KEGG.html";
print HTML <<XGL;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=gb2312" />
<meta NAME="Author" CONTENT="xugl\@biomarker.com.cn" /> 
<meta NAME="Version" CONTENT="2015821v1.0" /> 
<title>KEGG Enrichment Analysis</title>
<link href="src/base.css" type="text/css" rel="stylesheet">
</head>

<body>
<div>
<p><a name="home"><img class="normal" src="src/logo.png" /></a></p>
</div><br />
<h1>Different Expression Unigene KEGG Enrichment Analysis</h1>
<table class="gy">
<tr>
<th>Pahtway</th>
<th>Different expression Unigene number</th>
<th>Background Unigene number</th>
<th>P-value</th>
<th>Corrected P-value</th>
<th>Gene_id</th>
<th>KEGG Orthology</th>
</tr>
XGL
for my $pathway (sort {$KEGG{$a}{"Corrected_P-value"} <=> $KEGG{$b}{"Corrected_P-value"}} keys %KEGG) {
	print HTML "\t".'<tr>'."\n";
	print HTML "<td><a href=../DEG_Analysis/$group/pathway/kegg_map/".$KEGG{$pathway}{ko}.'.html target=_blank>'.$pathway.'</a></td><td>'.$KEGG{$pathway}{deg}.'</td><td>'.$KEGG{$pathway}{gene}.'</td><td>'.$KEGG{$pathway}{"P-value"}.'</td><td>'.$KEGG{$pathway}{"Corrected_P-value"}.'</td><td>'.$KEGG{$pathway}{geneName}.'</td><td>'.$KEGG{$pathway}{KO}.'</td>'."\n";
	print HTML "\t</tr>\n";
}
print HTML <<XGL;
</table>
</body>
</html>
XGL
close HTML;
######GO_Enrichment######
my %GO;
my @statFile=glob("$g/go_enrichment/$group.GO\.*.stat");
foreach my $file (@statFile) {
	open GO,"$file";
	while(<GO>){
	next if (m/^\#/ or m/^\s*$/);
	my @tmp=split/\t/;
	my ($root,$go)=split/:/,$tmp[0],2;
	$GO{$root}{$go}{"P-value"}=$tmp[3];
	$GO{$root}{$go}{"Corrected_P-value"}=$tmp[4];
}
close GO;
}
my @xlsFile=glob("$g/go_enrichment/$group.GO\.*.xls");
foreach my $file (@xlsFile) {
	open GO,"$file";
	while(<GO>){
	next if (m/^\#/ or m/^\s*$/);
	my @tmp=split/\t/;
	my ($root,$go)=split/:/,$tmp[0],2;
	$GO{$root}{$go}{deg}=$tmp[1];
	$GO{$root}{$go}{gene}=$tmp[2];
	$tmp[5]=~s/;/ ;/g;
	$GO{$root}{$go}{geneName}=$tmp[5];
}
close GO;
}

open HTML,">$outdir/$group.DEG_GO.html";
print HTML <<XGL;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=gb2312" />
<meta NAME="Author" CONTENT="xugl\@biomarker.com.cn" /> 
<meta NAME="Version" CONTENT="2015821v1.0" /> 
<title>GO Enrichment Analysis</title>
<link href="src/base.css" type="text/css" rel="stylesheet">
</head>

<body>
<div>
<p><a name="home"><img class="normal" src="src/logo.png" /></a></p>
</div><br />
<h1>Different Expression Unigene GO Enrichment Analysis</h1>
<table class="gy">
<tr>
<th>GO Root Term</th>
<th>GO Term</th>
<th>Different expression Unigene number</th>
<th>Background Unigene number</th>
<th>P-value</th>
<th>Corrected P-value</th>
<th>Gene_id</th>
</tr>
XGL
for my $root (keys %GO) {
	foreach my $go (sort {$GO{$root}{$a}{"Corrected_P-value"} <=> $KEGG{$root}{$b}{"Corrected_P-value"}} keys %{$GO{$root}}){
		print HTML "\t".'<tr>'."\n";
		print HTML '<td>'.$root.'</td><td>'.$go.'</td><td>'.$GO{$root}{$go}{deg}.'</td><td>'.$GO{$root}{$go}{gene}.'</td><td>'.$GO{$root}{$go}{"P-value"}.'</td><td>'.$GO{$root}{$go}{"Corrected_P-value"}.'</td><td>'.$GO{$root}{$go}{geneName}.'</td>'."\n";
		print HTML "\t</tr>\n";
	}
}
print HTML <<XGL;
</table>
</body>
</html>
XGL
close HTML;
}

open DEG,"$inputdir/DEG_Analysis/All_DEG/All.DEG_final.xls";

open HTML,">$outdir/All.DEG.html";
print HTML <<XGL;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=gb2312" />
<meta NAME="Author" CONTENT="xugl\@biomarker.com.cn" /> 
<meta NAME="Version" CONTENT="2015821v1.0" /> 
<title>DEG</title>
<link href="src/base.css" type="text/css" rel="stylesheet">
</head>

<body>
<div>
<p><a name="home"><img class="normal" src="src/logo.png" /></a></p>
</div><br />
<h1>Different Expression Unigenes Analysis</h1>
<table class="gy">
XGL
while(<DEG>){
	chomp;
	next if (m/^\s*$/);
	my @tmp=split/\s+/;
	if($.==1){
		print HTML "<tr>\n";
		for (my $i=0;$i<=$#tmp;$i++) {
			print HTML '<th>'.$tmp[$i].'</th>'."\n";
		}
		print HTML "</tr>\n";
	}
	else{
		print HTML "<tr>\n";
		for (my $i=0;$i<=$#tmp;$i++) {
			print HTML '<td>'.$tmp[$i].'</td>'."\n";
		}
		print HTML "</tr>\n";
	}
}
print HTML <<XGL;
</table>
</body>
</html>
XGL




