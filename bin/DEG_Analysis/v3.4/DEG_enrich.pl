#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Algorithm::Combinatorics qw(combinations permutations);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($id,$enrich,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$od,
                                "enrichment:s"=>\$enrich
				) or &USAGE;
&USAGE unless ($enrich and $od);

mkdir $od unless -d $od;
mkdir "$od/work_sh" unless -d "$od/work_sh";
$od=&ABSOLUTE_DIR($od);
my $notename=`hostname`;chomp $notename;

my $Rscript =${&readconf("$Bin/../../../config/sys.cfg")}{Rscript};


open (ANNO,">$od/work_sh/enrich.sh") or die $!;

foreach my $file (glob "$od/*_vs_*/*DEG*final.xls") {
	my $dir=dirname $file;
	my $group_name=basename $dir;
        mkdir "$dir/Anno_enrichment" unless ( -d "$dir/Anno_enrichment");
	print ANNO "perl $Bin/bin/Enrichment/draw_DEG_GO_graph.pl -i $enrich -deg $dir/$group_name.DEG_final.xls -k $group_name -od $dir/Anno_enrichment\n";
	print ANNO "perl $Bin/bin/Enrichment/draw_top_GO_graph.pl -i $enrich -deg $dir/$group_name.DEG_final.xls -k $group_name -od $dir/Anno_enrichment/\n";
	print ANNO "perl $Bin/bin/Enrichment/draw_DEG_KEGG_graph.pl -i $enrich -deg $dir/$group_name.DEG_final.xls -k $group_name -od $dir/Anno_enrichment\n";
	print ANNO "perl $Bin/bin/Enrichment/draw_DEG_COG_graph.pl -i $enrich -deg $dir/$group_name.DEG_final.xls -k $group_name -od $dir/Anno_enrichment\n\n";
        print ANNO "perl $Bin/bin/Enrichment/draw_DEG_KOG_graph.pl -i $enrich -deg $dir/$group_name.DEG_final.xls -k $group_name -od $dir/Anno_enrichment\n\n";
        print ANNO "perl $Bin/bin/Enrichment/draw_DEG_eggNOG_graph.pl -i $enrich -deg $dir/$group_name.DEG_final.xls -k $group_name -od $dir/Anno_enrichment\n\n";
	print ANNO "perl $Bin/bin/Enrichment/select_DEG_Anno.pl -i $enrich -deg $dir/$group_name.DEG_final.xls -o $dir/Anno_enrichment/$group_name.annotation.xls\n";
	}
close ANNO;


my @deg_zero=glob("$od/*/*.DEG_final.xls.zero");
if ($#deg_zero>=0) {
        for(my $i=0;$i<=$#deg_zero;$i++){
                my $filename=basename $deg_zero[$i];
                my ($sample)=$filename=~/^(\S+)\.DEG\_final\.xls\.zero/;
                system("sed -i /$sample/d $od/work_sh/enrich.sh");
        }
}
&timeLog("Enriching Annotation:\n$od/work_sh/enrich.sh");
&runOrDie("$od/work_sh/enrich.sh");


 

#######################################################################################
&timeLog("$Script Done");
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
		warn "Warning just for file and dir $in\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub max{#&max(lists or arry);
	#锟斤拷锟叫憋拷锟叫碉拷锟斤拷锟街?	
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
	#锟斤拷锟叫憋拷锟叫碉拷锟斤拷小值
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
	#锟斤拷取锟街凤拷锟斤拷锟斤拷锟叫的凤拷锟津互诧拷锟斤拷锟叫ｏ拷锟斤拷锟街凤拷锟斤拷锟斤拷式锟斤拷锟截★拷ATTCCC->GGGAAT
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

sub Shell_run
{
	my $sh = shift;
	my $t=GetTime();
	print "\[$t\] $sh starts:\n";
	unless(-e "$sh.finish"){
		my $cmd_res=system"sh $sh >$sh.log 2>&1";
		if($cmd_res){
			die"Err:$sh run wrong,$!\n";
		}
		else{
			system"touch $sh.finish";
		}
	}
	my $e=GetTime();
	print "\[$e\] $sh ends.\n";
}



sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Program Date:   2013.10.18
Usage:
  Options:

  -enrichment    All_Database_annotation.xls                 choice

  -od            Out dir                                     must be given;

  -h                help document

USAGE
	print $usage;
	exit;
}

