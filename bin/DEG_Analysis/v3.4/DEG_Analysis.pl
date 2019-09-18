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
my ($id,$cfg,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$id,
				"cfg:s"=>\$cfg,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($id and $cfg and $od );

mkdir $od unless -d $od;
mkdir "$od/work_sh" unless -d "$od/work_sh";
$id=&ABSOLUTE_DIR($id);
$cfg=&ABSOLUTE_DIR($cfg);
$od=&ABSOLUTE_DIR($od);
my $notename=`hostname`;chomp $notename;

my $Rscript =${&readconf("$Bin/../../../config/sys.cfg")}{Rscript};
############FDR FC

&timeLog("Calculate All Gene FDR and FC\nperl $Bin/bin/gene_counts_to_FDR_FC/gene_counts_to_FDR_FC.pl  -i $id -cfg $cfg -od $od");
system "perl $Bin/bin/gene_counts_to_FDR_FC/gene_counts_to_FDR_FC.pl  -i $id -cfg $cfg -od $od/ ";

############
my ($FC,$FDR,%deg);
open (IN,$cfg) or die $!;
while (<IN>) {
	chomp;
	if (/^fold/) {
		$FC=(split/\s+/,$_)[1];
		if ($FC<=0) {
			print "Fold Change Value Should Greater than zero!";die;
		}
	}
	if (/^FDR/) {
		$FDR=(split/\s+/,$_)[1];
	}
	if(/^Com/){
		my @tmp=split/,/,$_;
		for(my $i=0;$i<$#tmp;$i++){
			for(my $j=$i+1;$j<=$#tmp;$j++){
				$deg{$tmp[$i]}=$tmp[$j];
			}
		}
	}
	elsif(/^Sep/){
		my @tmp=split/;/,$_;
		for(my $i=0;$i<$#tmp;$i++){
			for(my $j=$i+1;$j<=$#tmp;$j++){
				$deg{$tmp[$i]}=$tmp[$j];
			}
		}
	}
}
close IN;

my $de_files;
my %groups;

open (SH1,">$od/work_sh/select.sh") or die $!;
open (SH2,">$od/work_sh/draw.sh") or die $!;


foreach my $dir (glob "$od/*_vs_*") {
	next unless -d $dir;
	my $group_name=basename $dir;
	my $cmd="perl $Bin/bin/filter_by_value/v1.0/filter_by_value.pl -i $dir/$group_name.final.xls -o $dir/$group_name.DEG_final.xls -f $od/All_gene_expression.list ";
#	$de_files.="$dir/$group_name.DEG_final.xls ";
#        if (defined $FDR) {
#		$cmd.="-FDR $FDR "
#	}
#	if (defined $FC) {
#		$cmd.="-FC $FC "
#	}
#	print SH1 "$cmd\n";
 #       &runOrDie("$cmd");      
	my ($control,$treated)=split/_vs_/,$group_name;
	$groups{$control}=1;
	$groups{$treated}=1;
	my @sample1=split/_/,$control;
	my @sample2=split/_/,$treated;
	open GRP,">$dir/$group_name.group";
	map {print GRP "$_\tgroup1\n"} @sample1;
	map {print GRP "$_\tgroup2\n"} @sample2;
	close GRP;
	if(-e "$dir/$group_name.DEG_final.xls" ){
        my $DEG_line =(split/\s/,`wc -l $dir/$group_name.DEG_final.xls`)[0];
         if($DEG_line > 2){
        $de_files.="$dir/$group_name.DEG_final.xls ";
        #########
	print SH2 "perl $Bin/bin/draw_DEG_cor_plot/v1.0/draw_DE_cor_plot.pl -treated $treated -control $control -all $dir/$group_name.final.xls -de $dir/$group_name.DEG_final.xls -o $dir/$group_name.DEG_cor.png \n";
	#########
	print SH2 "perl $Bin/bin/draw_FC_count_plot/v1.0/draw_FC_count_plot.pl -treated $treated -control $control -all $dir/$group_name.final.xls -de $dir/$group_name.DEG_final.xls -o $dir/$group_name.FC_count.png \n";
	#########
	print SH2 "perl $Bin/bin/draw_FC_FDR_plot/v1.0/draw_FC_FDR_plot.pl -all $dir/$group_name.final.xls -de $dir/$group_name.DEG_final.xls -o $dir/$group_name.FC_FDR.png -th $FC,$FDR\n";
	#########
	#print SH2 "perl /share/nas2/genome/bmksoft/tool/cluster_draw/v2.0/cluster_draw.pl -fn 3 -i $dir/$group_name.DEG_final.xls -od $dir/DEG_Cluster \n";
	
	print SH2 "$Rscript $Bin/bin/heatmap.R --infile $dir/$group_name.DEG_final.xls.forheatmap --outfile $dir/$group_name.DEG.cluster --groupfile $dir/$group_name.group -H 5000 -W 3000 --rowname F\n";

	}else{

            #########
                        print SH2 "perl $Bin/bin/draw_FC_count_plot/v1.0/draw_FC_count_plot.pl -treated $treated -control $control -all $dir/$group_name.final.xls -de $dir/$group_name.DEG_final.xls -o $dir/$group_name.FC_count.png \n";
            #########
                        print SH2 "perl $Bin/bin/draw_FC_FDR_plot/v1.0/draw_FC_FDR_plot.pl -all $dir/$group_name.final.xls -de $dir/$group_name.DEG_final.xls -o $dir/$group_name.FC_FDR.png -th $FC,$FDR\n";
            #########
#                        print SH2  " `touch $dir/$group_name.DEG_final.xls.zero`";
                        print SH2  " `rm $dir/$group_name.DEG_final.xls` "; 
                                                                                      }

}
 }
close SH1;close SH2;
#&timeLog("Selecting DEG:\n$od/work_sh/select.sh");
#&runOrDie("$od/work_sh/select.sh");

&timeLog("Drawing Graphs:\n$od/work_sh/draw.sh");
&runOrDie("$od/work_sh/draw.sh");



my @DEG_file = glob "$od/*_vs_*/*DEG*final.xls";
if(@DEG_file){
############veen 
&timeLog("draw_diff_group_veen\nperl $Bin/bin/draw_diff_group_veen/v1.0/draw_diff_group_veen.pl -id $od -od $od ");
system("perl $Bin/bin/draw_diff_group_veen/v1.0/draw_diff_group_veen.pl -id $od -od $od ");

&timeLog("draw_diff_group_veen\nperl $Bin/bin/draw_diff_group_veen/v2.0/draw_diff_group_veen.pl -id $od -od $od/Venn ");
system("perl $Bin/bin/draw_diff_group_veen/v2.0/draw_diff_group_veen.pl -id $od -od $od/Venn");


############
mkdir "$od/All_DEG" unless -d "$od/All_DEG";
&timeLog("Get All DEG\nperl $Bin/bin/filter_by_names/v1.0/filter_by_names.pl -i $od/All_gene_fpkm.list -s $de_files -o $od/All_DEG/All.DEG_final.xls");
system("perl $Bin/bin/filter_by_names/v1.0/filter_by_names.pl -i $od/All_gene_fpkm.list -s $de_files -o $od/All_DEG/All.DEG_final.xls");

###########

&timeLog("$Rscript $Bin/bin/heatmap.R  --infile $od/All_DEG/All.DEG_final.xls --outfile $od/All_DEG/All.DEG.cluster  -H 5000 -W 3000 --rowname F");
system "$Rscript $Bin/bin/heatmap.R --infile $od/All_DEG/All.DEG_final.xls --outfile $od/All_DEG/All.DEG.cluster -H 5000 -W 3000 --rowname F";
######################增加热图的网页形式2016-08-25

my $DEG_file="$od/All_DEG/All.DEG_final.xls ";
my  $sample_num=`head -n 1 $DEG_file|awk \'{print NF}\'`;
#print "sss$sample_num\n";
mkdir "$od/All_DEG/HTML" unless -d "$od/All_DEG/HTML";
my $cmd_sh="perl $Bin/bin/heatmap_html//heatmap_draw.pl  -infile $od/All_DEG/All.DEG_final.xls  --outfile  $od/All_DEG/HTML/All.DEG.cluster  -cluster row -color_type 1 -is_log  -div -anno   $od/../Uni_Anno/Result/  ";
$cmd_sh.=" -scale none" if($sample_num<4);
$cmd_sh.=" -scale row" if ($sample_num>=4) ;
&runOrDie("$cmd_sh");
}
#print "Cluster_draw\nperl /share/nas2/genome/bmksoft/tool/cluster_draw/v2.0/cluster_draw.pl -i $od/All_DEG/All.DEG_final.xls -od $od/All_DEG/DEG_Cluster \n\n";
#`perl /share/nas2/genome/bmksoft/tool/cluster_draw/v2.0/cluster_draw.pl -i $od/All_DEG/All.DEG_final.xls -od $od/All_DEG/DEG_Cluster `;



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
  -i             All_gene_counts.list                       must be given;
  
  -cfg           soft parameter to DE mRNA Analysis         must be given;

  -od            Out dir                                     must be given;

  -h                help document

USAGE
	print $usage;
	exit;
}
