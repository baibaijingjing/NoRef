#!/usr/bin/perl -w
use strict;
use warnings;
use Config::General;
use Encode;
use Cwd qw(abs_path getcwd);
use Getopt::Long;
use Data::Dumper;
use Encode qw(decode);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Algorithm::Combinatorics qw(combinations permutations);
use Spreadsheet::WriteExcel;
use Spreadsheet::ParseExcel;
use newPerlBase;
#use Array::Diff;
my $BEGIN_TIME=time();
my $version="1.0.0";

#######################################################################################
sub USAGE {#
	my $usage=<<"_USAGE_";
#-------------------------------------------------------------------------------------------------------------------
    Program: $Script
    Version: $version
    Contact: 
Description:

     Usage:
       Options:
            --id             <dir>  Web_report dir, required
            --od             <dir>  output dir,required
            --cleandata        <dir>  unforced , cleandata fastq file or dir,use ":" to separate files or dir; 
                                    example:  -cleandata T1.fq:T2.fq:T3.fq:dir_1:dir_2
            --inf            <str>  project basic information configration file.                [\$Bin/pro_inf.cfg]
            --oneSam                only one sample, no DEG and SNP Analysis
            --config_path     <str> path config file

            -h                      help
Example:
perl $Script --id ./ --od ./ --data_cfg data.cfg --detail_cfg detail.cfg
perl $Script --id ./ --od ./ --config_path ../Config/
#-------------------------------------------------------------------------------------------------------------------
_USAGE_
	print $usage;
	exit;
}
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($id,$cleandata_dir,$oneSample,$pro_inf,$od,$xls_inf,$data_cfg,$detail_cfg,$config_path,$noconfig);


GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$id,
                "inf:s"=>\$pro_inf,
				"cleandata:s"=>\$cleandata_dir,
				"oneSam"=>\$oneSample,
				"od:s"=>\$od,
				"config_path:s"=>\$config_path,
				) or &USAGE;
&USAGE unless ($id and $od);
if(-d $od){}else{mkdir $od;}

$od||=getcwd;
$od=abs_path($od);
$id=abs_path($id);
my $detail;

if($config_path){
	$config_path=abs_path($config_path);
    $data_cfg = (glob "$config_path/*data.cfg")[0];

	$detail=&load_cfg($config_path,$od);
	&load_data_cfg($data_cfg,"$od/data.cfg") if ($data_cfg);
}else{
	die "Please specify -data_cfg and -detail_cfg or -config_path";
}

#my$Only_One_Sam;




#$pro_inf ||= "$Bin/pro_inf.cfg";
#my %pro_inf;
#&LOAD_INF($pro_inf,\%pro_inf);
#print Dumper($detail);
my %pro_inf;
if (defined($pro_inf)){
	%pro_inf=&LOAD_INF1($pro_inf);
}else{
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	my $date=sprintf("%4d-%02d-%02d", $year+1900, $mon+1, $day);
	%pro_inf=(
		'project_name' =>(exists $$detail{projectInfo}{Project_name})?$$detail{projectInfo}{Project_name}:"default",
		'project_id'=> (exists $$detail{projectInfo}{Contract_NO})?$$detail{projectInfo}{Contract_NO}:"default",
		'user_name'=> (exists $$detail{projectInfo}{Customer_name})?$$detail{projectInfo}{Customer_name}:"default",
		'user_unit'=>(exists $$detail{projectInfo}{Customer_info})?$$detail{projectInfo}{Customer_info}:"default", 
		'pro_sam_receive_date'=>(exists $$detail{projectInfo}{First_time})?$$detail{projectInfo}{First_time}:"default",
		'pro_sam_pass_date'=>(exists $$detail{projectInfo}{Second_time})?$$detail{projectInfo}{Second_time}:"default",
		'pro_launch_date'=>(exists $$detail{projectInfo}{Third_time})?$$detail{projectInfo}{Third_time}:"default",
		'pro_finish_date'=>$date,
	);
}

#print Dumper(\%pro_inf);


mkdir "$od/biocloud/" unless (-d "$od/biocloud/");
my $cloud_dir = "$od/biocloud";
my $abstract;
&timeLog("$Script start.");
#######################################################################################
if (defined $cleandata_dir) {
	###############################################获取绝对路径
	$cleandata_dir = abs_path($cleandata_dir);

	###############################################原始数据链接
	mkdir "$cloud_dir/display_cleandata" unless (-d "$cloud_dir/display_cleandata");
	my @cleandata_dir = split ":",$cleandata_dir;
	foreach my $raw_dir (@cleandata_dir) {
		if (-d $raw_dir) {
			foreach (glob "$raw_dir/*.fq") {
				my $basename = basename($_);
				`ln -s $_ $cloud_dir/display_cleandata` unless (-f "$cloud_dir/display_cleandata/$basename");
			}
		}
		elsif (-f $cleandata_dir) {
			my $basename = basename($cleandata_dir);
			`ln -s $cleandata_dir $cloud_dir/display_cleandata` unless (-f "$cloud_dir/display_cleandata/$basename");
		}
	}
}



###############################################################图片预览文件夹
mkdir "$cloud_dir/images" unless (-d "$cloud_dir/images");


###############################################样品、分组、物种名
my $T_num=0;
my @samples = glob "$id/geneExpression/*.geneExpression.xls";
my $samples_num = @samples;
for (my $i=0;$i<$samples_num;$i++) {
	$samples[$i] =~ s/.*geneExpression\/(.*)\.geneExpression\.xls/$1/;
}
my $bio_rep = 0;
my @dir;
foreach my $dir (glob "$id/DEG_Analysis/*_vs_*") {
	$dir=basename $dir;
	push @dir,$dir;
	my $count_ = $dir =~ s/\_/\_/g; 
	$bio_rep = 1 if ($count_ > 2);
}

my $deg;
my ($treated,$control);
if (@dir) {
	($treated,$control)=(split/_vs_/,$dir[0])[0,1];
	if ($treated =~ /\_/) {
		$deg .= "分组$treated";
	}
	else {
		$deg .= "样品$treated";
	}
	$deg .= "和";
	if ($control =~ /\_/) {
		$deg .= "分组${control}间";
	}
	else {
		$deg .= "样品${control}间";
	}
}
my$Old=0;
if(-d "$id/Unigene/Unigene_CDS/"){
	$Old=1;
}elsif(-d "$id/Unigene/Unigene_Orf/"){
	$Old=2;
}else{
	die "Can't find dir '$id/Unigene/Unigene_CDS/' or '$id/Unigene/Unigene_Orf/'";

}
my$name=$$detail{projectInfo}{Project};
my@name;
if($Old==1){
	#@name = glob "$id/Unigene/Unigene_CDS/*.Unigene.cds.distribution.png";#Maize.Unigene.cds.distribution.png
	@name = glob "$id/Unigene/Unigene_CDS/*.distribution.png";#Maize.Unigene.cds.distribution.png

}elsif($Old==2){
#Old version
	@name = glob "$id/Unigene/Unigene_Orf/Graph/*.orf.distribution.png";
}else{
	
}

########################################################################示例文件
mkdir "$cloud_dir/template" unless (-d "$cloud_dir/template");
my @cp = ("Trinity_workflow.png","base_erro_ratio.txt","RNA-Seq_experimental_workflow.png","RNA-Seq_analysis_workflow.png","cleandata_FASTQ_format.png","FPKM_formula.png","enrichment_factor_formula.png","P26_pp_network.png");
foreach (@cp) {
`cp $Bin/template/$_ $cloud_dir/template`;
}





###############################################样品测序数据评估统计表
open ALLSAMP,"$id/cleandata/AllSample_GC_Q.stat" or die "$!:$id/cleandata/AllSample_GC_Q.stat";
#open ALLSAMPNEW,">$id/cleandata/AllSample_GC_Q_New.stat" or die "$!:$id/cleandata/AllSample_GC_Q_New.stat";
open ALLSAMPNEW,">$od/cleandata/AllSample_GC_Q_New.stat" or die "$!:$od/cleandata/AllSample_GC_Q_New.stat";

print ALLSAMPNEW "Samples\tBMK-ID\tRead Number\tBase Number\tGC Content\t",'%≥Q30',"\n";
#print ALLSAMPNEW "Samples\tBMK-ID\tTotal reads\tTotal bases\tGC%\tQ30%\n";
my $clean_data = 0;
my $clean_Q30 = 100;
my @samp_all = @samples;
while (<ALLSAMP>) {
	next if (/^#/ || /^\s*$/ || /^SampleID/);
	chomp;
	my @allsampnew = split /\s+/;
	my $len = 0;
	my $samp ="";
	foreach (@samp_all) {
		if ($allsampnew[0] =~ /$_/) {
			if (length $_ > $len) {
				$len = length $_;
				$samp = $_;
			}
		}
	}
	for (my $i=0;$i<@samp_all;$i++) {
		if ($samp_all[$i] eq $samp) {
			splice @samp_all,$i,1;
			$clean_data += $allsampnew[2];
			$clean_Q30 = $allsampnew[7] if ($allsampnew[7] < $clean_Q30);
		}
	}
	print ALLSAMPNEW "$$detail{projectInfo}{$samp}\t$samp\t$allsampnew[1]\t$allsampnew[2]\t$allsampnew[3]%\t$allsampnew[7]%\n" if ($samp);
	last unless (@samp_all);
}
close ALLSAMP;
close ALLSAMPNEW;
$clean_data = sprintf "%.2f",$clean_data/1000000000;




###############################################组装结果统计表
my %assembly;
my $assembly_num =0;
my $unigene_1kb;
# -------------------------------------------------
if (-d "$id/Assembly/All_Combination") { #合并组装
# -------------------------------------------------
	open ASSEMBLYC1,"$id/Assembly/All_Combination/contigs/All_Combination.contigs.stat.info.xls";
	while (<ASSEMBLYC1>) {
		next if (/^\D+/);
		chomp;
		$assembly_num = 200 if (/^200\~/);
		$assembly_num = 300 if (/^300\~/);
		$assembly_num = 500 if (/^500\~/);
		$assembly_num = 1000 if (/^1000\~/);
		$assembly_num = 2000 if (/^2000\~/);
		$assembly{'contigs'}{$assembly_num} += (split /\t/)[1];
		$assembly{'contigs'}{'total'} += (split /\t/)[1];
	}
	close ASSEMBLYC1;
	foreach my $num (keys %{$assembly{contigs}}) {
		next unless $num=~/^\d+$/;
		&assem_per($num);
	}
#	&assem_per(200);
#	&assem_per(300);
#	&assem_per(500);
#	&assem_per(1000);
#	&assem_per(2000);
	open ASSEMBLYC2,"$id/Assembly/All_Combination/contigs/All_Combination.contigs.stat.xls";
	while (<ASSEMBLYC2>) {
		chomp;
		$assembly{'contigs'}{'T_n'} = (split /\t/)[1] if (/^Total number/i);
		$assembly{'contigs'}{'T_l'} = (split /\t/)[1] if (/^Total length/i);
		$assembly{'contigs'}{'N50'} = (split /\t/)[1] if (/^N50/);
		$assembly{'contigs'}{'Mean'} = (split /\t/)[1] if (/^Mean/);
		$assembly{'contigs'}{'Mean'} = sprintf "%.2f",$assembly{'contigs'}{'Mean'} if (/^Mean/);
	}
	open ASSEMBLYT,"$id/Assembly/All_Combination/Transcripts/All_Combination.Transcripts.stat.xls";
	while (<ASSEMBLYT>) {
		chomp;
		&assem_t_u('Trans',$_);
	}
	close ASSEMBLYT;
	open ASSEMBLYU,"$id/Assembly/All_Combination/Unigenes/All_Combination.Unigenes.stat.xls";
	while (<ASSEMBLYU>) {
		chomp;
		&assem_t_u('Uni',$_);
	}
	close ASSEMBLYU;
	$unigene_1kb += (split '\(',$assembly{'Uni'}{1000})[0];
	$unigene_1kb += (split '\(',$assembly{'Uni'}{2000})[0];

	open ASSEMOUT,">$id/Assembly/All_Combination.stat.xls";
	print ASSEMOUT "Length range\tContigs\tTranscripts\tUnigenes\n";
	print ASSEMOUT "200\-300\t$assembly{'contigs'}{200}\t$assembly{'Trans'}{200}\t$assembly{'Uni'}{200}\n" if defined $assembly{contigs}{200};
	print ASSEMOUT "300\-500\t$assembly{'contigs'}{300}\t$assembly{'Trans'}{300}\t$assembly{'Uni'}{300}\n" if defined $assembly{contigs}{300};
	print ASSEMOUT "500\-1000\t$assembly{'contigs'}{500}\t$assembly{'Trans'}{500}\t$assembly{'Uni'}{500}\n";
	print ASSEMOUT "1000\-2000\t$assembly{'contigs'}{1000}\t$assembly{'Trans'}{1000}\t$assembly{'Uni'}{1000}\n";
	print ASSEMOUT "2000\+\t$assembly{'contigs'}{2000}\t$assembly{'Trans'}{2000}\t$assembly{'Uni'}{2000}\n";
	print ASSEMOUT "Total number\t$assembly{'contigs'}{'T_n'}\t$assembly{'Trans'}{'T_n'}\t$assembly{'Uni'}{'T_n'}\n";
	print ASSEMOUT "Total length\t$assembly{'contigs'}{'T_l'}\t$assembly{'Trans'}{'T_l'}\t$assembly{'Uni'}{'T_l'}\n";
	print ASSEMOUT "N50 length\t$assembly{'contigs'}{'N50'}\t$assembly{'Trans'}{'N50'}\t$assembly{'Uni'}{'N50'}\n";
	print ASSEMOUT "Mean length\t$assembly{'contigs'}{'Mean'}\t$assembly{'Trans'}{'Mean'}\t$assembly{'Uni'}{'Mean'}\n";
	close ASSEMOUT;
# ------------------------------------------------------------
} #合并组装
# ------------------------------------------------------------
else {  #分开组装
	my ($Uni_file1) = glob "$id/Unigene/*.stat.xls";
	open ASSEMBLYUA,"$Uni_file1";
	<ASSEMBLYUA>;
	while (<ASSEMBLYUA>) {
		chomp;
		&assem_t_u('Uni',$_);
	}
	close ASSEMBLYUA;

	my @Uni_group_file = glob "$id/Assembly/*/Unigenes/*.Unigenes.stat.xls";

	my $group_name_1 = basename($Uni_group_file[0]);
	$group_name_1 =~ s/\.Unigenes\.stat\.xls//;

	my $group_name_2 = basename($Uni_group_file[1]);
	$group_name_2 =~ s/\.Unigenes\.stat\.xls//;

	open ASSEMBLYU1,"$Uni_group_file[0]";
	<ASSEMBLYU1>;
	while (<ASSEMBLYU1>) {
		chomp;
		&assem_t_u('Uni_1',$_);
	}
	close ASSEMBLYU1;

	open ASSEMBLYU2,"$Uni_group_file[1]";
	<ASSEMBLYU2>;
	while (<ASSEMBLYU2>) {
		chomp;
		&assem_t_u('Uni_2',$_);
	}
	close ASSEMBLYU2;

	$unigene_1kb += (split '\(',$assembly{'Uni'}{1000})[0];
	$unigene_1kb += (split '\(',$assembly{'Uni'}{2000})[0];

	open ASSEMOUT,">$id/Assembly/All_Combination.stat.xls";
	
	print ASSEMOUT "Length range\t${group_name_1}_Unigenes\t${group_name_2}_Unigenes\tAll_Unigenes\n";
	print ASSEMOUT "200\-300\t$assembly{'Uni_1'}{200}\t$assembly{'Uni_2'}{200}\t$assembly{'Uni'}{200}\n" if defined $assembly{Uni_1}{200};
	print ASSEMOUT "300\-500\t$assembly{'Uni_1'}{300}\t$assembly{'Uni_2'}{300}\t$assembly{'Uni'}{300}\n" if defined $assembly{Uni_1}{300};
	print ASSEMOUT "500\-1000\t$assembly{'Uni_1'}{500}\t$assembly{'Uni_2'}{500}\t$assembly{'Uni'}{500}\n";
	print ASSEMOUT "1000\-2000\t$assembly{'Uni_1'}{1000}\t$assembly{'Uni_2'}{1000}\t$assembly{'Uni'}{1000}\n";
	print ASSEMOUT "2000\+\t$assembly{'Uni_1'}{2000}\t$assembly{'Uni_2'}{2000}\t$assembly{'Uni'}{2000}\n";
	print ASSEMOUT "Total number\t$assembly{'Uni_1'}{'T_n'}\t$assembly{'Uni_2'}{'T_n'}\t$assembly{'Uni'}{'T_n'}\n";
	print ASSEMOUT "Total length\t$assembly{'Uni_1'}{'T_l'}\t$assembly{'Uni_2'}{'T_l'}\t$assembly{'Uni'}{'T_l'}\n";
	print ASSEMOUT "N50 length\t$assembly{'Uni_1'}{'N50'}\t$assembly{'Uni_2'}{'N50'}\t$assembly{'Uni'}{'N50'}\n";
	print ASSEMOUT "Mean length\t$assembly{'Uni_1'}{'Mean'}\t$assembly{'Uni_2'}{'Mean'}\t$assembly{'Uni'}{'Mean'}\n";
	close ASSEMOUT;

}

###############################################################相关性系数表
my $de_limit=0;
my $de=(glob "$id/DEG_Analysis/*_vs_*")[0];
if ($de && -d $de) {
#	($de1,$de2)=$de=~/Diff_Analysis\/([^\/]+)_vs_([^\/]+)/;
	$de_limit=1;
}
if (-d "$id/DEG_Analysis/All_DEG") {
	$de_limit=2;
}

my %table_info;
my $ftable;

push @{$table_info{"生物学重复相关性统计"}{"table"}},["Sample 1","Sample 2",'R^2'];
if ($de_limit!=0) {
	$ftable = "$id/geneExpression/free_com.cor";
	open (TAB,$ftable) or die $!;
	$/="\n";
	my @Sam;
	my %COR;
	foreach my $key6 (split/\s+/,(scalar <TAB>)) {
		if ($key6=~/decor\./) {
			my $name6=$key6;
			$name6=~s/decor\.//;
			push @Sam,$name6;
		}
	}
	my $limit=@Sam;
	while (<TAB>) {
		chomp;
		next if(/^$/);
		my @A = split /\s+/,$_;
		for (my $i6=1;$i6<=$limit ;$i6++) {
			next if $A[0] eq $Sam[$i6-1];
			next if exists $COR{$Sam[$i6-1]}{$A[0]};
			$COR{$A[0]}{$Sam[$i6-1]}=$A[$i6];
		}
	}
	foreach my $key1 (keys %COR) {
		foreach my $key2 (keys %{$COR{$key1}}) {
			push @{$table_info{"生物学重复相关性统计"}{"table"}},[$key1,$key2,$COR{$key1}{$key2}];
		}
	}
	close (TAB) ;
}

open SWXCF,">$id/geneExpression/free_com.stat";
foreach my $a (@{$table_info{"生物学重复相关性统计"}{"table"}}) {
	foreach  (@{$a}) {
		print SWXCF $_,"\t";
	}
	print SWXCF "\n";
}


################################################比对结果统计表


my $mapstat = "#Sample\tTotal reads\tTotal mapping reads\tUniquely mapping reads\tMultiply mapping reads\n";
foreach my $sample_map (@samples) {
	open MAPSTAT,"$id/geneExpression/$sample_map.Mapped.stat.xls";
	$mapstat .= $sample_map;
	while (<MAPSTAT>) {
		next if (/^\s*$/);
		my @mapstat = split /\s+/,$_;
		if ($mapstat[0] eq 'Total') {
			$mapstat .= "\t$mapstat[-2]";
		}
		elsif ($mapstat[0] eq 'Mapped' || $mapstat[0] eq 'Uniq' || $mapstat[0] eq 'Multi') {
			$mapstat .= "\t$mapstat[-2]($mapstat[-1])";
		}
	}
	$mapstat .= "\n";
}
close MAPSTAT;
open MAPSTATOUT,">$id/geneExpression/Total.Mapped.stat";
print MAPSTATOUT $mapstat;
close MAPSTATOUT;

###############################SSR标记数量
#open SSR,"$id/Unigene/Unigene_SSR/$name.Unigene.1000.fa.stat.xls" or print "$!:$id/Unigene/Unigene_SSR/$name.Unigene.1000.fa.stat.xls\n";
#my $ssr_num;
#while (<SSR>) {
#	next unless (/^Total/);
#	chomp;
#	$ssr_num = (split "\t")[-1];
#}
#close SSR;

my $ssr_num;
my $ssr_stat_file = (glob "$id/Unigene/Unigene_SSR/*.statistics")[0];
open (SSR, $ssr_stat_file) or die $!;

while (<SSR>) {
    next unless (/^Total number of identified SSRs/);
    $ssr_num = (split /:/)[1];
    $ssr_num =~ s/\s+//g;
}

close SSR;

###############################Unigene注释数量
my $unigene_anno_num = `wc -l $id/Unigene/Unigene_Anno/Integrated_Function.annotation.xls`;
$unigene_anno_num = (split /\s+/,$unigene_anno_num)[0];
$unigene_anno_num --; 


######################################################################################
my $pic_num = 1;
my $tab_num = 1;
my $references;
my $references_num = 0;
######################################################################################
my $report_pic_path;

my $report_cfg .='<report_config>'."\n";
$report_cfg .=&info("report_version","3");
$report_cfg .=&info("report_name",$pro_inf{project_name});
$report_cfg .=&info("report_customer",'用户：'.$pro_inf{user_name}.' | 单位：'.$pro_inf{user_unit});
$report_cfg .=&info("report_code",$pro_inf{project_id});
$report_cfg .=&info("report_createtime",'样品到位 '.$pro_inf{pro_sam_receive_date}.' |样品检测合格 '.$pro_inf{pro_sam_pass_date}.' | 项目启动 '.$pro_inf{pro_launch_date}.' | 项目完成 '.$pro_inf{pro_finish_date});


$report_cfg .=&info("report_abstract","合同关键指标：&lt;br/&gt;完成${samples_num}个样品的转录组测序，每个样品产生不少于4Gb数据，保证Q30达到85%。完成组装转录本和构建Unigene库的分析。完成Unigenes的表达量分析和差异表达基因分析。完成SSR分析和SNP分析。完成Unigene功能注释和差异表达基因功能注释分析。&lt;br/&gt;分析结果概述：&lt;br/&gt;完成${samples_num}个样品的转录组测序，共获得${clean_data}Gb Clean Data，各样品测序数据量均达到4G，Q30碱基百分达到$clean_Q30\%及以上。De novo组装后共获得$assembly{'Uni'}{'T_n'}条Unigenes，其中长度在1kb以上Unigenes有${unigene_1kb}条。进行基于Unigene库的基因结构分析，包括CDS预测、SSR分析以及样品间SNP分析，获得SSR标记${ssr_num}个。进行Unigenes的生物信息学注释，包括与NR、Swiss-Prot、KEGG、COG、GO数据库的比对，共获得${unigene_anno_num}条Unigenes的注释结果。进行各样品中基因表达量和差异表达基因的分析。对差异表达基因进行模式聚类、功能注释以及富集性分析。");

my $out .= '<results>'."\n";

$abstract="完成${samples_num}个样品的转录组测序，共获得${clean_data}Gb Clean Data，各样品测序数据量均达到4G，Q30碱基百分达到$clean_Q30\%及以上。De novo组装后共获得$assembly{'Uni'}{'T_n'}条Unigenes，其中长度在1kb以上Unigenes有${unigene_1kb}条。进行基于Unigene库的基因结构分析，包括CDS预测、SSR分析以及样品间SNP分析，获得SSR标记${ssr_num}个。进行Unigenes的生物信息学注释，包括与NR、Swiss-Prot、KEGG、COG、GO数据库的比对，共获得${unigene_anno_num}条Unigenes的注释结果。进行各样品中基因表达量和差异表达基因的分析。对差异表达基因进行模式聚类、功能注释以及富集性分析。";

#############################################测序数据统计与评估
#现有流程生成的web_report中只有一个数据评估结果

$out .='<menu-1>'."\n";
$out .=&info("name","测序数据统计与评估");

$out .='<menu-2>'."\n";
$out .=&info("name","概述");
$out .=&info("desc","基于边合成边测序（Sequencing By Synthesis，SBS）技术，使用Illumina HiSeq2500高通量测序平台对cDNA文库进行测序，能够产出大量的高质量Reads，测序平台产出的这些Reads或碱基称为原始数据（Raw Data），其大部分碱基质量打分能达到或超过Q30。&lt;br/&gt;Raw Data以FASTQ格式存储，每个测序样品的Raw Data包括两个FASTQ文件，分别包含所有cDNA片段两端测定的Reads。");
$out .=&file("图$pic_num FASTQ格式文件示意图","","biocloud/template/cleandata_FASTQ_format.png","FASTQ格式文件示意图如下：","注：FASTQ文件中通常每4行对应一个序列单元：第一行以\@开头，后面接着序列标识（ID）以及其它可选的描述信息；第二行为碱基序列，即Reads；第三行以“+”开头，后面接着可选的描述信息；第四行为Reads每个碱基对应的质量值编码，长度必须和Reads的序列长度相同。","80k","png","png");
$pic_num++;
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","测序碱基质量值");
$out .=&info("desc","碱基质量值（Quality Score或Q-score）是碱基识别（Base Calling）出错的概率的整数映射。通常使用的Phred质量评估公式为：Q-score=-10*log10P，公式中，P为碱基识别出错的概率。下表给出了碱基质量值与碱基识别出错的概率的对应关系：");
$out .=&examplefile("表$tab_num 碱基质量值与碱基识别出错的概率的对应关系表","","biocloud/template/base_erro_ratio.txt","","","100k","xls","detail");
$tab_num++;
my $raw_qua_pic = (glob "$id/cleandata/PNG/*$samples[0]*.quality.png")[0];
$raw_qua_pic = basename($raw_qua_pic);
$out .=&file("图$pic_num 碱基质量值分布图","","cleandata/PNG/$raw_qua_pic","","注：横坐标为Reads的碱基位置， 纵坐标为单碱基错误率","100k","png","png");
$pic_num++;
my $atcg_map=(glob "$id/cleandata/PNG/*.acgtn.png")[0];
$atcg_map=basename($atcg_map);
$out .=&file("图$pic_num ATCG含量分布图","","cleandata/PNG/$atcg_map","","注：横坐标为Reads的碱基位置，纵坐标为单碱基所占比例。","100k","png","png");
$pic_num++;
$out .='</menu-2>'."\n";

my $data_stat=(glob "$id/cleandata/PNG/*.rawDataStat.png")[0];
if(defined $data_stat){
$out .='<menu-2>'."\n";
$out .=&info("name","测序质量控制");
$out .=&info("desc","在进行后续分析之前，首先需要确保这些reads有足够高的质量，以保证后续分析的准确。另外，一般Raw Data中会有极少部分的Reads带有测序引物、接头等人工序列，但它们并不是基因转录的产物，需要将其从Reads中截除。因此一系列的测序质量控制如下：&lt;br/&gt;(1) 去除测序接头以及引物序列；&lt;br/&gt;(2) 过滤低质量值数据，确保数据质量。&lt;br/&gt;经过上述一系列的质量控制之后得到的高质量reads，称之为Clean Data。Clean Data同样以FASTQ格式提供给客户。");

$data_stat=basename($data_stat);
$out .=&file("图$pic_num RawData数据分布统计图","","cleandata/PNG/$data_stat","","注：Adapter related：过滤掉的含有接头Reads数占总Raw Reads数的比例。Low quality：过滤掉的低质量Reads数占总Raw Reads数的比例。Clean Reads：经过以上过滤得到的Clean Reads 数占总Raw Reads 数的比例。","100k","png","png");
$pic_num++;
$out .='</menu-2>'."\n";
}
$out .='<menu-2>'."\n";
$out .=&info("name","测序数据产出统计");
$out .=&info("desc","本项目各样品数据产出统计见表${tab_num}：");
$out .=&examplefile("表$tab_num 样品测序数据评估统计表","暂无","cleandata/AllSample_GC_Q.stat",'','注：Samples：样品信息单样品名称；BMK-ID：百迈客样品分析编号；ReadSum：Clean Data中pair-end Reads总数；BaseSum：Clean Data总碱基数；GC(%)：Clean Data GC含量，即Clean Data中G和C两种碱基占总碱基的百分比；Q30(%)：Clean Data质量值大于或等于30的碱基所占的百分比。',"100k","xls","stat");
$tab_num++;
$out .='</menu-2>'."\n";

if (defined $cleandata_dir) {
	$out .='<ori_data>'."\n";
	foreach  (glob "$cloud_dir/display_cleandata/*.fq") {
		$_ =~ s/$cloud_dir\///;
		$out .=&examplefile("原始数据文件","","$_","测序原始数据——fastq格式","注：FASTQ文件中通常每4行对应一个序列单元：第一行以\@开头，后面接着序列标识（ID）以及其它可选的描述信息；第二行为碱基序列，即reads；第三行以&quot;+&quot;开头，后面接着可选的描述信息；第四行为reads每个碱基对应的质量值编码，长度必须和reads的序列长度相同。","100k","txt","fastq");
	}
	$out .='</ori_data>'."\n";
}

$out .='</menu-1>'."\n";

#############################################组装结果统计
$out .='<menu-1>'."\n";
$out .=&info("name","转录组测序数据组装");

$out .='<menu-2>'."\n";
$out .=&info("name","概述");
$out .=&info("desc","获得高质量的测序数据之后，需要对其进行序列组装。Trinity&references_and_link{(Grabherr MG, Haas BJ, Yassour M,, et al. Full length transcriptome assembly from RNA Seq data without a reference genome. Nature Biotechnology. 2011.(29): 644 -652.),(http://www.nature.com/nbt/journal/v29/n7/full/nbt.1883.html)}是一款专门为高通量转录组测序设计的组装软件。转录本测序深度除了受测序数据量等影响，还与该转录本的表达丰度有关。测序深度会直接影响组装的好坏。为了使各样品中表达丰度较低的转录本组装得更完整，对于同物种的测序样品推荐合并组装可以间接增加测序深度，从而使转录结果更完整，同时也有利于后续的数据分析；而对于不同物种的样品，由于基因组间存在差异，推荐采用分别组装或分开分析。");
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","组装原理");
$out .=&info("desc","Trinity软件首先将测序Reads打断为较短的片段（K-mer），然后将这些小片段延伸成较长的片段（Contig），并利用这些片段之间的重叠，得到片段集合（Component），最后利用De Bruijn图的方法和测序Read信息，在各个片段集合中分别识别转录本序列。&lt;br/&gt;Trinity软件具体组装过程：&lt;br/&gt;(1) 将测序Reads按照指定K-mer打断来构建K-mer库，去除可能包含错误的K-mer；&lt;br/&gt;(2) 选择频率最高的K-mer作为种子向两端进行贪婪延伸（以K-1个碱基的Overlap为标准，低复杂度或只出现一次的K-mer不能作为种子），不断循环此过程直至耗光K-mer库；&lt;br/&gt;(3) 对(2)中得到的Contig进行聚簇，得到Component（Contig之间包含K-1个碱基的Overlap，并且有一定数目K-mer分别有一半比对在两条Contig上，这样的Contig会聚为一个Component）；&lt;br/&gt;(4) 对每个Component中的Contig构建De Bruijn图；&lt;br/&gt;(5) 对(4)中得到的De Bruijn图进行简化（合并节点，修剪边沿）；&lt;br/&gt;(6) 以真实的Read来解开De Bruijn图，获得转录本序列。&lt;br/&gt;组装的原理图如图${pic_num}：");
$out .=&examplefile("图$pic_num Trinity组装程序原理图","","biocloud/template/Trinity_workflow.png","","","100k","png","png");
$pic_num++;
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","组装结果统计");
my $unigen_distribution_png = (glob "$id/Unigene/*.distribution.png");
$unigen_distribution_png = basename($unigen_distribution_png);

if (-d "$id/Assembly/All_Combination") {
	$out .=&info("desc","组装共得到".&format_figure($assembly{'Trans'}{'T_n'})."条转录本和".&format_figure($assembly{'Uni'}{'T_n'})."条Unigene，转录本与Unigene的N50分别为".&format_figure($assembly{'Trans'}{'N50'})."和".&format_figure($assembly{'Uni'}{'N50'})."，具体的统计信息见图${pic_num}与表${tab_num}。");
	$out .=&file("图$pic_num Unigene长度分布图","暂无","Unigene/$unigen_distribution_png",'','注：横坐标表示Unigene的不同长度区间；纵坐标表示某一区间内Unigene数量在总转录组的比例。',"100k","png","png");
	$out .=&examplefile("表$tab_num 组装结果统计表","暂无","Assembly/All_Combination.stat.xls",'','注：Length Range：表示Contig/Transcript/Unigene的不同长度区间；表格中的数字表示相应区间内Contig/Transcript/Unigene的数量，括号内的百分比表示相应长度区间内Contig/Transcript/Unigene所占的比例；Total Number：表示组装得到的Contig/Transcript/Unigene的总数；Total Length：表示组装得到的Contig/Transcript/Unigene的总长度；N50 Length：表示Contig/Transcript/Unigene的N50的长度；Mean Length：表示Contig/Transcript/Unigene的平均长度。',"100k","xls","xls");
}
else {
	$out .=&info("desc","组装共得到".&format_figure($assembly{'Uni'}{'T_n'})."条Unigene，N50为".&format_figure($assembly{'Uni'}{'N50'})."，具体的统计信息见图${pic_num}与表${tab_num}。");
	$out .=&file("图$pic_num Unigene长度分布图","暂无","Unigene/$unigen_distribution_png",'','注：横坐标表示Unigene的不同长度区间；纵坐标表示某一区间内Unigene数量在总转录组的比例。',"100k","png","png");
	$out .=&examplefile("表$tab_num 组装结果统计表","暂无","Assembly/All_Combination.stat.xls",'','注：Length Range：表示Contig/Transcript/Unigene的不同长度区间；表格中的数字表示相应区间内Contig/Transcript/Unigene的数量，括号内的百分比表示相应长度区间内Contig/Transcript/Unigene所占的比例；Total Number：表示组装得到的Contig/Transcript/Unigene的总数；Total Length：表示组装得到的Contig/Transcript/Unigene的总长度；N50 Length：表示Contig/Transcript/Unigene的N50的长度；Mean Length：表示Contig/Transcript/Unigene的平均长度。',"100k","xls","xls");
}

$tab_num++;$pic_num++;
$out .='</menu-2>'."\n";

$out .='<ori_data>'."\n";
if (-d "$id/Assembly/All_Combination") {
	$out .=&examplefile("Contigs序列文件","","Assembly/All_Combination/contigs/All_Combination.contigs.fa","","组装出的contig序列——fasta格式","100k","txt","fasta");
	$out .=&examplefile("Transcripts序列文件","","Assembly/All_Combination/Transcripts/All_Combination.Transcripts.fa","","组装出的transcript序列——fasta格式","100k","txt","fasta");
	$out .=&examplefile("Unigenes序列文件","","Assembly/All_Combination/Unigenes/All_Combination.Unigenes.fa","","组装出的Unigene序列——fasta格式","100k","txt","fasta");
}
else {
	foreach  (glob "$id/Assembly/*/Unigenes/*.Unigenes.fa") {
		my $name = basename($_);
		my $dir = $_;
		$dir =~ s/$id//;
		$name =~ s/\.Unigenes\.fa//;
		$out .=&examplefile("${name} Unigenes序列文件","","$dir","","组装出的Unigene序列——fasta格式","100k","txt","fasta");
	}
}
$out .='</ori_data>'."\n";

$out .='</menu-1>'."\n";

#############################################转录组文库质量评估
$out .='<menu-1>'."\n";
$out .=&info("name","转录组文库质量评估");

$out .='<menu-2>'."\n";
$out .=&info("name","概述");
$out .=&info("desc","合格的转录组测序文库是转录组数据分析结果可靠的必要条件，为确保测序文库的质量，从以下3个不同角度对转录组测序文库进行质量评估：&lt;br/&gt;(1) 通过检验插入片段在Unigene上的分布，评估mRNA片段化的随机性、mRNA的降解情况；&lt;br/&gt;(2) 通过绘制插入片段的长度分布图，评估插入片段长度的离散程度；&lt;br/&gt;(3) 通过绘制饱和度图，评估文库容量和比对到Unigene库的Reads（Mapped Reads）是否充足。");
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","mRNA片段化随机性检验");
$out .=&info("desc","mRNA片段化后的插入片段大小选择，可以理解为从mRNA序列中独立随机地抽取子序列，如果样本量（mRNA数目）越大、打断方式和时间控制得越合适，那么目的RNA每个部分被抽取到的可能性就越接近，即mRNA片段化随机性越高，mRNA上覆盖的reads越均匀。&lt;br/&gt;通过比对到Unigene的reads（mapped reads）在各mRNA转录本上的位置分布，模拟mRNA片段化结果，检验mRNA片段化的随机程度。各样品mapped reads在mRNA转录本上的位置分布图如图${pic_num}：");
if($Old==1){
	$out .=&file("图$pic_num Mapped Reads在mRNA上的位置分布图","","geneExpression/$samples[0].randcheck.png",'',"注：横坐标为mRNA位置，纵坐标为对应位置区间内Reads在Mapped Reads中所占百分比。由于参考的mRNA长度不同，作图时将每个mRNA按照长度划分成100个区间，进而统计每一区间内的Mapped Reads数目及所占的比例，图中反映的是所有mRNA各个区间内的Mapped Reads比例的汇总。","100k","png","png");
}elsif($Old==2){
	#Old version
	$out .=&file("图$pic_num 转录组测序数据饱和度模拟图","暂无","geneExpression/$samples[0].express.gene_tag.png",'',"注：通过将Reads等量地分成100份，逐渐增加数据查看检测到的基因数量来绘制饱和度曲线。横坐标为reads数目（以10^6为单位），纵坐标为检测到的基因数量（以10^3为单位）。表达量FPKM不小于0.5的基因为表达的基因。","100k","png","png");
}

$out .='</menu-2>'."\n";
$pic_num++;

$out .='<menu-2>'."\n";
$out .=&info("name","插入片段长度检验");
$out .=&info("desc","插入片段长度的离散程度能直接反映出文库制备过程中切胶或磁珠纯化的效果。样品$samples[0]的插入片段长度模拟分布图如图${pic_num}：");
$out .=&file("图$pic_num 插入片段长度模拟分布图","","geneExpression/$samples[0].insertSize.png","","注：横坐标为双端Reads在Unigene库中比对起止点之间的距离，范围为0到800bp；纵坐标为比对起止点之间不同距离的双端Reads或插入片段数量。","100k","","");
$out .='</menu-2>'."\n";
$pic_num++;

$out .='<menu-2>'."\n";
$out .=&info("name","转录组测序数据饱和度检验");
$out .=&info("desc","充足的有效数据量是信息分析准确的必要条件。相比传统的基因表达检测方法，转录组测序拥有较高的灵敏度，不仅能检测到高表达的基因，还能检测到低表达的基因。转录组测序检测到的基因数目与测序数据量成正相关性，即测序数据量越大，检测到的基因数目越多。但一个物种的基因数目是有限的，而且基因转录具有时间特异性和空间特异性，所以随着测序量的增加，检测到的基因数目会趋于饱和。&lt;br/&gt;为了评估数据是否充足，需要查看随着测序数据量的增加，新检测到的基因是否越来越少或没有，即检测到的基因数目是否趋于饱和。&lt;br/&gt;使用各样品的Mapped Reads对检测到的基因数目的饱和情况进行模拟，绘制曲线图如图${pic_num}：");
$out .=&file("图$pic_num 转录组测序数据饱和度模拟图","暂无","geneExpression/Total.gene_tag.png",'',"注：通过将Reads等量地分成100份，逐渐增加数据查看检测到的基因数量来绘制饱和度曲线。横坐标为reads数目（以10^6为单位），纵坐标为检测到的基因数量（以10^3为单位）。表达量FPKM不小于0.5的基因为表达的基因。","100k","png","png");
$out .='</menu-2>'."\n";
$pic_num++;
$out .='</menu-1>'."\n";

#############################################Unigenes注释
$out .='<menu-1>'."\n";
$out .=&info("name","Unigenes注释");

$out .='<menu-2>'."\n";
$out .=&info("name","概述");
$out .=&info("desc","使用BLAST&references_and_link{(Altschul SF, Madden TL, Schäffer AA, et al. Gapped BLAST and PSI BLAST: A New Generation of Protein Database Search Programs. Nucleic Acids Research. 1997. 25(17): 3389 -3402.),(http://blast.ncbi.nlm.nih.gov/Blast.cgi)}软件(version 2.2.26)将Unigene序列与nr&references_and_link{(Deng YY, Li JQ, Wu SF, et al. Integrated nr Database in Protein Annotation System and Its Localization. Computer Engineering. 2006. 32(5):71 -74.),(ftp://ftp.ncbi.nih.gov/blast/db/)}、Swiss-Prot&references_and_link{(Apweiler R, Bairoch A, Wu CH, et al. UniProt: the Universal Protein knowledgebase. Nucleic Acids Research. 2004. 32(Database issue):D115-9.),(http://www.uniprot.org/)}、GO&references_and_link{(Ashburner M, Ball C A, Blake J A, et al. Gene ontology: tool for the unification of biology. Nature genetics. 2000. 25(1): 25-29.),(http://www.geneontology.org/)}、COG&references_and_link{(Tatusov R L, Galperin M Y, Natale D A. The COG database: a tool for genome scale analysis of protein functions and evolution. Nucleic Acids Research. 2000. 28(1):33-36.),(http://www.ncbi.nlm.nih.gov/COG/)}、KEGG&references_and_link{(Kanehisa M, Goto S, Kawashima S, et al. The KEGG resource for deciphering the genome. Nucleic Acids Research. 2004. 32(Database issue):D277 -D280.),(http://www.genome.jp/kegg/)}数据库比对，使用KOBAS2.0&references_and_link{(Xie, C., Mao, X., Huang, J., Ding, Y., Wu, J., Dong, S., Kong, L., Gao, G., Li, C. and Wei, L. (2011) KOBAS 2.0: a web server for annotation and identification of enriched pathways and diseases. Nucleic Acids Res, 39, W316-322.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3125809/)}得到Unigene在KEGG中的KEGG Orthology结果，预测完Unigene的氨基酸序列之后使用HMMER&references_and_link{(Eddy S.R. Profile hidden Markov models (1998) [Bioinformatics Italic], 14 (9), pp. 755-763.),(http://bioinformatics.oxfordjournals.org/content/14/9/755.short)}软件与Pfam&references_and_link{(Finn RD, Bateman A, Clements J, et al. Pfam: the protein families database. [Nucleic Acids Research Italic], 2013: gkt1223.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965110/)}数据库比对，获得Unigene的注释信息。&lt;br/&gt;NR数据库是NCBI中的非冗余蛋白质数据库，包含了Swiss-Prot、PIR(Protein Information Resource)、PRF(Protein Research Foundation)、PDB(Protein Data Bank)蛋白质数据库及从GenBank和RefSeq的CDS数据翻译过来的蛋白质数据信息。&lt;br/&gt;Swiss-Prot数据库是由EBI（欧洲生物信息学研究所）负责维护的数据库，包含了有相关参考文献且经过校对的蛋白质注释信息数据库，可信度很高。&lt;br/&gt;COG (Clusters of Orthologous Groups) 数据库是对基因产物进行同源分类的数据库，是一个较早的识别直系同源基因的数据库，通过对完整的原核生物的蛋白质序列大量比较而来的，现在已经扩展到包含630个完整的基因组。&lt;br/&gt;GO (Gene Ontology) 数据库是一个国际标准化的基因功能分类体系，提供了一套动态更新的标准词汇表来全面描述生物体中基因和基因产物的功能属性。该数据库总共有三大类，分别是分子功能 (molecular function)，细胞组分 (cellular component) 和生物学过程 (biological process)，各自描述了基因产物可能行使的分子功能，以及所处的细胞环境和参与的生物学过程。GO数据库中最基本的概念是Term，每个条目都有一个Term名，比如&quot;cell&quot;、&quot;fibroblast growth factor receptor binding&quot;或者&quot;signal transduction&quot;，同时有一个唯一的编号，形如GO:nnnnnnn。&lt;br/&gt;KEGG(Kyoto Encyclopedia of Genes and Genomes)数据库是系统分析基因产物在细胞中的代谢途径以及这些基因产物功能的数据库。它整合了基因组、化学分子和生化系统等方面的数据，包括代谢通路(PATHWAY)、药物(DRUG)、疾病(DISEASE)、基因序列 (GENES) 及基因组 (GENOME) 等。利用该数据库有助于把基因及表达信息作为一个整体的网络进行研究。");
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","Unigenes注释统计");
$out .=&info("desc","本项目最终获得".&format_figure($unigene_anno_num)."个有注释信息的Unigene。基因注释的统计结果见表${tab_num}：");
$out .=&examplefile("表$tab_num Unigenes注释统计表","","Unigene/Unigene_Anno/Function_Annotation.stat.xls","","注：Anno_Database：表示各功能数据库；Annotated_Number：注释到数据库的Unigenes数；300&lt;=length&lt;1000 ：注释到数据库的长度大于300bp的Unigenes数；length&gt;=1000 ：注释到数据库的长度大于1000bp的Unigenes数。","100k","xls","stat");
$out .='</menu-2>'."\n";
$tab_num++;

$out .='</menu-1>'."\n";


#############################################基因结构分析
$out .='<menu-1>'."\n";
$out .=&info("name","基因结构分析");

$out .='<menu-2>'."\n";
$out .=&info("name","CDS预测");
$out .=&info("desc","TransDecoder软件基于开放阅读框（Open Reading Frame，ORF）长度、对数似然函数值（Log-likelihood Score）、氨基酸序列与Pfam数据库蛋白质结构域序列的比对等信息，能够从转录本序列中识别可靠的潜在编码区序列（Coding Sequence，CDS），是Trinity和Cuffinks等软件官方推荐的CDS预测软件。");
$out .='</menu-2>'."\n";
#$tab_num++;

$out .='<menu-2>'."\n";
$out .=&info("name","简单重复序列分析");
$out .=&info("desc","MISA（MIcroSAtellite identification tool）是一款鉴定简单重复序列（Simple Sequence Repeat，SSR）的软件，其参考网址见附表。它可以通过对Unigene序列的分析，鉴定出7种类型的SSR：单碱基（Mono-nucleotide）重复SSR、双碱基（Di-nucleotide）重复SSR、三碱基（Tri-nucleotide）重复SSR、四碱基（Tetra-nucleotide）重复SSR、五碱基（Penta-nucleotide）重复SSR和六碱基（Hexa-nucleotide）重复SSR。&lt;br/&gt;利用MISA软件对筛选得到的1kb以上的Unigene做SSR分析，统计结果见表${tab_num}：");
$out .=&examplefile("表$tab_num SSR分析结果统计表","暂无","Unigene/Unigene_SSR/$name.Unigene.1000.fa.statistics",'',"注：Total number of sequences examined：评估的序列数目；Total size of examined sequences (bp)：评估的序列总碱基数目；Total number of identified SSRs：识别的SSR总数；Number of SSR containing sequences：包含SSR的序列数目；Number of sequences containing more than 1 SSR ：包含1个以上SSR的序列数目；Number of SSRs present in compound formation：以复合物形式存在的SSR数目；1：单碱基重复SSR；2：双碱基重复SSR；3：三碱基重复SSR；4：四碱基重复SSR；5：五碱基重复SSR；6：六碱基重复SSR。","100k","xls","xls");
$tab_num++;
$out .='</menu-2>'."\n";

unless (defined $oneSample) {
	$out .='<menu-2>'."\n";
	$out .=&info("name","SNP分析");
#	$out .=&info("desc","SOAPsnp&references_and_link{(Li R, Li Y, Fang X, et al. SNP detection for massively parallel whole genome resequencing. Genome Research. 2009. (19):1124-1132.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/)}（version 1.00）是分析样本间SNP的常用软件软件。该软件利用基于贝叶斯理论而建立的一套方法，综合考虑了碱基质量、比对情况、测序错误率等因素，得到一致性序列质量值作为SNP可靠性的标准。&lt;br/&gt;将每个样品转录组测序得到的reads与组装得到的Unigene比对，可以观察到部分基因序列中存在多态性位点。进而可以分析这些SNP位点是否影响了基因的表达水平或者蛋白产物的种类。&lt;br/&gt;利用SOAPsnp软件进行样品间的SNP分析，两两样品间SNP数量统计见表${tab_num}。");
#	$out .=&examplefile("表$tab_num SNP数量统计表","暂无","SNP_Analysis/SNP.stat.xls","","注：Type：样品组合，S1和S2分别对应组合中的前一个和后一个样品，homo表示纯合基因型，hete表示杂合基因型。Total：所有类型的统计信息。","100k","xls","xls");

    if (-f "$id/SNP_Analysis/SNP.stat.xls") {
        # no_ref_trans pipeline v1.7 & previous versions
        $out .=&info("desc","SOAPsnp&references_and_link{(Li R, Li Y, Fang X, et al. SNP detection for massively parallel whole genome resequencing. Genome Research. 2009. (19):1124-1132.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/)}（version 1.00）是分析样本间SNP的常用软件软件。该软件利用基于贝叶斯理论而建立的一套方法，综合考虑了碱基质量、比对情况、测序错误率等因素，得到一致性序列质量值作为SNP可靠性的标准。&lt;br/&gt;将每个样品转录组测序得到的reads与组装得到的Unigene比对，可以观察到部分基因序列中存在多态性位点。进而可以分析这些SNP位点是否影响了基因的表达水平或者蛋白产物的种类。&lt;br/&gt;利用SOAPsnp软件进行样品间的SNP分析，两两样品间SNP数量统计见表${tab_num}。");
        $out .=&examplefile("表$tab_num SNP数量统计表","暂无","SNP_Analysis/SNP.stat.xls","","注：Type：样品组合，S1和S2分别对应组合中的前一个和后一个样品，homo表示纯合基因型，hete表示杂合基因型。Total：所有类型的统计信息。","100k","xls","xls");
    } elsif (-f "$id/SNP_Analysis/AllSample.snp.stat") {
        # no_ref_trans pipeline v1.7.2 +
        $out .=&info("desc","利用针对RNA-Seq的比对软件STAR&references_and_link{(Dobin A, Davis C A, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner[J]. Bioinformatics, 2013, 29(1): 15-21.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)}对每个样本的Reads与Unigene序列进行比对，并通过GATK&references_and_link{(McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data[J]. Genome Research. 2010, 20(9): 1297-1303.),(https://www.broadinstitute.org/gatk/)}针对RNA-Seq的SNP识别（SNP Calling）流程，识别单核苷酸多态性（Single Nucleotide Polymorphism，SNP）位点。进而可以分析这些SNP位点是否影响了基因的表达水平或者蛋白产物的种类。识别标准如下：&lt;br/&gt;(1) 35bp范围内连续出现的单碱基错配不超过3个；&lt;br/&gt;(2) 经过序列深度标准化的SNP质量值大于2.0。&lt;br/&gt;按照以上条件筛选，最终获得各样本SNP位点信息。&lt;br/&gt;根据SNP位点的等位（Allele）数目，即测序Reads支持的不同碱基的数目，可以将SNP位点分为纯合型SNP位点（只有一个等位）和杂合型SNP位点（两个或多个等位）。不同物种杂合型SNP所占的比例存在差异。各样品SNP位点数目统计见表${tab_num}。");
        $out .=&examplefile("表$tab_num SNP数量统计表","暂无","SNP_Analysis/AllSample.snp.stat","","注：Samples：样品编号；HomoSNP：纯合型SNP数目；HeteSNP：杂合型SNP数目；AllSNP：纯合型和杂合型SNP总数目。","100k","xls","xls");
    }

	$tab_num++;
	$out .='</menu-2>'."\n";
	$tab_num++;
}

$out .='<ori_data>'."\n";
if($Old==1){
	$out .=&examplefile("ORF预测编码区序列文件","","Unigene/Unigene_CDS/$name.Unigene.cds.fa","","Unigene编码区碱基序列——fasta格式","100k","txt","fa");
	$out .=&examplefile("ORF预测蛋白序列文件","","Unigene/Unigene_CDS/$name.Unigene.pep.fa","","Unigene编码的氨基酸序列——fasta格式","100k","txt","fa");
}elsif($Old==2){
	#Old version
	$out .=&examplefile("ORF预测编码区序列文件","","Unigene/Unigene_Orf/$name.Unigene.cds.fa","","Unigene编码区碱基序列——fasta格式","100k","txt","fa");
	$out .=&examplefile("ORF预测蛋白序列文件","","Unigene/Unigene_Orf/$name.Unigene.pep.fa","","Unigene编码的氨基酸序列——fasta格式","100k","txt","fa");
}


$out .=&examplefile("序列文件","","Unigene/Unigene_SSR/$name.Unigene.1000.fa","","碱基数大于1000的Unigene序列——fasta格式","100k","txt","fasta");
#$out .=&examplefile("序列及信息","","Unigene/Unigene_SSR/$name.Unigene.1000.fa.detail.xls","","MISA原始输出结果，ID：Unigene名；Length：Unigene长度；SSR_type：SSR类型；SSR：具体的重复碱基及重复次数；Start：SSR的起始位置；End：SSR的终止位置；Sequence：SSR的碱基","100k","txt","list");
$out .=&examplefile("统计信息","","Unigene/Unigene_SSR/$name.Unigene.1000.fa.SSR.result.xls","","每一行为一个SSR标记及其引物设计结果，每个SSR标记最多设计三对引物;Gene_ID：Unigene编号；SSR_nr：同一Unigene上的SSR序号；SSR_type：SSR类型，包括完美单碱基重复（p1）、完美双碱基重复（p2）、完美三碱基重复（p3）、完美四碱基重复（p4）、完美五碱基重复（p5）、完美六碱基重复（p6）和混合SSR（c，即包含至少两个完美SSR，且之间距离小于100bp）；SSR：SSR序列，括号内为重复单元，括号外数字表示重复次数；Size：SSR的长度；SSR_Start：SSR在Unigene上的开始位置；SSR_End：SSR在Unigene上的结束位置;FPr1(5'-3')：第一条正向引物序列；Tm：第一条正向引物序列的退火温度，单位为°C；Size：第一条正向引物序列的长度；RPr1(5'-3')：第一条反向引物序列；Tm：第一条反向引物序列的退火温度，单位为°C；Size：第一条反向引物序列的长度；Psize：产物的长度；PStart：产物在基因上的开始位置；PEnd：产物在基因上的结束位置。","100k","txt","list");
$out .=&examplefile("统计信息","","Unigene/Unigene_SSR/$name.Unigene.1000.fa.statistics","","MISA结果统计文件，内含每种SSR的个数","100k","txt","list");

if (-f "$id/SNP_Analysis/SNP.stat.xls") {
    # no_ref_trans pipeline v1.7 & previous versions
    foreach (glob "$id/SNP_Analysis/*_vs_*") {
        $_ = basename($_);
        my ($snp_1,$snp_2) = split "_vs_";
        $out .=&examplefile("${snp_1}杂合${snp_2}纯合","","SNP_Analysis/$_/${snp_1}.hete_${snp_2}.homo.snp.xls","","GeneID：基因名；pos：SNP在Unigene上的位置；reGenotype：Unigene的基因型；${snp_1}.Genotype：样品${snp_1}的基因型；MajorGenotype：首要基因型；Major_Depth：首要基因型深度；MinorGenotype：次要基因型；Minor_Depth：次要基因型深度；${snp_1}.Depth：样品${snp_1}在此位点的总深度；${snp_2}.Genotype：样品${snp_2}的基因型；${snp_2}.Depth：样品${snp_2}在此位点的总深度；Score：此位点打分","100k","xls","list");
        $out .=&examplefile("${snp_1}纯合${snp_2}杂合","","SNP_Analysis/$_/${snp_1}.homo_${snp_2}.hete.snp.xls","","GeneID：基因名；pos：SNP在Unigene上的位置；reGenotype：Unigene的基因型；${snp_1}.Genotype：样品${snp_1}的基因型；${snp_1}.Depth：样品${snp_1}在此位点的总深度；${snp_2}.Genotype：样品${snp_2}的基因型；MajorGenotype：首要基因型；Major_Depth：首要基因型深度；MinorGenotype：次要基因型；Minor_Depth：次要基因型深度；${snp_2}.Depth：样品${snp_2}在此位点的总深度；Score：此位点打分","100k","xls","list");
        $out .=&examplefile("${snp_1}、${snp_2}杂合","","SNP_Analysis/$_/${snp_1}.${snp_2}.hete.snp.xls","","GeneID：基因名；pos：SNP在Unigene上的位置；reGenotype：Unigene的基因型；${snp_1}.Genotype：样品${snp_1}的基因型；MajorGenotype：首要基因型；Major_Depth：首要基因型深度；MinorGenotype：次要基因型；Minor_Depth：次要基因型深度；${snp_1}.Depth：样品${snp_1}在此位点的总深度；${snp_2}.Genotype：样品${snp_2}的基因型；MajorGenotype：首要基因型；Major_Depth：首要基因型深度；MinorGenotype：次要基因型；Minor_Depth：次要基因型深度；${snp_2}.Depth：样品${snp_2}在此位点的总深度；Score：此位点打分","100k","xls","list");
        $out .=&examplefile("${snp_1}、${snp_2}纯合","","SNP_Analysis/$_/${snp_1}.${snp_2}.homo.snp.xls","","GeneID：基因名；pos：SNP在Unigene上的位置；reGenotype：Unigene的基因型；${snp_1}.Genotype：样品${snp_1}的基因型；${snp_1}.Depth：样品${snp_1}在此位点的总深度；${snp_2}.Genotype：样品${snp_2}的基因型；${snp_2}.Depth：样品${snp_2}在此位点的总深度；Score：此位点打分","100k","xls","list");
    }
} elsif (-f "$id/SNP_Analysis/AllSample.snp.stat") {
    # no_ref_trans pipeline v1.7.2 +
    foreach (glob "$id/SNP_Analysis/*.snp.list") {
        $_ = basename($_);
        my $snp_1 = (split '\.')[0];
        next if ($snp_1 eq 'final');
        $out .=&examplefile("样品${snp_1}SNP位点信息表","","SNP_Analysis/${snp_1}.snp.list","","注：GeneID：Unigene编号；Pos：SNP位点在Unigene上的位置；Ref：Unigene上的SNP等位；Alt：测序样品中识别到的其他的SNP等位；${snp_1}：样品${snp_1}该SNP位点的分型；Depth：样品${snp_1}该SNP位点的测序深度；AlleDepth：样品${snp_1}该SNP位点的各等位测序深度。","100k","xls","list");
    }
}

$out .='</ori_data>'."\n";
$out .='</menu-1>'."\n";

#############################################基因表达量分析
$out .='<menu-1>'."\n";
$out .=&info("name","基因表达量分析");

$out .='<menu-2>'."\n";
$out .=&info("name","Unigenes表达量估计");
#$out .=&info("desc","采用Bowtie&references_and_link{(Langmead B, Trapnell C, Pop M, et al. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology. 2009. 10(3): R25.),(http://bowtie-bio.sourceforge.net/index.shtml)}将各样品测序得到的reads与Unigene库进行比对，根据比对结果，结合RSEM&references_and_link{(Li B, Colin ND. RSEM: accurate transcript quantification from RNA Seq data with or without a reference genome. BMC Bioinformatics. 2011. (12):323.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3163565/)}进行表达量水平估计。利用FPKM值表示对应Unigene的表达丰度。&lt;br/&gt;FPKM&references_and_link{(Trapnell C, Williams B A, Pertea G, Mortazavi A, et al. Transcript assembly and quantification by RNA Seq reveals unannotated transcripts and isoform switching during cell differentiation. Nature Biotechnology 2010, 28(5):511 515.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3146043/)}（Fragments Per Kilobase of transcript per Million mapped reads）是每百万reads中来自某一基因每千碱基长度的reads数目，是转录组测序数据分析中常用的基因表达水平估算方法。FPKM能消除基因长度和测序量差异对计算基因表达的影响，计算得到的基因表达量可直接用于比较不同样品间的基因表达差异。FPKM计算公式如下：&lt;br/&gt;&lt;img src=&quot;/report/report/showjobimgByPath?filepath=$id/biocloud/template/FPKM_formula.png&quot; /&gt;&lt;br/&gt;公式中，cDNA Fragments表示比对到某一转录本上的片段数目，即双端Reads数目；Mapped Reads(Millions)表示Mapped Reads总数，以10^6为单位；Transcript Length(kb)：转录本长度，以10^3个碱基为单位。&lt;br/&gt;对每个基因的信息进行统计，结果文件见原始数据及更多操作模块。");

if (-d "$id/Assembly/All_Combination") {
    $out .=&info("desc","采用Bowtie&references_and_link{(Langmead B, Trapnell C, Pop M, et al. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology. 2009. 10(3): R25.),(http://bowtie-bio.sourceforge.net/index.shtml)}将各样品测序得到的reads与Unigene库进行比对，根据比对结果，结合RSEM&references_and_link{(Li B, Colin ND. RSEM: accurate transcript quantification from RNA Seq data with or without a reference genome. BMC Bioinformatics. 2011. (12):323.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3163565/)}进行表达量水平估计。利用FPKM值表示对应Unigene的表达丰度。&lt;br/&gt;FPKM&references_and_link{(Trapnell C, Williams B A, Pertea G, Mortazavi A, et al. Transcript assembly and quantification by RNA Seq reveals unannotated transcripts and isoform switching during cell differentiation. Nature Biotechnology 2010, 28(5):511 515.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3146043/)}（Fragments Per Kilobase of transcript per Million mapped reads）是每百万reads中来自某一基因每千碱基长度的reads数目，是转录组测序数据分析中常用的基因表达水平估算方法。FPKM能消除基因长度和测序量差异对计算基因表达的影响，计算得到的基因表达量可直接用于比较不同样品间的基因表达差异。FPKM计算公式如下：&lt;br/&gt;&lt;img src=&quot;/report/report/showjobimgByPath?filepath=$od/biocloud/template/FPKM_formula.png&quot; /&gt;&lt;br/&gt;公式中，cDNA Fragments表示比对到某一转录本上的片段数目，即双端Reads数目；Mapped Reads(Millions)表示Mapped Reads总数，以10^6为单位；Transcript Length(kb)：转录本长度，以10^3个碱基为单位。&lt;br/&gt;对每个基因的信息进行统计，结果文件见原始数据及更多操作模块。");
} else {
    $out .=&info("desc","采用BLAT&references_and_link{(Kent WJ. BLAT - the BLAST-like alignment tool. Genome Research. 2002 Apr;12(4):656-64.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC187518/)}将各样品测序得到的Reads与Unigene库进行比对，根据比对结果进行表达量水平估计。利用FPKM值表示对应Unigene的表达丰度。&lt;br/&gt;FPKM&references_and_link{(Trapnell C, Williams B A, Pertea G, Mortazavi A, et al. Transcript assembly and quantification by RNA Seq reveals unannotated transcripts and isoform switching during cell differentiation. Nature Biotechnology 2010, 28(5):511 515.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3146043/)}（Fragments Per Kilobase of transcript per Million mapped reads）是每百万reads中来自某一基因每千碱基长度的reads数目，是转录组测序数据分析中常用的基因表达水平估算方法。FPKM能消除基因长度和测序量差异对计算基因表达的影响，计算得到的基因表达量可直接用于比较不同样品间的基因表达差异。FPKM计算公式如下：&lt;br/&gt;&lt;img src=&quot;/report/report/showjobimgByPath?filepath=$od/biocloud/template/FPKM_formula.png&quot; /&gt;&lt;br/&gt;公式中，cDNA Fragments表示比对到某一转录本上的片段数目，即双端Reads数目；Mapped Reads(Millions)表示Mapped Reads总数，以10^6为单位；Transcript Length(kb)：转录本长度，以10^3个碱基为单位。&lt;br/&gt;对每个基因的信息进行统计，结果文件见原始数据及更多操作模块。");
}

$out .='</menu-2>'."\n";

if (@samples > 1) {
	my ($density_file) = glob "$id/geneExpression/all.*pkm_density.png";
	my ($box_file) = glob "$id/geneExpression/all.*pkm_box.png";
	$density_file =~ s/$id\///;
	$box_file =~ s/$id\///;
	$out .='<menu-2>'."\n";
	$out .=&info("name","样品基因表达量总体分布");
	$out .=&info("desc","利用转录组数据检测基因表达具有较高的灵敏度。通常情况下，能够测序到的蛋白质编码基因表达水平FPKM值横跨10-2到104六个数量级&references_and_link{(Djebali, Sarah and Mortazavi, et al. Landscape of transcription in human cells. Nature 2012, 489 (7414). pp. 101-108. ISSN 0028-0836.),(http://www.nature.com/nature/journal/v489/n7414/full/nature11233.html)}。");
	$out .=&file("图$pic_num 各样品FPKM分布密度图","暂无","$density_file",'',"注：图中不同颜色的曲线代表不同的样品，曲线上点的横坐标表示对应样品FPKM的对数值，点的纵坐标表示概率密度。","100k","png","png");
	$pic_num++;
	$out .=&file("图$pic_num 各样品FPKM箱线图","暂无","$box_file","从箱线图中不仅可以查看单个样品基因表达水平分布的离散程度，还可以直观的比较不同样品的整体基因表达水平。该项目各样品的FPKM分布箱线图如图${pic_num}：","注：图中横坐标代表不同的样品；纵坐标表示样品表达量FPKM的对数值。该图从表达量的总体离散角度来衡量各样品表达水平。","100k","png","png");
	$pic_num++;
	$out .='</menu-2>'."\n";
}

$out .='<ori_data>'."\n";
foreach (@samples) {
	$out .=&examplefile("${_}样品基因表达量结果文件","","geneExpression/$_.geneExpression.xls","","注：Gene_ID：Unigene编号；Effective_Length：Unigene有效长度，即该基因不同转录本的平均长度；Length：Unigene的长度；TPM：TPM方法标准化后的基因表达丰度值；FPKM：FPKM方法标准化后的基因表达丰度值；Transcript_ID(s)：转录本的编号；Expected_Count：标准化后的片段数。","100k","xls","list");
}
$out .='</ori_data>'."\n";
$out .='</menu-1>'."\n";


#################################只有一个样品时无差异表达分析及差异表达基因功能注释和富集分析两模块
unless (defined $oneSample) {
#################################只有一个样品时无差异表达分析及差异表达基因功能注释和富集分析两模块

#############################################差异表达分析
$out .='<menu-1>'."\n";
$out .=&info("name","差异表达分析");

$out .='<menu-2>'."\n";
$out .=&info("name","概述");
$out .=&info("desc","基因表达具有时间和空间特异性，外界刺激和内部环境都会影响基因的表达。在不同条件（如对照与处理、野生型和突变型、不同时间点、不同组织等）下，表达水平存在显著差异的基因，称之为差异表达基因（Differentially Expressed Gene，DEG）。同样地，表达水平存在显著差异的转录本，称之为差异表达转录本（Differentially Expressed Transcript，DET）。生物信息学中，寻找差异表达转录本或差异表达基因的过程叫做差异表达分析（Differential Expression Analysis）。&lt;br/&gt;最近的研究表明，基因的表达在不同的个体间存在生物学可变性&references_and_link{(Elowitz MB, Levine AJ, Siggia ED, Swain PS. Stochastic gene expression in a single cell. Science 2002; 297:1183–1186.),(http://www.sciencemag.org/content/297/5584/1183)} &references_and_link{(Kasper D. Hansen, Zhijin Wu, et al. Sequencing technology does not eliminate biological variability. Nat Biotech 2011, pp. 572-573, doi:10.1038/nbt.1910),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137276/)}（Biological Variability），而转录组测序技术（甚至qPCR、生物芯片等技术）不能消除这种可变性。可变性是基因表达的基本特征之一，不同的基因之间表达的可变程度存在差异。为了寻找真正感兴趣的差异表达基因，需要考虑和处理因生物学可变性造成的表达差异&references_and_link{(Robasky, K., Lewis, N. E. Church, G. M. The role of replicates for error mitigation in next-generation sequencing. Nature Reviews Genetics. 1–7 2013. doi:10.1038/nrg3655),(http://www.nature.com/nrg/journal/v15/n1/full/nrg3655.html)}。目前最常用和最有效的方法是在实验设计中设立生物学重复（Biological Replicates），即在同一条件下（最好能取材于相同的个体）制备多个生物学样品。重复条件限制越严格，重复样品数目越多，寻找到的差异表达基因越可靠。&lt;br/&gt;差异表达分析寻找到的基因集合叫做差异表达基因集。下文和分析结果中，使用&quot;A_vs_B&quot;的方式命名差异表达基因集，如T1_vs_T2或T1_T2_vs_T3_T4等。通常情况下，对于两个样品之间的差异表达基因集，A表示对照样品、野生型样品或前一个时间点样品；而B表示对应的处理样品、突变型样品或后一个时间点样品。相应地，对于两个条件（即两组样品）之间的差异表达基因集，A表达含有多个重复样品（Duplicates）的对照组、野生型组或前一个时间点样品组；B表示对应的处理组、突变型组、后一个时间点样品组。根据两（组）样品之间表达水平的相对高低，差异表达基因可以划分为上调基因（Up-regulated Gene）和下调基因（Down-regulated Gene）。上调基因在样品（组）B中的表达水平高于样品（组）A中的表达水平；反之为下调基因。因此，上调和下调是相对的，由所给A和B的顺序决定，更换A和B的顺序之后会完全反过来，但这不会对分析结果产生实质性的影响。");

#$out .=&info("desc","基因表达具有时间和空间特异性，外界刺激和内部环境都会影响基因的表达。在不同条件（如对照与处理、野生型和突变型、不同时间点、不同组织等）下，表达水平存在显著差异的基因，称之为差异表达基因（Differentially Expressed Gene，DEG）。同样地，表达水平存在显著差异的转录本，称之为差异表达转录本（Differentially Expressed Transcript，DET）。生物信息学中，寻找差异表达转录本或差异表达基因的过程叫做差异表达分析（Differential Expression Analysis）。&lt;br/&gt;最近的研究表明，基因的表达在不同的个体间存在生物学可变性&references_and_link{(Elowitz MB, Levine AJ, Siggia ED, Swain PS. Stochastic gene expression in a single cell. Science 2002; 297:1183–1186.),(http://www.sciencemag.org/content/297/5584/1183)} &references_and_link{(Kasper D. Hansen, Zhijin Wu, et al. Sequencing technology does not eliminate biological variability. Nat Biotech 2011, pp. 572-573, doi:10.1038/nbt.1910),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137276/)}（Biological Variability），而转录组测序技术（甚至qPCR、生物芯片等技术）不能消除这种可变性。可变性是基因表达的基本特征之一，不同的基因之间表达的可变程度存在差异。为了寻找真正感兴趣的差异表达基因，需要考虑和处理因生物学可变性造成的表达差异&references_and_link{(The role of replicates for error mitigation in next-generation sequencing.),(http://www.nature.com/nrg/journal/v15/n1/full/nrg3655.html)}。目前最常用和最有效的方法是在实验设计中设立生物学重复（Biological Replicates），即在同一条件下（最好能取材于相同的个体）制备多个生物学样品。重复条件限制越严格，重复样品数目越多，寻找到的差异表达基因越可靠。&lt;br/&gt;差异表达分析寻找到的基因集合叫做差异表达基因集。下文和分析结果中，使用&quot;A_vs_B&quot;的方式命名差异表达基因集，如T1_vs_T2或T1_T2_vs_T3_T4等。通常情况下，对于两个样品之间的差异表达基因集，A表示对照样品、野生型样品或前一个时间点样品；而B表示对应的处理样品、突变型样品或后一个时间点样品。相应地，对于两个条件（即两组样品）之间的差异表达基因集，A表达含有多个重复样品（Duplicates）的对照组、野生型组或前一个时间点样品组；B表示对应的处理组、突变型组、后一个时间点样品组。根据两（组）样品之间表达水平的相对高低，差异表达基因可以划分为上调基因（Up-regulated Gene）和下调基因（Down-regulated Gene）。上调基因在样品（组）B中的表达水平高于样品（组）A中的表达水平；反之为下调基因。因此，上调和下调是相对的，由所给A和B的顺序决定，更换A和B的顺序之后会完全反过来，但这不会对分析结果产生实质性的影响。");

$out .='</menu-2>'."\n";

if ($bio_rep == 1) {
	$out .='<menu-2>'."\n";
	$out .=&info("name","重复相关性评估");
	$out .=&info("desc","对于设立生物学重复的项目，评估生物学重复的相关性对于分析转录组测序数据非常重要。首先，生物学重复的相关性可以检验生物学实验操作的可重复性；其次，生物学重复的相关性可以评估差异表达基因的可靠性。最后，生物学重复的相关性可以辅助异常样品的筛查。&lt;br/&gt;将皮尔逊相关系数r（Pearson&apos;s Correlation Coefficient）作为生物学重复相关性的评估指标&references_and_link{(Schulze S K, Kanwar R, G?lzenleuchter M, et al. SERE: Single-parameter quality control and sample comparison for RNA-Seq. BMC genomics, 2012, 13(1): 524.),(http://www.biomedcentral.com/1471-2164/13/524/)}。r^2越接近1，说明两个重复样品相关性越强。&lt;br/&gt;为了使分析更加准确可信，原则上，同一条件的生物学重复样品数目不得少于3个。这是为了能够筛查异常样品，并且在剔除一个异常样品后，保证每个条件至少还有2个生物学重复样品，进而直接进行后续差异表达分析，减小整个实验失败的风险，提高效率。同时，百迈客也保证对同一条件的所有生物学重复样品进行同人同批样品提取、建库，同Run同Lane测序。对异常样品进行详细分析，并根据分析结果与沟通共识决定重新进行实验，还是剔除异常样品进行后续分析。&lt;br/&gt;该项目同一条件任意一对生物学重复样品的r^2统计如表${tab_num}：");
	$out .=&examplefile("表$tab_num 生物学重复相关性统计表","暂无","geneExpression/free_com.stat",'',"注：Sample 1表示样品1的编号；Sample 2表示样品2的编号；R^2表示皮尔逊相关系数的平方。","100k","xls","detail");
	$out .=&file("图$pic_num 样品相关性热图","暂无","geneExpression/sample_cluster.png",'',"注：该图反映的是两两样品间基因表达量相关性强弱，图中左边与底部标识为样品名，右边与上部为样品聚类关系图。图中第个方框的颜色代表对应的两个样品的皮尔逊相关系数的平方值的大小，即从红色到绿色渐变的过程中，R2是越来越大的，即相关程度越来越高。","100k","png","png");
    $pic_num++;
	if (-d "$id/DEG_Analysis/corr") {
		my $corr = (glob "$id/DEG_Analysis/corr/corr_result*cor.png")[0];
		$corr =~ s/.*corr_result_(.*)_cor.png/$1/;
		my ($yp1,$yp2) = split '_',$corr;
		$out .=&file("图$pic_num ${yp1}、${yp2}两样品的基因表达量散点图","暂无","DEG_Analysis/corr/corr_result_${corr}_cor.png","对同一条件的每一对生物学重复样品的基因表达量做相关性散点图，样品${yp1}和${yp2}的相关性散点图见图${pic_num}：","注：基因表达量散点图中每个点代表一个基因，横纵坐标分别是该基因在两个样品中的表达量（FPKM）的对数。点越偏离对角线，说明对应基因在两个样品间的表达水平差异越大。另外，偏离对角线的点越多，说明两样品表达量的相关性越低，表达量差异越大；反之亦然。","100k","png","png");
		$pic_num++;
	}
	$out .='</menu-2>'."\n";
	$tab_num++;
}

$out .='<menu-2>'."\n";
$out .=&info("name","差异表达筛选");
$out .=&info("desc","检测差异表达基因时，需要根据实际情况选取合适的差异表达分析软件。对于有生物学重复的实验，采用DESeq&references_and_link{(Anders S, Huber W. Differential expression analysis for sequence count data. Genome Biology. 2010. 11:R106.),(http://genomebiology.com/2010/11/10/R106)}进行样品组间的差异表达分析，获得两个条件之间的差异表达基因集；对于没有生物学重复的实验，则使用EBSeq&references_and_link{(Leng N, Dawson JA, Thomson JA, et al. EBSeq: An empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics. 2013. 29:1035-43.),(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624807/)}进行差异表达分析，获得两个样品之间的差异表达基因集。&lt;br/&gt;在差异表达分析过程中采用了公认有效的Benjamini-Hochberg方法对原有假设检验得到的显著性p值（p-value）进行校正，并最终采用校正后的p值，即FDR（False Discovery Rate）作为差异表达基因筛选的关键指标，以降低对大量基因的表达值进行独立的统计假设检验带来的假阳性。&lt;br/&gt;在筛选过程中，默认将FDR&lt;0.01且差异倍数（Fold Change）≥2作为筛选标准。其中，Fold Change表示两样品（组）间表达量的比值。&lt;br/&gt;差异表达基因部分结果见原始数据文件及更多操作。");
$out .=&examplefile("表$tab_num 差异表达基因数目统计表","暂无","DEG_Analysis/DEG.stat",'',"注：DEG Set：差异表达基因集名称；All DEG：差异表达基因数目；up-regulated：上调基因的数目；down-regulated：下调基因的数目。","100k","xls","stat");
$tab_num++;
$out .=&file("图$pic_num ${deg}差异表达基因火山图","暂无","DEG_Analysis/$dir[0]/$dir[0].FC_FDR.png","火山图（Volcano Plot）可以直观地展现所有基因的FDR与Fold Change之间的关系，以便快速查看基因在两组样品间的表达水平差异程度及其统计学显著性。&lt;br/&gt;${deg}差异表达火山图见图${pic_num}：","注：差异表达火山图中的每一个点表示一个基因，横坐标表示某一个基因在两样品中表达量差异倍数的对数值，其绝对值越大，说明表达量在两样品间的表达量倍数差异越大；纵坐标表示错误发现率的负对数值，其值越大，表明差异表达越显著，筛选得到的差异表达基因越可靠。图中绿色的点代表有显著性表达差异的基因，红色的点代表无显著性表达差异的基因。","100k","png","png");
$pic_num++;
$out .=&file("图$pic_num ${deg}差异表达基因MA图","暂无","DEG_Analysis/$dir[0]/$dir[0].FC_count.png","通过MA图可以直观地查看基因的两组样品的表达丰度和差异倍数的整体分布。&lt;br/&gt;基因的差异表达MA图见图${pic_num}：","差异表达基因MA图中每一个点代表一个基因。横坐标为A值：log2(FPKM)，即两样品中表达量均值的对数值；纵坐标为M值：log2(FC)，即两样品间基因表达量差异倍数的对数值，用于衡量表达量差异的大小。图中绿色的点代表显著差异表达的基因，红色的点代表表达差异不显著的基因。","100k","png","png");
$pic_num++;
$out .='</menu-2>'."\n";

if (-f "$id/DEG_Analysis/All_DEG_veen.png") {
	my $deg_num = @dir;
	$out .='<menu-2>'."\n";
	$out .=&info("name","差异表达基因集维恩图");
	$out .=&info("desc","当差异表达基因集在2个到5个之间时，可以对各基因集进行统计，绘制维恩图，直观展现出各个基因集共有的差异表达基因，及特有的差异表达基因。");
	$out .=&file("图$pic_num 差异表达基因集维恩图","venn_$deg_num\_class_draw","DEG_Analysis/All_DEG_veen.png",'',"每个圆形区域代表一个差异表达基因集，重叠区域中的数字即为不同集合共有元素个数。","100k","png","png");
	$out .='</menu-2>'."\n";
	$pic_num++;
}

$out .='<menu-2>'."\n";
$out .=&info("name","差异表达基因聚类分析");
$out .=&info("desc","对筛选出的差异表达基因做层次聚类分析，将具有相同或相似表达行为的基因进行聚类，用于展示不同实验条件下基因集的差异表达模式。差异表达基因聚类结果如图${pic_num}：");
#$out .=&file("图$pic_num ${deg}差异表达基因聚类图","expression_cluster_heatmap_draw","DEG_Analysis/$dir[0]/DEG_Cluster/hierarchical/$dir[0].DEG.cluster.png",'',"注：图中不同的列代表不同的样品，不同的行代表不同的基因。颜色代表了基因在样品中的表达量FPKM以2为底的对数值。","100k","png","png");
$out .=&file("图$pic_num ${deg}差异表达基因聚类图","expression_cluster_heatmap_draw","DEG_Analysis/$dir[0]/$dir[0].DEG.cluster.png",'',"注：图中不同的列代表不同的样品，不同的行代表不同的基因。颜色代表了基因在样品中的表达量FPKM以2为底的对数值。","100k","png","png");
$pic_num++;
$out .='</menu-2>'."\n";

$out .='<ori_data>'."\n";
foreach (@dir) {
	$out .=&examplefile("差异表达基因结果文件","","DEG_Analysis/$_/$_.DEG_final.xls","","注：GeneID：基因编号；FDR：错误发现率；log2FC：表达量差异倍数的对数值；regulated：上调基因（up）还是下调基因（down）；其它列为对应样品中基因的表达量FPKM值。","100k","xls","list");
#	$out .=&examplefile("差异表达分析聚类分析atr文件","","DEG_Analysis/$_/DEG_Cluster/hierarchical/$_.DEG.cluster.atr","","层次聚类样品聚类树文件，给出了每个样品间的相关系数。","100k","txt","list");
#	$out .=&examplefile("差异表达分析聚类分析cdt文件","","DEG_Analysis/$_/DEG_Cluster/hierarchical/$_.DEG.cluster.cdt","","层次聚类主体文件，内含每个基因在每个样品中的表达值","100k","txt","list");
#	$out .=&examplefile("差异表达分析聚类分析gtr文件","","DEG_Analysis/$_/DEG_Cluster/hierarchical/$_.DEG.cluster.gtr","","层次聚类基因聚类树文件，给出了基因聚类的情况。","100k","txt","list");
}
if (-d "$id/DEG_Analysis/All_DEG") {
#	$out .=&examplefile("差异表达分析聚类分析atr文件","","DEG_Analysis/All_DEG/DEG_Cluster/hierarchical/All_DEG.DEG.cluster.atr","","层次聚类样品聚类树文件，给出了每个样品间的相关系数。","100k","txt","list");
#	$out .=&examplefile("差异表达分析聚类分析cdt文件","","DEG_Analysis/All_DEG/DEG_Cluster/hierarchical/All_DEG.DEG.cluster.cdt","","层次聚类主体文件，内含每个基因在每个样品中的表达值","100k","txt","list");
	#$out .=&examplefile("差异表达分析聚类分析gtr文件","","DEG_Analysis/All_DEG/DEG_Cluster/hierarchical/All_DEG.DEG.cluster.gtr","","层次聚类基因聚类树文件，给出了基因聚类的情况。","100k","txt","list");
#	$out .=&examplefile("差异表达分析聚类分析gtr文件","","DEG_Analysis/All_DEG/All_DEG.DEG.cluster.gtr","","层次聚类基因聚类树文件，给出了基因聚类的情况。","100k","txt","list");
}
$out .='</ori_data>'."\n";

$out .='</menu-1>'."\n";




#############################################差异表达基因功能注释和富集分析
$out .='<menu-1>'."\n";
$out .=&info("name","差异表达基因功能注释和富集分析");

$out .='<menu-2>'."\n";
$out .=&info("name","差异表达基因注释统计");
$out .=&info("desc","对差异表达基因进行功能注释，各差异表达基因集注释到的基因数量统计见表${tab_num}：");
$out .=&examplefile("表$tab_num 注释的差异表达基因数量统计表","暂无","DEG_Analysis/DEG.anno.stat","","注：DEG Set：差异表达基因集名称；Annotated：注释到的差异表达基因数目；第三列到最后一列表示各功能数据库注释到的差异表达基因数目。","100k","xls","stat");
$tab_num++;
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","差异表达基因GO功能富集");
$out .=&info("desc","GO数据库是一个结构化的标准生物学注释系统，建立了基因及其产物知识的标准词汇体系，其信息适用于各物种。该数据库结构分为多个层级，层级越低，Term所代表的功能越具体。&lt;br/&gt;${deg}差异表达基因以及所有检测的基因在GO二级Term的注释结果见图${pic_num}：");
$out .=&file("图$pic_num ${deg}差异表达基因及所有基因的GO二级Term注释图","暂无","DEG_Analysis/$dir[0]/go_enrichment/$dir[0].GO.png",'',"注：横坐标为GO三大分类下的二级Term。纵坐标表示注释到该Term的基因数目及占所有基因数目的百分比，其中红色柱体表示所有基因的注释情况，蓝色柱体表示差异表达基因的主视情况。从图${pic_num}可以看出，GO各功能在差异表达基因和所有基因两个背景下的地位，其中红色柱体与蓝色柱体具有明显差异的Term可能与差异有关。","100k","png","png");
$pic_num++;
$out .=&file("图$pic_num ${deg}topGO有向无环图","暂无","DEG_Analysis/$dir[0]/Graph/$dir[0].topGO_MF.png",'然后，利用topGO软件对注释到GO数据库的样品组间差异表达基因进行富集分析，并对显著富集的节点在GO体系中的层级关系以有向无环图的形式进行直观展示。在有向无环图中，箭头代表包含关系，即该节点的所有基因同样注释到其上级节点中。&lt;br/&gt;差异表达基因利用topGO进行功能富集的分子功能的有向无环图如下图：',"注：对每个GO节点进行富集，最显著的10个节点在图中用方框表示，图中还包含其各层对应关系。每个方框（或椭圆）内给出了该GO节点的内容描述和富集显著性值。不同颜色代表不同的富集显著性，颜色越深，显著性越高。","100k","png","png");
$pic_num++;
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","差异表达基因COG分类");
$out .=&info("desc","COG（Cluster of Orthologous Groups of proteins）数据库是基于细菌、藻类、真核生物的系统进化关系构建得到的，利用COG数据库可以对基因产物进行直系同源分类。&lt;br/&gt;${deg}差异表达基因COG分类统计结果见图${pic_num}：");
$out .=&file("图$pic_num ${deg}差异表达基因COG注释分类统计图","暂无","DEG_Analysis/$dir[0]/Cog_Anno/$dir[0].Cog.classfy.png",'',"注：横坐标为COG各分类内容，纵坐标为基因数目。在不同的功能类中，基因所占多少反映对应时期和环境下代谢或者生理偏向等内容，可以结合研究对象在各个功能类的分布做出科学的解释。","100k","png","png");
$pic_num++;
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","差异表达基因KEGG注释");
$out .=&info("desc","在生物体内，不同的基因产物相互协调来行使生物学功能，进行差异表达基因的Pathway注释分析有助于进一步解读基因的功能。KEGG（Kyoto Encyclopedia of Genes and Genomes）数据库是关于Pathway的主要公共数据库。");

$out .=&file("图$pic_num ${deg}差异表达基因KEGG分类图","暂无","DEG_Analysis/$dir[0]/pathway/kegg_enrichment/$dir[0].KEGG.png","对差异表达基因KEGG的注释结果按照KEGG中通路类型进行分类，分类图见图${pic_num}：","注：横坐标为COG各分类内容，纵坐标为基因数目。在不同的功能类中，基因所占多少反映对应时期和环境下代谢或者生理偏向等内容，可以结合研究对象在各个功能类的分布做出科学的解释。","100k","png","png");
$pic_num++;

#$out .='<diagram>'."\n";
foreach (@dir){
	my @kegg_map = (glob "$id/DEG_Analysis/$_/pathway/kegg_map/*.png");
	my $kegg_map_num = int(@kegg_map/2);
	my $kegg_map =$kegg_map[$kegg_map_num];
	$kegg_map =~ s/.*\/(DEG_Analysis.*)/$1/;
	my $deg_kegg = &deg($_);
	$out .=&file("图$pic_num ${deg_kegg}间差异表达基因的KEGG通路注释图","暂无","$kegg_map",'',"注：相对于对照组来说，红色框标记的酶与上调基因有关，绿色框标记的酶与下调基因有关。蓝色框标记的酶与上调和下调基因均有关，框内的数字代表酶的编号（EC number），而整个通路由多种酶催化的复杂生化反应构成，此通路图中与差异表达基因相关的酶均用不同的颜色标出，研究人员可以根据自己的研究对象间的差异，重点研究某些代谢通路相关基因的差异表达情况，通过通路解释表型差异的根源。","100k","png","png");
	$pic_num++;
    last;
}
#$out .='</diagram>'."\n";
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","差异表达基因KEGG通路富集分析");
$out .=&info("desc","分析差异表达基因在某一通路上是否过出现（over-presentation）即为差异表达基因的Pathway富集分析。利用富集因子（Enrichment Factor）分析Pathway的富集程度，并利用Fisher精确检验方法计算富集显著性。其中富集因子的计算公式如下：&lt;br/&gt;&lt;img src=&quot;/report/report/showjobimgByPath?filepath=$od/biocloud/template/enrichment_factor_formula.png&quot; /&gt;&lt;br/&gt;&lt;br/&gt;差异表达基因的KEGG通路富集分析结果见图${pic_num}：");
$out .=&file("图$pic_num ${deg}差异表达基因KEGG通路富集散点图","暂无","DEG_Analysis/$dir[0]/Graph/$dir[0].KEGG.Phase.png",'',"注：图中每一个图形表示一个KEGG通路，通路名称见右侧图例。横坐标为富集因子（Enrichment Factor），表示注释到某通路的差异表达基因数目在所有注释到该通路的基因总数目中所占的比例。富集因子越大，表示差异表达基因在该通路中的富集水平越显著。纵坐标为Q值的对数值，其中Q值为多重假设检验较正之后的P值。因此，纵坐标越大，表示差异表达基因在该通路中的富集显著性越可靠。","100k","png","png");
$out .='</menu-2>'."\n";
$pic_num++;

$out .='<ori_data>'."\n";
$out .=&examplefile("Unigenes注释结果文件","","Unigene/Unigene_Anno/Integrated_Function.annotation.xls","","GeneID：Unigene名；COG_class：COG数据库中的蛋白功能分类编码；COG_class_annotation：COG数据库中具体的蛋白分类注释；GO_annotation：GO数据库具体功能注释；KEGG_annotation：KEGG数据库具体功能注释；Swissprot_annotation：Swiss-Prot具体功能注释；nr_annotation：nr具体注释结果","100k","xls","list");
foreach  (@dir) {
	$out .=&examplefile("${_}差异表达基因注释结果文件","","DEG_Analysis/$_/$_.annotation.xls","","ID：UnigeneID信息；FDR：false discovery rate，即假阳性比率；log2FC：表达量差异倍数的对数值；regulated：基因在实验组样品中的上调、下调信息；COG_class：COG数据库中的蛋白功能分类编码；COG_class_annotation：COG数据库中具体的蛋白分类注释；GO_annotation：GO数据库具体功能注释；KEGG_annotation：KEGG数据库具体功能注释；Swissprot_annotation：Swiss-Prot具体功能注释；nr_annotation：nr具体注释结果。","100k","xls","list");
	$out .=&examplefile("${_}topGO生物学过程富集结果文件","","DEG_Analysis/$_/Graph/$_.topGO_BP.xls","","Biological Process（生物学过程） topGO富集结果文件。GO.ID：GO功能注释编号；Term：具体功能注释；Annotated：所有基因中有此功能的基因数；Significant：差异基因中有此功能的基因数；Expected：差异基因中有此功能的期望基因数；KS：富集显著性KS值","100k","xls","list");
	$out .=&examplefile("${_}topGO细胞组分富集结果文件","","DEG_Analysis/$_/Graph/$_.topGO_CC.xls","","Cellular Component（细胞组分） topGO富集结果文件。GO.ID：GO功能注释编号；Term：具体功能注释；Annotated：所有基因中有此功能的基因数；Significant：差异基因中有此功能的基因数；Expected：差异基因中有此功能的期望基因数；KS：富集显著性KS值","100k","xls","list");
	$out .=&examplefile("${_}topGO分子功能富集结果文件","","DEG_Analysis/$_/Graph/$_.topGO_MF.xls","","Molecular Function（分子功能） topGO富集结果文件。GO.ID：GO功能注释编号；Term：具体功能注释；Annotated：所有基因中有此功能的基因数；Significant：差异基因中有此功能的基因数；Expected：差异基因中有此功能的期望基因数；KS：富集显著性KS值","100k","xls","list");
	$out .=&examplefile("${_}GO列表文件","","DEG_Analysis/$_/go_enrichment/$_.GO.list.txt","","第一列为Unigene名，其后每列都为此Unigene注释上的GO号","100k","txt","list");
	if (-f "$id/DEG_Analysis/$_/Cog_Anno/$_.Cog_class.txt") {
		$out .=&examplefile("${_}差异基因cog注释","","DEG_Analysis/$_/Cog_Anno/$_.Cog_class.txt","","Gene name：Unigene名；Portein name in COG：COG数据库中的蛋白名；E_value：blast比对E值；Identity：blast比对相似性值；Score：blast比对打分；Organism：物种名缩写；COG id：COG数据库的ID号；COG class defination：COG数据库注释分类；Function code：COG数据库功能编号；Functional categories：COG数据库注释功能类别；Function class defination：COG数据库功能分类","100k","xls","list");
	}
}

$out .='</ori_data>'."\n";

$out .='</menu-1>'."\n";

################差异基因蛋白互作
$out .='<menu-1>'."\n";
$out .=&info("name","差异表达基因蛋白互作网络");
$out .=&info("desc","STRING&references_and_link{(Franceschini A, Szklarczyk D, Frankild S, et al. STRING v9. 1: protein-protein interaction networks, with increased coverage and integration. Nucleic acids research. 2013, 41: D808-D815.),(http://nar.oxfordjournals.org/content/41/D1/D808.long)}是收录多个物种预测的和实验验证的蛋白质-蛋白质互作的数据库，包括直接的物理互作和间接的功能相关。结合差异表达分析结果和数据库收录的互作关系对，构建差异表达基因互作网络。&lt;br/&gt;对于数据库中包含的物种，可直接从数据库中提取出目标基因集的互作关系对构建互作网络；对于数据库中未收录信息的物种，使用BLAST软件，将目的基因与数据库中的蛋白质进行序列比对，寻找同源蛋白，根据同源蛋白的互作关系对构建互作网络。构建完成的蛋白质互作网络可导入Cytoscape&references_and_link{(Shannon P, Markiel A, Ozier O, et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome research. 2003, 13(11): 2498-2504.),(http://genome.cshlp.org/content/13/11/2498.long)}软件进行可视化。&lt;br/&gt;Cytoscape可视化的差异表达基因蛋白质互作网络如图${pic_num}：：");
$out .=&file("图$pic_num 差异表达基因蛋白质互作网络图","暂无","biocloud/template/P26_pp_network.png",'',"注：图中的节点为蛋白质，边为互作关系。互作网络中节点(node)的大小与此节点的度(degree)成正比，即与此节点相连的边越多，它的度越大，节点也就越大。节点的颜色与此节点的聚集系数(clustering coefficient)相关，颜色梯度由绿到红对应聚集系数的值由低到高；聚集系数表示此节点的邻接点之间的连通性好坏，聚集系数值越高表示此节点的邻接点之间的连通性越好。边(edge)的宽度表示此边连接的两个节点间的互相作用的关系强弱，互相作用的关系越强，边越宽。没有的组合代表没有互作关系。","100k","png","png");
$out .='</menu-1>'."\n";
#################################只有一个样品时无差异表达分析及差异表达基因功能注释和富集分析两模块
}
#################################只有一个样品时无差异表达分析及差异表达基因功能注释和富集分析两模块

$out .= '</results>'."\n";


$out .= '<report_method>'."\n";
$out .='<menu-1>'."\n";
$out .=&info("name","实验流程");

$out .=&info("desc","转录组测序实验流程包括样品检测、文库构建及其质量控制和上机测序。实验流程见下图：");
$out .='<menu-2>'."\n";  #################
$out .=&examplefile("图$pic_num 转录组测序实验流程图","暂无","biocloud/template/RNA-Seq_experimental_workflow.png",'',"","100k","png","png");
$out .='</menu-2>'."\n";################
$pic_num++;

$out .='<menu-2>'."\n";
$out .=&info("name","RNA样品检测");
$out .=&info("desc","分别采用Nanodrop、Qubit 2.0、Aglient 2100方法检测RNA样品的纯度、浓度和完整性等，以保证使用合格的样品进行转录组测序。");
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","RNA文库构建");
$out .=&info("desc","样品检测合格后，则启动文库构建，主要流程如下：&lt;br/&gt;(1) 用带有Oligo（dT）的磁珠富集真核生物mRNA（若为原核生物，则通过试剂盒去除rRNA富集mRNA）。&lt;br/&gt;(2) 加入fragmentation buffer将mRNA进行随机打断。&lt;br/&gt;(3) 以mRNA为模板，用六碱基随机引物（random hexamers）合成第一条cDNA链，然后加入缓冲液、dNTPs、RNase H和DNA polymerase I合成第二条cDNA链，利用AMPure XP beads纯化cDNA。&lt;br/&gt;(4) 纯化的双链cDNA再进行末端修复、加A尾并连接测序接头，然后用AMPure XP beads进行片段大小选择，&lt;br/&gt;(5) 最后通过PCR富集得到cDNA文库。");
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","文库质控");
$out .=&info("desc","文库构建完成后，分别使用Qubit2.0和Agilent 2100对文库的浓度和插入片段大小（Insert Size）进行检测，使用Q-PCR方法对文库的有效浓度进行准确定量，以保证文库质量。");
$out .='</menu-2>'."\n";

$out .='<menu-2>'."\n";
$out .=&info("name","上机测序");
$out .=&info("desc","库检合格后，用HiSeq2500进行高通量测序，测序读长为PE125。");
$out .='</menu-2>'."\n";
$out .='</menu-1>'."\n";

$out .='<menu-1>'."\n";
$out .=&info("name","生物信息学分析");
$out .=&info("desc","对Raw Data进行测序质量控制，去除其中的低质量reads和rRNA reads获得高质量的Clean Data数据。&lt;br/&gt;将Clean Data进行序列组装，获得该物种的Unigene库。基于此，可以进行随机性检验、饱和度检验等测序文库质量评估。文库质量评估合格后，进行表达量分析、基因结构分析，并根据基因在不同样品或不同样品组中的表达量进行差异表达分析、差异表达基因功能注释和功能富集等分析。&lt;br/&gt;无参考基因组的转录组生物信息分析流程见图${pic_num}：");
$out .='<menu-2>'."\n";  #################
$out .=&examplefile("图$pic_num 转录组测序生物信息分析流程图","暂无","biocloud/template/RNA-Seq_analysis_workflow.png",'',"","100k","png","png");
$out .='</menu-2>'."\n";################

$out .='</menu-1>'."\n";
$pic_num++;

$out .= '</report_method>'."\n";



$out .= '<report_attention>'."\n";
$out .= "&lt;SVG文件格式的查看&gt;结果文件中含有SVG格式的图片文件，SVG是矢量化的图片文件，可以随意放大而不失真。要查看SVG格式的文件，请先安装SVG插件。";
$out .= '</report_attention>'."\n";


$out .= '<report_ref>'."\n";

$out .= $references;

$out .= '</report_ref>'."\n";

$out .= '</report_config>'."\n";

=c
$out .=&info("name","");
$out .=&info("desc","");
$out .=&file("","","","","","","","");
=cut


#######################################################################################
open (OUT,">$od/configtest.xml") or die $!;
print OUT <<"_END_";
<?xml version="1.0" encoding="UTF-8"?>

$report_cfg
$report_pic_path
$out

_END_
close OUT;
open (OUT1,">$od/../abstract.txt") or die;
print OUT1 $abstract,"\n";
close(OUT1);

#######################################################################################
#system(qq(perl $Bin/get_gene_id.pl -inputdir $od/Web_Report -out $od/Web_Report));
if(defined($oneSample)){
	system(qq(perl $Bin/get_gene_id.pl -inputdir $id -out $od -oneSam));
}else{
	system(qq(perl $Bin/get_gene_id.pl -inputdir $id -out $od));
}

#######################################################################################
&timeLog("$Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s");
#######################################################################################

sub assem_per {
	my $contigs_len = shift;
	my $assem_per = 100*$assembly{'contigs'}{$contigs_len}/$assembly{'contigs'}{'total'};
	$assem_per = sprintf "%.2f",$assem_per;
	$assembly{'contigs'}{$contigs_len} .= "\($assem_per\%\)";
}


sub assem_t_u {
	my $type = shift;
	$_ = shift;
	if (/^(\d+)/) {
		my $num = $1;
		my @line = split /\t/;
		$line[2] =~ s/%//;
		$line[2] = sprintf "%.2f",$line[2];
		$line[2] .= "\%";
		$assembly{$type}{$num} = "$line[1]\($line[2]\)";
	}
	$assembly{$type}{'T_n'} = (split /\t/)[1] if (/^Total number/i);
	$assembly{$type}{'T_l'} = (split /\t/)[1] if (/^Total length/i);
	$assembly{$type}{'N50'} = (split /\t/)[1] if (/^N50/);
	$assembly{$type}{'Mean'} = (split /\t/)[1] if (/^Mean/);
	$assembly{$type}{'Mean'} = sprintf "%.2f",$assembly{$type}{'Mean'} if (/^Mean/);
}


sub deg {
	my $info = shift;
	my $deg1;
	my ($treated1,$control1)=(split/_vs_/,$info)[0,1];
	if ($treated1 =~ /\_/) {
		$deg1 .= "分组$treated1";
	}
	else {
		$deg1 .= "样品$treated1";
	}
	$deg1 .= "和";
	if ($control1 =~ /\_/) {
		$deg1 .= "分组$control1";
	}
	else {
		$deg1 .= "样品$control1";
	}
	return $deg1;
}


sub rpp {#&rpp("$path","$anno");
	my $min_path = shift;
	my $big_path = shift;
	my $anno = shift;
	return &info("report_pic_path",&info("min_pic_path",$min_path).&info("big_pic_path",$big_path).&info("anno",$anno));
}

sub references {
	my $in = shift;
	$references_num++;
	my $num = 1000 + $references_num;
	my $ref = "&lt;a num=&quot;testnum${num}&quot; href=&quot;#testnum${num}&quot;&gt;[$references_num]&lt;/a&gt;";
	$in =~ /&references\{\((.*?)\)\}/;
	$references .= "&lt;span id=&quot;testnum${num}&quot;&gt;[$references_num]&lt;/span&gt; $1&lt;br/&gt;";
	$in =~ s/&references\{\(.*?\)\}/$ref/;
	return $in;
}

#http://www.biomarker.com.cn

sub references_and_link {
	my $in = shift;
	$references_num++;
	my $num = 2000 + $references_num;
	my $ref = "&lt;a num=&quot;testnum${num}&quot; href=&quot;#testnum${num}&quot;&gt;[$references_num]&lt;/a&gt;";
              #<a num="testnum${num}"; href="#testnum${num}">[$references_num]</a>
	$in =~ /&references_and_link\{\((.*?)\),\((.*?)\)\}/;
	$references .= "[$references_num]&lt;/span&gt; &lt;a href=&quot;".$2."&quot; target=&quot;blank&quot;&gt; "."&lt;span id=&quot;testnum${num}&quot;&gt;$1&lt;br/&gt;"."&lt;/a&gt;";
                   #[$references_num]</span> <a href="$2" target="blank"> <span id="testnum2001">$1<br/></a>
	$in =~ s/&references_and_link\{\(.*?\),\(.*?\)\}/$ref/;
	return $in;
}

sub info{ #&info("type","12");
	my $type=shift;
	my $in=shift;
	print $type."\n" if(!defined $in);
	while ($in =~ /&references_and_link/) {
		$in = &references_and_link($in);
	}
	while ($in =~ /&references\{/) {
		$in = &references($in);
	}
	my $now='<'.$type.'>'.$in.'</'.$type.'>'."\n";
	return $now;
}

####################################################################################################图片预览有关文件
sub file{ #&file('$name','$software','$path','$desc','anno','$size','$type','$actiontype');
	my $name=shift;
	my $software=shift;
	my $path=shift;
	my $desc=shift;
	my $anno=shift;
	my $size=shift;
	my $type=shift;
	my $actiontype=shift;
	$path=~s/\.r// unless (-f "$id/$path");
	
	print "$id/$path file not exists\n" unless (-f "$id/$path");
	my $filename = basename($path);
	`convert -resize 540x180 $id/$path $cloud_dir/images/$filename` if (-f "$id/$path");
	$T_num++;
	$report_pic_path .= &rpp ("biocloud/images/$filename","$path","$name&lt;br/&gt;$anno");
	return &info("file",&info("id",$T_num).&info("name",$name).&info("software",$software).&info("path",$path).&info("desc",$desc).&info("anno",$anno).&info("size",$size).&info("type",$type).&info("actiontype",$actiontype));
}

####################################################################################################不预览文件
sub examplefile{ #&file('$name','$software','$path','$desc','anno','$size','$type','$actiontype');
	my $name=shift;
	my $software=shift;
	my $path=shift;
	my $desc=shift;
	my $anno=shift;
	my $size=shift;
	my $type=shift;
	my $actiontype=shift;
	$path=~s/_final/\.final/ unless (-f "$id/$path");
	$path=~s/All_Combination\.contigs\.fa/All_Combination\.contigs\.fa\.gz/ unless (-f "$id/$path");
	print "$id/$path file not exists\n" unless (-f "$id/$path");
	$T_num++;
	return &info("file",&info("id",$T_num).&info("name",$name).&info("software",$software).&info("path",$path).&info("desc",$desc).&info("anno",$anno).&info("size",$size).&info("type",$type).&info("actiontype",$actiontype));
}

####################################################################################################
sub LOAD_INF {
	my $para_file = shift;
	my $para = shift;
	my $error_status = 0;

	open (IN,$para_file) || die "fail open: $para_file";

	while (<IN>) {
		chomp;
		s/^\s+//; s/\s+$//; s/\r$//;
		next if (/^$/ or /^\#/);
		my ($para_key,$para_value) = split(/\s+/,$_);
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub LOAD_INF1 {
	my $para_file = shift;
	my $para = shift;
	my $error_status = 0;

    my $config=Config::General->new("$para_file");
    my%config=$config->getall;
    return %config;
	
}

####################################################################################################
#整数格式：3位一个逗号
sub Integer_Three_Digit{#
	my $interger = shift;
	$interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $interger;
}

####################################################################################################
#整数格式：3位一个逗号
#小数格式：小数点后两位
sub format_figure{#
	my $figure = shift;
	if (!defined $figure) {
		die;
	}
	if ($figure=~/\./) {
        if ($figure == 100) {
            $figure = 100;
        } else {
            $figure = sprintf("%.2f",$figure);
        }
	}else{
		$figure = Integer_Three_Digit($figure);
	}
	return $figure;
}

#sub get_inf_from_xls{
#		my $xls=shift;
#		my $pid=shift;
#	    my $parser   = Spreadsheet::ParseExcel->new();
#        my $workbook = $parser->parse($xls);
#		my %proInfo;
#        if ( !defined $workbook ) {
#            die $parser->error(), ".\n";
#        }
#
#        for my $worksheet ( $workbook->worksheets() ) {
#
#            my ( $row_min, $row_max ) = $worksheet->row_range();
#            my ( $col_min, $col_max ) = $worksheet->col_range();
#
#            for my $row ( $row_min .. $row_max ) {
#                for my $col ( $col_min .. $col_max ) {
#
#                    my $cell = $worksheet->get_cell( $row, $col );
#                    next unless $cell;
#                    chomp $cell;
#                    if($cell eq $pid){
#                    	
#                    	
#                    }
#                }
#            }
#        }
#	
#}
sub load_detail_cfg{
	my$detail_cfg=shift;
	my$out=shift;
	my%detail;
	`sed  -i 's/para_K-cov/para_K_cov/g'  $detail_cfg`;
	my@conItems=qw(projectInfo basicAnalysis geneAnn DEGAnalysis SNPAnalysis);
	if(&checkCfg($detail_cfg)){
		`cp $detail_cfg $out` if($detail_cfg ne $out);
		&checkDetailConfig($out);
		return &readConfig($detail_cfg);
	}
	$/="\n>>>>>";
	my $COM_limit=0;
	my $Check_num=0;
	my $key_word;
	open (IN,$detail_cfg) or die $!;
	open OUT,">$out" or die "$!";
	my@COM;
	while (<IN>) {
	    chomp;
	    #print STDOUT "$_\n";
	    $Check_num=$.;
		
			print OUT "<$conItems[$Check_num-1]>\n";
			print OUT "$_\n";
			print OUT "</$conItems[$Check_num-1]>\n";
		
	}

	    
	$/="\n";
	close IN;
	if ($Check_num!=5) {
	    print "Check Your detail_cfg,it should have 4 '>>>>>' \n";die;
	}
	
#	&writeConfig(\%detail,$out);
	close OUT;
	&checkDetailConfig($out);
	return &readConfig($detail_cfg);
}


sub load_data_cfg{
	my$data_cfg=shift;
	my$out=shift;
	my$Sample_name;
	my%Data;	
	if(&checkCfg($data_cfg)){
		`cp $data_cfg $out` if ($data_cfg ne $out);
		return;
	}
	open (IN,$data_cfg) or die $!;
	while (<IN>) {
	    if (/^Sample/) {
	        $Sample_name=(split/\s+/,$_)[1];
	    }
	    if ($_=~/^fq1/ || $_=~/^fq2/ || $_=~/^raw_fq1/ || $_=~/^raw_fq2/) {
	        my $file=(split/\s+/,$_)[1];
	        unless (-e $file) {
	            print "$file is not exist!\n";
	        }
	        $Data{$Sample_name}{fq1}=$file if $_=~/^fq1/;
	        $Data{$Sample_name}{fq2}=$file if $_=~/^fq2/;
	        $Data{$Sample_name}{raw_fq1}=$file if $_=~/^raw_fq1/;
	        $Data{$Sample_name}{raw_fq2}=$file if $_=~/^raw_fq2/;
	    }
	}
	close IN;
	&writeConfig(\%Data,$out);
	
	
}
sub checkCfg{
	my$file=shift;
	open CFG,"$file" or die "$!";
	while(<CFG>){
		if(/<.*>/){
			return 1;
		}
		if(/Sample/ || />>>>>/){
			return 0;
		}
		
	}
	
}
sub load_cfg{
	my$config_path=shift;
	my$od=shift;
	my%detail_cfg=(
	##'projectInfo'=>"",
    #'basicAnalysis'=>"",
    #'geneAnn'=>"",
    #'DEGAnalysis'=>"",
    #'SNPAnalysis'=>""
	);
	my%data_cfg;
	`sed  -i 's/para_K-cov/para_K_cov/g'  $config_path/basicAnalysis.cfg`;
	$detail_cfg{basicAnalysis}=&readConfig("$config_path/basicAnalysis.cfg");
	$detail_cfg{geneAnn}=&readConfig("$config_path/geneAnn.cfg");
	$detail_cfg{DEGAnalysis}=&readConfig("$config_path/DEGAnalysis.cfg");

    # deal with DEGAnalysis config contain 'Com    T1,T2,T3'
    if (defined $detail_cfg{DEGAnalysis}{Com}) {
        my %alt_com;
        if (ref($detail_cfg{DEGAnalysis}{Com}) eq 'ARRAY') {
            for my $com (@{$detail_cfg{DEGAnalysis}{Com}}) {
                my @sample = split ",",$com;

                if (@sample >2) {
                    my $iter = combinations(\@sample,2);
                    while (my $c = $iter->next) {
                        my ($boy, $girl) = @$c;
                        $alt_com{"$boy,$girl"} = 1;
                    }
                } else {
                    $alt_com{$com} = 1;
                }
            }
        } else {
            my $com = $detail_cfg{DEGAnalysis}{Com};
            my @sample = split ",",$com;

            if (@sample >2) {
                my $iter = combinations(\@sample,2);
                while (my $c = $iter->next) {
                    my ($boy, $girl) = @$c;
                    $alt_com{"$boy,$girl"} = 1;
                }
            } else {
                $alt_com{$com} = 1;
            }
        }

        delete $detail_cfg{DEGAnalysis}{Com};
        $detail_cfg{DEGAnalysis}{Com} = [keys %alt_com];
    }

	open IN,"$config_path/projectInfo.cfg" or die "$!";
	while(<IN>){
		chomp;
		next if (/#/ || /^\s+$/);
		my @tmp=split/\s+/;
		my @tmp2=map {$_=~s/\s+//g;$_} @tmp;
		$detail_cfg{'projectInfo'}{$tmp2[0]} = $tmp2[1];
#		$detail_cfg{'projectInfo'}{'Project_name'} = (split/\s+/,$_)[1] if /Project_name/;
#		$detail_cfg{'projectInfo'}{'Project_info'} = (split/\s+/,$_)[1] if /Project_info/;
#		$detail_cfg{'projectInfo'}{'Contract_NO'} = (split/\s+/,$_)[1] if /Contract_NO/;
#		$detail_cfg{'projectInfo'}{'Customer_info'} = (split/\s+/,$_)[1] if /Customer_info/;
#		$detail_cfg{'projectInfo'}{'Project'} = (split/\s+/,$_)[1] if /Project/;
	}
	close IN;
	
	my %para;
	my %sample;
	my@COM;

    if (-f "$config_path/SNPAnalysis.cfg") {
        open IN,"$config_path/SNPAnalysis.cfg";
        while (<IN>) {
            chomp;
            s/\r$//;s/^\s+//;s/\s+$//;
            next if (/^\#/ || /^$/ || /genome/ ||/Sample/||/fq/);
            my @tmp=split /\s+/,$_;
           # if ($tmp[0]=~m/KEY/) {
           #     my $fq1=<IN>;chomp $fq1;
           #     my $fq2=<IN>;chomp $fq2;
           #     my @fq_1=split /\s+/,$fq1;
           #     $sample{$tmp[1]}{fq1}=$fq_1[1];
           #     my @fq_2=split /\s+/,$fq2;
           #     $sample{$tmp[1]}{fq2}=$fq_2[1];
           # }
           # $para{$tmp[0]}=$tmp[1] unless /^COM/ or /^KEY/;
		   $para{$tmp[0]}=$tmp[1];
            #if (/^COM/) {
            #    push @COM,$tmp[1];
            #}
        }
        #$para{COM}=\@COM;
        close IN;
    }

	$detail_cfg{SNPAnalysis}=\%para;
#	&writeConfig(\%sample,"$od/data.cfg");
	&writeConfig(\%detail_cfg,"$od/detail.cfg");
	&checkDetailConfig("$od/detail.cfg");
	return \%detail_cfg;
}
sub transferDetail{
    my$configFile=shift @_;
    my$o=shift@_;
    my@conItems=qw(projectInfo basicAnalysis geneAnn DEGAnalysis SNPAnalysis);
    my $detail=Config::General->new("$configFile");
    my %detailCfg=$detail->getall;
    #print Dumper(\%detailCfg);
    #print Dumper(\%Data);
    my @projectInfo=qw(
                        Project_name
                        Project_info
                        Contract_NO
                        Customer_info
                        Project);
    my@basicAnalysis=qw(
                        Combine
                        Com_Data
                        Normalize
                        com_JM
                        max_cov
                        Separate
                        Sep_Data
                        SamG
                        Group1
                        Group1_Data
                        Group2
                        Group2_Data
                        Group3
                        Group3_Data
                        para_K_cov
                        para_K-cov
                        para_JM
                        para_thread
                        para_min_contig
                        para_seqType
                        para_lib_type
                        para_PEinsert
                        clu_identity
                        cpu
                        identity
                        overlap
                        max_len
                        );
    my @geneAnn=qw(
                        blast_cpu
                        hmmscan_cpu
                        blast_e
                        blast_cut
                        nr
                        nt
                        TrEMBL
                        Kegg
                        Swissprot
                        Pfam
                        Cog
                        Kog
                        );
    my@DEGAnalysis=qw(
                        Com
                        Sep
                        fold
                        FDR);
    my @SNPAnalysis=qw(
                        COM
                        min_ins
                        max_ins
                        max_read_len
                        score
                        ratio
                        depth_reg);
    my%orderHash=(
    'projectInfo'=>\@projectInfo,
    'basicAnalysis'=>\@basicAnalysis,
    'geneAnn'=>\@geneAnn,
    'DEGAnalysis'=>\@DEGAnalysis,
    'SNPAnalysis'=>\@SNPAnalysis);
    
    open OUT,">$o" or die "$!";
    my $i=0;
    for my $key (@conItems){
        my@aa=@{$orderHash{$key}};
        for my $conKey(@aa){
            if(exists $detailCfg{$key}{$conKey}){
                if(ref($detailCfg{$key}{$conKey}) eq 'ARRAY'){
                    for my $k(@{$detailCfg{$key}{$conKey}}){
                        print OUT "$conKey\t$k\n";
                    }
                
                }else{
                    print OUT "$conKey\t$detailCfg{$key}{$conKey}\n";
                }
                
            }
        }
         print OUT ">>>>>\n" if($i<4);
         $i++;
    }
    
    close OUT;
    return $o,\%detailCfg;
}

sub transferData{
    my$configFile=shift @_;
    my$o=shift@_;
    open IN,"$configFile" or die "$!";
    my $detail=Config::General->new("$configFile");
    my %Data=$detail->getall;
    #print Dumper(\%Data);

    open OUT,">$o" or die "$!";
    for my $sample(keys %Data){
        print OUT "Sample\t$sample\n";
        
            print OUT "fq1\t$Data{$sample}{fq1}\n" if exists $Data{$sample}{'fq1'};
            print OUT "fq2\t$Data{$sample}{fq2}\n" if exists $Data{$sample}{'fq2'};
            if (exists $Data{$sample}{'raw_fq1'}){
                print OUT "raw_fq1\t$Data{$sample}{raw_fq1}\n" ;
            
            }elsif(exists $Data{$sample}{fq1}){
                print OUT "raw_fq1\t$Data{$sample}{fq1}\n" ;
            }
            if (exists $Data{$sample}{'raw_fq2'}){
                print OUT "raw_fq2\t$Data{$sample}{raw_fq2}\n" ;
            }elsif(exists $Data{$sample}{fq2}){
                print OUT "raw_fq2\t$Data{$sample}{fq2}\n" ;
            }
            
    }
    
    
    close OUT;
    return $o,\%Data;
}

sub writeConfig{
    my($config,$outFile)=@_;
#    if(-e $outFile){
#        system(qq(mv $outFile  $outFile.bak));
#    }
    my $con=Config::General->new(-SplitPolicy=> 'custom' , -SplitDelimiter => '\s*=\s*',-StoreDelimiter=>'  ',-SaveSorted => 1);
     $con->save_file($outFile, \%$config);
}
sub readConfig{
	my$configFile=shift;
    my $d=Config::General->new(-ConfigFile => "$configFile");
    my %config=$d->getall;	
    return \%config;
}
sub writeOldDataConfig{
    my($config,$outFile)=@_;
    open OUT ,">$outFile" or die "$!";
    for my $sample (keys %$config){
        print OUT "Sample\t$sample\n";
        for my $faInfo (keys %{$config->{$sample}}){
            print OUT "$faInfo\t$$config{$sample}{$faInfo}\n";
        }
    }
    close OUT;
}

sub checkDataConfig{
	my $dataConfig=shift;
	my $analysisDir=shift;
	
	
	
}
sub checkDetailConfig{
	my $detailConfig=shift;
	#my $analysisDir=shift;
	my $Config=&readConfig($detailConfig);
	my@class=qw(
		Bacteria
		Archaea
		Plants
		Invertebrates
		Vertebrates
		Primates
		Human
		Rodents
		Viruses
		Mammals
		Fungi
		All
	);
	my@base=qw(
			nr
			nt
			Kegg
			Cog
			Kog
			Pfam
			Swissprot
			TrEMBL
	);
	my@annDataBase;
	for my $i ( @base){
		if(exists $$Config{geneAnn}{$i}){
			push @annDataBase,basename($$Config{geneAnn}{$i});
		}
	}
	#print STDOUT Dumper($dataBaseArray);
	my$dataBase=&readConfig("$Bin/dataBase.txt");
	my$dataBaseArray=&readConfig("$Bin/dataBaseArray.txt");
	for my $i (@class){
		my$diff=&compare(\@{$$dataBaseArray{$i}},\@annDataBase);
		#print STDOUT Dumper($dataBaseArray);
		#print $diff->count();
		#print $i,scalar @{$diff->added},$diff->added->[0],"\n";
		if ($diff){
			for my $j (@base){
					$$Config{geneAnn}{$j}=$$dataBase{$i}{$j} ;#unless(exists $$Config{geneAnn}{$j});
			}
			
			last;
		}
	}
	&writeConfig($Config,$detailConfig);

}
sub compare{
	my $a=shift;
	my $b=shift;
	
	my $flag=0;
	for(my $i=0;$i<@{$b};$i++){
	    for(my $j=0;$j<@{$a};$j++){
	         if($a->[$j] eq $b->[$i]){
	              $flag++;
	              last;
	          }
	    }
	}
	if($flag==@{$b}){
	   return 1;
	}else{
		return 0;
	}
}
####################################################################################################
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

####################################################################################################


