#!/usr/bin/perl -w
use strict;
no strict 'refs';
use warnings;
use Getopt::Long;
use Data::Dumper;
use autodie;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname fileparse);
my $BEGIN_TIME=time();
my $version="1.0.0";
use Encode qw(decode);

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($dir_template,$dir_data,$detail_cfg,$ptop,$oneSample);

GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$dir_data,
                "cfg:s"=>\$detail_cfg,
                "pp:s"=>\$ptop,
                "Only_One_Sam"=>\$oneSample,
				) or &USAGE;
&USAGE unless ($dir_data and $detail_cfg);

system"mkdir $dir_data/template" unless (-d "$dir_data/template");
$dir_data=AbsolutePath("dir",$dir_data);

$dir_template="$Bin/template";
`cp  -f $dir_template/* $dir_data/template`;
#`cp -r $Bin/images $dir_data `;
#`cp -r $Bin/js $dir_data `;
#`cp -r $Bin/stylesheets $dir_data `;

`rm $dir_data/template/*.txt`;
my $figure_num=1;
my $report_xml;
my %data_among_words;  ##文字部分项目数据信息
my %table_info;
my %table_file_info;
my %file_info;
my %file_all_info;
my %pic_info; 
my %pic_all_info; 
my %reference_info;
my %title_info;
my %image_info;
my %DEG;


my $ftemplate;
my $ftable;
my $sample_num;
my $abstract;
my $assembly ||="single";

my ($DEGcondition1,$DEGcondition2,$DEG_vs_name) = ('','','');  ## 差异基因目录名称
my ($br_replicate1,$br_replicate2);
my ($pic_num,$table_num,$first_title,$second_title,$third_title)=(1,1,0,0,0);

#unless (-d "$dir_data/Assembly/All_Combination") {
#    $assembly = 'multiple';
#}
my (%DEG_br_vs,%DEG_nbr_vs,$num_group);
if (-d "$dir_data/DEG_Analysis/") {
    my @DEGresult = (glob "$dir_data/DEG_Analysis/*_vs_*/*.DEG_final.xls");
	$num_group=@DEGresult;
    for my $i (0..$#DEGresult) {
        $DEGresult[$i] =~ /\/(([^\/]+)_vs_([^\/]+))\//;
        my ($vs,$A,$B) = ($1,$2,$3);
        if ($2=~/_/ or $3=~/_/) {
            $DEG_br_vs{$vs}{'A'} = $A;
            $DEG_br_vs{$vs}{'B'} = $B;
            ($br_replicate1,$br_replicate2) = (split /_/,$A)[0,1] unless $br_replicate1;
        } else {
            $DEG_nbr_vs{$vs}{'A'} = $A;
            $DEG_nbr_vs{$vs}{'B'} = $B;
            $br_replicate1 = $A unless $br_replicate1;
        }
    }

    if ((keys %DEG_br_vs) != 0) {
        $DEG_vs_name = (keys %DEG_br_vs)[0];
        $DEGcondition1 = $DEG_br_vs{$DEG_vs_name}{'A'};
        $DEGcondition2 = $DEG_br_vs{$DEG_vs_name}{'B'};
    } elsif ((keys %DEG_nbr_vs) != 0) {
        $DEG_vs_name = (keys %DEG_nbr_vs)[0];
        $DEGcondition1 = $DEG_nbr_vs{$DEG_vs_name}{'A'};
        $DEGcondition2 = $DEG_nbr_vs{$DEG_vs_name}{'B'};
    }
}
my $name=(glob "$dir_data/BMK_1_rawData//PNG/*.quality.png")[0];
my $oneSampleName=fileparse($name,".quality.png");
if($oneSample){
	$ftemplate="$dir_template/Trans_template_single.txt";
	$br_replicate1=$oneSampleName;
}
elsif($assembly eq 'single') {
    if(-f "$dir_data/BMK_1_rawData//PNG/$oneSampleName.rawDataStat.png"){
        $ftemplate = ($DEGcondition1 =~/_/ or $DEGcondition2 =~/_/) ? "$dir_template/Trans_template_sa_br.txt" : "$dir_template/Trans_template_sa_nbr.txt";
    }else{
        $ftemplate = ($DEGcondition1 =~/_/ or $DEGcondition2 =~/_/) ? "$dir_template/Trans_template_sa_br_cloud.txt" : "$dir_template/Trans_template_sa_nbr_cloud.txt";
    }
} else {
    if(-f "$dir_data/BMK_1_rawData//PNG/$oneSampleName.rawDataStat.png"){
        $ftemplate = ($DEGcondition1 =~/_/ or $DEGcondition2 =~/_/) ? "$dir_template/Trans_template_ma_br.txt" : "$dir_template/Trans_template_ma_nbr.txt";
    }else{
        $ftemplate = ($DEGcondition1 =~/_/ or $DEGcondition2 =~/_/) ? "$dir_template/Trans_template_ma_br_cloud.txt" : "$dir_template/Trans_template_ma_nbr_cloud.txt";
    }
}



print "\n[".&GetTime($BEGIN_TIME)."] $Script start ...\n";

# ------------------------------------------------------------------
# read in
# ------------------------------------------------------------------
# read detail config
my %detail_cfg;
&detail_cfg_read($detail_cfg,\%detail_cfg);
my $q30=(exists $detail_cfg{'Q30'})?$detail_cfg{'Q30'}:'85%';
my $demand=(exists $detail_cfg{'Contract_data'})?$detail_cfg{'Contract_data'}:'4';
	$demand=~s/G//ig;
my $first=(exists $detail_cfg{'First_time'})?$detail_cfg{'First_time'}:'XXXX/XX/XX';
my $second=(exists $detail_cfg{'Second_time'})?$detail_cfg{'Second_time'}:'XXXX/XX/XX';
my $third=(exists $detail_cfg{'Third_time'})?$detail_cfg{'Third_time'}:'XXXX/XX/XX';
my $seqPlatom=(exists $detail_cfg{SeqPlatom})?$detail_cfg{SeqPlatom}:"XXXX";

## picture correspond
&get_picture();

# read in table files
&get_table();

## 获得文字部分项目数据信息
&get_data_among_words();


# read in template file
$report_xml='<?xml version="1.0" encoding="utf-8"?>'."\n".'<report>'."\n\t".'<report_version value="v1.6" />'."\n\t".'<report_name value="'.$detail_cfg{Project_name}.'" />'."\n\t".'<report_code value="'.$detail_cfg{Contract_NO}.'" />'."\n\t".'<report_user value="'.$detail_cfg{Customer_name}.'" />'."\n\t".'<report_user_addr value="'.$detail_cfg{Customer_info}.'" />'."\n\t".'<report_time value="'."$first;$second;$third;".&GetDate.'" />'."\n\t".'<report_abstract value="';
#$report_xml='<?xml version="1.0" encoding="utf-8"?>'."\n".'<report>'."\n\t".'<report_version value="v1.6" />'."\n\t".'<report_name value="'.$detail_cfg{Project_name}.'" />'."\n\t".'<report_code value="'.$detail_cfg{Contract_NO}.'" />'."\n\t".'<report_user value="'.$detail_cfg{Customer_name}.'" />'."\n\t".'<report_user_addr value="'.$detail_cfg{Customer_info}.'" />'."\n\t".'<report_time value="'."$first;$second;$third;".&GetDate.'" />'."\n\t".'<report_abstract value="';

########html dir

########Fasta html
my $AssResult=(glob "$dir_data/BMK_2_Unigene_Assembly/Final_Unigene/*.Unigene.fa")[0];
my $tmpFile=$AssResult;
$tmpFile=~s/$dir_data\///;
$file_info{"组装结果序列"}=$tmpFile;

########SSR html

my $SSResult=(glob "$dir_data/BMK_5_Unigene_Structure/BMK_2_Unigene_SSR/*SSR.result.xls")[0];
$tmpFile=$SSResult;
$tmpFile=~s/$dir_data\///;
$file_info{"SSR结果和引物设计"}="$tmpFile";

######SNP html###
unless($oneSample){

	$file_info{"SNP分析结果文件"}="BMK_5_Unigene_Structure/BMK_3_SNP_Analysis/BMK_1_All_SNP/final.snp.xls";
}

######Exrepssion html
my @expression;
my $file="$dir_data/BMK_3_geneExpression/All_gene_expression.xls";
$tmpFile=$file;
$tmpFile=~s/$dir_data\///;
$file_info{"基因表达量结果文件"}=$tmpFile;

######DEG and DEG_KEGG  and DEG_GO html
unless($oneSample){
	chdir($dir_data);
	my @tmp=glob "BMK_6_DEG_Analysis/*_vs_*/BMK_3_GO_Enrichment/*/*.topGO_*.xls";
	my @tmp2;
	for my $file (@tmp){
		push @tmp2,$file if $file!~/gene/;
	}
	$file_all_info{"差异表达基因的KEGG富集结果"}=[glob "BMK_6_DEG_Analysis/*_vs_*/BMK_4_Pathway_Enrichment/KEGG_enrichment/*.KEGG.xls"];
	$file_all_info{"差异表达基因的GO富集结果"}=\@tmp2;
	$file_all_info{"差异表达分析结果"}=[glob "BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_1_Statistics_Visualization/*DEG.xls"];
}
##########open template
print "template:$ftemplate\n";
open (TEM, $ftemplate) or die $!;

$/="》";
while (<TEM>) {
	chomp;
	next if(/^$/ || /^\#/);
	if (!defined $ptop and $_=~/互作网络/){
		next;next;next;
	} 
	my (undef, $format, @context)=split /\n/,$_;
	&xml_write($format, \@context);
	
}
close(TEM);
$report_xml.='</report>'."\n";


my $picture_name;
foreach  my $pic_name (sort {$a<=>$b} keys %image_info ) {
	print "$image_info{$pic_name}->[2]\t$image_info{$pic_name}->[0]\n";
	if(defined $image_info{$pic_name}->[1]){
		$picture_name.="\t\t".'<pic name="图'.$pic_name.'" desc="'.$image_info{$pic_name}->[1].'" path="'.$image_info{$pic_name}->[0].'" />'."\n";
	}
	else{
		$picture_name.="\t\t".'<pic name="图'.$pic_name.'" desc="" path="'.$image_info{$pic_name}->[0].'" />'."\n";
	}
}

$picture_name=~s/delete//g;
$report_xml=~s/images_preview/\n$picture_name\t/g;

open (OUT1,">$dir_data/abstract.txt") or die;
print OUT1 $abstract,"\n";
close(OUT1);

open (OUT,">$dir_data/configtest_raw.xml") or die;
$report_xml=~s/$dir_data\///g;
print OUT $report_xml;
close(OUT);

#my $tmp1 = `iconv -f "UTF-8" -t "UTF-8" $dir_data/configtest_raw.xml -o $dir_data/configtest.xml`;
#my $tmp2 = `java -cp /share/nas1/cloud/test/cloud/WEB-INF/classes/:/share/nas1/cloud/test/cloud/WEB-INF/lib/*: com/bmkit/util/Xml2htmlDynamic $dir_data `;
#`rm $dir_data/configtest_raw.xml `;

#######################################################################################
print "\n[".&GetTime(time())."] $Script done.\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

#获取文字中项目数据信息
sub get_data_among_words {
	#$table_info{"组装结果统计表"}{"table"}->[5][0]="2000+";
    ########
    my $finish_date = &GetDate;

	######## 样品测序数据评估统计表
	# sample_num	base_num	min_Q30
	my $base_num=0;
	my $min_Q30=1000;
	my $min_data=100;   #Gb
	for (my $i=1; $i<@{$table_info{"样品测序数据评估统计表"}{"table"}}; $i++) {
		# base_num
		my $dnum = $table_info{"样品测序数据评估统计表"}{"table"}->[$i][3];
		$dnum =~s/,//g;
		$base_num+= $dnum/1000000000;   # convert to Gb
		$sample_num++;
		## min_Q30
		my $Q30 = $table_info{"样品测序数据评估统计表"}{"table"}->[$i][5];
		my $data = $table_info{"样品测序数据评估统计表"}{"table"}->[$i][3];
		$data =~s/,//g;
		$data = $data/1000000000;   # convert to Gb
		$Q30=~s/%$//;
		$min_Q30 = $Q30 if ($Q30 < $min_Q30);
		$min_data = $data if ($data < $min_data);
	}

	###SSR分析结果统计
#	my $SSR_num = $table_info{"SSR分析结果统计表"}{"table"}->[-1][1];
	my $SSR_num = $table_info{"SSR分析结果统计表"}{"number"};
	$SSR_num =~s/,//g;


	##### 组装结果统计表
	#All_Unigenes
	my $All_Unigenes = $table_info{"组装结果统计表"}{"table"}->[6][-1];
	$All_Unigenes=~s/,//g;
	#COMBIN
	my $com_trans="XXX";
	my $com_trans_N50="XXX";
	my $com_Uni_N50="XXX";
	if (-d "$dir_data/Assembly/All_Combination") {
		$com_trans=$table_info{"组装结果统计表"}{"table"}->[6][-2];
		$com_trans_N50=$table_info{"组装结果统计表"}{"table"}->[8][-2];
		$com_Uni_N50=$table_info{"组装结果统计表"}{"table"}->[8][-1];
	} else {
		$com_Uni_N50=$table_info{"组装结果统计表"}{"table"}->[8][-1];
    }


	#1k_Unigenes
	my $onek_twok_Unigene_str = $table_info{"组装结果统计表"}{"table"}->[4][-1];
	my $twok_more_Unigene_str = $table_info{"组装结果统计表"}{"table"}->[5][-1];
	my ($onek_twok_Unigene) = $onek_twok_Unigene_str=~/(\S+)\(\S+%\)/;
	my ($twok_more_Unigene) = $twok_more_Unigene_str=~/(\S+)\(\S+%\)/;
	$onek_twok_Unigene =~s/,//g;
	$twok_more_Unigene =~s/,//g;

	my $onek_Unigenes = $onek_twok_Unigene + $twok_more_Unigene;
	my $onek_percent = sprintf("%.2f",$onek_Unigenes/$All_Unigenes*100)."%";

	## $All_Anno
	my $All_Anno = $table_info{"Unigene注释统计表"}{"table"}->[-1][1];
	$All_Anno =~s/,//g;



	%data_among_words = (
		"\\\$first_time"=>$first,
		"\\\$second_time"=>$second,
		"\\\$third_time"=>$third,
        "\\\$finish_date"=>$finish_date,

		"\\\$sample_num"=>format_figure($sample_num),
		"\\\$base_num"=>format_figure($base_num)."Gb",
        "\\\$min_data"=>format_figure($min_data)."Gb",
		"\\\$min_Q30"=>$min_Q30."%",

		"\\\$SSR_num"=>format_figure($SSR_num),


		"\\\$com_trans"=>$com_trans,
		"\\\$com_trans_N50"=>$com_trans_N50,
		"\\\$com_Uni_N50"=>$com_Uni_N50,


		"\\\$All_Unigenes"=>format_figure($All_Unigenes), #"\\\$All_Unigenes" : \$All_Unigenes
		"\\\$1k_Unigenes"=>format_figure($onek_Unigenes),
		"\\\$1k_percent"=>$onek_percent,
		"\\\$All_Anno"=>format_figure($All_Anno),
		"\\\$seqPlatom"=>$seqPlatom,
		);

}

sub get_picture {#
	## get info
	my @qm_map=glob "$dir_data/BMK_1_rawData/PNG/*.quality.png";
	my @atcg_map=glob "$dir_data/BMK_1_rawData//PNG/*.acgtn.png";
	my @data_stat=glob "$dir_data/BMK_1_rawData//PNG/*.rawDataStat.png";
	my @saturation=glob "$dir_data/BMK_3_geneExpression/BMK_3_saturation/*.gene_tag.png";
	my $ssr_density_map=(glob "$dir_data/BMK_5_Unigene_Structure/BMK_2_Unigene_SSR/*.ssr.density.png")[0];
	%pic_all_info=(
		"碱基测序错误率分布图"=>[@qm_map],
		"ATCG含量分布图"=>[@atcg_map],
		"转录组测序数据饱和度模拟图"=>[@saturation],
		"Mapped Reads在mRNA上的位置分布图"=>[glob("$dir_data/BMK_3_geneExpression/BMK_1_randcheck/*.randcheck.png")],
#		"Raw Data数据分布统计图"=>[@data_stat],
	);
	if(@data_stat>0){
		$pic_all_info{"Raw Data数据分布统计图"}=[@data_stat];
	}
	#$pic_all_info{"相关性图"}=[glob("$dir_data/BMK_3_geneExpression/BMK_5_correlation/*.cor.png")] unless $oneSample;
	%pic_info=(
		"转录组测序实验流程图"=>"$dir_data/template/P01_RNA-Seq_experimental_workflow.png",
		"转录组测序生物信息分析流程图"=>"$dir_data/template/P02_RNA-Seq_analysis_workflow.png",
		"质量值计算公式"=>"$dir_data/template/F01_Qscore_formula.png",
		"FASTQ格式文件示意图"=>"$dir_data/template/P03_FASTQ_format.png",
		"Trinity组装程序原理图"=>"$dir_data/template/P05_Trinity_workflow.png",
		"Unigene长度分布图"=>(glob("$dir_data/BMK_2_Unigene_Assembly/Final_Unigene/*.Unigene.distribution.png"))[0],
        "IGV浏览器界面"=>"$dir_data/template/P07_IGV_interface.png",
		"Mapped Reads在mRNA上的位置分布图"=>"$dir_data/BMK_3_geneExpression/BMK_1_randcheck/Total.randcheck.png",
		"CDS分析结果文件示意图"=>"$dir_data/template/P10_CDS_example.png",
		"SSR密度分布图"=>"$ssr_density_map",
        'FPKM计算公式'=>"$dir_data/template/F02_FPKM_formula.png",
		"各样品FPKM密度分布对比图"=>"$dir_data/BMK_3_geneExpression/BMK_4_density/all.fpkm_density.png",
		"各样品FPKM箱线图"=>"$dir_data/BMK_3_geneExpression/BMK_4_density/all.fpkm_box.png",
		"SSR密度图"=>(glob("$dir_data/BMK_5_Unigene_Structure/BMK_2_Unigene_SSR/*.Unigene.1000.fa.ssr.density.png"))[0],
		"Unigenes的nr注释图"=>(glob("$dir_data/BMK_4_Unigene_Anno/BMK_2_statistic/*Unigene.fa.NR.lib.png"))[0],
		"Unigenes的cog注释图"=>(glob("$dir_data/BMK_4_Unigene_Anno/BMK_2_statistic/*.Unigene.fa.COG.cluster.png"))[0],
		"Unigenes的kog注释图"=>(glob("$dir_data/BMK_4_Unigene_Anno/BMK_2_statistic/*.Unigene.fa.KOG.cluster.png"))[0],
		"富集因子计算公式"=>"$dir_data/template/F03_enrichment_factor_formula.png",
		"差异表达基因的KEGG通路注释图示意图"=>"$dir_data/template/P20_pathway_example.png",
		"差异表达基因蛋白质互作网络图"=>"$dir_data/template/P26_pp_network.png",
	);
	unless($oneSample){
		$pic_info{"相关性热图"}="$dir_data/BMK_4_Unigene_Anno/BMK_2_statistic/cor.cluster.png";
		$pic_info{"差异表达基因维恩图"}="$dir_data/BMK_6_DEG_Analysis/BMK_1_All_DEG/All_DEG_veen.png" if -e "$dir_data/BMK_6_DEG_Analysis/BMK_1_All_DEG/All_DEG_veen.png";
	}
	if (-d "$dir_data/BMK_6_DEG_Analysis/" ) {
		$pic_all_info{"差异表达基因火山图"}=[glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_1_Statistics_Visualization/*.FC_FDR.png")];
		$pic_all_info{"差异表达基因MA图"}=[glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_1_Statistics_Visualization/*.FC_count.png")];
		$pic_all_info{"差异表达基因topGO富集有向无环图"}=[glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_3_GO_Enrichment/*/*.topGO_*.png")];
#		$pic_all_info{"差异表达基因topGO富集有向无环图（生物过程）"}=[glob("$dir_data/DEG_Analysis/*/Graph/*.topGO_BP.png")];
#		$pic_all_info{"差异表达基因topGO富集有向无环图（细胞组分）"}=[glob("$dir_data/DEG_Analysis/*/Graph/*.topGO_CC.png")];
		$pic_all_info{"差异表达基因COG注释分类统计图"}=[glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_2_DEG_Annotation/*.COG.classification.png")];
		$pic_all_info{"差异表达基因KOG注释分类统计图"}=[glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_2_DEG_Annotation/*.KOG.classification.png")];


		#$pic_info{"差异表达基因蛋白质互作网络示意图"}="$dir_data/template/P26_pp_network.png";

   
		$pic_all_info{"差异表达基因表达模式聚类图"} = [glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_1_Statistics_Visualization/*DEG.cluster.png")];
		$pic_all_info{"差异表达基因表达模式聚类图"} = ["$dir_data/template/none.jpg","$dir_data/template/none.jpg"] unless (exists $pic_all_info{"差异表达基因表达模式聚类图"} ) ;
		

		## 差异表达基因GO二级节点注释统计图
		$pic_all_info{"差异表达基因GO二级节点注释统计图"}=[glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_3_GO_Enrichment/*.GO.png")];
		$pic_all_info{"差异表达基因GO二级节点注释统计图"}=["$dir_data/template/none.jpg","$dir_data/template/none.jpg"] unless (exists $pic_all_info{"差异表达基因GO二级节点注释统计图"}) ;
	


		 ## 差异表达基因KEGG通路富集散点图
		$pic_all_info{"差异表达基因KEGG通路富集散点图"}=[glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_4_Pathway_Enrichment/KEGG_enrichment/*.KEGG.enrichment_factor.png")];
		$pic_all_info{"差异表达基因KEGG通路富集散点图"}=["$dir_data/template/none.jpg","$dir_data/template/none.jpg"] unless (exists $pic_all_info{"差异表达基因KEGG通路富集散点图"});


		## 差异表达基因KEGG分类图
		$pic_all_info{"差异表达基因KEGG分类图"}=[glob("$dir_data/BMK_6_DEG_Analysis/BMK_*_vs_*/BMK_4_Pathway_Enrichment/KEGG_enrichment/*.KEGG.png")];
		$pic_all_info{"差异表达基因KEGG分类图"}=["$dir_data/template/none.jpg","$dir_data/template/none.jpg"] unless (exists $pic_all_info{"差异表达基因KEGG分类图"}) ;

		#$pic_all_info{"差异表达基因的KEGG通路注释图"} = [glob("$dir_data/DEG")]
		#$pic_all_info{"差异表达基因的KEGG通路注释图"} = ["$dir_data/template/none.jpg","$dir_data/template/none.jpg"] unless (exists $pic_all_info{"差异表达基因的KEGG通路注释图"}); 
	}
		
	 ## 两样品的基因表达量散点图
		if ((keys %DEG_br_vs) != 0) {
            $pic_all_info{"两样品的基因相关性散点图"} = [glob("$dir_data/BMK_3_geneExpression/BMK_5_correlation/*.cor.png")];
			$pic_all_info{"两样品的基因相关性散点图"} =["$dir_data/template/none.jpg","$dir_data/template/none.jpg"] unless (exists $pic_all_info{"两样品的基因相关性散点图"} ) ;
		}
		 $pic_info{"样品的相关性热图"}="$dir_data/BMK_3_geneExpression/BMK_5_correlation/cor.cluster.png";
		
    ## 插入片段长度模拟分布图
    $pic_all_info{"插入片段长度模拟分布图"} = [glob("$dir_data/BMK_3_geneExpression/BMK_2_insertSize/*.insertSize.png")];
	$pic_all_info{"插入片段长度模拟分布图"} = ["$dir_data/template/none.jpg","$dir_data/template/none.jpg"] unless (exists $pic_all_info{"插入片段长度模拟分布图"});


    # SNP密度分布图
    $pic_info{"SNP密度分布图"} = (glob("$dir_data/BMK_5_Unigene_Structure/BMK_3_SNP_Analysis/BMK_1_All_SNP/AllSample.SNP_density.png"))[0];
	$pic_info{"SNP密度分布图"} = ["$dir_data/template/none.jpg","$dir_data/template/none.jpg"] unless (exists $pic_info{"SNP密度分布图"});


}

# ------------------------------------------------------------------
# get_table
# ------------------------------------------------------------------

sub get_table {#
    ## 表1 碱基质量值与碱基识别出错的概率的对应关系表
	$table_file_info{"碱基质量值与碱基识别出错的概率的对应关系表"}="$dir_data/template/base_error.xls";
    @{$table_info{"碱基质量值与碱基识别出错的概率的对应关系表"}{"width"}}= (2600,3608,2600);
    @{$table_info{"碱基质量值与碱基识别出错的概率的对应关系表"}{"align"}}= ('c','c','c');
    push @{$table_info{"碱基质量值与碱基识别出错的概率的对应关系表"}{"table"}},[("Phred Quality Score","Probability of Incorrect Base Call","Base Call Accuracy")];
    push @{$table_info{"碱基质量值与碱基识别出错的概率的对应关系表"}{"table"}},[("Q10","1/10","90%")];
    push @{$table_info{"碱基质量值与碱基识别出错的概率的对应关系表"}{"table"}},[("Q20","1/100","99%")];
    push @{$table_info{"碱基质量值与碱基识别出错的概率的对应关系表"}{"table"}},[("Q30","1/1000","99.9%")];
    push @{$table_info{"碱基质量值与碱基识别出错的概率的对应关系表"}{"table"}},[("Q40","1/10000","99.99%")];

	## 表2 样品测序数据评估统计表
	# width
	$ftable = "$dir_data/BMK_1_rawData/AllSample_GC_Q.stat";
	open (TAB,$ftable) or die $!;
	open (OUT,">$dir_data/template/AllSample_GC_Q.xls") or die $!;
	$/="\n";
	my $l1=<TAB>;
	$l1=~s/\s+$//;
	my @l1=split/\s+/,$l1;
	my $limit_1=@l1;

	if ($limit_1!=8) {
		print "Note your file: $ftable !";die;
	}
	if ($limit_1==8) {
		#@{$table_info{"样品测序数据评估统计表"}{"width"}}= (1350,1300,1600,1658,1450,1450);
		#@{$table_info{"样品测序数据评估统计表"}{"align"}}= ('c','c','c','c','c','c');

		push @{$table_info{"样品测序数据评估统计表"}{"table"}},[("Samples","BMK-ID","Read Number","Base Number","GC Content","%≥Q30")];
		print OUT "Samples\tBMK-ID\tRead Number\tBase Number\tGC Content\t",'%≥Q30',"\n";
		# info
        #SampleID       ReadSum BaseSum GC(%)   N(%)    Q20(%)  CycleQ20(%)     Q30(%)
        #T1      18723047        3781754973      52.90   0.00    97.03   100.00  92.64
		while (<TAB>) {
			chomp;
			next if(/^$/ or /^#/);
			my @info = split /\s+/,$_;
            die "Note your file: $ftable !\n" if (@info != 8);
			my @table = ($info[0],format_figure($info[1]),format_figure($info[2]),format_figure($info[3])."%",format_figure($info[7])."%");
			my $realname=(exists $detail_cfg{$info[0]})?$detail_cfg{$info[0]}:'XXX';
			unshift @table,$realname;
			print OUT join ("\t",@table),"\n";

			push @{$table_info{"样品测序数据评估统计表"}{"table"}},[@table];
		}
		close (TAB) ;
		close (OUT) ;
	}
	# my $change = `iconv -f "UTF-8" -t "UTF-8" $dir_data/template/AllSample_GC_Q.xls -o $dir_data/template/AllSample_GC_Q_final.xls`;
	#$change=` rm $dir_data/template/AllSample_GC_Q.xls `;
	$table_file_info{"样品测序数据评估统计表"}="$dir_data/template/AllSample_GC_Q.xls";

	## 表3 组装结果统计表
	if (-d "$dir_data/BMK_2_Unigene_Assembly/Final_Transcript") {

		for (my $i=0; $i<2; $i++) {
			if ($i==1) { ## Unigene
				$ftable = get_only_one_file("$dir_data/BMK_2_Unigene_Assembly/Final_Unigene/*.stat.xls");

			}
			elsif($i==0){  ## Trans
				$ftable = get_only_one_file("$dir_data/BMK_2_Unigene_Assembly/Final_Transcript/*.stat.xls");
			}
			#elsif($i==0){  ## contig
			#	$ftable = get_only_one_file("$dir_data/Assembly/All_Combination/contigs/*.stat.xls");
			#}
			open (TAB,$ftable) or die $!;
			$/="\n";
			my $line=0;
			my $mark=1;
			my $cor_limit=0;
			while (<TAB>) {
				chomp;
				next if(/^$/);
				if ($line==0) {
					$table_info{"组装结果统计表"}{"table"}->[0][0]="Length Range";
					if ($i==1) {
						$table_info{"组装结果统计表"}{"table"}->[0][$i+1]="Unigene";
					}
					elsif($i==0){
						$table_info{"组装结果统计表"}{"table"}->[0][$i+1]="Transcript";
					}
					$line++;
					next;
				}
				my @info2 = split /\t/,$_;
				$table_info{"组装结果统计表"}{"table"}->[$line][0]=$info2[0];
				if ($mark) {
					$mark=0 if($info2[0]=~/^\d+\+/);
					$info2[2] =~s/\%//;  
					$table_info{"组装结果统计表"}{"table"}->[$line][$i+1]= format_figure($info2[1])."(".format_figure($info2[2])."%)";
				}else{

					$table_info{"组装结果统计表"}{"table"}->[$line][$i+1]= format_figure($info2[1]);
				}

				$line++;
			}
			close (TAB) ;
		}

#        $table_info{"组装结果统计表"}{"table"}->[1][1] .= "*";
		open O,">$dir_data/template/Assembly.stat.xls" or die "$!";
		for(my $i=0;$i<@{$table_info{"组装结果统计表"}{"table"}};$i++){
			my $line="";
			for (my $j=0;$j<@{$table_info{"组装结果统计表"}{"table"}[$i]};$j++){
				$line.=$table_info{"组装结果统计表"}{"table"}[$i][$j]."\t";
			}
			$line=~s/\t$//;
			print O  "$line\n";
		}
		close O;
		$table_file_info{"组装结果统计表"}="$dir_data/template/Assembly.stat.xls";
	}
	else {
		## 所有样品的统计
		$ftable = get_only_one_file("$dir_data/BMK_2_Unigene_Assembly/Final_Unigene/*.stat.xls");
		## 单样品统计
		my @Assembly_Groups=glob("$dir_data/BMK_2_Unigene_Assembly/*_Assembly/Unigene/*.stat.xls");
		unshift @Assembly_Groups,$ftable;
		my $i2=0;
		for my $file (@Assembly_Groups){
			open (TAB,$file) or die $!;
			$/="\n";
			my $line2=0;
			my $name2;
			my $mark=1;
			while (<TAB>) {
				chomp;
				next if(/^$/);
				if ($line2==0) {
					if ($i2==0) {
						$table_info{"组装结果统计表"}{"table"}->[0][$i2]="All Unigenes";
					}
					else{
						($name2)=$file=~/BMK_2_Unigene_Assembly\/(.+)_Assembly/;
						$table_info{"组装结果统计表"}{"table"}->[0][$i2]=$name2." Unigenes";
					}
					$line2++;
					next;
				}
				my @info2 = split /\t/,$_;
#				if ($info3[0] ne $table_info{"组装结果统计表"}{"table"}->[$line2][0]) {
#					print "wrong assemlby statistic file: $ftable!\n";
#					die;
#				}

				$table_info{"组装结果统计表"}{"table"}->[$line2][0]=$info2[0];
				if ($mark) {
					$mark=0 if($info2[0]=~/^\d+\+/);
					$info2[2] =~s/\%//;  
					$table_info{"组装结果统计表"}{"table"}->[$line2][$i2]= format_figure($info2[1])."(".format_figure($info2[2])."%)";
				}else{

					$table_info{"组装结果统计表"}{"table"}->[$line2][$i2]= format_figure($info2[1]);
				}

				$line2++;
			}
			close (TAB) ;
			$i2++;
		}
		open O,">$dir_data/template/Assembly.stat.xls" or die "$!";
		for(my $i=0;$i<@{$table_info{"组装结果统计表"}{"table"}};$i++){
			my $line="";
			for (my $j=0;$j<@{$table_info{"组装结果统计表"}{"table"}[$i]};$j++){
				$line.=$table_info{"组装结果统计表"}{"table"}[$i][$j]."\t";
			}
			$line=~s/\t$//;
			print O  "$line\n";
		}
		close O;
		$table_file_info{"组装结果统计表"}="$dir_data/template/Assembly.stat.xls";
	}


	## 表4 测序数据与组装结果的比对统计表
	@{$table_info{"测序数据与组装结果的比对统计表"}{"width"}}=(1802,2402,2402,2202);
	@{$table_info{"测序数据与组装结果的比对统计表"}{"align"}}=('c','c','c','c');
    push @{$table_info{"测序数据与组装结果的比对统计表"}{"table"}}, [("BMK-ID","Clean Reads","Mapped Reads","Mapped Ratio")];
	my @ftable = glob "$dir_data/BMK_2_Unigene_Assembly/Map_stat/*.Mapped.stat.xls";

    for $ftable (@ftable) {
        my ($sample_id) = $ftable =~/\/([^\/]+).Mapped.stat.xls$/;
        my ($clean_num, $mapped_num, $mapped_ratio) = (0,0,'');
        open (TAB,$ftable) or die $!;
        $/="\n";
        while (<TAB>) {
            chomp;
            my @col = split /\t/;
            if ($col[0] =~/Total Reads/) {
                $clean_num = format_figure($col[1]);
            }
            if ($col[0] =~/Mapped Reads/) {
                $mapped_num = format_figure($col[1]);
                ($mapped_ratio) = $col[2] =~/([0-9\.]+)%$/;
                $mapped_ratio = format_figure($mapped_ratio)."%";
            }
        }
        push @{$table_info{"测序数据与组装结果的比对统计表"}{"table"}}, [($sample_id, $clean_num, $mapped_num, $mapped_ratio)];
        close (TAB) ;
    }
	

	open O,">$dir_data/template/Mapped.stat.xls";
	for(my $i=0;$i<@{$table_info{"测序数据与组装结果的比对统计表"}{"table"}};$i++){
			my $line="";
			for (my $j=0;$j<@{$table_info{"测序数据与组装结果的比对统计表"}{"table"}[$i]};$j++){
				$line.=$table_info{"测序数据与组装结果的比对统计表"}{"table"}[$i][$j]."\t";
			}
			$line=~s/\t$//;
			print O  "$line\n";
		}
		close O;
		$table_file_info{"测序数据与组装结果的比对统计表"}="$dir_data/template/Mapped.stat.xls";


	## 表5 Unigene注释统计表
	@{$table_info{"Unigene注释统计表"}{"width"}}=(2202,2202,2202,2202);
	@{$table_info{"Unigene注释统计表"}{"align"}}=('c','c','c','c');
	$ftable = "$dir_data/BMK_4_Unigene_Anno/Function_Annotation.stat.xls";
	open (TAB,$ftable) or die $!;
	$/="\n";
	while (<TAB>) {
		chomp;
		next if(/^$/);
		my ($head9,@info9) = split /\s+/,$_;
		if ($head9 =~/^\#/) {
			push @{$table_info{"Unigene注释统计表"}{"table"}}, [("Annotated databases","Unigene",T("≥300nt"),T("≥1000nt"))];
		}else{
			$head9 =~s/_Annotation//;
			$head9 =~s/_Annotated//;
            $head9 =~s/Swissprot/Swiss-Prot/;
			my $new_more=$info9[1]+$info9[2];
			@info9 = (format_figure($info9[0]), format_figure($new_more), format_figure($info9[2]));
			push @{$table_info{"Unigene注释统计表"}{"table"}},[$head9,@info9];
		}
	}
	close (TAB) ;
	$table_file_info{"Unigene注释统计表"}=$ftable;

    ###############################################################################
	## 表6 SSR分析结果统计表
	#width
	@{$table_info{"SSR分析结果统计表"}{"width"}}=(6306,2502);
	@{$table_info{"SSR分析结果统计表"}{"align"}}=('c','c');
	#title
	push @{$table_info{"SSR分析结果统计表"}{"table"}},["Searching Item","Number"];
	$ftable = get_only_one_file("$dir_data/BMK_5_Unigene_Structure/BMK_2_Unigene_SSR/*.stat.xls");
	$table_info{"SSR分析结果统计表"}{number}=`grep Total $ftable|cut -f 2`;
	$table_file_info{"SSR分析结果统计表"}=$ftable;
		
	## 表10 SNP数量统计表
		$table_file_info{"SNP数量统计表"}="$dir_data/BMK_5_Unigene_Structure/BMK_3_SNP_Analysis/BMK_1_All_SNP/AllSample.snp.stat";
	
	## 表12 生物学重复相关性统计表
	unless($oneSample){
		@{$table_info{"生物学重复相关性统计表"}{"width"}}=(2500,2500,3808);
		@{$table_info{"生物学重复相关性统计表"}{"align"}}=('c','c','c');
		push @{$table_info{"生物学重复相关性统计表"}{"table"}},["Sample 1","Sample 2",'r^2'];
		$ftable = "$dir_data/BMK_3_geneExpression/BMK_5_correlation/correlation.txt";
		$table_file_info{"生物学重复相关性统计表"}=$ftable;
		
		## 表13 差异表达分析结果示意表
	
		## 表14 差异表达基因数目统计表
		if ($DEGcondition1 =~ /_/) {
			@{$table_info{"差异表达基因数目统计表"}{"width"}}=(3408,1800,1800,1800);
		} else {
			@{$table_info{"差异表达基因数目统计表"}{"width"}}=(2408,2000,2200,2200);
		}
		@{$table_info{"差异表达基因数目统计表"}{"align"}}=('c','c','c','c');
		$ftable = "$dir_data/BMK_6_DEG_Analysis/BMK_1_All_DEG/All.DEG.stat.xls";
		$table_file_info{"差异表达基因数目统计表"}=$ftable;
	
	
		## 表15 注释的差异表达基因数量统计表
		$ftable = "$dir_data/BMK_6_DEG_Analysis/BMK_1_All_DEG/All.DEG.anno.stat.xls";
		$table_file_info{"注释的差异表达基因数量统计表"}=$ftable;
		my $DE_Anno_limit=0;
	
		## 表17 差异表达基因的KEGG富集部分结果

	}

    ## 附表1 软件列表
    # width
    @{$table_info{"软件列表"}{"width"}}= (1201,3002,4605);
    @{$table_info{"软件列表"}{"align"}}= ('c','w','w');
    #info
    push @{$table_info{"软件列表"}{"table"}},['Tools','Description','Linkages'];
    push @{$table_info{"软件列表"}{"table"}},['Trinity',T('转录组序列组装软件'),'http://trinityrnaseq.sourceforge.net/'];
    push @{$table_info{"软件列表"}{"table"}},['TransDecoder',T('寻找编码区序列的软件'),'http://sourceforge.net/projects/transdecoder/'];
    push @{$table_info{"软件列表"}{"table"}},['MISA',T('识别微卫星（SSR）的软件'),'http://pgrc.ipk-gatersleben.de/misa/misa.html'];
    push @{$table_info{"软件列表"}{"table"}},['STAR',T('RNA-Seq比对软件'),'http://code.google.com/p/rna-star/'];
    push @{$table_info{"软件列表"}{"table"}},['GATK',T('变体发现和分型软件'),'https://www.broadinstitute.org/gatk/'];
    if ($assembly eq 'single') {
        push @{$table_info{"软件列表"}{"table"}},['Bowtie',T('比对软件'),'http://bowtie-bio.sourceforge.net/index.shtml'];
    } else {
        push @{$table_info{"软件列表"}{"table"}},['BLAT',T('比对软件'),'http://genome.ucsc.edu/cgi-bin/hgBlat'];
    }
    push @{$table_info{"软件列表"}{"table"}},['BLAST',T('序列比对软件'),'http://blast.ncbi.nlm.nih.gov/Blast.cgi'];
    push @{$table_info{"软件列表"}{"table"}},['HMMER',T('蛋白结构域比对软件'),'http://hmmer.janelia.org/'];
    push @{$table_info{"软件列表"}{"table"}},['DESeq',T('差异表达基因软件，适用于有生物学重复的情况'),'http://www.bioconductor.org/packages/release/bioc/html/DESeq.html'];
    push @{$table_info{"软件列表"}{"table"}},['EBSeq',T('差异表达基因软件，适用于无生物学重复的情况'),'https://www.biostat.wisc.edu/~kendzior/EBSEQ/'];
    push @{$table_info{"软件列表"}{"table"}},['topGO',T('基于R的GO功能富集的软件'),'https://www.bioconductor.org/packages/2.12/bioc/html/topGO.html'];

	#close (TAB);

    ## 附表2 数据库列表
    # width
    @{$table_info{"数据库列表"}{"width"}}= (1401,4200,3207);
    @{$table_info{"数据库列表"}{"align"}}= ('c','w','w');
    #info
    push @{$table_info{"数据库列表"}{"table"}},['Database','Description','Homepage'];
    push @{$table_info{"数据库列表"}{"table"}},['NR',T('非冗余蛋白序列数据库'),'ftp://ftp.ncbi.nih.gov/blast/db/'];
    push @{$table_info{"数据库列表"}{"table"}},['Swiss-Prot',T('人工注释的非冗余蛋白序列数据库'),'http://www.uniprot.org/'];
    push @{$table_info{"数据库列表"}{"table"}},['GO',T('基因本体（Gene Ontology）数据库'),'http://www.geneontology.org/'];
    push @{$table_info{"数据库列表"}{"table"}},['COG',T('Clusters of Orthologous Groups'),'http://www.ncbi.nlm.nih.gov/COG/'];
    push @{$table_info{"数据库列表"}{"table"}},['KOG',T('euKaryotic Orthologous Groups'),'http://www.ncbi.nlm.nih.gov/COG/'];
    push @{$table_info{"数据库列表"}{"table"}},['KEGG',T('京都基因与基因组百科'),'http://www.genome.jp/kegg/'];
    push @{$table_info{"数据库列表"}{"table"}},['Pfam',T('蛋白家族数据库'),'http://pfam.xfam.org/'];

	#close (TAB);

    ## 附表3 核酸编码表
    # width
    @{$table_info{"核酸编码表"}{"width"}}= (2202,3303,3303);
    @{$table_info{"核酸编码表"}{"align"}}= ('c','c','c');
    #info
    push @{$table_info{"核酸编码表"}{"table"}},['Nucleic Acid Code','Meaning','Mnemonic'];
    push @{$table_info{"核酸编码表"}{"table"}},['A','A','Adenine'];
    push @{$table_info{"核酸编码表"}{"table"}},['C','C','Cytosine'];
    push @{$table_info{"核酸编码表"}{"table"}},['G','G','Guanine'];
    push @{$table_info{"核酸编码表"}{"table"}},['T','T','Thymine'];
    push @{$table_info{"核酸编码表"}{"table"}},['U','U','Uracil'];
    push @{$table_info{"核酸编码表"}{"table"}},['R','A or G','puRine'];
    push @{$table_info{"核酸编码表"}{"table"}},['Y','C, T or U','pYrimidines'];
    push @{$table_info{"核酸编码表"}{"table"}},['K','G, T or U','bases which are Ketones'];
    push @{$table_info{"核酸编码表"}{"table"}},['M','A or C','bases with aMino groups'];
    push @{$table_info{"核酸编码表"}{"table"}},['S','C or G','Strong interaction'];
    push @{$table_info{"核酸编码表"}{"table"}},['W','A, T or U','Weak interaction'];
    push @{$table_info{"核酸编码表"}{"table"}},['B','not A (i.e. C, G, T or U)','B comes after A'];
    push @{$table_info{"核酸编码表"}{"table"}},['D','not C (i.e. A, G, T or U)','D comes after C'];
    push @{$table_info{"核酸编码表"}{"table"}},['H','not G (i.e., A, C, T or U)','H comes after G'];
    push @{$table_info{"核酸编码表"}{"table"}},['V','neither T nor U (i.e. A, C or G)','V comes after U'];
    push @{$table_info{"核酸编码表"}{"table"}},['N','A C G T U','Nucleic acid'];

	#close (TAB);
}

# ------------------------------------------------------------------
# xml_write
# ------------------------------------------------------------------

sub xml_write {	#&xml_write($format, \@context);
	my ($format, $context)=@_;
#	print $context,"\n";
	my $line= scalar @{$context};
	if ($line > 1 && $format !~ "正文" && $format !~ "表格" && $format !~ "项目结果概述" && $format !~ "参考文献" && $format !~/图片/ && $format !~/级标题/ && $format !~/文件/ ) {
		printf("wrong format:$format \t\"$context->[0]\"!\n");
		die;
	}


	for (my $i=0; $i<$line;$i++) {

		## some data info
		my (@data_among_word) = $context->[$i] =~/(\$[a-z,A-Z,_,0-9]+)/g;

		if (@data_among_word !=0) {
			for (my $j=0;$j<@data_among_word;$j++) {
				$data_among_word[$j] =~s/\$/\\\$/;
				#print "$data_among_word[$j]\n" unless (exists $data_among_words{$data_among_word[$j]});
				$context->[$i] =~s/$data_among_word[$j]/$data_among_words{$data_among_word[$j]}/;
			}
		}
	}
		

	if ($format eq "项目结果概述") {

		$report_xml.='&lt;p class=&quot; p-abstract&quot; &gt;合同关键指标:&lt;/p&gt;&lt;p class=&quot; p-abstract&quot; &gt;';
		$report_xml.='完成'.$sample_num.'个样品的转录组测序，每个样品测序产出不少于'.$demand.'Gb Clean Data，Q30碱基百分比达到'.$q30.'。完成Unigene的拼接和功能注释，并进行Unigene的CDS预测、SSR分析、SNP分析、表达定量、差异分析和功能富集分析。';
		$report_xml.='&lt;/p&gt;&lt;p class=&quot; p-abstract&quot; &gt;分析结果概述：&lt;/p&gt;&lt;p class=&quot; p-abstract&quot; &gt;';
		$report_xml.=join("",@{$context}).'&lt;/p&gt;'.'" />'."\n";
		$abstract=join("",@{$context});

	}
	elsif($context->[0] =~ /差异表达基因维恩图/){}
	elsif($format eq "表格"){## 表格
		&table_write($format,$context);
    }elsif($format=~/级标题/){## 附表
		&title_write($format,$context);
    }elsif($format=~/文件集合/) {## 公式
        &file_all_write($format,$context);
	}elsif($format=~/文件/) {## 公式
        &file_write($format,$context);
	}elsif($format =~/图片集合/) {
		&picture_all_write($format,$context);
	}elsif($format =~/图片/){## 图片
		&picture_write($format,$context);
	}elsif($format =~/正文/){## 图片
		&text_write($format,$context);
	}elsif($format =~/参考文献/){## 图片
		&ref_all_write($format,$context);
	}

}


sub picture_write{# 
	my ($format,$context) = @_;
	my $line= scalar @{$context};
	my ($picture, $name)=split /\s+/,$context->[0],2;
	$context->[0]=~s/图\d*/图$pic_num/;
	my $pic_n=$pic_num;
	#$pic_n=~s/图//g;
	# my @pic_info = @{$pic_info{$name}};
	print "$name\n" unless(exists $pic_info{$name});
	my $pic_info = $pic_info{$name};
	if ($line==1) {
		if($name=~/公式|流程图/){
			$report_xml.="\t".'<pic name="'.$context->[0].'" type="img-width-normal" desc="'.' '.'" path="'.$pic_info.'" /> '."\n";
		}
		else{
			$report_xml.="\t".'<pic name="'.$context->[0].'" type="'.'type1'.'" desc="'.' '.'" path="'.$pic_info.'" /> '."\n";
		}

		$image_info{$pic_n}=[$pic_info,'delete',$context->[0]] if ($picture=~/图/); 
	}
	if ($line==2) {
		if($name=~/公式|流程图/){
			$report_xml.="\t".'<pic name="'.$context->[0].'" type="img-width-normal" desc="'.$context->[1].'" path="'.$pic_info.'" /> '."\n";
		}
		else{
			$report_xml.="\t".'<pic name="'.$context->[0].'" type="'.'type1'.'" desc="'.$context->[1].'" path="'.$pic_info.'" /> '."\n";
		}
		$image_info{$pic_n}=[$pic_info,$context->[1],$context->[0]] if ($picture=~/图/);
	}
	$pic_num++;
}

sub picture_all_write{# 
	my ($format,$context) = @_;
	my $line= scalar @{$context};
	my ($picture,$name)=split /\s+/,$context->[0],2;
	return(1) if($name =~/Data数据分布统计图/ and !exists $pic_all_info{"Raw Data数据分布统计图"});
	$context->[0]=~s/图\d*/图$pic_num/;
	my $pic_n=$pic_num;
	#$pic_n=~s/图//g;
	my $number=0;
	my @pic_info = @{$pic_all_info{$name}};
	if ($line==1) {
		$report_xml.="\t".'<pic_list name="'.$context->[0].'" type="'.'type1'.'" desc="">'."\n";
		$image_info{$pic_n}=[$pic_info[0],'delete',$context->[0]] if ($picture=~/图/);

		foreach my $path (@pic_info) {
			my $title=basename$path;
			$title=~s/(^\w+).*/$1/;
			$report_xml.="\t\t".'<pic name="'.$title.'" desc="" path="'.$path.'" />'."\n";
			$number++;
			last if ($number>=100) ;
		}
		$report_xml.="\t".'</pic_list>'."\n";
	}
	else {
		$report_xml.="\t".'<pic_list name="'.$context->[0].'" type="'.'type1'.'" desc="'.$context->[1].'">'."\n";
		$image_info{$pic_n}=[$pic_info[0],$context->[1],$context->[0]] if ($picture=~/图/);

		foreach my $path (@pic_info) {
			my $title=basename$path;
			#$title=~s/(^\w+).*/$1/;
			$report_xml.="\t\t".'<pic name="'.$title.'" desc="" path="'.$path.'" />'."\n";
			$number++;
			last if ($number>=100) ;	
		}
		$report_xml.="\t".'</pic_list>'."\n";
	}
	$pic_num++;
}
sub text_write{# 
	my ($format,$context) = @_;
	my $line= scalar @{$context};
	my $AllText="";
	foreach my $text (@{$context}) {
		$text=~s/\[(\d+)\]/\&lt\;a href=\&quot\;\#ref$1\&quot\;\&gt\;\[$1\]\&lt\;\/a\&gt\;/g;
		#$AllText.=$text."</br>";
		#$AllText.=$text."\&lt;br\/\&gt;\&nbsp;\&nbsp;\&nbsp;\&nbsp;\&nbsp;\&nbsp;\&nbsp;\&nbsp;";
		$report_xml.="\t".'<p type="'.'type1'.'" desc="'.$text.'" />'."\n";
	}
	#$AllText=~s/\&nbsp;\&nbsp;\&nbsp;\&nbsp;\&nbsp;\&nbsp;\&nbsp;\&nbsp;$//;
	#$report_xml.="\t".'<p type="'.'type1'.'" desc="'.$AllText.'" />'."\n";
}

sub file_write{# 
	my ($format,$context) = @_;
	my $line= scalar @{$context};
	my ($name)= split /\s+/,$context->[0],2;
	my $path=$file_info{$name};
	#print "$name\n" unless exists $table_file_info{$name};
	if ($name=~/组装结果序列/) {
		$report_xml.="\t".'<file name="'.$name.'" type="'.'type1'.'" desc="'.''.'" path="'.$path.'" action="'.'fasta'.'" />'."\n";		
	}
	else{
		if (defined $context->[1]) {
			$report_xml.="\t".'<file name="'.$name.'" type="'.'xls'.'" desc="'.$context->[1].'" path="'.$path.'" action="'.'xls'.'" />'."\n";
		}else {
			$report_xml.="\t".'<file name="'.$name.'" type="'.'xls'.'" desc="'.''.'" path="'.$path.'" action="'.'xls'.'" />'."\n";
		}
	}
}

sub table_write{# 
	my ($format,$context) = @_;
	my $line= scalar @{$context};
	my (undef,$name)= split /\s+/,$context->[0],2;
	$context->[0]=~s/表\d*/表$table_num/;
	#print "$name\n" unless exists $table_file_info{$name};
	my $path=$table_file_info{$name};
	if ($line==1) {
		$report_xml.="\t".'<table name="'.$context->[0].'" type="'.'full'.'" desc="'.''.'" path="'.$path.'" /> '."\n";
	}
	if ($line==2) {
		$report_xml.="\t".'<table name="'.$context->[0].'" type="'.'full'.'" desc="'.$context->[1].'" path="'.$path.'" /> '."\n";
	}
	$table_num++;
}
sub file_all_write{# 
	my ($format,$context) = @_;
	my $line= scalar @{$context};
	my ($name)= split /\s+/,$context->[0],2;
	
	my @file_info = @{$file_all_info{$name}};
	if ($line==1) {
		$report_xml.="\t".'<file_list name="'.$context->[0].'" type="'.'xls'.'" desc="">'."\n";
	}
	else{
		$report_xml.="\t".'<file_list name="'.$context->[0].'" type="'.'xls'.'" desc="'.$context->[1].'">'."\n";
	}
	foreach my $path (@file_info) {
		my $title=basename$path;
		#$title=~s/(^\w+).*/$1/;
		$report_xml.="\t\t".'<file name="'.$title.'" desc="" type="xls" path="'.$path.'" action="'.'xls'.'" />'."\n";
	}
	$report_xml.="\t".'</file_list>'."\n";

}

sub ref_all_write{# 
	my ($format,$context) = @_;
	$first_title++;
	$format=~s/^[\d\.\s]*/$first_title /;
	my $line= scalar @{$context};

	$report_xml.="\t".'<ref_list name="'.$format.'" type="'.'type1'.'" desc="'.''.'">'."\n";
	foreach my $path (@{$context}) {

		my ($num,$name)= split /\]\s*/,$path,2;
		#my $ref_info=$reference_info{$name};
		$num=~s/\[//g;
#		print $num,"\t$ref_info\n";
		next if (!defined $ptop and $name=~/Cytoscape|STRING/);
		$report_xml.="\t\t".'<ref id="'.$num.'" name="'.$name.'" link="" />'."\n" ;
		#$report_xml.="\t\t".'<ref id="'.$num.'" name="'.$name.'" link="'.$ref_info.'" />'."\n" ;
		
	}
	$report_xml.="\t".'</ref_list>'."\n";

}
sub table_html(){
	my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
	print "title:$title\n";##for
	open HH,">$outHtml" or die "$!";
	print HH <<HTML;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
<meta charset="UTF-8"></meta>
<!--[if lt IE 9]>
<script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
<![endif]-->
<meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
<meta content="width=device-width, initial-scale=1" name="viewport"></meta>
<link href="$srcPath/css/bootstrap.min.css" type="text/css" rel="stylesheet" />
<link href="$srcPath/css/index.css" type="text/css" rel="stylesheet" />
<link href="$srcPath/js/fancyBox/jquery.fancybox.css" type="text/css" rel="stylesheet" />
<link href="$srcPath/css/nav.css" type="text/css" rel="stylesheet" />
<link href="$srcPath/css/raxus.css" type="text/css" rel="stylesheet" />
<script src="$srcPath/js/jquery-1.11.3.min.js" type="text/javascript"></script>
<script src="$srcPath/js/nav.js" type="text/javascript"></script>
<script src="$srcPath/js/raxus-slider.min.js" type="text/javascript"></script>
<script src="$srcPath/js/fancyBox/jquery.fancybox.pack.js" type="text/javascript"></script>
<script src="$srcPath/js/fancyBox/jquery.mousewheel-3.0.6.pack.js" type="text/javascript"></script>
<script src="$srcPath/js/bootstrap.min.js" type="text/javascript"></script>
<script src="$srcPath/js/ready.js" type="text/javascript"></script>
<script src="$srcPath/js/scrolltop.js" type="text/javascript"></script>
</head>
<body style=\"overflow:hidden\"> 
<div class="container shadow"  id="box"><header><img src="$srcPath/images/logo.jpg" class="pull-right" />
<div role="main" ><header><h2 id="title" class="text-center">$title</h2>
</header>
</div>
HTML
	if($text){
		print HH "<div class=\"table-responsive\" id=\"textbox\"><p>$text</p></div>\n";
	}
	print HH "<div class=\"table-responsive\" id=\"box2\"><table class=\"table table-bordered table-hover table-striped\"><thead><tr class=\"bg-info\">\n";
	if($text){
		print HH <<HTML;
<script type="text/javascript">
    var textbox=\$("#textbox").height();
    var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8-textbox;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
        var height2 = height*0.8-textbox;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML
	}
	else{
	print HH <<HTML;
<script type="text/javascript">
	var windowHeight = \$(window).height();
    var height = windowHeight*0.98;
    var height2 = height*0.8;
    \$("#box").css('height',height);
    \$("#box2").css('height',height2);
    \$(window).resize(function() {
        var windowHeight = \$(window).height();
        var height = windowHeight*0.98;
        var height2 = height*0.8;
        \$("#box").css('height',height);
        \$("#box2").css('height',height2);
    });
</script>
HTML

	}
	
	for (my $i=0;$i<=$#{$$input[0]};$i++){
		print HH "<th>$$input[0][$i]</th>\n";	
	}
	print HH "</tr></thead>\n<tbody>\n";
	for (my $k=1;$k<=$#{$input};$k++) {
		print HH "<tr>";
		for (my $i=0;$i<=$#{$$input[$k]};$i++){
			if($linkColNum){
				my $j=$i+1;
				if($linkColNum=~/,$j,/){
					#print "out:$outHtml\n$k  $i $$input[$k][$i]\n"  if(!exists $$linkHash{$$input[$k][$i]});exit;
					print HH "<td><a href=\"$$linkHash{$$input[$k][$i]}\" target=\"_blank\">$$input[$k][$i]</a></td>";
				}
				else{
					print HH "<td>$$input[$k][$i]</td>";
				}
			}
			else{
				print HH "<td>$$input[$k][$i]</td>";
			}
		}
		print HH "</tr>\n";
	}	
print HH <<XGL;
</tbody>
</table>
</div>
</body>
</html>
XGL
close HH;
}

sub text_html(){
	my ($input,$outHtml,$title,$linkColNum,$linkHash,$srcPath,$text)=@_;
	open H,">$outHtml";
	print H <<HTML;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html lang="zh_CN" xmlns="http://www.w3.org/1999/xhtml"><head><title>$title</title>
<meta charset="UTF-8"></meta>
<!--[if lt IE 9]>
<script src="$srcPath/js/html5shiv.min.js"></script><script src="$srcPath/js/respond.min.js"></script>
<![endif]-->
<meta content="IE=edge" http-equiv="X-UA-Compatible"></meta>
<meta content="width=device-width, initial-scale=1" name="viewport"></meta>
<link href="$srcPath/css/bootstrap.min.css" type="text/css" rel="stylesheet" />
<link href="$srcPath/css/index.css" type="text/css" rel="stylesheet" />
<link href="$srcPath/js/fancyBox/jquery.fancybox.css" type="text/css" rel="stylesheet" />
<link href="$srcPath/css/nav.css" type="text/css" rel="stylesheet" />
<link href="$srcPath/css/raxus.css" type="text/css" rel="stylesheet" />
<script src="$srcPath/js/jquery-1.11.3.min.js" type="text/javascript"></script>
<script src="$srcPath/js/nav.js" type="text/javascript"></script>
<script src="$srcPath/js/raxus-slider.min.js" type="text/javascript"></script>
<script src="$srcPath/js/fancyBox/jquery.fancybox.pack.js" type="text/javascript"></script>
<script src="$srcPath/js/fancyBox/jquery.mousewheel-3.0.6.pack.js" type="text/javascript"></script>
<script src="$srcPath/js/bootstrap.min.js" type="text/javascript"></script>
<script src="$srcPath/js/ready.js" type="text/javascript"></script>
<script src="$srcPath/js/scrolltop.js" type="text/javascript"></script>
</head>
<body>
<div class="container shadow"><header><img src="$srcPath/images/logo.jpg" class="pull-right" />
<div role="main" ><header><h2 id="title" class="text-center">$title</h2>
</header>
</div>
<div style="word-wrap:break-word;word-break:break-all">
HTML
	if($text){
		print H "<p>$text</p>\n";
	}
	for my $key(keys %{$input}){
		print H "$key<br/>\n$$input{$key}<br/>\n";	
	}
print H <<XGL;
</div>
</body>
</html>
XGL
}


sub title_write{# 
	my ($format,$context) = @_;
	if ($format=~/一级标题/) {
		$first_title++;
		($second_title,$third_title)=(0,0);
		$context->[0]=~s/^[\d\.\s]*/$first_title /;
		$report_xml.="\t".'<h1 name="'.$context->[0].'" type="'.'type1'.'" desc="'.$context->[0].'" />'."\n";
		$report_xml.="\t".'<pic_list name="" type="type1" desc="">images_preview</pic_list>'."\n" if ($context->[0]=~/图片预览/) ;
	}elsif ($format=~/二级标题/) {
		$second_title++;
		$third_title=0;
		$context->[0]=~s/^[\d\.\s]*/$first_title.$second_title /;
		$report_xml.="\t".'<h2 name="'.$context->[0].'" type="'.'type1'.'" desc="'.$context->[0].'" />'."\n";
	}elsif ($format=~/三级标题/) {
		$third_title++;
		$context->[0]=~s/^[\d\.\s]*/$first_title.$second_title.$third_title /;
		$report_xml.="\t".'<h3 name="'.$context->[0].'" type="'.'type1'.'" desc="'.$context->[0].'" />'."\n";
	}
}
#整数格式：3位一个逗号
sub Integer_Three_Digit{#
	my $interger = shift;
	$interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $interger;
}

#整数格式：3位一个逗号
#小数格式：小数点后两位
sub format_figure{#
	my $figure = shift;
	#print "figure$figure_num:$figure\n";
	$figure_num++;
	if (!defined $figure) {
		print "Err:$figure\n";
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

sub AbsolutePath {		#获取指定目录或文件的决定路径
    my ($type,$input) = @_;
    my $return;
    $/="\n";

    if ($type eq 'dir') {
        my $pwd = `pwd`;
        chomp $pwd;
        chdir($input);
        $return = `pwd`;
        chomp $return;
        chdir($pwd);
    }
    elsif($type eq 'file') {
        my $pwde = `pwd`;
        chomp $pwde;

        my $dir=dirname($input);
        my $file=basename($input);
        chdir($dir);
        $return = `pwd`;
        chomp $return;
        $return .="\/".$file;
        chdir($pwde);
    }
    return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub GetDate {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d\/%02d\/%02d", $year+1900, $mon+1, $day);
}

=c
sub T {
   my $text = shift;
   return decode( 'UTF-8', $text);
}
=cut

sub detail_cfg_read {
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        my ($key, $value) = (split /\s+/,$_,2)[0,1];
            $detail_cfg->{$key} = $value;
    }
    close CFG;
}
### get_name($path_name_star)
# $path_name_star: contain *; only one file in practical, when more or no file, then die;
sub get_only_one_file{# get_name("$dir_data/Unigene/*.Unigene.distribution.png")
	my $path_name_star = shift;
	my @name = glob("$path_name_star");

	if (@name >1) {
		print "Wrong: more file of $path_name_star\n";
		die;
	}elsif (@name == 0) {
		print "Wrong: no file of $path_name_star\n";
		die;
	}
	return $name[0];
}

sub T {
   my $text = shift;
   return decode( 'UTF-8', $text);
}

sub USAGE {
        my $usage=<<"USAGE";
#-----------------------------------------------------------------------------------------
  Program: $Script
  Version: $version
  Contact: zeng huaping<zenghp\@biomarker.com.cn>
     Date: 
 Modifier: WangYajing <wangyj\@biomarker.com.cn>
     Date: 2015-04-22
    Usage:
      Options:
      --id  <dir>   directory of Web report.
      --cfg <dir>   html config
      --pp  <str>   protein to protein analysis have done
      --Only_One_Sam Only one sample

#-----------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}