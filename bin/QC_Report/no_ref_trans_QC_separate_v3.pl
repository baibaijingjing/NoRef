#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.1.6";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($cfg,$od,$data_type,$sep,$cut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$od,
				"a:s"=>\$sep,
				"b:s"=>\$cut,
				"type:s"=>\$data_type,
				"cfg:s"=>\$cfg,
				) or &USAGE;
&USAGE unless ($cfg and $od);
open (CFG,$cfg) or die $!;
$cut||=15;
$sep||=2;
$data_type||='pe';
if(!-d $od){mkdir $od};  #如果不存在目录，则建一个目录

my %cfg;my %stat;my @stat;my %lib;
while (<CFG>) {
	chomp;
	s/\r//;
	next if (/^$/);
	next if (/^\#/);
	if (/^Project ID/) {
		$cfg{'project_id'}=<CFG>;
		$cfg{'project_id'}=~s/\r//;
		chomp $cfg{'project_id'};
	}
	if (/^Library path/) {
		$cfg{'Library_path'}=<CFG>;
		$cfg{'Library_path'}=~s/\r//;
		chomp $cfg{'Library_path'};
	}
	if (/^Analysis path/) {
		$cfg{'Analysis_path'}=<CFG>;
		$cfg{'Analysis_path'}=~s/\r//;
		chomp $cfg{'Analysis_path'};
	}
	if (/^DataAssess path/) {
		$cfg{'DataAssess_path'}=<CFG>;
		$cfg{'DataAssess_path'}=~s/\r//;
		chomp $cfg{'DataAssess_path'};
	}
	if (/^DEG_Analysis path/) {
		$cfg{'DEG_Analysis_path'}=<CFG>;
		$cfg{'DEG_Analysis_path'}=~s/\r//;
		chomp $cfg{'DEG_Analysis_path'};
	}
}
close CFG;
################################文库
my $num=0;
if ($cfg{'Library_path'} ne '') {
	my @lib=glob "$cfg{'Library_path'}/*.txt";
	foreach my $lib (@lib) {
		open (LIB,$lib) or die $!;
		while (<LIB>) {
			chomp;
			s/\r//;
			if (/$cfg{'project_id'}/ and /Trans/) {
				my @in=split /\t/,$_;
				$in[3]=~/^[^-]+-([^_]+)-.*/;
				my $sample=$1;
				$lib{$sample}{'lib_id'}=$in[3];
				$lib{$sample}{'qubit'} = &format_figure($in[14]);
				$lib{$sample}{'quality'}=$in[19];
			}
		}
	close LIB;
	}
      if($num>0){
	open OUT,">$od/Library.stat" or die $!;

	print OUT "#SampleID\tLibraryID\tQbit\tLibraryQuality\n";
	foreach my $sample (sort keys %lib) {
		print OUT "$sample\t$lib{$sample}{'lib_id'}\t$lib{$sample}{'qubit'}\t$lib{$sample}{'quality'}\n";
	}
	close OUT;
    }
}
###############################数据评估
#Data_Assess/AllSample_GC_Q.stat
#SampleID        ReadSum BaseSum GC(%)   N(%)    Q20(%)  CycleQ20(%)     Q30(%)
#T1      50244012        10148671105     42.68   0.05    92.22   100.00  88.08
#Data_Assess/AllSample.data.stat
#Run_ID Sample_ID       Required_Data   Obtained_Data   Obtained_Reads  GC_content      Q20%    Q30%    Adapter%        rRNA%   Inferior%
#Run136  T1      3000000000      2811459597      13919233        52.70   99.03   94.88   0.53    1.45    18.18

my $GC_div;
if ($cfg{'DataAssess_path'} ne '') {
	open DATA,"$cfg{'DataAssess_path'}/AllSample_GC_Q.stat" or die $!;
	my @GC;
	while (<DATA>) {
		chomp;
        next if ($.==1 || /^#/);
		my @in=split /\t/,$_;
		my $sample=$in[0];
		$stat{$sample}{'GC'} = &format_figure($in[3]);
		push @GC,$in[3];
		$stat{$sample}{'Data'} = &format_figure($in[2]);
		$stat{$sample}{'Q30'}= &format_figure($in[-1]).'%';
	}
	$GC_div = &max(@GC)-&min(@GC);
	$GC_div = &percent($GC_div/100);
	close DATA;

	my @GC_Assess = glob "$cfg{'DataAssess_path'}/*.acgtn";
	foreach my $GC_Assess (@GC_Assess) {
#		$GC_Assess=~/.*\/([a-zA-Z][a-zA-Z0-9]*).acgtn/;
		$GC_Assess=~/.*\/([^\/]+)\.acgtn/;
		my $sample=$1;
		($stat{$sample}{'GC_iso'},$stat{$sample}{'GC_wave'})=&GC_Assess($GC_Assess,$data_type,$sep,$cut);
		$stat{$sample}{'GC_iso'}=&percent($stat{$sample}{'GC_iso'}/100);
		$stat{$sample}{'GC_wave'}=&percent($stat{$sample}{'GC_wave'});
	}
	open OUT,">$od/DataAssess.stat" or die $!;
	push @stat,"$od/DataAssess.stat";
	print OUT "#SampleID\tObtainedData\tObtainedQ30%\tGCcontent\tGCcontentRange\tBaseIso\tBaseWave\n";
	foreach my $sample (sort keys %stat) {
		print OUT "$sample\t$stat{$sample}{'Data'}\t$stat{$sample}{'Q30'}\t";
		print OUT "$stat{$sample}{'GC'}%\t$GC_div\t$stat{$sample}{'GC_iso'}\t$stat{$sample}{'GC_wave'}\n";
	}
	close OUT;
   if(-e "$cfg{'DataAssess_path'}/AllSample.data.stat"){
	open ASS,"$cfg{'DataAssess_path'}/AllSample.data.stat" or die $!;
	my ($rRNA_col,$ReqData_col);
	while (<ASS>) {
		chomp;
		if($.==1){
			my @tmp=split/\t/;
			for(my $i=0;$i<=$#tmp;$i++){
				$rRNA_col=$i if $tmp[$i]=~/^rRNA/;
				$ReqData_col=$i if $tmp[$i]=~/^Required_data/;
			}
			next;
		}
        next if (/^\s+$/|| /^#/);
        my @col = (split /\t/);
        #die "Note your file: $cfg{'DataAssess_path'}/AllSample.data.stat \n" if (@col != 11);##by xugl 2016-1-25
		my $sample = $col[1];
		$stat{$sample}{'rRNA'} = &format_figure($col[$rRNA_col]).'%';
		$stat{$sample}{'Required_data'}= &format_figure($col[$ReqData_col]);
	}
	close ASS;
   }
}

#########################################主流程分析
my $gene_div = 0;
my $rmap_div = 0;

if ($cfg{'Analysis_path'} ne '') {
	#my @rmap = glob "$cfg{'Analysis_path'}/Remap/rmap/*/*.statMapped.xls";
	my @rmap = glob "$cfg{'Analysis_path'}/Remap/rmap/*/*.Mapped.stat.xls";
	my @rmap_div;
	foreach my $rmap (@rmap) {
		open RMAP,$rmap or die $!;
		$rmap=~/Remap\/rmap\/([^\/]+)\//;
		my $sample=$1;
		while (<RMAP>) {
			chomp;
			if (/^Mapped Reads/) {
				$stat{$sample}{'rmap'}=(split /\t/,$_)[2];
				$stat{$sample}{'rmap'}=~s/\%//;
				push @rmap_div,$stat{$sample}{'rmap'};
			}
		}
		close RMAP;
	}
	$rmap_div = &max(@rmap_div)-&min(@rmap_div);
	$rmap_div = &percent($rmap_div/100);
	my @assembly = glob "$cfg{'Analysis_path'}/Assembly/Trinity_assembly/*/*_Result/Unigenes/*.Unigenes.stat.xls";
	my @gene_div;
	foreach my $assembly (@assembly) {
		open ASSEMBLY,$assembly or die $!;
#		$assembly=~/(T\d+).Unigenes.stat.xls/;
		#$assembly=~/([a-zA-Z][a-zA-Z0-9]*).Unigenes.stat.xls/;
		$assembly=~/([^\/]+).Unigenes.stat.xls/;
		my $sample=$1;
        die "ERROR: sample id: $sample doesn't exist, please check!\n" unless (exists $stat{$sample});
		while (<ASSEMBLY>) {
			chomp;
			if (/^N50/) {
				$stat{$sample}{'N50'} = &format_figure((split /\t/,$_)[-1]);
			}
			if (/^1000/) {
				$stat{$sample}{'1K+'}=(split /\t/,$_)[1];
				$stat{$sample}{'1K+'}+=(split /\t/,<ASSEMBLY>)[1];
				push @gene_div,$stat{$sample}{'1K+'};
                $stat{$sample}{'1K+'} = &format_figure($stat{$sample}{'1K+'});
			}
			if (/^Total Length/) {#####################################
				$stat{$sample}{'length'} = &format_figure((split /\t/,$_)[-1]);
			}
			if (/^200\-/) {
				$stat{$sample}{'300'}=(split /\t/,$_)[-1];
				$stat{$sample}{'300'} = &percent($stat{$sample}{'300'});
			}
			if (/^Total Number/) {
				$stat{$sample}{'number'} = &format_figure((split /\t/,$_)[-1]);
			}
		}
		close ASSEMBLY;
	}


	$gene_div=(&max(@gene_div)-&min(@gene_div))/&min(@gene_div);
	$gene_div=&percent($gene_div);
	my @insert = glob "$cfg{'Analysis_path'}/Remap/rmap/*/*.insertSize";
	foreach my $insert (@insert) {
		open INSERT,$insert or die $!;
		$insert=~/Remap\/rmap\/(.*)\//;
		my $sample =$1; 
		my %insert;
		while (<INSERT>) {
			chomp;
			my @in = split /\t/,$_;
			$insert{$in[0]}=$in[1];
		}
		close INSERT;
		my @len = sort {$insert{$a}<=>$insert{$b}} (keys %insert);
		$stat{$sample}{'insert_max'}=$len[-1];
		my $count=0;my $sum=0;
		for (my $i=$len[-1]-75;$i<=$len[-1]+75;$i++) {
			$insert{$i}||=0;
			$count+=$insert{$i};
		}
		foreach my $len (keys %insert) {
			$sum+=$insert{$len};
		}
		$stat{$sample}{'insert_ratio'} = &percent($count/$sum);
	}
	my @randcheck = glob "$cfg{'Analysis_path'}/Remap/rmap/*/*.randcheck.list";

	foreach my $randcheck (@randcheck) {
		open RAND,$randcheck or die $!;
		$randcheck=~/Remap\/rmap\/(.*)\//;
		my $sample=$1;
		$stat{$sample}{'rand'} = 0;

		while (<RAND>) {
			chomp;
			if (/^[\.\d]+\:([\.\d]+)/) {
#				$stat{$sample}{'rand'}+=($1-1)**2;  #DepthVariance
                $stat{$sample}{'rand'} = ($stat{$sample}{'rand'}>$1) ? $stat{$sample}{'rand'} : $1; #MaxDepthPer
			}
		}
		close RAND;
        $stat{$sample}{'rand'} = &format_figure($stat{$sample}{'rand'}).'%';
	}

	open OUT,">$od/Analysis.stat" or die $!;
	push @stat,"$od/Analysis.stat";
	print OUT "#SampleID\tMapped%\tMapped%Range\tMaxDepthPer\tInsertMode\tInsertAround%\tN50Length\t<300nt%\tUnigenLength\n";
	foreach my $sample (sort keys %stat) {
		$stat{$sample}{'rmap'}=&percent($stat{$sample}{'rmap'}/100);
		print OUT "$sample\t$stat{$sample}{'rmap'}\t$rmap_div\t$stat{$sample}{'rand'}\t$stat{$sample}{'insert_max'}\t$stat{$sample}{'insert_ratio'}\t";
		print OUT "$stat{$sample}{'N50'}\t$stat{$sample}{'300'}\t$stat{$sample}{'length'}\n";
	}
	close OUT;
}

######################################差异基因分析
my $DEG_max=0;
my $DEG_min=100000000;

if (exists $cfg{'DEG_Analysis_path'}) {
	my @deg = glob "$cfg{'DEG_Analysis_path'}/*_vs_*/*.DEG_final.xls";
	my @deg_num;
    my %DEGset;
	foreach my $deg (@deg) {
		my $deg_num=`grep -v "#" $deg |wc -l`; chomp($deg_num);
        my ($set_name) = &basename($deg) =~/^(\S+).DEG_final.xls$/;
        $DEGset{$set_name}{'up'} = `grep -w 'up' $deg |wc -l`;
        $DEGset{$set_name}{'down'} = `grep -w 'down' $deg |wc -l`;
        chomp($DEGset{$set_name}{'up'});
        chomp($DEGset{$set_name}{'down'});
        $DEGset{$set_name}{'all'} = $deg_num;

		push @deg_num,$deg_num;
        $DEGset{$set_name}{'all'} = &format_figure($DEGset{$set_name}{'all'});
        $DEGset{$set_name}{'up'} = &format_figure($DEGset{$set_name}{'up'});
        $DEGset{$set_name}{'down'} = &format_figure($DEGset{$set_name}{'down'});
	}
	$DEG_max=&max(@deg_num);
	$DEG_min=&min(@deg_num);
	open OUT,">$od/DEG.stat" or die $!;
	push @stat,"$od/DEG.stat";
	print OUT "#DEG_Set\tDEG_All\tUp-regulated\tDown-regulated\n";
	foreach my $set_name (sort keys %DEGset) {
		print OUT "$set_name\t$DEGset{$set_name}{'all'}\t$DEGset{$set_name}{'up'}\t$DEGset{$set_name}{'down'}\n";
	}
	close OUT;
}

########################################结果整理输出
#类型	项目期号	样品名称（BMK）	物种类别	有无参考基因组信息	合同要求数据量（Trans:bp /DEG:Reads）#	实际得到数据量（Trans:bp /DEG:Reads）#	合同要求质量(Q30)	分析用数据量（Trans:bp /DEG:Reads）#	文库Qubit值#	文库胶图质量	碱基分离值%	碱基波动值%	NT库检测污染（菌类）%	NT库检测污染（非近源）	核糖体比例	是否NEB建库	插入片段（目标范围内集中片段比例）%	插入片段主峰位置	GC含量%	GC含量差异%	随机性分布评估	比对效率%	样品间比对效率差异%	各样品组装结果（无 ref，N50）#/-	各样品组装结果（无ref，<300bp的序列比例）	组织部位（无ref）	组装差距(无ref）%/-	差异表达基因数量#/-	可变剪接数占ref基因数的（有ref）%/-	检测到的表达基因%（有ref和DEG）%/-	预测到的新基因数（有ref）#/-	总基因数（无Ref）	生物学重复（组内相关系数-最小值）	生物学重复（组间相关系数-最大值）
#Trans	BMK140122-A62	T1	动植物	 noRef 	 4,000,000,000 	4153200194	80.00%	 3,574,515,871 	4.77 	合格	0.52%	0.00%	0.00%	0.00%	0.00%	否	98.96%	221	37.17%	0.00%	5.3461	81.58%	0.00%	1846	25.48%	 same 	 - 	 - 	 - 	 - 	 - 	20,600	-	-	2014/7/23	2014/8/13

open OUT,">$od/$cfg{'project_id'}.QC.stat" or die $!;

print OUT join("\t",
    '#ProjectType','ProjectID','SampleID','SpeciesType','RefSeq?','RequiredData','ObtainedData','RequiredQuality','UseData','Qubit','LibraryQuality',
    'BaseIso','BaseWave','NT%(contaminate)','NT%(irrelative)','rRNA%','NEBorNOT','InsertAround%','InsertMode','GCcontent','GCcontentRange','MaxDepthPer',
    'Mapped%','Mapped%Range','N50Length','<300nt%','Tissue','assemble_div',
    'DEGnum','AlterSpliceGene%','Saturation','NovelGeneNum','TotalGene','IntraSq.Rmin','InterSq.Rmax','Deadline','FinishDate'
)."\n";

foreach my $key (sort keys %stat) {
    my $DEG_num = ($key eq (sort keys %stat)[0]) ? $DEG_max : $DEG_min;
    my $FinishDate = &GetDate;
    $DEG_num = &format_figure($DEG_num);
    $lib{$key}{'qubit'}="-" if(!defined $lib{$key}{'qubit'});
    $lib{$key}{'quality'}="-" if(!defined $lib{$key}{'quality'});
    $stat{$key}{'Required_data'}="-" if(!defined $stat{$key}{'Required_data'});
    $stat{$key}{'GC_iso'}="-" if(!defined $lib{$key}{'GC_iso'});
    $stat{$key}{'GC_wave'}="-" if(!defined $lib{$key}{'GC_wave'});
    $stat{$key}{'rRNA'}="-" if(!defined $lib{$key}{'rRNA'});
    $stat{$key}{'insert_ratio'}="-" if(!defined $lib{$key}{'insert_ratio'});
    $stat{$key}{'insert_max'}="-" if(!defined $lib{$key}{'insert_max'});
    $stat{$key}{'rRNA'}="-" if(!defined $lib{$key}{'rRNA'});
    $stat{$key}{'GC'}="-" if(!defined $lib{$key}{'GC'});
    $GC_div="-" if(!defined $GC_div);
    $stat{$key}{'rand'}="-" if(!defined $lib{$key}{'rand'});
    $stat{$key}{'rmap'}="-" if(!defined $lib{$key}{'rmap'});
    $rmap_div ||="-";
    $stat{$key}{'number'}||="-";
    $stat{$key}{'300'}||="-";
    $gene_div||="-";

    print OUT join("\t",
        'Trans',$cfg{'project_id'},
        $key,'?','noRef',
        $stat{$key}{'Required_data'},
        $stat{$key}{'Data'},'85.00%',
        $stat{$key}{'Data'},
        $lib{$key}{'qubit'},
        $lib{$key}{'quality'},

        $stat{$key}{'GC_iso'},
        $stat{$key}{'GC_wave'},'?','?',
        $stat{$key}{'rRNA'},'?',
        $stat{$key}{'insert_ratio'},
        $stat{$key}{'insert_max'},
        $stat{$key}{'GC'}."%",
        $GC_div,$stat{$key}{'rand'},

        $stat{$key}{'rmap'},
        $rmap_div,$stat{$key}{'N50'},
        $stat{$key}{'300'},'?',$gene_div,

        $DEG_num,'-','-','-',
        $stat{$key}{'number'},'?','?','?',
        $FinishDate
    )."\n";
}

close OUT;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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
sub GC_Assess {
	my $fIn=shift;
	my $data_type=shift;
	my $sep=shift;
	my $cut=shift;
	my (%GC,$GC_iso,$GC_sep);
	open IN,$fIn or die $!;
	<IN>;
	while (<IN>) {
		my @in=split/\t/,$_;
		$GC{$in[0]}{A}=$in[1];
		$GC{$in[0]}{C}=$in[2];
		$GC{$in[0]}{G}=$in[3];
		$GC{$in[0]}{T}=$in[4];
	}
	close IN;
	my $sum_AT=0;my $sum_GC;my $num=0;
	if ($data_type eq 'se') {
		foreach my $key (keys %GC) {
			if ($key>$cut) {
				$sum_AT+=abs ($GC{$key}{A}-$GC{$key}{T});
				$sum_GC+=abs ($GC{$key}{G}-$GC{$key}{C});
				if (abs ($GC{$key}{A}-$GC{$key}{T})>$sep or abs ($GC{$key}{G}-$GC{$key}{C})>$sep) {
					$num++;
				}
			}
		}
		$GC_iso=&max($sum_AT,$sum_GC)/(101-$cut);
		$GC_sep=$num/(101-$cut);
		return ($GC_iso,$GC_sep);
		}else{
		foreach my $key (keys %GC) {
			if ($key>$cut and $key <=101) {
				$sum_AT+=abs ($GC{$key}{A}-$GC{$key}{T});
				$sum_GC+=abs ($GC{$key}{G}-$GC{$key}{C});
				if ((abs ($GC{$key}{A}-$GC{$key}{T}))>=$sep or (abs ($GC{$key}{G}-$GC{$key}{C}))>=$sep) {
					$num++;
				}
			}elsif ($key>101+$cut) {
				$sum_AT+=abs ($GC{$key}{A}-$GC{$key}{T});
				$sum_GC+=abs ($GC{$key}{G}-$GC{$key}{C});
				if (abs ($GC{$key}{A}-$GC{$key}{T})>$sep or abs ($GC{$key}{G}-$GC{$key}{C})>$sep) {
					$num++;
				}
			}
		}
		$GC_iso=&max($sum_AT,$sum_GC)/(202-$cut*2);
		$GC_sep=$num/(202-$cut*2);
		return ($GC_iso,$GC_sep);
	}
}
################################################################################################################
sub percent{
	my $num =shift;
	if ($num=~s/\%//) {
		$num=sprintf "%.2f",$num;
		$num.='%';
	}else{
		$num=sprintf "%.2f",$num*100;
		$num.='%';
	}
}
################################################################################################################

sub max{#&max(lists or arry);
	#求列表中的最大值
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
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub GetDate {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d", $year+1900, $mon+1, $day);
}

################################################################################################################
#整数格式：3位一个逗号
sub Integer_Three_Digit{#
	my $interger = shift;
	$interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $interger;
}

################################################################################################################
#整数格式：3位一个逗号
#小数格式：小数点后两位
sub format_figure{#
	my $figure = shift;
#	if (!defined $figure) {
#		die;
#	}
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

################################################################################################################
sub USAGE {#
	my $usage=<<"USAGE";
#-----------------------------------------------------------------------------------------
  Program: $Script
  Version: $version
  Contact: Yuan ZhengWen <yuanzw\@biomarker.com.cn> 
     Date: 
 Modifier: Simon Young  <yangxh\@biomarker.com.cn>
     Date: 2014-09-12
 Modified: Based-v1.0.0, for No_Ref_Trans_Free.pl v1.6

    Usage:
           -cfg  <STR>   config file
           -o    <STR>   output dir

           Options:
           -type <STR>   data type, pair-ends or single-ends, pe|se         ['pe']
           -a    <INT>   the separation of AT or GC more than the number    [2]
           -b    <INT>   the number of removing bases per read              [15]
           -h            help document

  Example:
        perl $Script -cfg Config/QC.cfg -o QC_Report/

#-----------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}

