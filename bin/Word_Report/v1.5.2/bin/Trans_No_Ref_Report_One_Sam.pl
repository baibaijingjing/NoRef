#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="2.0.0";
use RTF::Writer;
use Encode qw(decode);

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($dir_template, $dir_data,$fkey,$outdir);

GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$dir_data,
				"pf:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($dir_data and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

$dir_template="$Bin/../template";
my $ftemplate = "$dir_template/Trans_template_single.txt";
my %data_among_words;  ##���ֲ�����Ŀ������Ϣ
my %table_info;
my %pic_info;  #"ת¼�����ʵ������ͼ"=>"$dir_template/experiment_flow_chart.jpg",
my %formula_info;
my $ftable;
my $sample_name;

print "\n[".&GetTime($BEGIN_TIME)."] $Script start ...\n";
# ------------------------------------------------------------------
# read in
# ------------------------------------------------------------------

## picture correspond
&get_picture();

## get formula picture
&get_formula();

# read in table files
&get_table();

## ������ֲ�����Ŀ������Ϣ
&get_data_among_words();

# ------------------------------------------------------------------
# output rtf
# ------------------------------------------------------------------

## open rtf
 my $rtf = RTF::Writer->new_to_file("$outdir/${fkey}_Report.rtf");
 $rtf->prolog(
    'fonts' => ["Times New Roman",T("����"),T("����"),T("������κ")],
	'title' => "abc",
	colors => [
				undef,        # color 0 == black
				[255,0,0],    # color 1 == red
				[0,128,0],    # color 2 == green
				[0,0,255],    # color 3 == blue
				[255,255,0],  # color 4 == yellow
				[255,255,255],# color 5 == white
				[200,200,200],# color 6 == light gray
				[51,102,153], # color 7 == table blue
				[144,172,202],# color 8 == table header blue
				[102,102,102],# color 9 == gray
			   ],
 );
#$rtf->number_pages;

my $limit=0;

## write cover page
&cover_page_write($rtf);

## blank for content
&content_page($rtf);

## write context
# read in template file
open (TEM, $ftemplate) or die $!;
$/="��";
my $n=0;
while (<TEM>) {
	chomp;
	next if(/^$/ || /^\#/);
	my (undef, $format, @context)=split /\n/,$_;

	&rtf_write($format, \@context);
}

$rtf->close;

system "sed -i 's/\\\\_/-/g' $outdir/${fkey}_Report.rtf";    #�����к��߻��ɶ��к���
system "zip $outdir/${fkey}_Report.rtf.zip $outdir/${fkey}_Report.rtf";
#######################################################################################
print "\n[".&GetTime(time())."] $Script done. Total elapsed time: ".(time()-$BEGIN_TIME)."s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

#��ȡ��������Ŀ������Ϣ
sub get_data_among_words {
	#$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[5][0]="2000+";
    ########
    my $finish_date = &GetDate;

	######## ��Ʒ������������ͳ�Ʊ�
	# sample_num	base_num	min_Q30
	my $base_num=0;
    my $sample_num = 0;
	my $min_Q30=1000;
	my $min_data=100;   #Gb
	for (my $i=1; $i<@{$table_info{"��Ʒ������������ͳ�Ʊ�"}{"table"}}; $i++) {
		# base_num
		my $dnum = $table_info{"��Ʒ������������ͳ�Ʊ�"}{"table"}->[$i][3];
		$dnum =~s/,//g;
		$base_num+= $dnum/1000000000;   # convert to Gb
        $sample_num++;

		## min_Q30
		my $Q30 = $table_info{"��Ʒ������������ͳ�Ʊ�"}{"table"}->[$i][5];
		my $data = $table_info{"��Ʒ������������ͳ�Ʊ�"}{"table"}->[$i][3];
        $data =~s/,//g;
        $data = $data/1000000000;   # convert to Gb
		$Q30=~s/%$//;

        $min_Q30 = $Q30 if ($Q30 < $min_Q30);
        $min_data = $data if ($data < $min_data);

	}

	###SSR�������ͳ��
#	my $SSR_num = $table_info{"SSR�������ͳ�Ʊ�"}{"table"}->[-1][1];
	my $SSR_num = $table_info{"SSR�������ͳ�Ʊ�"}{"table"}->[3][1];
	$SSR_num =~s/,//g;


	##### ��װ���ͳ�Ʊ�
	#All_Unigenes
	my $All_Unigenes = $table_info{"��װ���ͳ�Ʊ�"}{"table"}->[6][-1];
	$All_Unigenes=~s/,//g;
	#COMBIN
	my $com_trans="XXX";
	my $com_trans_N50="XXX";
	my $com_Uni_N50="XXX";
	if (-d "$dir_data/Assembly/All_Combination") {
		$com_trans=$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[6][-2];
		$com_trans_N50=$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[8][-2];
		$com_Uni_N50=$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[8][-1];
	}


	#1k_Unigenes
	my $onek_twok_Unigene_str = $table_info{"��װ���ͳ�Ʊ�"}{"table"}->[4][-1];
	my $twok_more_Unigene_str = $table_info{"��װ���ͳ�Ʊ�"}{"table"}->[5][-1];
	my ($onek_twok_Unigene) = $onek_twok_Unigene_str=~/(\S+)\(\S+%\)/;
	my ($twok_more_Unigene) = $twok_more_Unigene_str=~/(\S+)\(\S+%\)/;
	$onek_twok_Unigene =~s/,//g;
	$twok_more_Unigene =~s/,//g;

	my $onek_Unigenes = $onek_twok_Unigene + $twok_more_Unigene;
	my $onek_percent = sprintf("%.2f",$onek_Unigenes/$All_Unigenes*100)."%";

	## $All_Anno
	my $All_Anno = $table_info{"Unigeneע��ͳ�Ʊ�"}{"table"}->[-1][1];
	$All_Anno =~s/,//g;



	%data_among_words = (
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
		"\\\$All_Anno"=>format_figure($All_Anno)
		);

}
sub get_picture {#
	## get info
	my $qm_map=(glob "$dir_data/cleandata/PNG/*.quality.png")[0];
	my $unigene_map=(glob "$dir_data/Unigene/*.distribution.png")[0];
    my $insertSize_map = (glob "$dir_data/geneExpression/*.insertSize.png")[0];
    my $insertSize_map_r = (glob "$dir_data/geneExpression/*.insertSize.r.png")[0];
    my $ssr_density_map=(glob "$dir_data/Unigene/Unigene_SSR/*.ssr.density.png")[0];

	%pic_info=(
		"ת¼�����ʵ������ͼ"=>["$dir_template/P01_RNA-Seq_experimental_workflow.png",[ 100, 100]],
		"ת¼�����������Ϣ��������ͼ"=>["$dir_template/P02_RNA-Seq_analysis_workflow.png",[ 100, 100]],
        "FASTQ��ʽ�ļ�ʾ��ͼ"=>["$dir_template/P03_FASTQ_format.png",[ 100, 100]],
		"�������ֵ�ֲ�ͼ"=>["$qm_map", [45,45]],
        "Trinity��װ����ԭ��ͼ"=>["$dir_template/P05_Trinity_workflow.png",[ 100, 100]],
        "Unigene���ȷֲ�ͼ"=>[$unigene_map,[ 53, 53]],
		"Mapped Reads��mRNA�ϵ�λ�÷ֲ�ͼ"=>["$dir_data/geneExpression/Total.randcheck.png",[ 50, 50]],
#       "����Ƭ�γ���ģ��ֲ�ͼ"=>[$insertSize_map,[ 60, 60]],
        "ת¼��������ݱ��Ͷ�ģ��ͼ"=>["$dir_data/geneExpression/Total.gene_tag.png",[ 60, 60]],
		"CDS��������ļ�ʾ��ͼ"=>["$dir_template/P10_CDS_example.png",[ 100, 100]],
		"SSR�ܶȷֲ�ͼ"=>["$ssr_density_map",[ 60, 60]],

		"����ƷFPKM�ܶȷֲ��Ա�ͼ"=>["$dir_data/geneExpression/all.fpkm_density.png",[ 70, 70]],
		"����ƷFPKM����ͼ"=>["$dir_data/geneExpression/all.fpkm_box.png",[ 70, 70]],
#		"��Ʒ�������ͼ"=>["$dir_data/geneExpression/sample_cluster.png",[ 70, 70]],
	);

    if (defined $insertSize_map_r) {
        $pic_info{"����Ƭ�γ���ģ��ֲ�ͼ"} = [$insertSize_map_r,[60,60]];
    } else {
        $pic_info{"����Ƭ�γ���ģ��ֲ�ͼ"} = [$insertSize_map,[50,50]];
    }
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


sub get_formula {#
	## get formula
	## get info
	%formula_info=(
        '����ֵ���㹫ʽ'=>["$dir_template/F01_Qscore_formula.png",[ 100, 100]],
        'FPKM���㹫ʽ'=>["$dir_template/F02_FPKM_formula.png",[ 100, 100]],
        '�������Ӽ��㹫ʽ'=>["$dir_template/F03_enrichment_factor_formula.png",[ 100, 100]],
	);
}

# ------------------------------------------------------------------
# get_table
# ------------------------------------------------------------------

sub get_table {#
    ## ��1 �������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��
    @{$table_info{"�������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��"}{"width"}}= (2600,3608,2600);
    @{$table_info{"�������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��"}{"align"}}= ('c','c','c');
    push @{$table_info{"�������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��"}{"table"}},[("Phred Quality Score","Probability of Incorrect Base Call","Base Call Accuracy")];
    push @{$table_info{"�������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��"}{"table"}},[("Q10","1/10","90%")];
    push @{$table_info{"�������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��"}{"table"}},[("Q20","1/100","99%")];
    push @{$table_info{"�������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��"}{"table"}},[("Q30","1/1000","99.9%")];
    push @{$table_info{"�������ֵ����ʶ�����ĸ��ʵĶ�Ӧ��ϵ��"}{"table"}},[("Q40","1/10000","99.99%")];

	## ��2 ��Ʒ������������ͳ�Ʊ�
	# width
	$ftable = "$dir_data/cleandata/AllSample_GC_Q.stat";
	open (TAB,$ftable) or die $!;
	$/="\n";
	my $l1=<TAB>;
	$l1=~s/\s+$//;
	my @l1=split/\s+/,$l1;
	my $limit_1=@l1;

	if ($limit_1!=8) {
		print "Note your file: $ftable !";die;
	}
	if ($limit_1==8) {
		@{$table_info{"��Ʒ������������ͳ�Ʊ�"}{"width"}}= (1350,1300,1600,1600,1450,1450);
		@{$table_info{"��Ʒ������������ͳ�Ʊ�"}{"align"}}= ('c','c','c','c','c','c');

		push @{$table_info{"��Ʒ������������ͳ�Ʊ�"}{"table"}},[("Samples","BMK-ID","Read Number","Base Number","GC Content",T("%��Q30"))];
		# info
        #SampleID       ReadSum BaseSum GC(%)   N(%)    Q20(%)  CycleQ20(%)     Q30(%)
        #T1      18723047        3781754973      52.90   0.00    97.03   100.00  92.64
		while (<TAB>) {
			chomp;
			next if(/^$/ or /^#/);
			my @info = split /\s+/,$_;
			my @table = ($info[0],format_figure($info[1]),format_figure($info[2]),format_figure($info[3])."%",format_figure($info[7])."%");
			unshift @table,'XXX';
			push @{$table_info{"��Ʒ������������ͳ�Ʊ�"}{"table"}},[@table];
		}
		close (TAB) ;
	}

	## ��3 ��װ���ͳ�Ʊ�
	if (-d "$dir_data/Assembly/All_Combination") {
		@{$table_info{"��װ���ͳ�Ʊ�"}{"width"}}=(2150,2200,2200,2200);
		@{$table_info{"��װ���ͳ�Ʊ�"}{"align"}}=('c','c','c','c');

		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[0][0]="Length Range";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[1][0]="200-300";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[2][0]="300-500";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[3][0]="500-1000";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[4][0]="1000-2000";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[5][0]="2000+";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[6][0]="Total Number";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[7][0]="Total Length";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[8][0]="N50 Length";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[9][0]="Mean Length";

		for (my $i=0; $i<3; $i++) {
			if ($i==2) { ## Unigene
				$ftable = get_only_one_file("$dir_data/Unigene/*.stat.xls");

			}
			elsif($i==1){  ## Trans
				$ftable = get_only_one_file("$dir_data/Assembly/All_Combination/Transcripts/*.stat.xls");
			}
			elsif($i==0){  ## contig
				$ftable = get_only_one_file("$dir_data/Assembly/All_Combination/contigs/*.stat.xls");
			}
			open (TAB,$ftable) or die $!;
			$/="\n";
			my $line=0;
			my $cor_limit=0;
			while (<TAB>) {
				chomp;
				next if(/^$/);
				if ($line==0) {
					if ($i==2) {
						$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[0][$i+1]="Unigene";
					}
					elsif($i==1){
						$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[0][$i+1]="Transcript";
					}
					elsif($i==0){
						$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[0][$i+1]="Contig";
					}
					$line++;
					next;
				}
				my @info2 = split /\t/,$_;
				if ($info2[0] ne $table_info{"��װ���ͳ�Ʊ�"}{"table"}->[$line][0]) {
					$cor_limit++;
					if ($cor_limit>1) {
						print "wrong assemlby statistic file: $ftable!\n";
						die;
					}
				}

				if ($line <6) {
					$info2[2] =~s/\%//;  
					$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[$line][$i+1]= format_figure($info2[1])."(".format_figure($info2[2])."%)";
				}else{

					$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[$line][$i+1]= format_figure($info2[1]);
				}

				$line++;
			}
			close (TAB) ;
		}

        $table_info{"��װ���ͳ�Ʊ�"}{"table"}->[1][1] .= "*";
	}
	else {
		@{$table_info{"��װ���ͳ�Ʊ�"}{"width"}}=(2202,2202,2202,2202);
		@{$table_info{"��װ���ͳ�Ʊ�"}{"align"}}=('w','e','e','e');

		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[0][0]="Length Range";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[1][0]="200-300";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[2][0]="300-500";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[3][0]="500-1000";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[4][0]="1000-2000";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[5][0]="2000+";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[6][0]="Total Number";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[7][0]="Total Length";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[8][0]="N50 Length";
		$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[9][0]="Mean Length";
		for (my $i2=0; $i2<3; $i2++) {
			my @Assembly_Groups=glob "$dir_data/Assembly/*/Unigenes/*.Unigenes.stat.xls";
			if ($i2==2) { ## ������Ʒ��ͳ��
				$ftable = get_only_one_file("$dir_data/Unigene/*.stat.xls");

			}else{  ## ����Ʒͳ��
				$ftable = $Assembly_Groups[$i2];
			}
			open (TAB,$ftable) or die $!;
			$/="\n";
			my $line2=0;
			my $name2;
			while (<TAB>) {
				chomp;
				next if(/^$/);
				if ($line2==0) {
					if ($i2==2) {
						$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[0][$i2+1]="All Unigenes";
					}else{
						($name2) =split /\./,$_;
						$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[0][$i2+1]=$name2." Unigenes";
					}
					$line2++;
					next;
				}
				my @info3 = split /\t/,$_;
				if ($info3[0] ne $table_info{"��װ���ͳ�Ʊ�"}{"table"}->[$line2][0]) {
					print "wrong assemlby statistic file: $ftable!\n";
					die;
				}

				if ($line2 <6) {
					$info3[2] =~s/\%//;  
					$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[$line2][$i2+1]= format_figure($info3[1])."(".format_figure($info3[2])."%)";
				}else{

					$table_info{"��װ���ͳ�Ʊ�"}{"table"}->[$line2][$i2+1]= format_figure($info3[1]);
				}

				$line2++;
			}
			close (TAB) ;
		}
	}


	## ��4 ������������װ����ıȶ�ͳ�Ʊ�
	@{$table_info{"������������װ����ıȶ�ͳ�Ʊ�"}{"width"}}=(1802,2402,2402,2202);
	@{$table_info{"������������װ����ıȶ�ͳ�Ʊ�"}{"align"}}=('c','c','c','c');
    push @{$table_info{"������������װ����ıȶ�ͳ�Ʊ�"}{"table"}}, [("BMK-ID","Clean Reads","Mapped Reads","Mapped Ratio")];
	my @ftable = glob "$dir_data/geneExpression/*.Mapped.stat.xls";

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
        push @{$table_info{"������������װ����ıȶ�ͳ�Ʊ�"}{"table"}}, [($sample_id, $clean_num, $mapped_num, $mapped_ratio)];
        close (TAB) ;
    }


	## ��5 Unigeneע��ͳ�Ʊ�
	@{$table_info{"Unigeneע��ͳ�Ʊ�"}{"width"}}=(2202,2202,2202,2202);
	@{$table_info{"Unigeneע��ͳ�Ʊ�"}{"align"}}=('c','c','c','c');
	$ftable = "$dir_data/Unigene/Unigene_Anno/Function_Annotation.stat.xls";
	open (TAB,$ftable) or die $!;
	$/="\n";
	while (<TAB>) {
		chomp;
		next if(/^$/);
		my ($head9,@info9) = split /\s+/,$_;
		if ($head9 =~/^\#/) {
			push @{$table_info{"Unigeneע��ͳ�Ʊ�"}{"table"}}, [("Annotated databases","Unigene",T("��300nt"),T("��1000nt"))];
		}else{
			$head9 =~s/_Annotation//;
			$head9 =~s/_Annotated//;
            $head9 =~s/Swissprot/Swiss-Prot/;
			my $new_more=$info9[1]+$info9[2];
			@info9 = (format_figure($info9[0]), format_figure($new_more), format_figure($info9[2]));
			push @{$table_info{"Unigeneע��ͳ�Ʊ�"}{"table"}},[$head9,@info9];
		}
	}
	close (TAB) ;


#	## ��6 SSR�������ͳ�Ʊ�
#	#width
#	@{$table_info{"SSR�������ͳ�Ʊ�"}{"width"}}=(6306,2502);
#	@{$table_info{"SSR�������ͳ�Ʊ�"}{"align"}}=('c','c');
#	#title
#	push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},["Searching Item","Number"];
#	$ftable = get_only_one_file("$dir_data/Unigene/Unigene_SSR/*.stat.xls");
#	open (TAB,$ftable) or die $!;
#	$/="\n";
#		while (<TAB>) {
#			chomp;
#			next if(/^$/);
#			my @info5 = split /\t+/,$_;
#			my @table5 = ($info5[0],format_figure($info5[1]));
#			push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},[@table5];
#		}
#		close (TAB) ;
    ###############################################################################
	## ��6 SSR�������ͳ�Ʊ�
	#width
	@{$table_info{"SSR�������ͳ�Ʊ�"}{"width"}}=(6306,2502);
	@{$table_info{"SSR�������ͳ�Ʊ�"}{"align"}}=('c','c');
	#title
	push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},["Searching Item","Number"];
	$ftable = get_only_one_file("$dir_data/Unigene/Unigene_SSR/*.statistics");
	open (TAB,$ftable) or die $!;
	$/="\n";
    while (<TAB>) {
        chomp;
        next unless (/^Total/ or /^Number/ or /^\d/);
        my $number;
        if (/^Total number of sequences examined/) {
            $number = (split /:/,$_)[1];
            $number =~ s/\s+//g;
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Total number of sequences examined',format_figure($number)];
        } elsif (/^Total size of examined sequences/) {
            $number = (split /:/,$_)[1];
            $number =~ s/\s+//g;
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Total size of examined sequences (bp)',format_figure($number)];
        } elsif (/^Total number of identified SSRs/) {
            $number = (split /:/,$_)[1];
            $number =~ s/\s+//g;
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Total number of identified SSRs',format_figure($number)];
        } elsif (/^Number of SSR containing sequences/) {
            $number = (split /:/,$_)[1];
            $number =~ s/\s+//g;
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Number of SSR containing sequences',format_figure($number)];
        } elsif (/^Number of sequences containing more than 1 SSR/) {
            $number = (split /:/,$_)[1];
            $number =~ s/\s+//g;
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Number of sequences containing more than 1 SSR',format_figure($number)];
        } elsif (/^Number of SSRs present in compound formation/) {
            $number = (split /:/,$_)[1];
            $number =~ s/\s+//g;
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Number of SSRs present in compound formation',format_figure($number)];
        } elsif (/^1/) {
            $number = (split /\s+/,$_)[1];
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Mono nucleotide',format_figure($number)];
        } elsif (/^2/) {
            $number = (split /\s+/,$_)[1];
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Di nucleotide',format_figure($number)];
        } elsif (/^3/) {
            $number = (split /\s+/,$_)[1];
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Tri nucleotide',format_figure($number)];
        } elsif (/^4/) {
            $number = (split /\s+/,$_)[1];
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Tetra nucleotide',format_figure($number)];
        } elsif (/^5/) {
            $number = (split /\s+/,$_)[1];
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Penta nucleotide',format_figure($number)];
        } elsif (/^6/) {
            $number = (split /\s+/,$_)[1];
            push @{$table_info{"SSR�������ͳ�Ʊ�"}{"table"}},['Hexa nucleotide',format_figure($number)];
        }
    }
    close (TAB) ;

	## ��7 SSR��������ļ�ʾ���
    # width
    @{$table_info{"SSR��������ļ�ʾ���"}{"width"}}= (2308,1000,1208,1736,840,854,854);
    @{$table_info{"SSR��������ļ�ʾ���"}{"align"}}= ('c','c','c','c','c','c','c');
    #info
    $ftable = (glob "$dir_data/Unigene/Unigene_SSR/*.misa")[0];
    open (TAB,$ftable) or die $!;
	$/="\n";
	while (<TAB>) {
		chomp;
		next if(/^$/);
        my $cn = () = $_ =~/\)/g;
        next if ($cn>1);
		my @table6 = split /\t/,$_;
        $table6[0]=~s/^#//;
		push @{$table_info{"SSR��������ļ�ʾ���"}{"table"}},[@table6];
        last if ($.>8);
	}
	close (TAB);

	## ��8 SSR������ƽ���ļ�ʾ���
    # width
    @{$table_info{"SSR������ƽ���ļ�ʾ���"}{"width"}}= (1668,1750,600,500,1750,600,500,625,625,625);
    @{$table_info{"SSR������ƽ���ļ�ʾ���"}{"align"}}= ('c','c','c','c','c','c','c','c','c','c');
    #info
    $ftable = (glob "$dir_data/Unigene/Unigene_SSR/*.ssr.primer.xls")[0];
    open (TAB,$ftable) or die $!;
	$/="\n";
	while (<TAB>) {
		chomp;
		next if(/^$/);
        my ($gene_id,$pr,$ptm,$psize,$rp,$rtm,$rsize,$size,$start,$end) = split /\t/,$_;

        if ($.==1) {
            $gene_id =~ s/^#//;
        } else {
            $ptm = sprintf("%.2f",$ptm);
            $rtm = sprintf("%.2f",$rtm);
        }

        push @{$table_info{"SSR������ƽ���ļ�ʾ���"}{"table"}},[$gene_id,$pr,$ptm,$psize,$rp,$rtm,$rsize,$size,$start,$end];
        last if ($.>8);
	}
	close (TAB);

    ## ��9 ������������ļ�ʾ���
    # width
    @{$table_info{"������������ļ�ʾ���"}{"width"}}= (1776,1072,858,858,958,2258,1008);
    @{$table_info{"������������ļ�ʾ���"}{"align"}}= ('c','c','c','c','c','c','c');
    #info
    $ftable = (glob "$dir_data/geneExpression/*.geneExpression.xls")[0];
    open (TAB,$ftable) or die $!;
	$/="\n";
	while (<TAB>) {
		chomp;
		next if(/^$/);
		my @table8 = split /\s+/,$_;
        $table8[0]=~s/^#//;
		push @{$table_info{"������������ļ�ʾ���"}{"table"}},[@table8];
        last if ($.>8);
	}
	close (TAB);

    ## ����1 ����б�
    # width
    @{$table_info{"����б�"}{"width"}}= (1201,3002,4605);
    @{$table_info{"����б�"}{"align"}}= ('c','w','w');
    #info
    push @{$table_info{"����б�"}{"table"}},['Tools','Description','Linkages'];
    push @{$table_info{"����б�"}{"table"}},['Trinity',T('ת¼��������װ���'),'http://trinityrnaseq.sourceforge.net/'];
    push @{$table_info{"����б�"}{"table"}},['TransDecoder',T('Ѱ�ұ��������е����'),'http://sourceforge.net/projects/transdecoder/'];
    push @{$table_info{"����б�"}{"table"}},['MISA',T('ʶ��΢���ǣ�SSR�������'),'http://pgrc.ipk-gatersleben.de/misa/misa.html'];
    push @{$table_info{"����б�"}{"table"}},['BLAST',T('���бȶ����'),'http://blast.ncbi.nlm.nih.gov/Blast.cgi'];
    push @{$table_info{"����б�"}{"table"}},['HMMER',T('���׽ṹ��ȶ����'),'http://hmmer.janelia.org/'];
    push @{$table_info{"����б�"}{"table"}},['topGO',T('����R��GO���ܸ��������'),'https://www.bioconductor.org/packages/2.12/bioc/html/topGO.html'];

	close (TAB);

    ## ����2 ���ݿ��б�
    # width
    @{$table_info{"���ݿ��б�"}{"width"}}= (1401,4200,3207);
    @{$table_info{"���ݿ��б�"}{"align"}}= ('c','w','w');
    #info
    push @{$table_info{"���ݿ��б�"}{"table"}},['Database','Description','Homepage'];
    push @{$table_info{"���ݿ��б�"}{"table"}},['NR',T('�����൰���������ݿ�'),'ftp://ftp.ncbi.nih.gov/blast/db/'];
    push @{$table_info{"���ݿ��б�"}{"table"}},['Swiss-Prot',T('�˹�ע�͵ķ����൰���������ݿ�'),'http://www.uniprot.org/'];
    push @{$table_info{"���ݿ��б�"}{"table"}},['GO',T('�����壨Gene Ontology�����ݿ�'),'http://www.geneontology.org/'];
    push @{$table_info{"���ݿ��б�"}{"table"}},['COG',T('Clusters of Orthologous Groups'),'http://www.ncbi.nlm.nih.gov/COG/'];
    push @{$table_info{"���ݿ��б�"}{"table"}},['KOG',T('euKaryotic Orthologous Groups'),'http://www.ncbi.nlm.nih.gov/COG/'];
    push @{$table_info{"���ݿ��б�"}{"table"}},['KEGG',T('���������������ٿ�'),'http://www.genome.jp/kegg/'];
    push @{$table_info{"���ݿ��б�"}{"table"}},['Pfam',T('���׼������ݿ�'),'http://pfam.xfam.org/'];

	close (TAB);

    ## ����3 ��������
    # width
    @{$table_info{"��������"}{"width"}}= (2202,3303,3303);
    @{$table_info{"��������"}{"align"}}= ('c','c','c');
    #info
    push @{$table_info{"��������"}{"table"}},['Nucleic Acid Code','Meaning','Mnemonic'];
    push @{$table_info{"��������"}{"table"}},['A','A','Adenine'];
    push @{$table_info{"��������"}{"table"}},['C','C','Cytosine'];
    push @{$table_info{"��������"}{"table"}},['G','G','Guanine'];
    push @{$table_info{"��������"}{"table"}},['T','T','Thymine'];
    push @{$table_info{"��������"}{"table"}},['U','U','Uracil'];
    push @{$table_info{"��������"}{"table"}},['R','A or G','puRine'];
    push @{$table_info{"��������"}{"table"}},['Y','C, T or U','pYrimidines'];
    push @{$table_info{"��������"}{"table"}},['K','G, T or U','bases which are Ketones'];
    push @{$table_info{"��������"}{"table"}},['M','A or C','bases with aMino groups'];
    push @{$table_info{"��������"}{"table"}},['S','C or G','Strong interaction'];
    push @{$table_info{"��������"}{"table"}},['W','A, T or U','Weak interaction'];
    push @{$table_info{"��������"}{"table"}},['B','not A (i.e. C, G, T or U)','B comes after A'];
    push @{$table_info{"��������"}{"table"}},['D','not C (i.e. A, G, T or U)','D comes after C'];
    push @{$table_info{"��������"}{"table"}},['H','not G (i.e., A, C, T or U)','H comes after G'];
    push @{$table_info{"��������"}{"table"}},['V','neither T nor U (i.e. A, C or G)','V comes after U'];
    push @{$table_info{"��������"}{"table"}},['N','A C G T U','Nucleic acid'];

	close (TAB);
}

# ------------------------------------------------------------------
# rtf_write
# ------------------------------------------------------------------

sub rtf_write {	#&rtf_write($format, \@context);
	my ($format, $context)=@_;
#	print $context,"\n";
	my $line= scalar @{$context};
	if ($line > 1 && $format !~ /����/ && $format ne "ע��" && $format ne "�ο�����" && $format !~/ͼƬ/ && $format ne "��ʽ") {
		printf("wrong format:$format \t\"$context->[0]\"!\n");
		die;
	}

	## exchange
	for (my $i=0; $i<$line;$i++) {
		## some data info
		my (@data_among_word) = $context->[$i] =~/(\$[a-z,A-Z,_,0-9]+)/g;
		if (@data_among_word !=0) {
			for (my $j=0;$j<@data_among_word;$j++) {
				$data_among_word[$j] =~s/\$/\\\$/; ##�滻ʱ��"$All_unigenes"ƥ�䲻�ˣ�ֻ��ƥ�䣺"\$All_unigenes"
				$context->[$i] =~s/$data_among_word[$j]/$data_among_words{$data_among_word[$j]}/;
			}
		}
	}
		

	if($format eq "���"){## ���
		&table_write($context->[0]);
    }elsif($format eq "����"){## ����
		&table_appendix_write($context->[0]);
    }elsif($format eq "��ʽ") {## ��ʽ
        &formula_write($format,$context);
	}elsif($format =~/ͼƬ/){## ͼƬ
		&picture_write($format,$context);
	}else{ #�����ϱ꣬���word
		for (my $ie=0; $ie<$line;$ie++) { ## һ��һ�δ���
			## ���ϱ� ��ֶ���
			my @super = $context->[$ie]=~/\[(\S+)\sSuperscript\]/g;
			my @part = split /\[\S+\sSuperscript\]/,$context->[$ie];

			##paragraph
			$rtf->print(\'\pard');
			for (my $je=0; $je<@part-1; $je++) {
				&format_print_word($format, $part[$je],$je);
				$rtf->print(
					\'\super',
					$super[$je],
					\'\nosupersub',
				);
			}
			&format_print_word($format, $part[@part-1], @part-1);

			$rtf->print(\'\par');
		}
	}

}

### ����ģ���ʽ������index���word
### ģ���ʽ�����ġ�n�������
### ����index�� �����ڸö���ĵ�index�Σ���Ϊ0�����Ƕ��俪ͷ�����ĺʹ��༭������Ҫ���������Ϊ0��������
sub format_print_word {#
	my ($format, $context,$index)=@_;
	## ����
	if ($format eq "һ������") {
        if ((split /\s+/,$context)[0]==1) {
#            $rtf->print(
#                \'\f2\f0\fs30\i0\sl-800\slmult1\cf0\outlinelevel0',  ##1�����⣨����С���� #����ǰ�㣺\keepn\widctlpar\wrapdefault\faauto
#                T($context)
#            );
            $rtf->print(
                \'\f2\f0\fs32\i0\sl-440\slmult1\cf0\outlinelevel0\lisa100\lisb100',  ##1�����⣨�������ţ���ǰ���κ���Ϊ1�У��о�̶�ֵ22����
                T($context)
            );
        } else {
#            $rtf->print(
#                \'\f2\f0\page\fs30\i0\sl-800\slmult1\cf0\outlinelevel0',  ##1�����⣨����С���� #����ǰ�㣺\keepn\widctlpar\wrapdefault\faauto
#                T($context)
#            );
            $rtf->print(
                \'\f2\f0\page\fs32\i0\sl-440\slmult1\cf0\outlinelevel0\lisa100\lisb100',  ##1�����⣨�������ţ���ǰ���κ���Ϊ1�У��о�̶�ֵ22����
                T($context)
            );
        }
	}elsif( $format eq "��������"){
#        $rtf->print(
#            \'\f2\f0\fs28\i0\sl-600\slmult1\cf0\outlinelevel1',  ##2�����⣨�����ĺţ�
#            T($context)
#        );
        $rtf->print(
            \'\f2\f0\fs28\i0\sl-440\slmult1\cf0\outlinelevel1\lisa100\lisb100',  ##2�����⣨�����ĺţ���ǰ���κ���Ϊ1�У��о�̶�ֵ22����
            T($context)
        );
	}elsif($format eq "��������"){
#        $rtf->print(
#            \'\f2\f0\fs24\i0\sl-600\slmult1\cf0\outlinelevel2',  ##3�����⣨����С�ģ�
#            T($context)
#        );
        $rtf->print(
            \'\f2\f0\fs24\i0\sl-440\slmult1\cf0\outlinelevel2\lisa100\lisb100',  ##3�����⣨����С�ģ���ǰ���κ���Ϊ1�У��о�̶�ֵ22����
            T($context)
        );
	}elsif($format eq "�ļ�����"){
		if ($limit==1) {
			$rtf->print(
				\'\page\fs24\b\i0\sl-600\slmult1\cf0\outlinelevel3',  ##4�����⣨����С�ļӴ֣�
				T($context)
			);
		}
		if ($limit==0) {
			$rtf->print(
				\'\fs24\b\i0\sl-600\slmult1\cf0\outlinelevel3',  ##4�����⣨����С�ļӴ֣�
				T($context)
			);
		}
		$limit=0;
	}elsif($format eq "�弶����"){
		if ($limit==1) {
			$rtf->print(
				\'\page\fs24\b\i0\sl-600\slmult1\cf0\outlinelevel4',  ##5�����⣨����С�ļӴ֣�
				T($context)
			);
		}
		if ($limit==0) {
			$rtf->print(
				\'\fs24\b\i0\sl-600\slmult1\cf0\outlinelevel4',  ##5�����⣨����С�ļӴ֣�
				T($context)
			);
		}
		$limit=0;
	}elsif($format =~/����/){ 
		if ($format eq "����") {
            if ($context=~/��ʽ��/) {
                $rtf->print(
                    \'\f1\f0\fs24\b0\i0\sl-500\slmult0\cf0\qj',  ##���ģ�����С�ģ�
                    T($context)
                );
            } else {
                if ($context=~/\[[^\[]+\sItalic\]/) {   ## ����б��
                    my @italic = $context =~/\[([^\[]+)\sItalic\]/g;
                    my @part = split /\[[^\[]+\sItalic\]/,$context;

                    for (my $je=0; $je<@part-1; $je++) {
                        $rtf->print(
                            \'\f1\f0\fs24\b0\i0\sl-500\slmult0\cf0\cufi200\qj',
                            T($part[$je]),
                            \'\i',
                            T($italic[$je]),
                        );
                    }
                    $rtf->print(
                        \'\f1\f0\fs24\b0\i0\sl-500\slmult0\cf0\cufi200\qj',
                        T($part[@part-1]),
                    );
                } else {
                    $rtf->print(
                        \'\f1\f0\fs24\b0\i0\sl-500\slmult0\cf0\cufi200\qj',  ##���ģ�����С�ģ�
                        T($context)
                    );
                }
            }
		}elsif( $format eq "���༭����"){
            $rtf->print(
                \'\f1\f0\fs24\b0\i0\sl-500\slmult0\cf1\cufi200\qj',  ##���༭���ģ�����С�� ��ɫ�� #\cf Foreground color (the default is 0).
                T($context)
            );
		}
	}elsif($format eq "ע��"){
		$limit=1;
		$rtf->print(
			\'\f1\f0\fs18\b0\i0\sl-350\slmult0\cf0\qj',  ##ע�ͣ�����С�壩
			T($context)
		);
	}elsif($format eq "�ο����ױ���"){
		$rtf->print(
			\'\f2\f0\sect\sbkpage\fs28\qc\cf0',  ##�ο����ױ��⣨�����ĺž��У� #\sect\sbkpage ��ӷֽڷ�
			\'\par',T($context),\'\par'
		);
	}elsif($format eq "�ο�����"){
        if ($context=~/\[[^\[]+\sItalic\]/) {   ## ����б��
            my @italic = $context =~/\[([^\[]+)\sItalic\]/g;
            my @part = split /\[[^\[]+\sItalic\]/,$context;

            for (my $je=0; $je<@part-1; $je++) {
                $rtf->print(
                    \'\f1\f0\fs21\b0\i0\sl-500\slmult0\cf0\cufi-200\culi0\qj',
                    $part[$je],
                    \'\i',
                    $italic[$je],
                );
            }
            $rtf->print(
                \'\f1\f0\fs21\b0\i0\sl-500\slmult0\cf0\cufi-200\culi0\qj',
                $part[@part-1],
            );
        } else {
            $rtf->print(
                \'\f1\f0\fs21\b0\i0\sl-500\slmult0\cf0\cufi-200\culi0\qj',  ##ע�ͣ���ţ� #\cufi-200\culi0 ������������ #\qj ˫�˶���
                $context
            );
        }
	} elsif ($format eq "��¼����") {
		$rtf->print(
			\'\f2\f0\sect\sbkpage\fs28\ql\cf0',  ##�ο����ױ��⣨�����ĺž��У� #\sect\sbkpage ��ӷֽڷ�
			T($context)
		);
    }
}

sub table_write {#
	my $context = shift;
	my (undef, $name)=split /\s+/,$context,2;
	if (!exists $table_info{$name}{"table"} || !exists $table_info{$name}{"width"}) {
		print $name,"\n";
		die;
	}

	my @table = @{$table_info{$name}{"table"}};
	my @width = @{$table_info{$name}{"width"}};
	my @align = @{$table_info{$name}{"align"}};
	my $row_num = @table;
	my $column_num = @width;

	## title
	$rtf->paragraph(
		\'\f2\f0\fs21\qc\sl-500\slmult0', ##ͼ����⣨������ž��У�
		T($context)
	);
	## border
	my @decl=();
	# first row
	$decl[0] = RTF::Writer::TableRowDecl->new(
		widths => [@width],  ## ҳ��һ��8808
		borders=>["n-15-th s-15-s","n-15-th s-15-s","n-15-th s-15-s","n-15-th s-15-s","n-15-th s-15-s"],
		align=>['c','c','c','c','c','c','c','c','c','c'],
	);
	# middle row 
	for (my $i=1; $i< $row_num-1; $i++) {
		$decl[$i]= RTF::Writer::TableRowDecl->new(
			widths => [@width],  ## ҳ��һ��8808
			borders=>[],
		    align=>[@align],
		);
	}
	# last row
	$decl[$row_num-1] = RTF::Writer::TableRowDecl->new(
		widths => [@width],  ## ҳ��һ��8808
		borders=>["s-15-th","s-15-th","s-15-th","s-15-th","s-15-th"],
		align=>[@align],
	);
	
	## write

	for (my $io=0; $io<$row_num; $io++) {
		my @info =();
		for (my $j=0;$j<$column_num; $j++) {
			if ($io==0) {
				push @info, [\'\b\fs21\trrh580\sl800\slmult0',$table[$io][$j]] ;  #\rtlch\fcs1 \af50 \ltrch\fcs0 \insrsid12659029 \hich\af50\dbch\af31505\loch\f50
#				push @info, [\'\b\fs21\trrh580\sl800\slmult0\rtlch\fcs1 \af50 \ltrch\fcs0 \insrsid12659029 \hich\af50\dbch\af31505\loch\f50',$table[$io][$j]] ;  #\rtlch\fcs1 \af50 \ltrch\fcs0 \insrsid12659029 \hich\af50\dbch\af31505\loch\f50
			}else{
				push @info, [\'\fs21\slmult0\intbl',$table[$io][$j]] ;
#				push @info, [\'\fs21\sl800\slmult0\rtlch\fcs1 \af50 \ltrch\fcs0 \insrsid12659029 \hich\af50\dbch\af31505\loch\f50',$table[$io][$j]] ;
			}
		}

		$rtf->row_bmk($decl[$io], @info); 
	}
}

sub table_appendix_write {#
	my $context = shift;
	my (undef, $name)=split /\s+/,$context,2;
	if (!exists $table_info{$name}{"table"} || !exists $table_info{$name}{"width"}) {
		print $name,"\n";
		die;
	}

	my @table = @{$table_info{$name}{"table"}};
	my @width = @{$table_info{$name}{"width"}};
	my @align = @{$table_info{$name}{"align"}};
	my $row_num = @table;
	my $column_num = @width;

	## title
    $rtf->paragraph(
        \'\f2\f0\fs21\ql\sl-500\slmult0', ##ͼ����⣨������ž��У�
        T($context)
    );

	## border
	my @decl=();
	# first row
	$decl[0] = RTF::Writer::TableRowDecl->new(
		widths => [@width],  ## ҳ��һ��8808
		borders=>["n-15-th s-15-s","n-15-th s-15-s","n-15-th s-15-s","n-15-th s-15-s","n-15-th s-15-s"],
		align=>['c','c','c','c','c','c','c','c','c','c'],
	);
	# middle row 
	for (my $i=1; $i< $row_num-1; $i++) {
		$decl[$i]= RTF::Writer::TableRowDecl->new(
			widths => [@width],  ## ҳ��һ��8808
			borders=>[],
		    align=>[@align],
		);
	}
	# last row
	$decl[$row_num-1] = RTF::Writer::TableRowDecl->new(
		widths => [@width],  ## ҳ��һ��8808
		borders=>["s-15-th","s-15-th","s-15-th","s-15-th","s-15-th"],
		align=>[@align],
	);
	
	## write

	for (my $io=0; $io<$row_num; $io++) {
		my @info =();
		for (my $j=0;$j<$column_num; $j++) {
			if ($io==0) {
				push @info, [\'\f1\f0\b\fs21\slmult0',$table[$io][$j]] ;  #\rtlch\fcs1 \af50 \ltrch\fcs0 \insrsid12659029 \hich\af50\dbch\af31505\loch\f50
			}else{
				push @info, [\'\f1\f0\fs21\slmult0\intbl',$table[$io][$j]] ;
			}
		}

		$rtf->row($decl[$io], @info); 
	}
}

sub picture_write{# 
	my ($format,$context) = @_;

	my @format =split /\s+/,$format;

	my $line= scalar @{$context};
	my (undef, $name)=split /\s+/,$context->[$line-1],2;

	my @pic_info = @{$pic_info{$name}};
	
	## check
	if (@pic_info%2 !=0 ) {
		print "Wrong! pic_info/2 !=0\t$context\n";
		die;
	}

	## draw
	if (@format == 1 || $format[1] == 1) { ## 1��ͼƬ
		$rtf->paragraph(
			\'\qc',
			$rtf->image( 'filename' => $pic_info[0],'scalex'=>$pic_info[1][0], 'scaley'=>$pic_info[1][1]),
		);
	}else{
		if ($format[1] == 2) {
			$rtf->paragraph(
				\'\qc',
				$rtf->image( 'filename' => $pic_info[0],'scalex'=>$pic_info[1][0], 'scaley'=>$pic_info[1][1]),
				$rtf->image( 'filename' => $pic_info[2],'scalex'=>$pic_info[3][0], 'scaley'=>$pic_info[3][1]),
			);
		}elsif($format[1] == 3) {
			$rtf->paragraph(
				\'\qc',
				$rtf->image( 'filename' => $pic_info[0],'scalex'=>$pic_info[1][0], 'scaley'=>$pic_info[1][1]),
				$rtf->image( 'filename' => $pic_info[2],'scalex'=>$pic_info[3][0], 'scaley'=>$pic_info[3][1]),
				$rtf->image( 'filename' => $pic_info[4],'scalex'=>$pic_info[5][0], 'scaley'=>$pic_info[5][1]),
			);
		}elsif($format[1] == 4) {
			$rtf->paragraph(
				\'\qc',
				$rtf->image( 'filename' => $pic_info[0],'scalex'=>$pic_info[1][0], 'scaley'=>$pic_info[1][1]),
				$rtf->image( 'filename' => $pic_info[2],'scalex'=>$pic_info[3][0], 'scaley'=>$pic_info[3][1]),
				$rtf->image( 'filename' => $pic_info[4],'scalex'=>$pic_info[5][0], 'scaley'=>$pic_info[5][1]),
				$rtf->image( 'filename' => $pic_info[6],'scalex'=>$pic_info[7][0], 'scaley'=>$pic_info[7][1]),
			);
		}
	}
	for (my $i=0; $i<$line; $i++) {
		$rtf->paragraph(
			\'\f2\f0\fs21\qc\sl-500\slmult0', ##ͼ����⣨������ž��У�
			T($context->[$i])
		);
	}

}

sub formula_write {#
	my ($format,$context) = @_;

	my @format =split /\s+/,$format;

	my $line= scalar @{$context};
	my (undef, $name)=split /\s+/,$context->[$line-1],2;

	my @formula_info = @{$formula_info{$name}};
	
	## check
	if (@formula_info%2 !=0 ) {
		print "Wrong! formula_info/2 !=0\t$context\n";
		die;
	}

	## draw
	if (@format == 1 || $format[1] == 1) { ## 1��ͼƬ
		$rtf->paragraph(
			\'\qc',
			$rtf->image( 'filename' => $formula_info[0],'scalex'=>$formula_info[1][0], 'scaley'=>$formula_info[1][1]),
		);
	}
}

sub cover_page_write {#
	my $rtf = shift;

	my $report_name  = "XXX��Ŀ���ⱨ��";
	my $customer_unit= "                                      ";
	my $company_name = "��������������Ƽ����޹�˾            ";
	my $contact_name = "                                      ";
	my $phone =		   "010-570450XX                         ";
	my $fax =		   "010-57045001                          ";
	my $date         = &GetDate."                        ";
	my $in_charge =    "               ";
	my $verifier =     "               ";

	## logͼƬ
	$rtf->paragraph(
		$rtf->image( 'filename' => "$dir_template/logo.png",'scalex'=>45, 'scaley'=>45),
	);
	## ����
	$rtf->paragraph(
        \'\fs20\qc\par\par\par\par',
	);
	### �������
	$rtf->paragraph(
		\'\f3\f0\fs72\i0\qc',  # ������⣺����36����С����,�Ӵ�:bold, i0: ��б��, qc������
		T($report_name)
	);

	## ����
	$rtf->paragraph(
        \'\fs36\qc\par\par\par\par\par\par',
	);

	## ���浥λ����Ϣ
	$rtf->paragraph(
		\'\f2\f0\fs28\i0\ul0\sl-550\slmult4\li1000',    #����ʵ�־���
		T('�ͻ���λ��'),
		\'\f1\f0\fs28\b0\i0\ul',
		T($customer_unit),\'\par',

		\'\f2\f0\fs28\i0\ul0',
		T('���浥λ��'),
		\'\f1\f0\fs28\b0\i0\ul',
		T($company_name),\'\par',
		
		\'\f2\f0\fs28\i0\ul0\sl-550\slmult4',
		T('�� ϵ �ˣ�'),
		\'\f1\f0\fs28\b0\i0\ul',
		T($contact_name),\'\par',

		\'\f2\f0\fs28\i0\ul0',   
		T('��ϵ�绰��'),
		\'\f1\f0\fs28\b0\i0\ul',
		T($phone),\'\par',

		\'\f2\f0\fs28\i0\ul0',   
		T('��    �棺'),
		\'\f1\f0\fs28\b0\i0\ul',
		T($fax),\'\par',

		\'\f2\f0\fs28\i0\ul0',   
		T('�������ڣ�'),
		\'\f1\f0\fs28\b0\i0\ul',
		T($date),\'\par',
		
		\'\f2\f0\fs28\i0\ul0',   
		T('��Ŀ������:'),
		\'\f1\f0\fs28\b0\i0\ul',
		T($in_charge),

		\'\f2\f0\fs28\i0\ul0',   
		T('����ˣ�'),
		\'\f1\f0\fs28\b0\i0\ul',
		T($verifier),

	);

}

sub content_page {
	my $rtf = shift;
    $rtf->Page();

    $rtf->print(
        \'\f2\f0\fs32\qc\cf0\par',
        T("Ŀ  ¼"),\'\f1\f0\fs24\par\par',
        T("�����뷽��������->Ŀ¼->����Ŀ¼�������ֺ����ã���ʼ->���壨����+Times New Roman�����ֺţ�С�ģ����о�22������ǰ���κ�0.5�С���"),
        \'\qj\par\par\sect\sbkpage'
    );

    ## ҳü
    $rtf->paragraph(
        \'\header\fi-1800', #\fi-1800 ��ͷ
        $rtf->image( 'filename' => "$dir_template/header.png",'scalex'=>100, 'scaley'=>100,),
    );

    $rtf->paragraph(
        \'\footer\qc\chpgn',
    );

}

sub T {
   my $text = shift;
   return decode( 'gb2312', $text);
}

#������ʽ��3λһ������
sub Integer_Three_Digit{#
	my $interger = shift;
	$interger=~s/(?<=\d)(?=(\d\d\d)+$)/,/g;
	return $interger;
}

#������ʽ��3λһ������
#С����ʽ��С�������λ
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

sub AbsolutePath {		#��ȡָ��Ŀ¼���ļ��ľ���·��
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
	return sprintf("%4d��%02d��%02d��", $year+1900, $mon+1, $day);
}

sub USAGE {
        my $usage=<<"USAGE";
#-----------------------------------------------------------------------------------------
  Program: $Script
  Version: $version
  Contact: zeng huaping<zenghp\@biomarker.com.cn>
     Date: 
 Modifier: Simon Young  <yangxh\@biomarker.com.cn>
     Date: 2014-08-29
 Modified: Based-v1.0.0, optimizing layout and fonts:
           (1) add first-line indent, hanging indent, specify justified alignment;
           (2) add various fonts, title in "������κ", outlines in "����", and so on;
           (3) add header and footer, add section break at proper place;
           (4) can deal with italic; adjust table row height.
           (5) ...

    Usage:
      Options:
      --id <dir>   directory of noref trans anayslsis result.
      --pf <str>   prefix of output file.

      --od <dir>   output directory        [./]

#-----------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
