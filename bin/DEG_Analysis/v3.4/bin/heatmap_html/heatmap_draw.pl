#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0";

my($in,$out,$scale,$cluster,$color,$log,$line,$col,$id_file,$div,$Anno_dir,$sRNA);
GetOptions(
	"help|?"	=>\&USAGE,
	"infile:s"		=>\$in,
	"outfile:s"		=>\$out,
	"id:s"			=>\$id_file,
	"scale:s"		=>\$scale,
	"cluster:s"		=>\$cluster,
	"color_type:s"		=>\$color,
	"is_log"		=>\$log,
	"line:s"		=>\$line,
	"col:s"			=>\$col,
	"div"			=>\$div,
	"anno:s"		=>\$Anno_dir,
	"sRNA:s"		=>\$sRNA,
) or &USAGE;
&USAGE unless ($in and $out);
&USAGE unless (($div && -d $Anno_dir) || (!$div && !-d $Anno_dir));
$color ||= 1;
die "-color must be int 1-6." unless ($color == 1 || $color == 2 || $color == 3 || $color == 4 || $color == 5 || $color == 6);
die "-cluster must be row or column or both or none." unless ($cluster eq "row" || $cluster eq "column" || $cluster eq "both" || $cluster eq "none");
die "-scale must be row or column or none." unless ($scale eq "row" || $scale eq "column" || $scale eq "none");
$out = &OutFileCheck($out);



###########font for row: font=>rows
my %row=(
#1=>600,
#2=>300,
#3=>200,
4=>150,
5=>120,
6=>100,
7=>86,
8=>75,
9=>66,
10=>60,
11=>54,
12=>50,
13=>46,
14=>43,
15=>40,
#16=>37,
#17=>35,
#18=>33,
#19=>31,
#20=>30
);

###########font for column: font=>column
my %column=(
#1=>400,
#2=>200,
#3=>133,
4=>100,
5=>80,
6=>66,
7=>57,
8=>50,
9=>44,
10=>40,
11=>36,
12=>33,
13=>30,
14=>28,
15=>26,
#16=>25,
#17=>23,
#18=>22,
#19=>21,
#20=>20
);



################################################################################################################
#main code
################################################################################################################
#读入ID文件
my %id;
my $id_opt=0;
if ($id_file && -f $id_file) {
	open ID,"$id_file";
	while (<ID>) {
		chomp;
		next if (/^\s*$/ or /^#/);
		my @line = split "\t",$_;
		$id{$line[0]} = 1;
		$id_opt=1;
	}
	close ID;
}


#提取列数据
if ($col && $col =~ /^(\d+(\-\d+)?)(\,(\d+(\-\d+)?))*$/) {
	`cut -f 1,$col $in > $out.data.cutted`;
	$in = "$out.data.cutted";
}


#检查输入文件格式
my $title;
my $title_sep_num = 0;
my $data_col_num;
my %data;
my %uniq;
my %sep;
my $sep = "\t";
open IN,"$in";
while (<IN>) {
	chomp;
	s/^\s*//;#去掉行首空白符
	s/\s*$//;#去掉行尾空白符
	#记录表头，如存在\t，则将\t+换为\t,并计数。如没有\t，则将“ +”换为\t,并计数,将分隔符$sep换为空格。\t与空格混合的情况，等确定数据内容格式后在处理。
	if ($.==1) {
		if (/^\s*$/) {
			&LogErrorFile("Title is empty, please check input file.\n");
			&ERROR_AND_DIE("Input file error, empty title.");
		}
		$title = $_;
		$title_sep_num += $title =~ s/\t+/\t/g;
		if ($title_sep_num == 0) {
			$title_sep_num += $title =~ s/ +/\t/g;
			$sep = " ";
		}
		if ($title_sep_num == 0) {
			&LogErrorFile("No tab or space in title, please check input file.\n");
			&ERROR_AND_DIE("Input file error, illegal title.");
		}
		next;
	}
	#跳过除第一行表头外，剩下以#开头的注释行，舍弃这些行。
	next if (/^\s*$/ or /^#/);
	#去除重复行
	if (defined $uniq{$_}) {
		next;
	}
	$uniq{$_} = 1;
	#按分隔符$sep将每行数据切开。如果存在行id相同，但数据内容不同的情况，则报错。生产错误日志文件。因为无法确定哪条数据正确。
	my @data = split "$sep+",$_;
	if ($id_opt == 1) {#根据ID文件筛选数据，优先级低于按前XX行。
		my $sep2 = $sep;
		$sep2 =~ tr/\t / \t/;
		my @data2 = split "$sep+",$_;
		next unless (defined $id{$data[0]} or defined $id{$data2[0]});
	}
	if (defined $data{$data[0]}) {
		&LogErrorFile("There are same id with different following data:\n$data{$data[0]}\n$_\n");
		&ERROR_AND_DIE("Input file error, repetitive column id.");
	}
	$data{$data[0]} = $_;
	#记录每行中" +"和"\t+"的数量。
	my $space = $_ =~ / +/g;
	$space ||= 0;
	$sep{s}{$space}{$data[0]} = 1;
	my $tab = $_ =~ /\t+/g;
	$tab ||= 0;
	$sep{t}{$tab}{$data[0]} = 1;
	if ($line && $. > $line) {
		last;
	}
}
close IN;

##遍历%sep，判断能否用\t或空格作为分隔符去得到正确的矩阵数据（整数，浮点数，科学计数法），如不能，则报错。生成错误日志文件。
##优先级： \t > 空格
my @tab_type = keys %{$sep{t}};
my @space_type = keys %{$sep{s}};
my $data_column = 0;
if (@tab_type == 1 && $tab_type[0] != 0) {
	&CheckDataFormat("\t+",$tab_type[0]);
	if ($title_sep_num >= $tab_type[0]) {
		$title =~ s/ +/_/g;
		$title =~ s/^([^\t](\t[^\t]+){$tab_type[0]}).*$/$1/;
	}
	else {
		my $space_num = $title =~ m/ +\S/g;
		$space_num ||= 0;
		if ($space_num + $title_sep_num >= $tab_type[0]) {
			for my $i (1..$tab_type[0]-$title_sep_num) {
				$title =~ s/ +/\t/;
			}
			$title =~ s/ +/_/g;
		}
		else {
			my $times = $title =~ s/ +(\S)/\t$1/g;
			$times ||= 0;
			for my $i (1..$tab_type[0]-$title_sep_num-$times) {
				$title .= "\tNone";
			}
			$title =~ s/ +//g;
		}
	}
	$data_column = $tab_type[0];
}
elsif (@space_type == 1 && $space_type[0] != 0) {
	&CheckDataFormat(" +",$space_type[0]);
	if ($title_sep_num >= $tab_type[0]) {
		$title =~ s/ +/_/g;
		$title =~ s/^([^\t](\t[^\t]+){$tab_type[0]}).*$/$1/;
	}
	else {
		my $times = $title =~ s/ +(\S)/\t$1/g;
		$times ||= 0;
		for my $i (1..$tab_type[0]-$title_sep_num-$times) {
			$title .= "\tNone";
		}
		$title =~ s/ +//g;
	}
	$data_column = $space_type[0];
}
else {
	&LogErrorFile("The data can not be separated as a matrix by space or tab.");
	&ERROR_AND_DIE("Input file error, not a matrix data.");
}


#输出校验后的表头和数据矩阵
my $data_line=0;
open OUT,">$out.data";
print OUT "$title\n";
foreach my $key (sort keys %data) {
	print OUT "$data{$key}\n";
	$data_line++;
}
close OUT;



#画图
$scale||="none";#标准化，按行/列：		row/column/none
$cluster||="both";#聚类，按行/列：	row/column/both/none

my $cmd = "Rscript $Bin/pheatmap.r ";
$cmd .= " --infile $out.data --outfile $out --scale $scale --color.type $color --file.type 1 --legend ";
$cmd .= " --is.log " if ($log);

if ($cluster eq 'both') {
	$cmd .= " --cluster_rows --cluster_cols ";
}
elsif ($cluster eq 'row') {
	$cmd .= " --cluster_rows ";
}
elsif ($cluster eq 'column') {
	$cmd .= " --cluster_cols ";
}

if ($div) {
	$cmd .= " --div.clust ";
}

my $row_font=0;
foreach my $r_font (sort {$a<=>$b} keys %row) {
	if ($data_line < $row{$r_font}) {
		$row_font=$r_font;
	}
	else {
		last;
	}
}
unless ($row_font==0) {
	$cmd .= " --show_rownames --fontsize_row $row_font ";
}

my $col_font=0;
foreach my $c_font (sort {$a<=>$b} keys %column) {
	if ($data_column < $column{$c_font}) {
		$col_font=$c_font;
	}
	else {
		last;
	}
}
unless ($col_font==0) {
	$cmd .= " --show_colnames --fontsize_col $col_font ";
}
if ($data_column >= 10 && $data_line <= 60) {
	$cmd .= " --width 4000 ";
}
#print $cmd,"\n";
&cmd_call($cmd);
#`rm $out.data`;
if ($col && $col =~ /^(\d+(\-\d+)?)(\,(\d+(\-\d+)?))*$/) {
#	`rm $out.data.cutted`;
}
#


my $count = 0;
my %srna;
my %gene;
if ($div) {
	&cmd_call("convert -resize 14\%x14\% $out $out.png ");
	if ($sRNA) {
		$cmd = "perl $Bin/plot_expression_patterns2.pl $out.clustered.data_subcluster_\*.data --title miRNAs >/dev/null 2>&1 ";
	}
	else {
		$cmd = "perl $Bin/plot_expression_patterns2.pl $out.clustered.data_subcluster_\*.data >/dev/null 2>&1 ";
	}
	&cmd_call($cmd);
	
	my @div_data = glob "$out.clustered.data_subcluster_\*.data";
	foreach my $div_data (@div_data) {
		if ($sRNA) {
			die "$sRNA does no exists!" unless (-f $sRNA);
			if ($count == 0) {
				open SRNA,"$sRNA";
				while (<SRNA>) {
					chomp;
					next if (/^#/ or /^\s*$/);
					my @srna = split /\s+/,$_;
					$srna{$srna[0]} = $srna[1];
				}
				$count = 1;
				close SRNA;
			}
			%gene = ();
			open DATA,"$div_data";
			open GENE,">$div_data.gene.data";
			print GENE "#Gene_ID\n";
			while (<DATA>) {
				chomp;
				next if (/^\s*$/ or $. == 1);
				my @data = split /\s+/,$_;
				if (defined $srna{$data[0]}) {
					my @gene = split ";",$srna{$data[0]};
					foreach my $gene (@gene) {
						next if (defined $gene{$gene} && $gene{$gene}==1);
						print GENE "$gene\n";
						$gene{$gene} = 1;
					}
				}
			}
			close DATA;
			close GENE;
			$div_data = "$div_data.gene.data";
		}
		$cmd = "perl $Bin/go_enrichment/go_enrichment.pl -i $div_data -dir $Anno_dir -od $div_data.go_enrichment ";
		&cmd_call($cmd);
#		sleep(1);
		my $go_result = "$div_data.go_enrichment/go_enrichment/go_enrichment.list";
		if (!-f $go_result) {
			&cmd_call("echo \"There is no result of GO enrichment.\" > $div_data.go_enrichment.list");
		}
		else {
			my $go_num = (split " ",`wc -l $go_result`)[0];
			if (defined $go_num and $go_num > 1) {
				my %go_enrich = ();
				open INGO,"$go_result";
				while (<INGO>) {
					chomp;
					next if (/^\s*$/ or $. == 1);
					my @line = split /\t/;
					$go_enrich{$line[0]} = $line[3];
				}
				close INGO;
				
				open GO,">$div_data.go_enrichment.list";
				print GO "GO_Term\tP-value\n";
				my $go_count = 0;
				foreach my $go_term (sort {$go_enrich{$a}<=>$go_enrich{$b}} keys %go_enrich) {
					print GO "$go_term\t$go_enrich{$go_term}\n";
					$go_count ++;
					last if $go_count >= 5;
				}
				close GO;
			}
			else {
				&cmd_call("echo \"There is no result of GO enrichment.\" > $div_data.go_enrichment.list");
			}
		}
	}
}

if ($div) {
	my $out_dir = dirname($out);
	$cmd = "perl $Bin/make_html/make_heatmap_html.pl -dir $out_dir";
	&cmd_call($cmd);

	$cmd = "perl $Bin/make_html/make_heatmap_html_cloud.pl -dir $out_dir";
	&cmd_call($cmd);
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################


################################################################################################################
#sub functions
################################################################################################################

sub CheckDataFormat {#&CheckDataFormat($sep,$data_num) $sep为分隔符，$data_num为匹配到的分隔符的数量，也等于除第一列ID外，剩余数据的列数。
	my ($sep,$data_num) = @_;
	my $check = 1;
	my $err_content;
	my $zero_count = 0;
	$data_num+=1;
	foreach my $key (keys %data) {
		my @data = split "$sep",$data{$key};
		for my $i (1..$#data) {
			if ($sep eq '\t+') {
				$data[$i] =~ s/^ +//;
				$data[$i] =~ s/ +$//;
			}
			elsif ($sep eq ' +') {
				$data[$i] =~ s/^\t+//;
				$data[$i] =~ s/\t+$//;
			}
			
			unless ($data[$i] =~ /^\-?\d+(\.\d+)?$|^\d(\.\d+)?[Ee]\-?\d+$/) {
				$check = 0;
				$err_content .= $data{$key};
				last;
			}
			
			if ($data[$i] =~ /^\d(\.\d+)?[Ee]\-(\d+)$/) {
				if ($2 >= 8) {
					$data[$i] = 0;
				}
			}
			
			$zero_count++ if ($data[$i] == 0);
			
		}
		
		if ($zero_count > $#data) {
			die "Logical ERROR, \$zero_count:$zero_count should not be bigger than \$#data $#data.";
		}
		elsif ($zero_count == $#data) {
			print STDERR "All data <= 1e-8: $data{$key}\nIt will be deleted.\n";
			delete $data{$key};
		}
		$zero_count = 0;
		
		next unless defined ($data{$key});
		
		$data[0] =~ s/\t/ /g;
		$data{$key} = join "\t",@data;
	}
	if ($check == 0) {
		&LogErrorFile("The data countains non-numeric info:\n$err_content");
		&ERROR_AND_DIE("Input file error, not a numeric matrix data.");
	}
}


################################################################################################################

sub LogErrorFile {
	my $err_info = shift;
	open ERROUT,">$out.error";
	print ERROUT "$err_info";
	close ERROUT;
}

################################################################################################################

sub OutFileCheck {#检查输出文件路径是否为文件夹，如果是，则添加 _数字 后缀; 并会创建输出目录
	my $out_file = shift;
	my $out_file_check = $out_file;
	while (-d $out_file_check) {
		my $num ++;
		$out_file_check = $out_file."_$num";
	}
	my $out_dir = dirname($out_file_check);
	&MAKE_DIR($out_dir);
	return $out_file_check;
}

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
	-infile		input file, must be separated as a matrix by space or tab
	-outfile	output png file, must end with .png
	-id		file contains id to draw, separated by tab, firse column
	-scale		scale by column or row: column or row or none
	-cluster	cluster by column or row: column or row or both or none
	-color_type	color type: 1-6
	-is_log		log data
	-line		first int line to draw
	-col		column to draw, example: 2,3,4 or 2-3,4,5-7
	
	-div		output results of every cluster
	-anno		annotation directory of genes in "-infile", containing *.GO_tree.stat.xls and *.GO.list.txt
		must choose both of "-div" and "-anno", or none of them

	-sRNA		miRNA 2 targetGene file, *.mir2target.list, for sRNA. 
			If the first column of "-infile" is miRNA ID, you must give this.

	-h		Help

USAGE
	print $usage;
	exit;
}
