#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use threads;
use FindBin qw($Bin $Script);
use Config::General;
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
use newPerlBase;

#######################################################################################
my $BEGIN_TIME=time();
my $Title="no_Ref_Trans";
my $version="v2.4.0";
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($data_cfg,$detail_cfg,$step,$step_by_step,$od,$backup_od,$PL,$CSE,$Nt_TrEMBL,$min_kmer_cov,$biocloud,$bioinfo,$parallel,$normalize_cpu,$debug,$queue1,$queue2,$test);
GetOptions(
				"help|h" =>\&USAGE,
				"cfg1:s"=>\$data_cfg,
				"cfg2:s"=>\$detail_cfg,
				"sbs"  =>\$step_by_step,
				"step:s"=>\$step,
				"od:s"=>\$od,
				"bd:s"=>\$backup_od,
				"PL:s"=>\$PL,
				"CSE:s"=>\$CSE,
				"Nt_TrEMBL"=>\$Nt_TrEMBL,
				"ptd_pc:s"=>\$normalize_cpu,
                                "k_cov:i"=>\$min_kmer_cov,
				"bioinfo"=>\$bioinfo,
				"biocloud"=>\$biocloud,
				"parallel"=>\$parallel,
				"debug"=>\$debug,
				"test"=>\$test,
				"queue1:s"=>\$queue1,
				"queue2:s"=>\$queue2,
				) or &USAGE;
&USAGE unless ($data_cfg and $od and $detail_cfg and(($bioinfo and $PL and $CSE and $backup_od) or $biocloud) );

system "mkdir -p $backup_od" if ($bioinfo && !-d $backup_od);
system "mkdir -p $od" unless (-d $od);
$od=&ABSOLUTE_DIR($od);
#$backup_od=&ABSOLUTE_DIR($backup_od);
$data_cfg=&ABSOLUTE_DIR($data_cfg);
$detail_cfg=&ABSOLUTE_DIR($detail_cfg);
my ($Config_dir, $work_sh_dir) = ("$od/Config", "$od/work_sh");
system "mkdir $Config_dir" unless (-d $Config_dir);
system "mkdir $work_sh_dir" unless (-d $work_sh_dir);
my $Q;
$queue1||="general.q";
$queue2||="middle.q";
###
my $Contract_No;
my $Only_One_Sam;
#########################set log ###########
#####################set the step info#########
my %sys_cfg=%{&readconf("$Bin/config/sys.cfg")};
for my $key (keys %sys_cfg){
	next if "PATH headnode"=~/$key/;
	&logAndDie("sys.cfg $key does not exist") unless -e $sys_cfg{$key};
}
$step||='1';
my %step;
&steps_process($step,$step_by_step,\%step);


##################################################Check data_cfg
#&log_current_time("load data config:",\*STDOUT);
&timeLog("load data config:$data_cfg.");
my (%Data,$sampleStat);



if($biocloud){
	my $input_data = &data_cfg_convert($data_cfg,"$Config_dir/data.cfg");
	system "cp $data_cfg $Config_dir/noref_trans.data.cfg";

	print "------Check and Change Phred to 33-----\n";
        system "mkdir -p $od/Data_Assess" unless (-d "$od/Data_Assess"); 
        system "perl $Bin/bin/Data_Assess/bin/Phred_Change/Phred_Change.pl -cfg $Config_dir/data.cfg -od $od/Data_Assess";
my $Sample_name;
open (BIO, "$Config_dir/data.cfg" ) or die $!;
        while (<BIO>) {
                $_=~s/(^s*|\s*$)//;
                next if /^\s*$/;
                if (/^Qphred/) {
                $Q = (split /\s+/,$_)[1];
                   }

                if (/^Sample/) {
                        $Sample_name=(split/\s+/,$_)[1];
                }

                if ($_=~/^fq1/ || $_=~/^fq2/ || $_=~/^raw_fq1/ || $_=~/^raw_fq2/) {
                        my $file=(split/\s+/,$_)[1];
                        if (!-e $file and !$debug) {
                                print "$file is not exist!\n";die;
                        }
                        $Data{$Sample_name}{fq1}=$file if $_=~/^fq1/;
                        $Data{$Sample_name}{fq2}=$file if $_=~/^fq2/;
                }

                if(/^Basecall_stat/){
                        $sampleStat=(split/\s+/,$_)[1];
                        if(!-e $sampleStat and !$debug){
                                die "$sampleStat does not exist!\n";
                        }    
                }
               
   }
        close BIO;
#	%Data=%{$input_data};
}elsif($bioinfo){

        print "------Check and Change Phred to 33-----\n";
        system "mkdir -p $od/Data_Assess" unless (-d "$od/Data_Assess");   
        system "perl $Bin/bin/Data_Assess/bin/Phred_Change/Phred_Change.pl -cfg $data_cfg -od $od/Data_Assess";
       
	system"cp $data_cfg $Config_dir/data.cfg";
	my $Sample_name;
	open (IN,$data_cfg) or die $!;
	while (<IN>) {
		$_=~s/(^s*|\s*$)//;
		next if /^\s*$/;
		
		if (/^Qphred/) {
                $Q = (split /\s+/,$_)[1];
                }
		if (/^Sample/) {
			$Sample_name=(split/\s+/,$_)[1];
		}
		if ($_=~/^fq1/ || $_=~/^fq2/ || $_=~/^raw_fq1/ || $_=~/^raw_fq2/) {
			my $file=(split/\s+/,$_)[1];
			if (!-e $file and !$debug) {
				print "$file is not exist!\n";die;
			}
			$Data{$Sample_name}{fq1}=$file if $_=~/^fq1/;
			$Data{$Sample_name}{fq2}=$file if $_=~/^fq2/;
		}
		if(/^Basecall_stat/){
			$sampleStat=(split/\s+/,$_)[1];
			if(!-e $sampleStat and !$debug){
				die "$sampleStat does not exist!\n";
			}
		}
	}
	close IN;
}
my $sampleNumber=keys %Data;
if($sampleNumber==1){
	$Only_One_Sam=1;
}
else{
	$Only_One_Sam=0;
}
&timeLog("load data config done.");
if ($Only_One_Sam) {
	delete $step{3};
	}
#&log_current_time("load data config done.\n",\*STDOUT);
##################################################detail Configs
#&log_current_time("load parameters config:",\*STDOUT);
&timeLog("load parameters config:$detail_cfg");
my ($key_word,%detail_cfg,%parameters);
if($biocloud){
	my $para_object = Config::General->new($detail_cfg);
    %parameters = $para_object->getall;
	$key_word=$parameters{'projectInfo'}{'Project'};
	for my $type (keys %parameters){
		open (OUT,">$Config_dir/$type.cfg") or die $!;
		for my $label (keys %{$parameters{$type}}){
			if($parameters{$type}{$label}=~/ARRAY/){
				for my $value (@{$parameters{$type}{$label}}){
					print OUT "$label\t$value\n";
				}
			}
			else{
				print OUT "$label\t$parameters{$type}{$label}\n";
			}
		}
		close OUT;
	}
	for my $type (qw(projectInfo geneAnn SNPAnalysis basicAnalysis)){
		if($type eq "projectInfo"){
			open OUT,">>$Config_dir/$type.cfg";
#			print OUT "First_time\txxxx\n";
#			print OUT "Second_time\txxxx\n";
#			print OUT "Third_time\txxxx\n";
#			print OUT "Contract_data\t4G\n";
#			print OUT "Customer_name\txxxx\n";
#			print OUT "Q30\t85%\n";
#			foreach my $key (keys %Data) {
#				print OUT "$key\txxxx\n";
#			}
			print OUT "Rawdata\t$od/Data_Assess\n";
			print OUT "Basic\t$od/Basic_Analysis\n";
                        print OUT "DEG\t$od/DEG_Analysis\n";
#			print OUT "DEG\t$od/DEG_Analysis\n" unless $Only_One_Sam;
			print OUT "ANNO\t$od/Uni_Anno/Result\n";
			print OUT "SNP\t$od/SNP_Analysis\n" unless $Only_One_Sam;
			close OUT;
		}
		elsif($type eq "geneAnn"){
			open OUT,">>$Config_dir/$type.cfg";
			print OUT "mRNA\t$od/Basic_Analysis/Cluster/Unigene/$key_word.Unigene.fa\n\n";
			print OUT "queue $queue1\n";
			close OUT;
		}
		elsif($type eq "SNPAnalysis"){
			open OUT,">>$Config_dir/$type.cfg";
			print OUT "genome\t$od/Basic_Analysis/Cluster/Unigene/$key_word.Unigene.fa\n";
			foreach my $key (keys %Data) {
				print OUT "Sample\t$key\n" .
					"fq1\t$Data{$key}{fq1}\n" .
					"fq2\t$Data{$key}{fq2}\n";
			}
			print OUT "queue $queue2\n";
			close OUT;
		}
		elsif($type eq "basicAnalysis"){
			open OUT,">>$Config_dir/$type.cfg";
			print OUT "para_queue $queue1\n";
			close OUT;
		}
	}
}
elsif($bioinfo){
	$/="\n>>>>>";
	my $Check_num=0;
	open (IN,$detail_cfg) or die $!;
	#system  "iconv -f "UTF-8" -t "UTF-8" $dir_data/template/AllSample_GC_Q.xls -o $dir_data/template/AllSample_GC_Q_final.xls";
	while (<IN>) {
		chomp;
		$Check_num=$.;
		my @lines=split/\n+/;
		$_=~s/\r//g;
		if ($Check_num==1) {
			$_=~/\s*Project\s+(\S+)/;
			$key_word=$1;
			$_=~/\s*Contract_NO\s+(\S+)/;
			$Contract_No=$1;
			for my $line (@lines){
				next if ($line=~/^\#/ or $line=~/^\s*$/);
				my @tmp=split/\s+/,$line;
				$detail_cfg{$tmp[0]}=$tmp[1];
			}
			open (OUT,">$od/Config/projectInfo.cfg") or die $!;
			print OUT $_."\n";
			print OUT "Rawdata\t$od/Data_Assess\n";
			print OUT "Basic\t$od/Basic_Analysis\n";
                        print OUT "DEG\t$od/DEG_Analysis\n"; 
#			print OUT "DEG\t$od/DEG_Analysis\n" unless $Only_One_Sam;
			print OUT "ANNO\t$od/Uni_Anno/Result\n";
			print OUT "SNP\t$od/SNP_Analysis\n" unless  $Only_One_Sam;
			close OUT;
		}
		if ($Check_num==2) {
			open (OUT,">$od/Config/basicAnalysis.cfg") or die $!;
			print OUT $_;
			print OUT "\npara_queue    $queue1\n";
			print OUT "para_normalize_cpu $normalize_cpu\n" if $normalize_cpu;
			close OUT; 
			#&logAndDie("Err:$detail_cfg doesnot contain the parameter para_K_cov\n") if $_!~/para_K_cov\s+\d+/;
		}
		if ($Check_num==3) {
			open (OUT,">$od/Config/geneAnn.cfg") or die $!;
			print OUT "mRNA\t$od/Basic_Analysis/Cluster/Unigene/$key_word.Unigene.fa\n\n";
			print OUT $_;
			print OUT "\nqueue   $queue1\n";
			close OUT;
			if($_!~/kobas\/seq_pep\//){
				die "Error:KEGG database does not supply right\n";
			}
			my %database=("nr"=>1,"nt"=>1,"Kegg"=>1,"Swissprot"=>1,"Pfam"=>1,"Cog"=>1,"Kog"=>1,"TrEMBL"=>1,"eggNOG"=>1);
			my %existDB;
			foreach my $line (@lines) {
				next if ($line=~/^\#/ or $line=~/^\s*$/);
				my @tmp=split/\s+/,$line;
				if(exists $database{$tmp[0]}){
					$database{$tmp[0]}=$tmp[1];
					unless (-e $tmp[1]) {
						die "Error:Database $tmp[0] path is not right\n";
					}
				}
			}
		}
		if ($Check_num==4) {
			for my $line (@lines){
				next if ($line=~/^\#/ or $line=~/^\s*$/);
				my @tmp=split/\s+/,$line;
				if($line=~/^Com|Sep/){
					for my $s (split/,|;/,$tmp[1]){
						die "Error:Samples do not have $s!!!\n" unless exists $Data{$s};
					}
				}
			}
			open (OUT,">$od/Config/DEGAnalysis.cfg") or die $!;
			print OUT $_;
			close OUT;
		}
		if ($Check_num==5) {
			open (OUT,">$od/Config/SNPAnalysis.cfg") or die $!;
			print OUT "genome                     $od/Basic_Analysis/Cluster/Unigene/$key_word.Unigene.fa\n";
			foreach my $key (keys %Data) {
				print OUT "Sample                     $key\n" .
						"fq1                        $Data{$key}{fq1}\n" .
						"fq2                        $Data{$key}{fq2}\n";
			}
			print OUT "\n$_";
			print OUT "\nqueue $queue2\n";
			close OUT;
		}
	}
	close IN;
	$/="\n";
	if ($Check_num!=5) {
		print "Check Your detail_cfg,it should have 4 '>>>>>' \n";die;
	}
	for my $info (qw(Customer_name First_time Second_time Third_time)){
		die "Error:Detail_cfg does not have the $info infomation!!!\n" unless exists $detail_cfg{$info};
	}
	for my $s (keys %Data){
		unless (exists $detail_cfg{$s}){
			die"Error:Sample $s does not have ID with the detail_cfg\n";
		}
	}
}
open B,"$od/Config/basicAnalysis.cfg";
###check Assembly type right
my $assemblyTypeNum=0;
while(<B>){
	chomp;
	next if (/^\#/ or /^\s*$/);
	my @tmp=split/\s+/;
	$assemblyTypeNum++ if (/^(Combine|SamG|Separate)/ && $tmp[1]==1);
}
die "Error:the assembly type(Combine,Separate and SamG) has $assemblyTypeNum method\n" if $assemblyTypeNum!=1;
&timeLog("load parameters config done.");
 #&log_current_time("load parameters config done.\n",\*STDOUT);
 #######################################################################################
 if($bioinfo){&createLog($Title,$version,$$,"$od/",$test)};
################################################## MAIN

open (SS,">>$od/../bioinfor_pipeline.log");
&log_current_time("This project will be analysed in ".scalar(keys %step)." steps.",\*SS);
###################RNA-Seq_Data_Assess
if (exists $step{1} &&(!-e "$od/work_sh/step1.$step{1}.sh.finish")) {
	open (SH1,">$od/work_sh/step1.$step{1}.sh") or die $!;
	mkdir "$od/Data_Assess" unless -d "$od/Data_Assess";
	mkdir "$od/Data_Assess/PNG" unless -d "$od/Data_Assess/PNG";
	print SH1 "perl $Bin/bin/Data_Assess/rna_seq_data_assess.pl -Q $Q -config $Config_dir/data.cfg -byraw 0 -outdir $od/Data_Assess -queue $queue1"; #cmd
	print SH1 "&&cp $sampleStat $od/Data_Assess/AllSample.data.stat &&$sys_cfg{Rscript} $Bin/bin/Data_Assess/pie_plot.R --infile $sampleStat --outfile $od/Data_Assess/PNG " if (defined $sampleStat);  ##2015-7-3 add pie png
    print SH1 "&&echo \"[`date '+%Y-%m-%d %T'`] $Script Data Assess Step is done! \" " if $bioinfo;
    print SH1 "&&perl -e 'print \"Dear $PL,\\n\\t$Contract_No Data Assessment is done. The result is shown below:\\n\"; system \"cat $od/Data_Assess/AllSample_GC_Q.stat\"; print \"\\nResult Directory: $od/Data_Assess/\";' |mail -s \"$Contract_No Data Assessment Result\" $PL\@biomarker.com.cn \n" if $bioinfo;   #mail result to PL
	close SH1;
	&process_cmd("$od/work_sh/step1.$step{1}.sh");
}
######################for data sieze for assembly
open B,"$od/Config/basicAnalysis.cfg" or die "$!";
my %tmpHash;
while(<B>){
	chomp;
	next if (/^\#/ or /^\s*$/);
	my @tmp=split/\s+/;
	$tmpHash{$tmp[0]}=$tmp[1];
}
close B;
my %dataSize;
if(!$debug){
	open IN,"$od/Data_Assess/AllSample_GC_Q.stat" or die "$!";
	
	while(<IN>){
		chomp;
		next if (/^\#/ or /^\s*$/);
		my @tmp=split/\s+/;
		$dataSize{$tmp[0]}=$tmp[2];
	}
	close IN;
}
else{
	for my $sample (keys %Data){
		my $dataSize=(-s $Data{$sample}{fq1})+(-s $Data{$sample}{fq2});
		$dataSize{$sample}=$dataSize*0.4;
	}
}
####only one sample ,set assembly type to be combine
if($Only_One_Sam==1){
	$tmpHash{Combine}=1;
	$tmpHash{Separate}=0;
	$tmpHash{SamG}=0;
}
if($tmpHash{Combine}==1){
	my ($totalSize);
	for my $sample (keys %dataSize){
		$totalSize+=$dataSize{$sample};
	}
	$totalSize=int($totalSize/1000000000)."G";
	$tmpHash{Com_Data}="$totalSize";
}
elsif($tmpHash{Separate}==1){
	my ($totalSize);
	for my $sample (keys %dataSize){
		$totalSize+=$dataSize{$sample};
	}
	my $sampleNum=keys %dataSize;
	$tmpHash{Sep_Data}=int($totalSize/$sampleNum/1000000000)."G";
}
elsif($tmpHash{SamG}==1){
	my $groupNum=0;
	for my $group (keys %tmpHash){
		if($group=~/^(Group\d+)$/){
			$groupNum++;
			my $groupPrefix=$1;
			my ($totalSize);
			for my $sample (split/,/,$tmpHash{$groupPrefix}){
				die "Error:sample:$sample is not in the $od/Data_Assess/AllSample_GC_Q.stat file\n" unless (exists $dataSize{$sample});
				$totalSize+=$dataSize{$sample};
			}
			$tmpHash{$groupPrefix."_Data"}=int($totalSize/1000000000)."G";
		}
	}
	die "Error:samG assembly type has $groupNum groups which is illegal\n" if $groupNum<2;
}

open B,">$od/Config/basicAnalysis.cfg" or die "$!";
for  my $key (keys %tmpHash){
	print B "$key\t$tmpHash{$key}\n";
}
close B;
###################Basic_Analysis
if (exists $step{2} && (!-e "$od/work_sh/step2.$step{2}.sh.finish")) {
	open (SH2,">$od/work_sh/step2.$step{2}.sh") or die $!;
	mkdir "$od/Basic_Analysis" unless -d "$od/Basic_Analysis";
        if(defined $min_kmer_cov){
        print SH2 "perl $Bin/bin/Basic_Analysis/v1.4/Basic_Analysis.pl -basic_config $od/Config/basicAnalysis.cfg -data_config $Config_dir/data.cfg -i $key_word -qphred $Q -od $od/Basic_Analysis -min_kmer_cov $min_kmer_cov ";
        }else{
	print SH2 "perl $Bin/bin/Basic_Analysis/v1.4/Basic_Analysis.pl -basic_config $od/Config/basicAnalysis.cfg -data_config $Config_dir/data.cfg -i $key_word -qphred $Q -od $od/Basic_Analysis "; }
    print SH2 "&&echo \"[`date '+%Y-%m-%d %T'`] $Script Basic Analysis Step is done! \" " if $bioinfo;
    print SH2 "&&perl -e 'print \"Dear $PL,\\n\\t$Contract_No Basic Analysis is done. The result is shown below:\\n\"; system \"more $od/Basic_Analysis/Cluster/*/$key_word.Unigene.stat.xls $od/Basic_Analysis/Remap/geneExpression/*.Mapped.stat.xls\"; print \"\\nResult Directory: $od/Basic_Analysis/ \";' |mail -s \"$Contract_No Basic Analysis Result\" $PL\@biomarker.com.cn \n" if $bioinfo;   #mail result to PL
	close SH2;
	&process_cmd("$od/work_sh/step2.$step{2}.sh");
}

while(exists $step{2}){
	last if $debug;
	if(-e "$od/Basic_Analysis/Cluster/Unigene/$key_word.Unigene.fa"){
		last;
	}
	else{
		writeLog("Warning:can not find the trinity results,wait!\n");
		sleep(3600);
	}
}
###################SNP_Analysis

if (exists $step{3} &&(!-e "$od/work_sh/step3.$step{3}.sh.finish")) {
	open (SH,">$od/work_sh/step3.$step{3}.sh") or die $!;
	mkdir "$od/SNP_Analysis" unless -d "$od/SNP_Analysis";
	print SH "perl $Bin/bin/SNP_Analysis/v2.0/SNP_Trans_main_noRef.pl -cfg $od/Config/SNPAnalysis.cfg -od $od/SNP_Analysis -qphred $Q";  #cmd
    print SH "&&echo \"[`date '+%Y-%m-%d %T'`] $Script SNP Analysis Step is done! \"  " if $bioinfo;
    print SH "&&perl -e 'print \"Dear $PL,\\n\\t$Contract_No SNP Analysis is done. The result statistic is shown below:\\n\"; system \"cat $od/SNP_Analysis/SNP/stat/AllSample.snp.stat \"; print \"\\nResult Directory: $od/SNP_Analysis/ \";' |mail -s \"$Contract_No SNP Analysis Result\" $PL\@biomarker.com.cn \n" if $bioinfo;   #mail result to PL
	close SH;
}
##############Unigene_Function_Annotation
my $unigene_file = "$od/Basic_Analysis/Cluster/Unigene/$key_word.Unigene.fa";
my $geneExpression_dir = "$od/Basic_Analysis/Remap/geneExpression";
my $anno_result_file = "$od/Uni_Anno/Result/All_Database_annotation.xls";

if (exists $step{4} && (!-e "$od/work_sh/step4.$step{4}.sh.finish")) {
	open (SH,">$od/work_sh/step4.$step{4}.sh") or die $!;
	mkdir "$od/Uni_Anno" unless -d "$od/Uni_Anno";

	if (defined $Nt_TrEMBL) {
        print SH "perl $Bin/bin/Gene_Anno/v1.5/Gene_Func_Anno_Pipline.pl --query $unigene_file --all --cfg $od/Config/geneAnn.cfg --od $od/Uni_Anno "; #cmd
        print SH "&&cp -r $od/Uni_Anno/Unigene_CDS_Predict/ $od/Basic_Analysis/Cluster/Unigene/Unigene_CDS "; #cmd
    } else {
        print SH "perl $Bin/bin/Gene_Anno/v1.5/Gene_Func_Anno_Pipline.pl --query $unigene_file --nr --swissprot --cog --kog --kegg --pfam --GO --eggNOG --cfg $od/Config/geneAnn.cfg --od $od/Uni_Anno "; #cmd
        print SH "&&cp -r $od/Uni_Anno/Unigene_CDS_Predict/ $od/Basic_Analysis/Cluster/Unigene/Unigene_CDS"; #cmd
    }

    print SH "&&echo \"[`date '+%Y-%m-%d %T'`] $Script Gene Function Annotation Step is done! \"  " if $bioinfo;
    print SH "&&cat $od/Uni_Anno/Result/Function_Annotation.stat.xls ";
    print SH "&&perl -e 'print \"Dear $PL,\\n\\t$Contract_No Gene Function Annotation is done. The result statistic is shown below:\\n\"; system \"cat $od/Uni_Anno/Result/Function_Annotation.stat.xls \"; print \"\\nResult Directory: $od/Uni_Anno/Result/ \";' |mail -s \"$Contract_No Gene Function Annotation Result\" $PL\@biomarker.com.cn \n" if $bioinfo;   #mail result to PL
	close SH;
}

if(($bioinfo or ($biocloud &&$parallel)) && exists $step{3} &&exists $step{4}){
	&process_cmds_parallel("$od/work_sh/step3.$step{3}.sh", "$od/work_sh/step4.$step{4}.sh");
}
else{
	&process_cmd("$od/work_sh/step3.$step{3}.sh") if exists $step{3};
    &process_cmd("$od/work_sh/step4.$step{4}.sh") if exists $step{4};
}



###################DEG_Analysis
if (exists $step{5} && (!-e "$od/work_sh/step5.$step{5}.sh.finish")) {
	open (SH,">$od/work_sh/step5.$step{5}.sh") or die $!;
	mkdir "$od/DEG_Analysis" unless -d "$od/DEG_Analysis";
        print SH "perl $Bin/bin/DEG_Analysis/v3.4/ExtractExpression.pl -indir $od/Basic_Analysis/Remap/geneExpression -out $od/DEG_Analysis ";
        unless($Only_One_Sam){
	print SH "&& perl $Bin/bin/DEG_Analysis/v3.4/fpkm_plot.pl -i $od/DEG_Analysis/All_gene_counts.list -o $od/DEG_Analysis/density ";
        print SH "&& perl $Bin/bin/DEG_Analysis/v3.4/DEG_Analysis.pl -i $od/DEG_Analysis/All_gene_counts.list -cfg $od/Config/DEGAnalysis.cfg  -od $od/DEG_Analysis ";
        print SH "&& perl $Bin/bin/DEG_Analysis/v3.4/DEG_enrich.pl -enrichment $od/Uni_Anno/Result -od $od/DEG_Analysis";
        print SH "&& perl $Bin/bin/DEG_Analysis/v3.4/Protein_to_protein_deg.pl -id $od/DEG_Analysis -od $od/DEG_Analysis/DEG_PPI/  -ppi $od/Basic_Analysis/Cluster/Unigene/Unigene_PPI/PPI.txt ";
    print SH "&&echo \"[`date '+%Y-%m-%d %T'`] $Script DEG Analysis Step is done! \"" if $bioinfo;
    print SH "&&perl -e 'print \"Dear $PL,\\n\\t$Contract_No DEG Analysis is done. The result statistic is shown below:\"; system \"perl $Bin/bin/Tool/stat_deg.pl -id $od/DEG_Analysis \"; print \"Result Directory: $od/DEG_Analysis/ \";' |mail -s \"$Contract_No DEG Analysis Result\" $PL\@biomarker.com.cn \n" if $bioinfo;      
	                 }
        close SH;
	&process_cmd("$od/work_sh/step5.$step{5}.sh");
}

#####################6.Result_Extract_xml
if(exists $step{6}){
	open (SH,">$od/work_sh/step6.$step{6}.sh") or die $!;
		mkdir "$od/Web_Report" unless -d "$od/Web_Report";
		mkdir "$od/BMK_Results" unless -d "$od/BMK_Results";
                mkdir "$od/Needed_Data" unless -d "$od/Needed_Data";
                mkdir "$od/Needed_Data/Config" unless -d "$od/Needed_Data/Config";
		print SH "perl $Bin/bin/Result_Extract/biocloud/Result_extract_biocloud.pl -cfg $od/Config/projectInfo.cfg -od $od/Web_Report\n";
		if($biocloud){print SH "cp $data_cfg $od/Needed_Data/Config/data.cfg && cp $detail_cfg $od/Needed_Data/Config/detail.cfg \n";}elsif($bioinfo){
		print SH "perl $Bin/bin/Biocloud_Report/No_ref_bioCloud_cfg.pl --id $od/Needed_Data/Config --cfg1 $data_cfg --cfg2 $detail_cfg\n";
		}
                print SH "perl $Bin/bin/Result_Extract/local/Result_extract_Needed.pl -cfg $od/Config/projectInfo.cfg -od $od/Needed_Data ";
                if($Only_One_Sam){print SH "--oneSam\n";}else{print SH "\n";}
		print SH "perl $Bin/bin/Result_Extract/local/Result_extract_local.pl -cfg $od/Config/projectInfo.cfg -od $od/BMK_Results\n";
		print SH "perl $Bin/bin/Biocloud_Report/get_gene_id.pl -cfg $od/Config/projectInfo.cfg -inputdir $od/Web_Report -out $od/Web_Report";
		if ($Only_One_Sam){
			print SH " --oneSam\n"; 
		}
		else{print SH "\n"}
		print SH "perl $Bin/bin/Html_Report/local_xml.pl -cfg $od/Config/projectInfo.cfg -id $od/BMK_Results -pp";
		if ($Only_One_Sam){
			print SH " --Only_One_Sam\n";
		}
		else{print SH "\n"}
                print SH "cd $od/BMK_Results\n";
                print SH "$sys_cfg{python} $Bin/bin/Html_Report/bin/xml2HtmlConverter.py -i configtest_raw.xml\n";
		print SH "cp -r $od/BMK_Results/* $od/Web_Report\n";
                print SH "cp $od/Web_Report/no_ref_trans_all_table.xls $od/BMK_Results\n";
		print SH "cd $od/Web_Report/HTML/ &&for i in *html;do cp \$i \$i.cloud ;done&&sed -i 's/\.\.\\\/.*logo.jpg/logo_image_path/' *.cloud\n";
                unless($Only_One_Sam){
                print SH "cd $od/Web_Report/HTML/ &&for i in *DEG.html.cloud;do sed -i 's/<a href=\\\"\\.\\.\\\/\\(.*\\\)\\\" >.*<\\\/a>\/ \\1 \/' \$i;done\n";
                print SH "cd $od/Web_Report/HTML/ &&for i in *DEG_KEGG.html.cloud;do sed -i 's/<a href=\\\"\\.\\.\\\//<a href=\\\"project_Template_Path\\\/\\.\\.\\\//' \$i;done\n";
                print SH "sed -i 's/<a href=\\\"\\.\\.\\\/\\(.*\\\)\\\" >.*<\\\/a>\/ \\1 \/' SNP.html.cloud\n";
                for my $file (glob("$od/Web_Report/HTML/*DEG_KEGG.html.cloud")){
                     print SH "sed -i 's/<a href=\"\.\.\//<a href=\"project_Template_Path\/\.\.\//' $file\n";
                         }
                }
                 for my $file (qw(expression.html.cloud SSR.html.cloud)){
                         print SH "sed -i 's/<a href=\\\"\\.\\.\\\/\\(.*\\\)\\\" >.*<\\\/a>\/ \\1 \/' $file\n";
                 }
		print SH "rm -r $od/BMK_Results/configtest_raw.xml $od/BMK_Results/abstract.txt $od/BMK_Results/template\n";
		print SH "perl $Bin/bin/package/package.pl -in $od/BMK_Results/\n";
 	    if($bioinfo){
			print SH "perl -e 'print \"Dear $PL,\\n\\t$Contract_No Final Analysis Result and Web Report is prepared. Please check!\\n\\nResult Directory: $od/Web_Report/ \";' |mail -s \"$Contract_No Final Analysis Result and Web Report\" $PL\@biomarker.com.cn \n";  #mail result to PL
		}

	close SH;
	process_cmd("$od/work_sh/step6.$step{6}.sh");
}

###################Word Report
if (exists $step{word_report}) {
	unless (-d "$od/Web_Report"){print "You Must Have $od/Web_Report!\n";die;}
	if ($Only_One_Sam) {
        print SH4 "perl $Bin/bin/Word_Report/v1.5.2/bin/Trans_No_Ref_Report_One_Sam.pl --id $od/Web_Report --pf $Contract_No --od $od/ && ";
    } else {
        print SH4 "perl $Bin/bin/Word_Report/v1.5.2/bin/Trans_No_Ref_Report2.pl --id $od/Web_Report --pf $Contract_No --od $od/ && ";
    }

    print SH4 "echo \"[`date '+%Y-%m-%d %T'`] $Script Analysis is completely done! \" && ";
#   print SH4 "ssh cluster echo Your Project Is Completely Done,Please Check!GOOD LUCK!|mail -s S8_Word_Report_is_Done_$od $PL\@biomarker.com.cn \n";
    print SH4 "perl -e 'my \@stat=stat(\"$od/${Contract_No}_Report.rtf.zip\"); if (\$stat[7]<5242880) {print \"Dear $CSE and $PL,\\n\\t$Contract_No Data Analysis is Totally Finished.\\n\\tFinal Report is attached as ${Contract_No}_Report.rtf.zip. Please unzip it and open rtf file with Microsoft Office, and save it as word format. Appropriate adjustment is needed.\\n\"; system \"uuencode $od/${Contract_No}_Report.rtf.zip ${Contract_No}_Report.rtf.zip\";}else{print \"Dear $CSE and $PL,\\n\\t$Contract_No Data Analysis is Totally Finished.\\n\\tFinal Report is $od/${Contract_No}_Report.rtf. Please unzip it and open rtf file with Microsoft Office, and save it as word format. Appropriate adjustment is needed.\\n\";} ' |mail -s \"$Contract_No Data Analysis Final Report\" $CSE\@biomarker.com.cn -c $PL\@biomarker.com.cn \n";  #mail result to CSE and PL
}

##################QC report
if(exists $step{7} &&(!-e "$od/work_sh/step7.$step{7}.sh.finish") &&$bioinfo){
	open (SH,">$od/work_sh/step7.$step{7}.sh") or die $!;
	if (exists $step{7} && defined $Contract_No) {
		mkdir "$od/QC_Report" unless -d "$od/QC_Report";

		# get assembly groups
		my ($SampleGroup,$SamG,$COM_limit) = ("Samples Group\n",0,0);
		open (BASIC,"<$od/Config/basicAnalysis.cfg") or die $!;
		while (<BASIC>) {
			chomp;
			next if (/^\s+$/ || /^#/);
			$SamG = $1 if (/^SamG\s+(\d)/);
			$SampleGroup.= "$_\n" if (/^Group\d\s+/);
			$COM_limit=$1 if (/^\s*Combine\s+(\d+)/);
		}
		close BASIC;

		#make config
		open (CFG,">$od/Config/QC.cfg") or die $!;
			print CFG "Project ID\n$Contract_No\nLibrary path\n/share/nas1/litt/Config_dir/\nDataAssess path\n$od/Data_Assess\nAnalysis path\n$od/Basic_Analysis\n";
			print CFG "DEG_Analysis path\n$od/DEG_Analysis\n" unless $Only_One_Sam;
			print CFG "$SampleGroup" unless ($SamG == 0);
			close CFG;

			if ($COM_limit==1 && $SamG==0) {
			print SH "perl $Bin/bin/QC_Report/no_ref_trans_QC_combine_v3.pl -cfg $od/Config/QC.cfg -o $od/QC_Report "; #cmd
			print SH "&&echo \"[`date '+%Y-%m-%d %T'`] $Script Analysis Quality Contral Step is done! \"";
			print SH "&&perl -e 'print \"Dear $PL and $CSE,\\n\\t$Contract_No Analysis Quality Contral is done. The result is shown in $Contract_No.QC.stat attached.\\n\\n\"; print \"Analysis Directory: $od/\\n\"; system \"uuencode $od/QC_Report/$Contract_No.QC.stat $Contract_No.QC.stat\"; ' |mail -s \"$Contract_No Quality Contral Report\" $PL\@biomarker.com.cn -c $CSE\@biomarker.com.cn \n";   #mail result to PL and CSE
		}
		elsif ($COM_limit==0 && $SamG==0) {
			print SH "perl $Bin/bin/QC_Report/no_ref_trans_QC_separate_v3.pl -cfg $od/Config/QC.cfg -o $od/QC_Report && "; #cmd
			print SH "echo \"[`date '+%Y-%m-%d %T'`] $Script Analysis Quality Contral Step is done! \" && ";
			print SH "perl -e 'print \"Dear $PL and $CSE,\\n\\t$Contract_No Analysis Quality Contral is done. The result is shown in $Contract_No.QC.stat attached.\\n\\n\"; print \"Analysis Directory: $od/\\n\"; system \"uuencode $od/QC_Report/$Contract_No.QC.stat $Contract_No.QC.stat\"; ' |mail -s \"$Contract_No Quality Contral Report\" $PL\@biomarker.com.cn -c $CSE\@biomarker.com.cn \n";   #mail result to PL and CSE
		}
		elsif ($COM_limit==0 && $SamG!=0) {
			print SH "perl $Bin/bin/QC_Report/no_ref_trans_QC_group_v3.pl -cfg $od/Config/QC.cfg -o $od/QC_Report && "; #cmd
			print SH "echo \"[`date '+%Y-%m-%d %T'`] $Script Analysis Quality Contral Step is done! \" && ";
			print SH "perl -e 'print \"Dear $PL and $CSE,\\n\\t$Contract_No Analysis Quality Contral is done. The result is shown in $Contract_No.QC.stat attached.\\n\\n\"; print \"Analysis Directory: $od/\\n\"; system \"uuencode $od/QC_Report/$Contract_No.QC.stat $Contract_No.QC.stat\"; ' |mail -s \"$Contract_No Quality Contral Report\" $PL\@biomarker.com.cn -c $CSE\@biomarker.com.cn \n";   #mail result to PL and CSE
		}

	}
	close SH;
	process_cmd("$od/work_sh/step7.$step{7}.sh");
}
###################8.back up
if(exists $step{8} &&$bioinfo){
	open (SH,">$od/work_sh/step8.$step{8}.sh") or die $!;
	print SH "perl $Bin/bin/Backup/v1.0/no_ref_trans_backup.pl --ai $Contract_No --id $od --od $backup_od --bm --ck ";
	print SH "&&echo \"[`date '+%Y-%m-%d %T'`] $Script Backup Step is done! \" ";
	print SH "&&perl -e 'print \"Dear $CSE and $PL,\\n\\t$Contract_No Rawdata and Analysis Result for Client is prepared. \\n\\nBackup Directory: $backup_od/\";' |mail -s \"$Contract_No Rawdata and Analysis Result for Client\" $CSE\@biomarker.com.cn -c $PL\@biomarker.com.cn \n";  #mail result to CSE and PL
	close SH;
	process_cmd("$od/work_sh/step8.$step{8}.sh");
}
&log_current_time("Analysis of this project is completed.",\*SS);
close SS;
#######################################################################################
&totalTime();
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
		die;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################
sub data_cfg_convert{
    my ($input_cfg, $output_cfg) = @_;
    my $data_object = Config::General->new($input_cfg);
    my %data = $data_object->getall;
    print "Sample Number: ".(scalar keys %data)."\n";

    open (OUT, ">$output_cfg") or die $!;
    for my $sample (sort keys %data) {
#       die "ERROR: The data of $sample is invalid.\n" unless (exists $data{$sample}{fq1} and exists $data{$sample}{fq2} and -e $data{$sample}{fq1} and -e $data{$sample}{fq2});
        die "ERROR: The data of $sample is invalid.\n" unless (exists $data{$sample}{fq1} and exists $data{$sample}{fq2}); #tmp
        print OUT "Sample\t$sample\n";
        print OUT "fq1\t$data{$sample}{fq1}\n";
        print OUT "fq2\t$data{$sample}{fq2}\n";
        print OUT "raw_fq1\t$data{$sample}{raw_fq1}\n" if (exists $data{$sample}{raw_fq1});
        print OUT "raw_fq2\t$data{$sample}{raw_fq2}\n" if (exists $data{$sample}{raw_fq2});
        print OUT "";
        print "${sample}_fq1: $data{$sample}{fq1}\n${sample}_fq2: $data{$sample}{fq2}\n";
    }
    close OUT;
    return \%data;
}
##################
sub steps_process {
    &timeLog("step check:");
    my ($step_str, $step_by_step, $step_hash) = @_;
    my %_step_ = (
        '1' => 'RNA-Seq_Data_Assess',
        '2' => 'Basic_Analysis',
        '3' => 'SNP_Analysis',
        '4' => 'Unigene_Function_Annotation',
        '5' => 'DEG_Analysis',
        '6' => 'Result_Extract_xml',
		'7'=> 'QC_report',
		'8' => 'backup',
    );

    if ($step_by_step) {
        &writeLog("step-by-step: ON");
        my $total_step=scalar (keys %_step_);
        if ($step_str <=$total_step) {
            for my $s ($step_str..$total_step) {
                $step_hash->{$s} = $_step_{$s};
            }
        } else {
            print "ERROR: illegal steps specified by --step, or options specified by --step and --sbs clashed. \n";
            die;
        }
    }
	else {
		for my $s (split /,/,$step_str) {
            if (exists $_step_{$s}) {
                $step_hash->{$s} = $_step_{$s};
            } else {
                print "ERROR: illegal steps specified by --step.\n";
                die;
            }
        }
    }

#    print "steps_to_run: ".(join ", ",(map {sprintf "$_.$step_hash->{$_}"} sort keys %{$step_hash})).".\n";
    &timeLog("step check done.");
}
#######################

##################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
################################################################################################################
sub process_cmd {
    my ($cmd) = @_;
    my ($stepNum)=$cmd=~/work_sh\/step(\d+)/;
	my $start_time = time();
	log_current_time("step_$stepNum: $step{$stepNum} start.",\*SS);
	#log_current_time("step_$stepNum: $step{$stepNum} start.",\*STDOUT);
	&stepStart($stepNum,"$step{$stepNum}");
	unless($debug){
		unless(-e "$cmd.finish"){
			 my $ret = system"sh $cmd >$cmd.log 2>&1";
			if ($ret) {
				 die "Error, cmd: $cmd died with ret $ret";
			}
			system"touch $cmd.finish";
		}
	}
	log_current_time("step_$stepNum: $step{$stepNum} finished.",\*SS);
	#log_current_time("step_$stepNum: $step{$stepNum} finished.",\*STDOUT);
	&stepTime($stepNum,"$step{$stepNum}");
}
#####################
sub log_current_time {
     # get parameter
     my ($info,$fh) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print $fh "$curr_time:\t$info\n";

}
################################################################################################################
sub process_cmds_parallel {
    my @cmds = @_;
    my %threads;
    foreach my $cmd (@cmds) {
        my $thread = threads->create(\&process_cmd,$cmd);
		my $tid=$thread->tid();
		$threads{$tid}=basename($cmd);
		sleep(2) if $debug;
    }
    foreach my $thread (threads->list(threads::all)) {
        my $tid=$thread->tid();
		$thread->join();
		#&log_current_time("The thread $tid ($threads{$tid}) has joined.",\*STDOUT);
		&timeLog("The thread $tid ($threads{$tid}) has joined.");
        if (my $error = $thread->error()) {
            die"Error, $threads{$tid} exited with error $error\n";
        }
    }
}

################################################################################################################
sub USAGE {#
	my $usage=<<"USAGE";
---------------------------------------------------------------------------------------------------------------------------------------
  Program: $Script
  Version: $version
  Contact: Simon Young <yangxh\@biomarker.com.cn>
     Date: 2014-10-21

    Usage:
           -cfg1    <file>  data config, forced
           -cfg2   <file>  detail config, forced
           -od           <dir>   analysis output dir, forced
           -bd           <dir>   backup output dir, forced only bioinfo
           -sbs                  run  by step
           -step         <str>   steps to run or to start
                         1         RNA-Seq_Data_Assess, default to start
                         2         Basic_Analysis
                         3         SNP_Analysis
                         4         Unigene_Function_Annotation
                         5         DEG_Analysis
                         6         Result_Extract_xml
                         7         QC report (only for bioinfo)
                         8         backup    (only for bioinfo)
           -ptd_pc     <str>   the cpu name to run trinity normalize step or combine assembly
           -PL           <str>   abbr. of Project Leader\'s name, forced only bioinfo
           -CSE          <str>   abbr. of Customer Service Executive\'s name, forced only bioinfo
           -Nt_TrEMBL            Do Nt and TrEMBL Anno in Step 4
           -k_cov     <num>      min coverage of k-mer
           -bioinfo              Run the script for bioinfor
           -biocloud             Run the script for biocloud
           -parallel             Some Script can run parallelly(only for biocloud)
           -queue1    <str>      the queue for data_assess,quantifygraph,buterfly,geneAnno,
           -queue2    <str>      the queue for snp_analysis
           -debug                only print the run sh file
           -test                 do not create sys log file to stat the script run station
           -h                    Help

  Example:
      bioinfo:  perl $Script -cfg1 data.cfg -cfg2 detail.cfg -od Analysis/ -bd Backup/ -sbs -PL xugl -CSE xugl -bioinfo -queue1 general.q -queue2 middle.q
      bioinfo:  perl $Script  -cfg1 data.cfg -cfg2 detail.cfg -od Analysis/ -bd Backup/ -step 1,2,3 -PL xugl -CSE xugl -k_cov 2 -bioinfo -queue1 general.q -queue2 middle.q
      biocloud: perl $Script -cfg1 data.cfg -cfg2 detail.cfg -od Analysis -sbs  -biocloud -queue1 general.q -queue2 middl.q (-parallell )

---------------------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
