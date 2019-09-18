#!/usr/bin/perl -w
my $version="1.0";
my $BEGIN=time();

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);  
use newPerlBase;
my $programe_dir=basename($0);
my $path=dirname($0);

#########################
my ($unigene_prefix,$list,$trans,$step,$cfg,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"=>\$cfg,
				"i:s"=>\$unigene_prefix,
				"l:s"=>\$list,
				"s:n"=>\$step,
				"od=s"=>\$od,
				) or &USAGE;
&USAGE unless ($list and $cfg and $od and $unigene_prefix);

$cfg = &ABSOLUTE_DIR ("$cfg");
&MKDIR ("$od");
$od = &ABSOLUTE_DIR ("$od");
&MKDIR ("$od/work_sh");

my $notename=`hostname`;chomp $notename;
my %sys_cfg=%{&readconf("$Bin/../../../../../config/sys.cfg")};
my $cd_hit_est = $sys_cfg{"cd-hit-est"}; #2014-12-08 ~ 
my $TGICL_HOME = $sys_cfg{"tgicl_home"}; #2014-12-08 ~ 

#####################
my @trans;
if (-f $list && $list !~/fa$/) {
	open IN,"$list" || die $!;
	while (<IN>) {
		chomp;
		s/\s+$//; s/\r+$//;
		next if (/^$/);
		push (@trans, $_);
	}
	close IN;
}
else {
	@trans = split /,/, $list;
}

my %para;
&LOAD_PARA ($cfg,\%para);

$step = $step || 1;

###################### cat fa for cluster
if ($step == 1) {
	print "Step1: cat trans fasta to All trans......\n\n";
	open SH,">$od/work_sh/step1_cat_fasta.sh" || die $!;
	my $fasta;
	foreach (@trans) {
		$fasta .= "$_ ";
	}

	print SH "cat $fasta >$od/All_trans.fasta \n";
	close SH;
	system "sh $od/work_sh/step1_cat_fasta.sh";
	print "Step1: Done!\n\n";
	$step = 2;
}

##################### cd-hit remove redundancy

if ($step == 2) {
	print "Step2: cd-hit removing trascript redundancy......\n\n";
	open SH,">$od/work_sh/step2_cd-hit.sh" || die $!;
	print SH "$cd_hit_est -i $od/All_trans.fasta -c $para{clu_identity} -M 0 -o $od/All_trans.cd-hit_clu.fasta && \n";
	close SH;

	&Shell_qsub ("$od/work_sh/step2_cd-hit.sh","general.q","25G");##2015-11-25   50G
	print "Step2: Done!\n\n";
	$step = 3;
}

##################### tgicl clustering & cap3 assembly

if ($step == 3) {
	print "Step3: tgicl clustering trans & cap3 assembly cluster sequence......\n\n";
	open SH,">$od/work_sh/step3_tgicl.sh" || die $!;
	print SH "cp -f $TGICL_HOME/conf/tgicl.cfg $od &&export LD_LIBRARY_PATH=$Bin\/lib\/:\$LD_LIBRARY_PATH&&";
	print SH "cd $od && $TGICL_HOME/bin/tgicl -F $od/All_trans.cd-hit_clu.fasta -c $para{cpu} -p $para{identity} -l $para{overlap} && "; 
	print SH "rm $od/tgicl.cfg &&\n";
	close SH;

	&Shell_qsub ("$od/work_sh/step3_tgicl.sh","general.q","25G");##2015-11-25  60G
	print "Step3: Done!\n\n";
	$step = 4;
}

################### cat tgicl result

if ($step == 4) {
	print "Step4: cat tgicl assembly result......\n\n";
	open SH,">$od/work_sh/step4_cat_tgicl.sh" || die $!;
	print SH "cat $od/asm_*/contigs >$od/cap3.contigs.fa\n";
	print SH "$TGICL_HOME/bin/cdbyank $od/All_trans.cd-hit_clu.fasta.cidx < $od/All_trans.cd-hit_clu.fasta.singletons >$od/singletons.fa\n";
	print SH "cat $od/cap3.contigs.fa $od/singletons.fa >$od/trans_TGICL_results.fa\n";
	close SH;

	system "sh $od/work_sh/step4_cat_tgicl.sh";
        `perl $Bin/bin/Change_gene_id.pl -i $od/trans_TGICL_results.fa -o $od`;
	print "Step4: Done!\n\n";
	$step = 5;
}

#################### stat final Unigene & Check the Assembly result

if ($step == 5) {
	print "Step5: stat assembly result & exclude abnormal assembly sequence......\n\n";
	open SH,">$od/work_sh/step5_assess_result.sh" || die $!;
	print SH "cd $od && perl $Bin/bin/assembly_stat.pl -i $unigene_prefix -fa trans_TGICL_results_ChangeId.fa -s 300,500,1000,2000 -len $para{max_len} -od Unigene\n";
	close SH;

	system "sh $od/work_sh/step5_assess_result.sh";

	`perl $Bin/bin/Makegff.pl -i $od/Unigene/$unigene_prefix.Unigene.fa -o $od/Unigene/$unigene_prefix.Unigene.gff`;
#	`perl $Bin/bin/Assembly_stat_graph.pl -i $unigene_prefix.Unigene -fa $od/Unigene/$unigene_prefix.Unigene.fa -s 300,500,1000,2000 -od $od/Unigene`;
	`perl $Bin/bin/Getorf_Extract.pl -fa $od/Unigene/$unigene_prefix.Unigene.fa -od $od/Unigene/Unigene_Orf`;
	`perl $Bin/bin/assembly_stat.pl -i $unigene_prefix.orf -fa $od/Unigene/Unigene_Orf/$unigene_prefix.Unigene.cds.fa -s 300,500,1000,2000 -od $od/Unigene/Unigene_Orf/Graph`;
	`perl $Bin/bin/SSR_Analysis.pl  -fa $od/Unigene/$unigene_prefix.Unigene.fa -od $od/Unigene/Unigene_SSR`;
        `perl $Bin/../ppi_network/ppi_network.pl --qseq $od/Unigene/$unigene_prefix.Unigene.fa --odir $od/Unigene/Unigene_PPI`;
	print "cluster2unigene Process is Done!\n\n";
}


##############subs
####################subs
sub LOAD_PARA
{
	my $para_file= shift;
	my $para= shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
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

sub Shell_qsub
{ # &Shell_qsub($sh,$qeue,$vf);
	my $sh = shift;
	my $qeue = shift;
	my $vf = shift;

	if ($notename=~/cluster/)
	{
		`qsub-sge.pl --queue $qeue --resource vf=$vf --reqsub --independent $sh `;
	}
	else
	{
		`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $qeue --resource vf=$vf --reqsub --independent $sh `;
	}
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}
sub USAGE {#
	my $usage=<<"USAGE";

Program: Clustering transcripts or ESTs to Unigene
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:

Usage:
  -cfg     <str>       config file (cluster & assembly parameters)                                      must be given;
  -i       <str>       Index of clustered unigene                                                       must be given;
  -l       <str>       list of files for clustering (comma separated files or list file for dir)        must be given;
  -od      <str>       out files dir                                                                    must be given;
  -s       <int>          step of process
                             1 cat all list trans(ests) into All trans file        default;
                             2 use cd-hit remove redundancy
                             3 use tgicl cluster trans & cap3 assembly
                             4 state final assembly unigenes

USAGE
	print $usage;
	exit;
}
