#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.1.1";
my @Original_ARGV=@ARGV;

# ==============================================================
# Get Options
# ==============================================================
my ($aligndir,$genomeFile,$dOut,$step,$window,$cluster,$QD,$FS,$doINDEL,$auto,$ploidy,$doRecal,$cpu,$vf,$queue,$qphred);

GetOptions(
				"help|?" =>\&USAGE,
				"aligndir:s"=>\$aligndir,
				"ploidy:s"=>\$ploidy,
				"ref:s"=>\$genomeFile,
				"win:s"=>\$window,
				"clu:s"=>\$cluster,
				"QD:s"=>\$QD,
				"FS:s"=>\$FS,
				"od:s"=>\$dOut,
				"step:s"=>\$step,
				"doIndel:s"=>\$doINDEL,
				"doRecal:s"=>\$doRecal,
				"queue:s"=>\$queue,
				"vf:s"=>\$vf,
				"cpu:s"=>\$cpu,
				"auto:s"=>\$auto,
                                "qphred"=>\$qphred,
				) or &USAGE;
&USAGE unless ($aligndir and $genomeFile and $queue);

&timeLog("Starting SNP Calling");

#===============================================================
# Default optional value
#===============================================================
###### software dir ###########
my %sys_cfg=%{&readconf("$Bin/../../../../config/sys.cfg")};
my $picardpath    =  $sys_cfg{picard};
#my $picardpath	=	"/share/nas2/genome/biosoft/picard-tools/1.94";   # 2015-05-27
my $_GroupRead_	=	"$picardpath/AddOrReplaceReadGroups.jar";
my $_MarkDup_	=	"$picardpath/MarkDuplicates.jar";
my $_CreateDict_	=	"$picardpath/CreateSequenceDictionary.jar";
my $_samtools_    =   $sys_cfg{samtools};
#my $_samtools_	=	"/share/nas2/genome/biosoft/samtools/0.1.18/samtools";   # 2015-05-27
my $_GATK_    =   $sys_cfg{gatk};
#my $_GATK_	=	"/share/nas2/genome/biosoft/GenomeAnalysisTK/3.2-2/GenomeAnalysisTK.jar";   # 2015-05-27
my $_tools_	=	"$Bin/bin";
#################### resource request ################
my %sample;
my %para;
my $chr_num;
my $notename = `hostname`; chomp $notename;
my $Thread		=	5;
   #$queue||=	"general.q";
   $vf	||=	"20G";
   $cpu	||=	"50";
   $doRecal	||=	0;
   $step	||=	1;
   $doINDEL	||=	0;
   $auto	||=	1;
   $ploidy	||=	2;
#############  default filter threshold #################
$window	||=	35;
$cluster||=	3;
$FS		||=	30.0;
$QD		||=	2.0;
my $window_strict	=	20;
my $cluster_strict	=	2;
my $FS_strict		=	20.0;
my $QD_strict		=	5.0;
#############  default file path ################
$dOut||=abs_path."/SNP_Result";&MKDIR($dOut);
$dOut=abs_path($dOut);
$genomeFile=abs_path($genomeFile);
my $tmppath		=	"$dOut/tmp";
my $GenomeDir	=	"$dOut/Genome";
my $_java_="$sys_cfg{java} -XX:ParallelGCThreads=$Thread -Xmx30g -Djava.io.tmpdir=$tmppath";
my $_vcf2snplist_="$_tools_/vcf_to_snplist_v1.5.1.pl";
&MKDIR($tmppath);
&MKDIR($GenomeDir);
#if (-e "$dOut/work_sh") {
#	die("work shell exists!");
#}
#system("rm -r $dOut/work_sh/") if (-d "$dOut/work_sh");
&MKDIR("$dOut/work_sh") unless(-d "$dOut/work_sh");
my $resultdir="$dOut/Result";
&MKDIR($resultdir);
my $java_tmp;
#===============================================================
# Process
#===============================================================
my $newgenomeFile="$GenomeDir/".basename($genomeFile);
$aligndir=abs_path($aligndir);
my @x=glob("$aligndir/*.bam");
die "no alignment bam file!!!" if(@x<1);
`ln -s $genomeFile $newgenomeFile` unless (-e $newgenomeFile) ;
`ln -s $aligndir/*.bam $resultdir/`;
$genomeFile=$newgenomeFile;
$resultdir=abs_path($resultdir);

my @inBam=glob("$aligndir/*.bam");
for(my $i=0;$i<@inBam;$i++){
    $inBam[$i]="$resultdir/".basename($inBam[$i]);
}

my $stepNum=1;
############################################# create index and dictionary file ############################
open (SH0,">$dOut/work_sh/s$stepNum.CreatIndex.sh") or die $!;
if (!-f "$genomeFile.fai"){
	my $cmd = "$_samtools_ faidx $genomeFile";
	print SH0 $cmd."\n";
}
(my $dict_file = $genomeFile) =~ s/\.fa$|\.fasta$/.dict/;
if (!-f $dict_file){
	my $cmd = "$_java_ -jar $_CreateDict_ R=$genomeFile O=$dict_file";
	print SH0 $cmd."\n";
}
close(SH0);
$stepNum++;
#########################################  sort and mark duplicate #################################
#system "mv $dOut/work_sh/GroupRead.sh $dOut/work_sh/GroupRead.sh.bak" if(-e "$dOut/work_sh/*GroupRead.sh");
#system "mv $dOut/work_sh/MarkDup.sh $dOut/work_sh/MarkDup.sh.bak" if(-e "$dOut/work_sh/*MarkDup.sh");
#system "mv $dOut/work_sh/SplitNTrim.sh $dOut/work_sh/SplitNTrim.sh.bak" if(-e "$dOut/work_sh/*SplitNTrim.sh");
#system "mv $dOut/work_sh/RealignerTargetCreator.sh $dOut/work_sh/RealignerTargetCreator.sh.bak" if(-e "$dOut/work_sh/*RealignerTargetCreator.sh");
#system "mv $dOut/work_sh/IndelRealigner.sh $dOut/work_sh/IndelRealigner.sh.bak" if(-e "$dOut/work_sh/*IndelRealigner.sh");
#system "mv $dOut/work_sh/BaseRecalibrator.sh $dOut/work_sh/BaseRecalibrator.sh.bak" if(-e "$dOut/work_sh/*BaseRecalibrator.sh");
#system "mv $dOut/work_sh/PrintReads.sh $dOut/work_sh/PrintReads.sh.bak" if(-e "$dOut/work_sh/*PrintReads.sh");

my @inBam_split;
my @inBam_recal;
open (SH1,">$dOut/work_sh/s$stepNum.GroupRead.sh") or die $!;
&MKDIR("$tmppath/GroupRead");
my $GroupRead_tmp="$_java_/GroupRead";
$stepNum++;

open (SH2,">$dOut/work_sh/s$stepNum.MarkDup.sh") or die $!;
&MKDIR("$tmppath/MarkDup");
my $MarkDup_tmp="$_java_/MarkDup";
$stepNum++;

open (SH3,">$dOut/work_sh/s$stepNum.SplitNTrim.sh") or die $!;
&MKDIR("$tmppath/SplitNTrim");
my $SplitNTrim_tmp="$_java_/SplitNTrim";
$stepNum++;

foreach my $inBam (@inBam){
(my $samName = basename($inBam)) =~ s/\.bam$//;
#print $samName."\n";
my $cmd1=join " ",($GroupRead_tmp,'-jar',$_GroupRead_,
								"I=$inBam",
								"O=$inBam.sorted.bam",
								"SO=coordinate",
								"RGID=$samName",
								"RGLB=library",
								"RGPL=ILLUMINA",
								"RGPU=machine",
								"RGSM=$samName");print SH1 $cmd1."\n";

my $cmd2=join " ",($SplitNTrim_tmp,'-jar',$_MarkDup_, "I=$inBam.sorted.bam",
											"O=$inBam.dedupped.bam",
											"CREATE_INDEX=true",
											"VALIDATION_STRINGENCY=SILENT",
											"M=output.metrics");			print SH2 $cmd2."\n";

#########################################    split n trim
my $cmd3=join " ",($SplitNTrim_tmp,'-jar',$_GATK_, "-T SplitNCigarReads",  # 1h
											"-R $genomeFile",
											"-I $inBam.dedupped.bam",
											"-o $inBam.split.bam",
											"-rf ReassignOneMappingQuality",
											"-RMQF",255,     ##########   need attention
											"-RMQT",60,     ###########    need attention
											"-U ALLOW_N_CIGAR_READS");	

$cmd3=join " ",($cmd3," -fixMisencodedQuals") if (defined $qphred);
print SH3 $cmd3."\n";
push @inBam_split,"$inBam.split.bam";
} #end of @inBam
close(SH1);close(SH2);close(SH3);


if($doINDEL){
open (SH4,">$dOut/work_sh/s$stepNum.RealignerTargetCreator.sh") or die $!;
&MKDIR("$tmppath/RealignerTargetCreator");
my $RealignerTargetCreator_tmp="$_java_/RealignerTargetCreator";
$stepNum++;

open (SH5,">$dOut/work_sh/s$stepNum.IndelRealigner.sh") or die $!;
&MKDIR("$tmppath/Indel_Realigner");
my $Indel_Realigner_tmp="$_java_/Indel_Realigner";
$stepNum++;

foreach my $inBam (@inBam) {
#########################################   Indel Realigner  (optional) #######################
	my $cmd4=join " ",($RealignerTargetCreator_tmp,'-jar',$_GATK_, "-T RealignerTargetCreator",
												"-R $genomeFile",
												"-I $inBam.split.bam",
												"-o $inBam.target_intervals.list ",
												);			print SH4 $cmd4."\n";

	my $cmd5=join " ",($Indel_Realigner_tmp,'-jar',$_GATK_, "-T IndelRealigner",
												"-R $genomeFile",
												"-I $inBam.split.bam",
												"-o $inBam.realigned.bam",
												"-targetIntervals","$inBam.target_intervals.list",
												);			print SH5 $cmd5."\n";
}
}
close(SH4);close(SH5);
##########################################    variants Calling #############################
if($doRecal){
my ($cmd6,$cmd7_1,$cmd7_2,$cmd7_3,$cmd_filt1)=&variant_Calling("$resultdir/noRecal.vcf",1,\@inBam_split);
open (SH6,">$dOut/work_sh/s$stepNum.VariantCalling1.sh") or die $!;
print SH6 $cmd6;
close(SH6);
$stepNum++;
open (SH7_1,">$dOut/work_sh/s$stepNum.VariantCombine1_1.sh") or die $!;
print SH7_1 $cmd7_1;
close(SH7_1);
$stepNum++;
open (SH7_2,">$dOut/work_sh/s$stepNum.VariantCombine1_2.sh") or die $!;
print SH7_2 $cmd7_2;
close(SH7_2);
$stepNum++;
open (SH7_3,">$dOut/work_sh/s$stepNum.VariantConvert1.sh") or die $!;
print SH7_3 $cmd7_3;
close(SH7_3);
$stepNum++;
open (FILT1,">$dOut/work_sh/s$stepNum.VariantFiltration1.sh") or die $!;
print FILT1 $cmd_filt1;
close(FILT1);
$stepNum++;
}
#########################################   Base Recalibrator ##########################
if($doRecal){
open (SH8,">$dOut/work_sh/s$stepNum.BaseRecalibrator.sh") or die $!;
&MKDIR("$tmppath/BaseRecalibrator");
my $BaseRecalibrator_tmp="$_java_/BaseRecalibrator";
$stepNum++;
open (SH9,">$dOut/work_sh/s$stepNum.PrintReads.sh") or die $!;
&MKDIR("$tmppath/PrintReads");
my $PrintReads_tmp="$_java_/PrintReads";
$stepNum++;
foreach my $inBam (@inBam){
	my $cmd8=join " ",($BaseRecalibrator_tmp,'-jar',$_GATK_, "-T BaseRecalibrator",
												"-R $genomeFile",
												"-I $inBam.split.bam",
												"-o $inBam.recal.table",
												"-knownSites ","$resultdir/noRecal.vcf",
												);			print SH8 $cmd8."\n";


	my $cmd9=join " ",($PrintReads_tmp,'-jar',$_GATK_, "-T PrintReads",
												"-R $genomeFile",
												"-I $inBam.split.bam",
												"-o $inBam.recal.bam",
												"-BQSR $inBam.recal.table",
												);			print SH9 $cmd9."\n";

	push @inBam_recal,"$inBam.recal.bam";
}
}
else{
	@inBam_recal=@inBam_split;
}
close(SH8);close(SH9);
##########################################    variants Calling 2 ####################
my ($cmd10,$cmd11_1,$cmd11_2,$cmd11_3,$cmd_filt2)=&variant_Calling("$resultdir/final.vcf",0,\@inBam_recal);

open (SH10,">$dOut/work_sh/s$stepNum.VariantCalling2.sh") or die $!;
print SH10 $cmd10;
close(SH10);
$stepNum++;
open (SH11_1,">$dOut/work_sh/s$stepNum.VariantCombine2_1.sh") or die $!;
print SH11_1 $cmd11_1;
close(SH11_1);
$stepNum++;
open (SH11_2,">$dOut/work_sh/s$stepNum.VariantCombine2_2.sh") or die $!;
print SH11_2 $cmd11_2;
close(SH11_2);
$stepNum++;
open (SH11_3,">$dOut/work_sh/s$stepNum.VariantConvert2.sh") or die $!;
print SH11_3 $cmd11_3;
close(SH11_3);
$stepNum++;
open (FILT2,">$dOut/work_sh/s$stepNum.VariantFiltration2.sh") or die $!;
print FILT2 $cmd_filt2;
close(FILT2);
$stepNum++;
#######################  pipeline ##########################
if($step==1){
	&Cut_shell_qsub("$dOut/work_sh/*CreatIndex.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*CreatIndex.sh");
	$step++ if ($auto==1);
}if($step==2){
	&Cut_shell_qsub("$dOut/work_sh/*GroupRead.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*GroupRead.sh");
	$step++ if ($auto==1);
}if($step==3){
	&Cut_shell_qsub("$dOut/work_sh/*MarkDup.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*MarkDup.sh");
	$step++ if ($auto==1);
}if($step==4){
	&Cut_shell_qsub("$dOut/work_sh/*SplitNTrim.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*SplitNTrim.sh");
	$step++ if ($auto==1);
}if($step==5){
	if($doINDEL!=0){
		&Cut_shell_qsub("$dOut/work_sh/*RealignerTargetCreator.sh",$cpu,$vf,$queue);
		&Check_qsub_error("$dOut/work_sh/*RealignerTargetCreator.sh");
		&Cut_shell_qsub("$dOut/work_sh/*IndelRealigner.sh",$cpu,$vf,$queue);
		&Check_qsub_error("$dOut/work_sh/*IndelRealigner.sh");
		for my $inBam (@inBam){
			`mv $inBam.split.bam $inBam.split.bam.old`;
			`mv $inBam.split.bai $inBam.split.bai.old`;
			`ln -s $inBam.realigned.bam $inBam.split.bam`;
			`ln -s $inBam.realigned.bai $inBam.split.bai`;
		}
	}
	$step++ if ($auto==1);
}if($step==6){
	if ($doRecal==1) {
	&Cut_shell_qsub("$dOut/work_sh/*VariantCalling1.sh",$cpu,$vf,$queue);
	&check_nct_undo_reqsub("$dOut/work_sh/*VariantCalling1.sh");
	&Cut_shell_qsub("$dOut/work_sh/*VariantCombine1_1.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*VariantCombine1_1.sh");
	&Cut_shell_qsub("$dOut/work_sh/*VariantCombine1_2.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*VariantCombine1_2.sh");
	&Cut_shell_qsub("$dOut/work_sh/*VariantConvert1.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*VariantConvert1.sh");
	&Cut_shell_qsub("$dOut/work_sh/*VariantFiltration1.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*VariantFiltration1.sh");
	`grep -P "^\#" $resultdir/noRecal.vcf >$resultdir/tmp`;
	`grep -P "PASS" $resultdir/noRecal.vcf >>$resultdir/tmp`;
	`rm -rf $resultdir/noRecal.vcf*`;
	`mv $resultdir/tmp $resultdir/noRecal.vcf`;
	&run_or_die("$_java_ -jar $_GATK_ -T ValidateVariants -V $resultdir/noRecal.vcf -R $genomeFile");
}
	$step++ if ($auto==1);
}if($step==7){
	if ($doRecal==1) {
	&Cut_shell_qsub("$dOut/work_sh/*BaseRecalibrator.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*BaseRecalibrator.sh");
	&Cut_shell_qsub("$dOut/work_sh/*PrintReads.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*PrintReads.sh");
	}
	$step++ if ($auto==1);
}if($step==8){
	if ($doRecal==0) {
		for my $inBam (@inBam){
			`ln -s $inBam.split.bam $inBam.recal.bam`;
			`ln -s $inBam.split.bai $inBam.recal.bai`;
		}
	}
	&Cut_shell_qsub("$dOut/work_sh/*VariantCalling2.sh",$cpu,$vf,$queue);
	&check_nct_undo_reqsub("$dOut/work_sh/*VariantCalling2.sh");
	&Cut_shell_qsub("$dOut/work_sh/*VariantCombine2_1.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*VariantCombine2_1.sh");
	&Cut_shell_qsub("$dOut/work_sh/*VariantCombine2_2.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*VariantCombine2_2.sh");
	&Cut_shell_qsub("$dOut/work_sh/*VariantConvert2.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*VariantConvert2.sh");
	&Cut_shell_qsub("$dOut/work_sh/*VariantFiltration2.sh",$cpu,$vf,$queue);
	&Check_qsub_error("$dOut/work_sh/*VariantFiltration2.sh");
	`grep -P "^\#" $resultdir/final.vcf >$resultdir/tmp`;
	`grep -P "PASS" $resultdir/final.vcf >>$resultdir/tmp`;
	`rm -rf $resultdir/final.vcf*`;
	`mv $resultdir/tmp $resultdir/final.vcf`;
	&run_or_die("$_java_ -jar $_GATK_ -T ValidateVariants -V $resultdir/final.vcf -R $genomeFile");
	$step++ if ($auto==1);
}
if ($step==9) {
	my ($snp_vcf_file, $indel_vcf_file)=&extract_SNP_and_Indel("$resultdir/final.vcf");
	`perl $_vcf2snplist_ -i $snp_vcf_file -o $dOut/final.snp.list -ref 1`;
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

#################################################################################################
sub check_nct_undo_reqsub{
	my $sh=shift;
        chomp($sh=`ls $sh`);
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";
	my %checkfile;
	my $nctshell=$sh.".nct.redo.sh";
	if ($#sh_file!=$#Check_file) {
		foreach my $check (@Check_file) {
			$checkfile{basename($check)}=1;
		}
		foreach my $shfile (@sh_file) {
			my $name=basename($shfile);
			if (!exists $checkfile{"$name.Check"}){
				`cat $shfile|sed \"s/-nct *[0-9]*//g\" >> $nctshell`;
			}
		}
		&Cut_shell_qsub("$nctshell",$cpu,$vf,$queue);
		&Check_qsub_error("$nctshell");
	}
	else{
		&Check_qsub_error("$sh");
	}
}

sub variant_Calling{
	my ($vcfname,$strict,$in)=@_;
	my $cmd1="";
	my $cmd2="";
	my $cmd3="";
	my @list_files = &get_chr_list_file();
	my $dir=dirname($vcfname);
	my @vcf_files = ();
	if ($ploidy==2) {  # using HC
		my $cmd2_1="";
		my $cmd2_2="";
		my @samgvcf=();

		&MKDIR("$tmppath/HaplotypeCaller");
		my $HaplotypeCaller_tmp="$_java_/HaplotypeCaller";

		foreach my $inBam (@$in){
			my @gvcf_files=();
			for (my $i=1; $i<=@list_files; $i++){
				my $list_file = $list_files[$i-1] ;
				my $basename=basename($inBam).".$i.gvcf";
				my $outfile = "$dir/$basename";
				$cmd1.=join " ",($HaplotypeCaller_tmp,'-jar',$_GATK_, "-T HaplotypeCaller",
													"-R $genomeFile",
													"-recoverDanglingHeads -dontUseSoftClippedBases",
													"-stand_call_conf","20.0",    #   min Q value of base
													"-stand_emit_conf","20.0",
													"-nct",$Thread,"-mte",
													"--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000",
#													"-pairHMM VECTOR_LOGLESS_CACHING", # new feature in GATK 3.1+ for HaplotypeCaller ,but works with 2.6.30 kernel in the official suggestions
													"-o $outfile -L $list_file",
													" -I $inBam",
													);
				push @gvcf_files,$outfile;
				$cmd1.="\n";
			}
			&MKDIR("$tmppath/CombineGVCFs");
			my $CombineGVCFs_tmp="$_java_/CombineGVCFs";
			$cmd2_1.="$CombineGVCFs_tmp -jar $_GATK_ -T CombineGVCFs -R $genomeFile "; # combine  list gvcf file for one sample
			foreach my $gvcf (@gvcf_files) {
				$cmd2_1.="-V $gvcf ";
			}
			my $samfile="$inBam.gvcf";
			$cmd2_1.="-o $samfile\n";
			push @samgvcf,$samfile;
		}
		my $CombineGVCFs_tmp="$_java_/CombineGVCFs";
		$cmd2_2="$CombineGVCFs_tmp -jar $_GATK_ -T CombineGVCFs -R $genomeFile -o $vcfname.rawgvcf "; # combine  list gvcf file for all sample
		foreach my $gvcf (@samgvcf) {
			$cmd2_2.="-V $gvcf ";
		}
		$cmd2_2.="\n";

		&MKDIR("$tmppath/GenotypeGVCFs");
		my $GenotypeGVCFs_tmp="$_java_/GenotypeGVCFs";
		my $cmd2_3 = "$GenotypeGVCFs_tmp -jar $_GATK_ -T GenotypeGVCFs -R $genomeFile -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $vcfname.rawvcf --variant $vcfname.rawgvcf\n"; # convert gvcf to vcf

		&MKDIR("$tmppath/VariantFiltration");
		my $VariantFiltration_tmp="$_java_/VariantFiltration";
		$cmd3="$VariantFiltration_tmp -jar $_GATK_ -T VariantFiltration -R $genomeFile -V $vcfname.rawvcf -o $vcfname "; # filter vcf file
		if($strict==0){
			$cmd3.="-window $window -cluster $cluster -filterName FS -filter \"FS > $FS\" -filterName QD -filter \"QD < $QD\"";
		}else{
			$cmd3.="-window $window_strict -cluster $cluster_strict -filterName FS -filter \"FS > $FS_strict\" -filterName QD -filter \"QD < $QD_strict\"";
		}
		$cmd3.="\n";
		return ($cmd1,$cmd2_1,$cmd2_2,$cmd2_3,$cmd3);
	}
	else{ # Using UG
		die("ploidy more than 2 is NOT available!!\n\n");
		for (my $i=1; $i<=@list_files; $i++){
			my $list_file = $list_files[$i-1] ;
			(my $basename = basename($vcfname)) =~ s/\.vcf$/".$i.vcf"/e ;
			my $outfile = "$dir/$basename" ;
			$cmd1.=join " ",($_java_,'-jar',$_GATK_, "-T UnifiedGenotyper",
												"-R $genomeFile",
												"-ploidy $ploidy",
												"-glm BOTH",
												"-stand_call_conf","20.0",    #   min Q value of base
												"-stand_emit_conf","20.0",
												"-nct",$Thread,
												"-mte",
												"-L $list_file -o $outfile ",
												);
			foreach my $inBam (@$in){
				$cmd1.=" -I $inBam ";
			}
			$cmd1.="\n" ;
			push @vcf_files, $outfile ;
		}
	$cmd2 = "$_java_ -jar $_GATK_ -T CombineVariants -R $genomeFile " ;
	for my $vcf_file (@vcf_files){
		$cmd2 .= " -V $vcf_file" ;
	}
	$cmd2 .= " -o $vcfname.rawvcf\n" ;
	$cmd3=join " ",($_java_,'-jar',$_GATK_, "-T VariantFiltration",
											"-R $genomeFile",
											"-V $vcfname.rawvcf",
											"-o","$vcfname ",
											);
	if($strict==0){
		$cmd3.=join " ",("-window","$window",
							"-cluster","$cluster",
							"-filterName FS",
							"-filter \"FS > $FS\"",
							"-filterName QD",
							"-filter \"QD < $QD\"",
							);
	}else{
		$cmd3.=join " ",("-window","$window_strict",
					"-cluster","$cluster_strict",
					"-filterName FS",
					"-filter \"FS > $FS_strict\"",
					"-filterName QD",
					"-filter \"QD < $QD_strict\"",
					);
	}
	$cmd3.="\n";
	return ($cmd1,$cmd2,$cmd3);
	}
}
sub read_config {#
	my($cfg)=@_;
	open (IN,"$cfg") || die "$!";
	while (<IN>) {
		chomp;
		s/\r$//;s/\s+//;s/\s+$//;
		next if (/\#/ || /$/);
		my @tmp=split /\s+/,$_;
		if ($tmp[0]=~m/Sample/) {
			my $fq1=<IN>;chomp $fq1;
			my $fq2=<IN>;chomp $fq2;
			my @fq_1=split /\s+/,$fq1;
			$sample{$tmp[1]}{FQ1}=$fq_1[1];
			my @fq_2=split /\s+/,$fq2;
			$sample{$tmp[1]}{FQ2}=$fq_2[1];
		}
		$para{$tmp[0]}=$tmp[1];
	}
	close IN;
}
sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
#	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}
sub get_chr_list_file()
{
	my $dir_name;
	my $fa_len_file = "$genomeFile.len" ;
	$dir_name = dirname($genomeFile);
	if (!-f $fa_len_file){
		my $cmd = "perl $_tools_/ref_GC_len.pl -ref $genomeFile -od $dir_name ";
		&run_or_die($cmd);
	}
	$chr_num=`cat $dir_name/*.GC|wc -l`;
	$chr_num=$chr_num-2;
	$chr_num=30 if ($chr_num>30);   #    set 30 if chr > 30
	# get list file
	&MKDIR("$dOut/tmp");
	my $outfile = "$dOut/tmp/chr_allocation" ;
	if (!-e "$outfile.1.list") {
		my $cmd = "perl $_tools_/distribute.pl -i $fa_len_file -o $outfile -n $chr_num" ;
		&run_or_die($cmd) ;
	}
	my @list_files = ();
	for (my $i=1; $i<=$chr_num; $i++){
		push @list_files, "$outfile.$i.list";
	}
	return (@list_files);
}
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	chomp($shell=`ls $shell`);
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;
	die "$shell does not exist\n" unless -e $shell;
	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		&qsubOrDie($shell,$queue,$cpu,$vf);
		#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
	}
	if ($line>1000) {
		my @div=glob "$shell.div*sh";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div=();
		my $div_index=1;
		my $line_num=1;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index.sh" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=0;
				close OUT;
			}
		}
		if ($line_num!=0) {
			close OUT;
		}
		@div=glob "$shell.div*sh";
		foreach my $div_file (@div) {
			&qsubOrDie($div_file,$queue,$cpu,$vf);
			#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
}
sub check_file {#
	my ($file) = @_;
	if (-e $file) {
		return 1;
	}else{
		return 0;
	}
}
sub Check_qsub_error {#
	# Check The qsub process if error happend
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";

	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else {
		print "$sh qsub is Done!\n";
	}
}
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

## show log
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &GetTime();
     print "[$Time] $txt\n";
	return ($time) ;
}

#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}
sub extract_SNP_and_Indel()
{
	my ($variant_file) = @_ ;
	my $dir=$dOut;
	my $basename = basename($variant_file) ;
	(my $snp_vcf_file = "$dir/$basename") =~ s/vcf$/snp.vcf/ ;
	(my $indel_vcf_file = "$dir/$basename") =~ s/vcf$/indel.vcf/ ;
	my $cmd = "$_java_ -jar $_GATK_ -T SelectVariants -R $genomeFile -V $variant_file -selectType SNP -o $snp_vcf_file" ;
	&run_or_die($cmd) ;
	$cmd = "$_java_ -jar $_GATK_ -T SelectVariants -R $genomeFile -V $variant_file -selectType INDEL -o $indel_vcf_file" ;
	&run_or_die($cmd) ;

	return ($snp_vcf_file, $indel_vcf_file);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Shi Tongwei <shitw\@biomarker.com.cn>
Discription:
	2015/02/04 v1.1 update: adapt to STAR pipeline v1.1 (bam output)
                                  nct bug shell reqsub
	2015/03/12 v1.1.1 update: fix a "filename" bug while step>4 (bam output)
		----------------------------------------------------------
		step 1 : create index and dictionary file;
		step 2 : Add read groups and sort;
		step 3 : mark duplicates;
		step 4 : Split'N'Trim and reassign mapping qualities;
		step 5 : (optional)Indel Realignment,                       using -doIndel parameter to define;
		step 6 : (optional)Generate base recalibration requirement, using -doRecal parameter to define;
		step 7 : (optional)Base Recalibration,                      using -doRecal parameter to define;
		step 8 : Variant calling and filtering.
		step 9 : Select SNP and Indel.
		----------------------------------------------------------
Usage:
  Options:
  -aligndir	<file>	required	input alignment bam/sam dir
  -ref		<file>	required	reference fasta file
  -queue    <str>    required   the queue to run the cmd
  -ploidy	<num>	optional	sample ploidy,default[2]
  -win		<num>	optional	filter window size,default [35]
  -clu		<num>	optional	filter cluster of >=\$clu SNPs within \$win bp window,default [3]
  -QD		<num>	optional	filter SNPs according to Qual Depth value,default [2.0]
  -FS		<num>	optional	filter SNPs according to Fisher Strand value,default [30.0]
  -step		<num>	optional	run from this step,default [1]
  -auto		<num>	optional	auto run from the step,default[1]
  -doIndel	<num>	optional	whether to run Indel Realignment,1/0, default [0]
  -doRecal	<num>	optional	whether to run Base Recalibration,1/0, default [0]
  -od		<str>	optional	Directory where output file produced,default [./SNP_Result]
  -help		Help
USAGE
	print $usage;
	exit;
}