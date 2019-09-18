#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $version="1.0.0";
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
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($cfg,$od,$Assembly_dir,$min_kmer_cov,$Result_dir,$step,$worksh);
GetOptions(
				"help|?" =>\&USAGE,
				"cfg:s"=>\$cfg,
                                "min_kmer_cov:i"=>\$min_kmer_cov,
				"od:s"=>\$od,
                "s:s"=>\$step,
				) or &USAGE;
&USAGE unless ($cfg and $od) ;

#########################
my %para;
&LOAD_PARA ($cfg,\%para);

my ($Index,$fq1,$fq2,$seqType,$fragment_length,$min_contig_length,$max_number_of_paths_per_node,$p,$JM,$lib_type);

$Index = $para{Index};
$fragment_length = $para{PEinsert};
$min_contig_length = $para{min_contig};
$seqType = $para{seqType};
$p = $para{thread};
$lib_type = $para{lib_type};
my $queue = $para{queue};
my %sys_cfg=%{&readconf("$Bin/../../../../../../config/sys.cfg")};
$fq1 = $para{FQ1};
$fq2 = $para{FQ2};
####set the min_kmer_cov and JM values
my $file_size;
for my $fq (split/\s+/,$fq1){
	$file_size+=(-s $fq);
}
for my $fq (split/\s+/,$fq2){
	$file_size+=(-s $fq);
}
$file_size=$file_size/(1024**3);
unless(defined $min_kmer_cov){
    if($file_size<=30){
	$min_kmer_cov=1;
    }
    elsif($file_size<=150){
	$min_kmer_cov=2;
    }
    else{
	$min_kmer_cov=3;
    }
}
$JM=(int($file_size*0.4*2)+1)."G";  ###0.4 is the ration of file size and data size 
$fq1=~s/\s+/,/g;
$fq2=~s/\s+/,/g;

&MKDIR ($od);
$od = &ABSOLUTE_DIR ($od);
$worksh = "$od/work_sh";
&MKDIR ($worksh);
$Assembly_dir = "$od/$Index"."_Trinity";
&MKDIR ($Assembly_dir);
$Result_dir = "$od/$Index"."_Result";
&MKDIR ($Result_dir);

$step = $step || 1;

###############Time
&timeLog("$programe_dir Start.");

################################ Program
my $trinity_HOME=$sys_cfg{"trinity"};
########## step1 Trinity Assembly
my $cmd="";
START:
if ($step==1) {
	&stepStart(1,"Trinity phase 1: Clustering of RNA-Seq Reads Process");
	open OUT,">$worksh/$Index.Trinity.chrysalis.sh" || die $!;
	if ($lib_type =~ /none/) {
		if (exists $para{jaccard_clip}){
			$cmd="perl $trinity_HOME/Trinity --seqType $seqType --max_memory $JM --CPU $p --min_contig_length $min_contig_length --group_pairs_distance $fragment_length  --left $fq1 --right $fq2 --min_kmer_cov $min_kmer_cov --bflyHeapSpaceMax 20G --output $Assembly_dir --no_cleanup --jaccard_clip --bflyGCThreads $p --no_distributed_trinity_exec --no_version_check " ;
		}
		else{
			$cmd="perl $trinity_HOME/Trinity --seqType $seqType --max_memory $JM --CPU $p --min_contig_length $min_contig_length --group_pairs_distance $fragment_length  --left $fq1 --right $fq2 --min_kmer_cov $min_kmer_cov --bflyHeapSpaceMax 20G --output $Assembly_dir --no_cleanup --bflyGCThreads $p --no_distributed_trinity_exec --no_version_check " ;
		}
		print OUT $cmd;
		close OUT;
		&runOrDie($cmd);
	}
	elsif ($lib_type !~ /none/) {
		if (exists $para{jaccard_clip}){
			$cmd="perl $trinity_HOME/Trinity --seqType $seqType --SS_lib_type $lib_type --max_memory $JM --CPU $p --min_contig_length $min_contig_length --left $fq1 --right $fq2 --output $Assembly_dir --min_kmer_cov $min_kmer_cov --group_pairs_distance $fragment_length --bflyHeapSpaceMax 20G --no_cleanup --jaccard_clip --bflyGCThreads $p --no_distributed_trinity_exec --no_version_check " ;
		}
		else{
			$cmd="perl $trinity_HOME/Trinity --seqType $seqType --SS_lib_type $lib_type --max_memory $JM --CPU $p --min_contig_length $min_contig_length --left $fq1 --right $fq2 --output $Assembly_dir --min_kmer_cov $min_kmer_cov --group_pairs_distance $fragment_length --bflyHeapSpaceMax 20G --no_cleanup --bflyGCThreads $p --no_distributed_trinity_exec --no_version_check " ;
		}
		print OUT $cmd;
		close OUT;
		&runOrDie($cmd);
	}
	&stepTime(1);
	$step=2;
}

if ($step==2) {
	&stepStart(2,"Trinity phase 2: Parallel Assembly of Read Clusters");
        if(-e "$Assembly_dir/recursive_trinity.cmds.ok"){
        system "mv $Assembly_dir/recursive_trinity.cmds $Assembly_dir/chrysalis/Phase2TrinityCommands";}else{print "Warnings: Phase2 Trinity cmd not exists!\n";}
	open OUT,">$worksh/$Index.Trinity.ParalellAssembly.sh" || die $!;
	print OUT "perl $Bin/quantifyGraph_qsub.pl -id $Assembly_dir -worksh $worksh  -queue $queue";
	close OUT;
        system"sh $worksh/$Index.Trinity.ParalellAssembly.sh >$worksh/$Index.Trinity.ParalellAssembly.sh.log 2>&1";
	&stepTime(2);
	$step=3;
}

my $cmd_cat="find $Assembly_dir/read_partitions/  -name '*inity.fasta'  | $trinity_HOME/util/support_scripts/partitioned_trinity_aggregator.pl TRINITY_DN > $Assembly_dir/$Index.Trinity.fasta";
my $flag=system($cmd_cat);
if($flag!=0){print "Error : conmmand fail: $cmd_cat";die $!;}

#system" find $Assembly_dir/read_partitions/  -name '*inity.fasta'  | $trinity_HOME/util/support_scripts/partitioned_trinity_aggregator.pl TRINITY_DN > $Assembly_dir/    $Index.Trinity.fasta";

if (!-f "$Assembly_dir/$Index.Trinity.fasta" && $step == 3) {
	$step = 1;
	goto START;
}
elsif (-f "$Assembly_dir/$Index.Trinity.fasta" && $step == 3) {
	&stepStart(3,"Trinity Results to Unigenes");
	&MKDIR ("$Result_dir/contigs");
	&MKDIR ("$Result_dir/Transcripts");
	&MKDIR ("$Result_dir/Unigenes");
	`cp $Assembly_dir/$Index.Trinity.fasta $Result_dir/Transcripts/$Index.Transcripts.fa`;
	`cp $Assembly_dir/inchworm.*.fa $Result_dir/contigs/$Index.contigs.fa` unless $para{jaccard_clip};
	system"cp $Assembly_dir/inchworm.*.clipped.fa $Result_dir/contigs/$Index.contigs.fa" if $para{jaccard_clip};
	`perl $Bin/trinity_cluster_with_id.pl -i $Assembly_dir/$Index.Trinity.fasta -od $Result_dir/Unigenes -key $Index `;
	#system "rm -r $Assembly_dir/chrysalis/Component_bins&";
	$step = 4;
	&stepTime(3);
}

if ($step==4) {
	&stepStart(4,"Stat assembly results.");
	open OUT,">$worksh/$Index.result.stat.sh" || die $!;
	print OUT "perl $Bin/assembly_stat.pl -i $Index.contigs -fa $Result_dir/contigs/$Index.contigs.fa -od $Result_dir/contigs \n";
	print OUT "perl $Bin/assembly_stat.pl -i $Index.Transcripts -fa $Result_dir/Transcripts/$Index.Transcripts.fa -od $Result_dir/Transcripts \n";
	print OUT "perl $Bin/assembly_stat.pl -i $Index.Unigenes -fa $Result_dir/Unigenes/$Index.Unigenes.fa -od $Result_dir/Unigenes \n";
	close OUT;

	runOrDie("$worksh/$Index.result.stat.sh");
	&stepTime(4);
}

###############Time

&timeLog("End $programe_dir");
&totalTime();


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
		my ($para_key,$para_value) = split(/\s+/,$_,2);
		next if$para_key eq "jaccard_clip" && $para_value!=1;
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub parse_config
{ # load config file
	my $config_file= shift;
	my $DataBase= shift;
	
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/\s+//;s/\s+$//;s/\r$//;
		next if(/$/ or /\#/);
		my ($software_name,$software_address) = split(/\s+/,$_);
		$DataBase->{$software_name} = $software_address;
		if (! -e $software_address){
			warn "Non-exist:  $software_name  $software_address\n"; 
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}

sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
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
		warn "Warning just for file and dir in $in\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: Pair End Reads Assembly With Trinity
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:

Usage:
  -cfg                <str>                  Index of Output                                 must be given;
  -od                 <str>                  out directory                                   must be given;
  -s                  <int>                  step of process
                                                           1 from Trinity Cluster(that is, chrysalis)          default;
                                                           2 from Trinity Assembly Paralell
                                                           3 produce Unigene
							   4 stastistics the assemble results
USAGE
	print $usage;
	exit;
}
