#!/usr/bin/perl -w
#
#Copyright (c) BMK 2011
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011
my $ver="1.0.0";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;

##############################
my %opts;
GetOptions(\%opts,"cfg=s","i=s","qphred=s","od=s","h","queue=s");
if (!defined($opts{i}) ||!defined($opts{od}) || !defined($opts{cfg}) || defined($opts{h})) {
	&help;
	exit;
}

#######################

&timeLog("$Script Start.");

######################
my $notename=`hostname`;chomp $notename;
`mkdir $opts{od} ` if(!-d $opts{od});
my $i=&ABSOLUTE_DIR($opts{i});
my $od=&ABSOLUTE_DIR($opts{od});
my $cfg=&ABSOLUTE_DIR($opts{cfg});
my $qphred=$opts{qphred};

my $shdir="$od/work_sh"; `mkdir $shdir` if(!-d $shdir);
my $rmapdir="$od/rmap";`mkdir $rmapdir` if (!-d $rmapdir);
my $geneExpression="$od/geneExpression";`mkdir $geneExpression` if (!-d $geneExpression);
################### Load Parmeters
#
#Load the cnofig file
#
#############################################
open (IN,"$cfg") || die "$!";
my %sample;
while (<IN>) {
	chomp;
	s/\r$//;s/^\s+//;s/\s+$//;
	next if (/^\#/ || /^$/);
	my @tmp=split /\s+/,$_;
	if ($tmp[0]=~m/Sample/) {
		my $fq1=<IN>;
		$fq1=(split/\s+/,$fq1)[1];
		my $fq2=<IN>;
		$fq2=(split/\s+/,$fq2)[1];
		$sample{$tmp[1]}{FQ1}=$fq1;
		$sample{$tmp[1]}{FQ2}=$fq2;
	}
}
close IN;

#my $trinity_HOME = "/share/nas2/genome/biosoft/trinity/trinityrnaseq_r2013-02-25"; #~ 2014-12-09
#my $trinity_HOME = "/share/nas2/genome/biosoft/trinity/r20131110"; #2014-12-09 ~
my $trinity_HOME  =${&readconf("$Bin/../../../../../../config/sys.cfg")}{trinity};
my $UTIL_DIR = "$trinity_HOME/util";


#############################################
#Sample Transcripts TO Sample Unigene
#############################################

my $sam_number;
open SH1,">$od/work_sh/Rmap_pre.sh"|| die "$!";
open SH2,">$od/work_sh/Rmap.sh" || die "$!";
&MKDIR("$od/rmap/RefIndex");
system "ln -s $i $od/rmap/RefIndex/TRANS";
print SH1 "cd $od/rmap/RefIndex && export PATH=\$PATH:".${&readconf("$Bin/../../../../../../config/sys.cfg")}{PATH}."  && perl  $UTIL_DIR/align_and_estimate_abundance.pl --transcripts $od/rmap/RefIndex/TRANS --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference \n";
&timeLog("BuildIndex:\n$od/work_sh/Rmap_pre.sh");
&runOrDie("$od/work_sh/Rmap_pre.sh");

foreach my $sam (keys %sample) {
	&MKDIR("$od/rmap/$sam");
        print SH2 "ln -s $od/rmap/RefIndex/*  $od/rmap/$sam/ ";
	print SH2 "&& cd $od/rmap/$sam && export PATH=\$PATH:".${&readconf("$Bin/../../../../../../config/sys.cfg")}{PATH}." && perl  $UTIL_DIR/align_and_estimate_abundance.pl  --transcripts $od/rmap/$sam/TRANS --left $sample{$sam}{FQ1} --right $sample{$sam}{FQ2} --seqType fq  --est_method RSEM --aln_method bowtie --trinity_mode --output_dir ./  --output_prefix $sam";
        print SH2 " --qphred 64" if($qphred == 64);
        print SH2 "\n";
}

	&Shell_qsub ("$od/work_sh/Rmap.sh",$opts{queue},20,"20G");
	&qsubCheck("$od/work_sh/Rmap.sh");

foreach my $sam (keys %sample) {
	open (IN,"$od/rmap/$sam/$sam.genes.results") or die $!;
	<IN>;
	open (OUT,">$od/geneExpression/$sam.geneExpression.xls") or die $!;
	print OUT "#Gene_ID\tLength\tEffective_Length\tTPM\tFPKM\tTranscript_ID(s)\tExpected_Count\n";
	while (<IN>) {
		chomp;
		my @A=split/\s+/,$_;
		print OUT "$A[0]\t$A[2]\t$A[3]\t$A[5]\t$A[6]\t$A[1]\t$A[4]\n";
	}
	close IN;
	close OUT;
}

###############Time
&timeLog("$Script End.");
&Runtime($BEGIN);

###########subs
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
sub Shell_qsub
{ # &Shell_qsub($sh,$qeue,$cpu);
	my $sh = shift;
	my $queue = shift;
	my $cpu = shift;
	my $vf  = shift;
	&qsubOrDie($sh,$queue,$cpu,$vf);
	#`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --reqsub --maxproc $cpu --independent $sh  `;
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

sub help
{
	print <<"	Usage End.";
	Description:
		Version  : $ver
		Writer   : Zhang XueChuan <zhangxc\@biomarker.com.cn>
		Usage    :
		-cfg
			The cleandata and parameter config set;
		-i
			input file *.Transcripts.fa;
		-od
			outdir for Result;
	Usage End.

	exit;
}
