#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.1";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($fRef,$fF1,$fF2,$fK,$Verbose,$queue,$SPLIT_SIZE,$MAX_MISMATCH,$MAX_INDEL,$MAX_LEN_OF_INTRO,$MIN_LEN_OF_INTRO,$MAX_INSERT_SIZE,$MAX_Proc,$Step,$NOISE);
GetOptions(
				"help|?" =>\&USAGE,
				"r:s"=>\$fRef,
				"1:s"=>\$fF1,
				"2:s"=>\$fF2,
				"o:s"=>\$fK,
				"v"=>\$Verbose,
				"queue:s"=>\$queue,
				"SPLIT_SIZE:i"=>\$SPLIT_SIZE,
				"MAX_INSERT_SIZE:i"=>\$MAX_INSERT_SIZE,
				"MAX_MISMATCH:i"=>\$MAX_MISMATCH,
				"MAX_INDEL:i"=>\$MAX_INDEL,
				"MIN_LEN_OF_INTRO:i"=>\$MIN_LEN_OF_INTRO,
				"MAX_LEN_OF_INTRO:i"=>\$MAX_LEN_OF_INTRO,
				"MAX_Proc:i"=>\$MAX_Proc,
				"NOISE:i"=>\$NOISE,
				"Step:i"=>\$Step,
				) or &USAGE;
&USAGE unless ($fRef and $fF1 and $fF2 and $fK) ;

$SPLIT_SIZE||=30_000;
$MAX_INSERT_SIZE||=300;
$MAX_MISMATCH||=3;
$MAX_INDEL||=2;
$MIN_LEN_OF_INTRO||=30;
$MAX_LEN_OF_INTRO||=8_000;
$MAX_Proc||=40;
$Step||=1;
$NOISE||=10;

$fRef=&AbsolutePath("file",$fRef);
$fF1=&AbsolutePath("file",$fF1);
$fF2=&AbsolutePath("file",$fF2);
mkdir dirname($fK) unless (-d dirname($fK)) ;
$fK=&AbsolutePath("file",$fK);
my $PRD=`pwd`;chomp $PRD;$PRD.="/";
my $WD=dirname($fK);$WD.="/";
my $Key=basename($fK);

my $file1Format=&fileFormat($fF1);
my $file2Format=&fileFormat($fF2);
if (&fileFormat($fRef) ne "fasta" or $file1Format ne "fasta" and $file1Format ne "fastq" or $file2Format ne "fasta" and $file2Format ne "fastq") {
	print STDERR "Genome file must be fasta format and reads file must be fastq or fasta format!\n";
	exit;
}

print STDERR "Genome: $fRef\n";
print STDERR "read1: $fF1, $file1Format format\n";
print STDERR "read2: $fF2, $file2Format format\n";
print STDERR "outKey: $fK\n";
print STDERR "Program Run Direction: $PRD\n";
print STDERR "Working Direction: $WD\n";


# ------------------------------------------------------------------
# Split Reads into "splitedReads" direction
# ------------------------------------------------------------------
if ($Step==1) {
	&stepStart(1,"Split Reads.");
	mkdir "$WD/splitedReads" unless (-d "$WD/splitedReads") ;
	my @subFileList1;
	my @subFileList2;
	my $i=1;
	open (READ1,"<",$fF1) or die $!;
	open (READ2,"<",$fF2) or die $!;
	open (OUT1,">",$WD."splitedReads/".basename($fF1)."_$i.fa") or die $!;
	open (OUT2,">",$WD."splitedReads/".basename($fF2)."_$i.fa") or die $!;
	push @subFileList1,$WD."splitedReads/".basename($fF1)."_$i.fa";
	push @subFileList2,$WD."splitedReads/".basename($fF2)."_$i.fa";
	while (my $id1=<READ1>) {
		my $seq1=<READ1>;
		my $id2=<READ2>;
		my $seq2=<READ2>;
		$id1=~s/^@/>/,<READ1>,<READ1> if ($file1Format eq "fastq") ;
		$id2=~s/^@/>/,<READ2>,<READ2> if ($file2Format eq "fastq") ;
		print OUT1 $id1,$seq1;
		print OUT2 $id2,$seq2;
		$i++;
		if ($i%$SPLIT_SIZE==0) {
			close (OUT1) ;
			close (OUT2) ;
			open (OUT1,">",$WD."splitedReads/".basename($fF1)."_$i.fa") or die $!;
			open (OUT2,">",$WD."splitedReads/".basename($fF2)."_$i.fa") or die $!;
			push @subFileList1,$WD."splitedReads/".basename($fF1)."_$i.fa";
			push @subFileList2,$WD."splitedReads/".basename($fF2)."_$i.fa";
		}
	}
	close (OUT1) ;
	close (OUT2) ;
	close (READ1) ;
	close (READ2) ;
	my $Total_Reads=$i-1;
	$Step++;
	&stepTime(1);
}



# ------------------------------------------------------------------
# Blat
# ------------------------------------------------------------------
if ($Step==2) {
	&stepStart(2,"Blat.");

	my @subFileList1=glob("$WD/splitedReads/".basename($fF1)."_*.fa");
	my @subFileList2=glob("$WD/splitedReads/".basename($fF2)."_*.fa");

	my $fSH=$WD."blat.sh";
	open (SH,">",$fSH) or die $!;
	for (my $i=0;$i<@subFileList1 ;$i++) {

		if ($Step<=2) {
			print SH "blat $fRef $subFileList1[$i] $subFileList1[$i].psl -tileSize=10 -stepSize=6 -noHead   &&  ";
			print SH "blat $fRef $subFileList2[$i] $subFileList2[$i].psl -tileSize=10 -stepSize=6 -noHead   &&  ";
		}

		print SH "perl $Bin/subClassify.pl $subFileList1[$i].psl $subFileList2[$i].psl $subFileList1[$i].rmap $subFileList1[$i].Classify.stat $MAX_MISMATCH $MAX_INDEL $MAX_LEN_OF_INTRO $MIN_LEN_OF_INTRO $MAX_INSERT_SIZE $NOISE&&  ";
		print SH "\n";
	}
	close (SH) ;
	&Cut_shell_qsub($fSH,$MAX_Proc,"3G",$queue);
	#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $MAX_Proc --resource vf=3G --reqsub --independent $fSH ";
	&stepTime(2);
	$Step++;
}


# ------------------------------------------------------------------
# Get all the result
# ------------------------------------------------------------------
if ($Step==3) {
	&stepStart(3,"Get all sub classified reads");
	my $dPSL=$WD."splitedReads";
	`cat $dPSL/*.rmap >$fK.rmap`;
	&stepTime(3);
	$Step++;
}


# ------------------------------------------------------------------
# get insertSize
# ------------------------------------------------------------------
if ($Step==4) {
	&stepStart(4,"get data insertSize");
	print STDERR "perl $Bin/draw_insertlen_by_rmap.pl -i $fK.rmap -od $WD\n";
	`perl $Bin/draw_insertlen_by_rmap.pl -i $fK.rmap -od $WD `;
	&stepTime(4);
	$Step++;
}

# ------------------------------------------------------------------
# Statistics maped info.
# ------------------------------------------------------------------
if ($Step==5) {
	&stepStart(5,"Stat mapped information");
	print STDERR "perl $Bin/rmap_file_stat.pl -i $fK.rmap -file1 $fF1 -file1Format $file1Format -od $WD\n";
	`perl $Bin/rmap_file_stat.pl -i $fK.rmap -file1 $fF1 -file1Format $file1Format -od $WD`;
	&stepTime(5);
	$Step++;
}




#######################################################################################
&timeLog("$Script Done. Total elapsed time : ".time()-$BEGIN_TIME."s");
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub fileFormat()
{	#��ѯ�ļ���ʽ��Ŀǰ֧��fasta��fastq��GenBank,unkonw
	my $file =shift;
	open (IN,"<",$file) or die $!;
	my $type=<IN>;
	close (IN) ;
	chomp $type;
	if    ($type=~/^>/){
		$type="fasta";
	}elsif($type=~/^@/){
		$type="fastq";
	}elsif($type=~/LOCUS/){
		$type="GenBank";
	}else{
		$type="nuknow";
	}
	return $type;
}

sub AbsolutePath
{		#��ȡָ��Ŀ¼���ļ��ľ���·��
        my ($type,$input) = @_;

        my $return;

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	chomp(my $line = `wc -l $shell|cut -f 1 -d" "`);
	if ($line<=1000) {
		&qsubOrDie($shell,$queue,$cpu,$vf);
		#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $cpu --resource vf=$vf --independent --reqsub $shell ";
	}
	if ($line>1000) {
		my @div=glob "$shell.div*sh";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
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
			#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $cpu --resource vf=$vf --reqsub $div_file ";
		}
	}
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: Rmap.pl (map mRNA reads to reference, generate related files.)
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>

Description:
	This program is a memeber of Transcript Analysis Kit, It can map mRNA reads to reference
	and generate PE, MP, SO, ST reads info ,junction info and InsertSize info which can be
	used by other program for Transcript Analysis.
Usage:
  Forced Options:
  -r  <file>   Genome(genome in fasta format)
  -1  <file>   mRNA-seq read1 file, Both fasta and fastq format are accepted
  -2  <file>   mRNA-seq read2 file, Both fasta and fastq format are accepted
  -o  <Str>    Key of output files.  (key.ramp, key.insertSize, key.Mapped.stat Info)

  Optional Options:
  -Step               <int>    0,From begining; 1,From Split; 2,From Blat; 3,From Classify; 4,From Cat; 5,From InsertSize; 6,From STAT; default 0;

  -SPLIT_SIZE         <int>    Split size, split reads file into small pieces to svae align time. default is 30,000;
  -MAX_INSERT_SIZE    <int>    Max InsertSize, default 250;
  -MAX_MISMATCH       <int>    Max MisMatch, default 3;
  -MAX_INDEL          <int>    Max InDel, default 2;
  -MIN_LEN_OF_INTRO   <int>    Min len of intro, default 30;
  -MAX_LEN_OF_INTRO   <int>    Max len of intro, default 8,000;
  -MAX_Proc           <int>    Max proc number, used for qsub, default 40;
  -NOISE              <int>    Noise, This used in subClassify.pl or GAP_COVER

Output:
  The program will output these files:
	outKey.rmap               Classified Reads:PE, MP, SO, ST.
	outKey.Mapped.stat Info     Reads mapped information statistics file.
	outKey.insertSize         Distribution of insertSize of Reads.
	outKey.insertSize.svg     SVG file of the Distribution of insertSize.

Advice:
	The parameter -MAX_INSERT_SIZE is very important for classify PE, MP, SO, ST READS, you should give the Biggest InsertSize
	but average InsertSize.

USAGE
	print $usage;
	exit;
}
