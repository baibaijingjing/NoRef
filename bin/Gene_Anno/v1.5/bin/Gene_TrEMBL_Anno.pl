#!/usr/bin/perl -w
#
#Copyright (c) BMK 2012
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2012
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2012
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
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;

#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info
# ------------------------------------------------------------------
my ($id,$Blast_cpu,$Blast_e,$Database,$odir,$HELP,$queue);
my @anno;

GetOptions(
		"id:s"=>\$id,
		"Database:s"=>\$Database,
		"od:s"=>\$odir,
		"Blast_cpu:s"=>\$Blast_cpu,
		"Blast_e:s"=>\$Blast_e,
		"queue:s"=>\$queue,
		"help"=>\$HELP
	) or &USAGE;

&USAGE if (!defined $id || !defined $odir || !defined $Database || $HELP) ;


###------------------���·��----------------------------###
# all the program are in this dir
chomp $Bin;
#my $blastall = "/share/nas2/genome/biosoft/blast/2.2.26/bin/blastall";  # 2015-02-04
#my $blastall = "/share/nas2/genome/biosoft/ncbi-blast/2.2.31/bin/blastx";  # lium: 2015-8-21
my %sys_cfg=%{&readconf("$Bin/../../../../config/sys.cfg")};
my $blastx=$sys_cfg{blast}."/blastx";

#=================== һЩ���� ================================
&MKDIR($odir);
$id=&ABSOLUTE_DIR($id);
$odir=&ABSOLUTE_DIR($odir);
$Database=&ABSOLUTE_DIR($Database);


$Blast_cpu||=50;
$Blast_e||=1e-5;
my $notename=`hostname`;chomp $notename;

###############Time
&timeLog("$programe_dir Start.");


my $Q_name=(glob "$id/*.fa")[0];
$Q_name=basename $Q_name;
if ($Q_name=~/(.+)\.\d+\.fa$/) {
	$Q_name=$1;
}
else {print "Your file name is Wrong!\n";die;}

&MKDIR("$odir/Result");
my $Result_dir="$odir/Result";
&MKDIR("$odir/TrEMBL_Dir");
&MKDIR("$odir/work_sh");
&MKDIR("$odir/02.gene-annotation");
my $Tab_dir="$odir/02.gene-annotation";
&MKDIR("$odir/work_sh");
&MKDIR("$odir/work_sh/TrEMBL_sh");
my $sh_dir="$odir/work_sh/TrEMBL_sh";

my $blast_shell_file = "$odir/work_sh/TrEMBL_sh/TrEMBL.blast.sh";
my @subfiles;
@subfiles = glob("$id/*.fa");

############creat shell file
open OUT,">$blast_shell_file" || die "fail $blast_shell_file";
foreach my $subfile (@subfiles) {
	my $name=basename $subfile;
	#print OUT "$blastall -b 100 -v 100 -p blastx -e $Blast_e -F F -d $Database -i $subfile -a 2 -o $odir/TrEMBL_Dir/$name.TrEMBL.blast && \n";
	print OUT "$blastx -task blastx-fast -num_descriptions 100 -num_alignments 100 -evalue $Blast_e -db $Database -query $subfile -num_threads 2 -out $odir/TrEMBL_Dir/$name.TrEMBL.blast && \n";
}
close OUT;

####################run the shell file
&qsubOrDie("$blast_shell_file",$queue,$Blast_cpu,"6G");

`perl $Bin/bin/blast_format.pl -id $odir/TrEMBL_Dir -od $Tab_dir -shdir $sh_dir  -middir $odir/mid --trembl -queue $queue`;


####################
open (TrEMBL,"$Tab_dir/$Q_name.TrEMBL.blast.tab.best")or die "cant open file $Tab_dir/$Q_name.TrEMBL.blast.tab.best";
open (OUT,">$Tab_dir/$Q_name.TrEMBL.blast.tab.best.anno")or die "cant open file $Tab_dir/$Q_name.TrEMBL.blast.tab.best.anno";
print OUT "#NtGeneID\tDatabase_ID\tE_value\tIdentity\tScore\tAnnotation\n";
while (<TrEMBL>)
{
	chomp;
	&ann;
}
close OUT;
close TrEMBL;


###### annotation of each library ######
`cp $Tab_dir/$Q_name.TrEMBL.blast.tab.best.anno $Result_dir/$Q_name.TrEMBL.anno.txt`;



###############Time
&timeLog("$programe_dir End .")
&Runtime($BEGIN);

#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+

sub LOAD_PARA
{
	my $para_file= shift;
	my $para= shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
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

sub runcmd
{
	my ($program,$cmd)=@_;
	open (SH, ">>work.sh") ||  die "Can't write work: $!\n";
	print SH "$cmd \n";
	system($cmd) && LOGFILE(1,$program);
	close SH;
}

sub LOGFILE
{
	my $flog=shift;
	my $program=shift;
	my $Time= sub_format_datetime(localtime(time()));
	print LOG "[$Time +0800]\ttask\t0\tstart\t$program\tStart to analysis......\n";
	if($flog==0){
		print LOG "[$Time +0800]\ttask\t0\tend\t$program\tDone.\n";
	}else{
		print LOG "[$Time +0800]\ttask\t0\terror\t$program\tAt least one $program in this section is in error.\n";
		close LOG;
		exit;
	}
}
close LOG;
sub ann
{
	@anno = split /\t/, $_;
	print OUT "$anno[0]\t$anno[4]\t$anno[13]\t$anno[8]\t$anno[12]\t$anno[15]\n";
}

sub LOAD_SEQ
{
	my ($fa,$info) = @_;

	open IN,"$fa" || die $!;
	$/='>';
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r+$//;
		next if (/^$/ || /^\#/);
		my ($head,$seq)=split/\n+/,$_,2;
		my $id=(split/\s+/,$head)[0];
		$info->{$id}=$seq;
	}
	$/="\n";
	close IN;
}


sub CUTFA
{
	my ($fa,$od,$cut,$name) = @_;

	&MKDIR("$od/$name.div");
	my %seq=%$fa;
	my @aa=sort(keys %seq);
	my $index=0;
	LAB: for (my $i=1;;) {
		my $num=0;
		open OUT,">$od/$name.div/$name.$i.fa" || die $!;
		for ($index..$#aa) {
			$index++;
			if ($num<$cut) {
				print OUT ">$aa[$_]\n$seq{$aa[$_]}\n";
				$num++;
			}
			if ($num>=$cut) {
				$num=0;
				$i++;
				close OUT;
				if ($index==$#aa+1) {
					last;
				}
				else {
					next LAB;
				}
			}
		}
		if ($num) {
			close OUT;
		}
		last;
	}
}



sub cut_str
{
	my $string = shift;
	my @str = split /\s+/,$string;
	if (@str > 2) {
		return "$str[0] $str[1]"
	}else{
		return $string;
	}
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
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}

	chdir $cur_dir;
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
			#system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $queue --maxproc $cpu --resource vf=$vf --reqsub $div_file ";
		}
	}
}

sub USAGE
{
	print <<"	Usage End.";
	Program:Blast && interproscan with NT_DataBase Annotate New Gene
	Version: $version
	Contact: zhang xuechuan <zhangxc\@biomarker.com.cn>

	Description:

      -id                    mRNA dir,forced
      -Database              Nt Database file,forced
      -od                    OUT DIR,forced
      -Blast_cpu             cpu number,default 50
      -Blast_e               e-value of cutoff,default 1e-5


      -help

	Usage End.
	exit;
}
