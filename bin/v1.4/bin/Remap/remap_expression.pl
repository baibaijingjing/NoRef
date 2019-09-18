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
GetOptions(\%opts,"cfg=s","id=s","od=s","h","queue=s");
if (!defined($opts{id}) ||!defined($opts{od}) || !defined($opts{cfg}) || defined($opts{h})) {
	&help;
	exit;
}

#######################
&timeLog("$Script Start.");

######################
my $notename=`hostname`;chomp $notename;
`mkdir $opts{od} ` if(!-d $opts{od});
my $id=&ABSOLUTE_DIR($opts{id});
my $od=&ABSOLUTE_DIR($opts{od});
my $cfg=&ABSOLUTE_DIR($opts{cfg});
my $U_fa=(glob "$id/*.Unigene.fa")[0];
my $U_gff=(glob "$id/*.Unigene.gff")[0];
my $U_orf=(glob "$id/Unigene_Orf/*.Unigene.cds_pep.stat.xls")[0];

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
		<IN>;<IN>;
		my $fq1=<IN>;
		$fq1=(split/\s+/,$fq1)[1];
		my $fq2=<IN>;
		$fq2=(split/\s+/,$fq2)[1];
		$sample{$tmp[1]}{FQ1}=$fq1;
		$sample{$tmp[1]}{FQ2}=$fq2;
	}
}
close IN;


#############################################
#Sample Transcripts TO Sample Unigene
#############################################

my $sam_number;

open SH2,">$od/work_sh/Rmap.sh" || die "$!";
foreach my $sam (keys %sample) {
	&MKDIR("$od/rmap/$sam");
	print SH2 "cd $od/rmap/$sam && perl $Bin/bin/Rmap.pl -r $U_fa -1 $sample{$sam}{FQ1} -2 $sample{$sam}{FQ2} -o $sam -MAX_Proc 50 -queue $opts{queue}&& ";
	print SH2 "perl $Bin/bin/GeneExpression_v1.1.pl -g $U_gff -i $sam -orf $U_orf &&";
#	print SH2 "perl $Bin/bin/SVG_cover_length.pl -i $sam.geneExpression.xls -o $sam &&";
	print SH2 "perl $Bin/bin/gene_tagnum_graph.pl -i1 $sample{$sam}{FQ1} -i2 $sam.geneExpression.xls -value 0.1 -o $sam \n";
}

	&Shell_qsub ("$od/work_sh/Rmap.sh",$opts{queue},20,"10G");
	&qsubCheck("$od/work_sh/Rmap.sh");
	`cp $od/rmap/*/*.geneExpression.xls $od/geneExpression `;
	`cp $od/rmap/*/*.png $od/geneExpression `;
	`cp $od/rmap/*/*.Mapped.stat.xls $od/geneExpression `;


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
{ # &Shell_qsub($sh,$qeue,$cpu,$vf);
	my $sh = shift;
	my $queue = shift;
	my $cpu = shift;
	my $vf=shift;
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
		-id
			indir include *.Unigene.fa *.Unigene.gff and /Unigene_Orf/*.Unigene.cds_pep.stat.xls;
		-od
			outdir for Result;
	Usage End.

	exit;
}

