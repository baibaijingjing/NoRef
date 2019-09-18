#!/usr/bin/perl -w
#
# Copyright (c) BMK 2017
# Writer:         baij <baij@biomarker.com.cn>
# Program Date:   2017.1.21
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $program_name=basename($0);

my $ver="1.0";
############################################
my ($cfg, $odir,$oneSam);
GetOptions(

                                "help|?" =>\&USAGE,
                                "cfg:s" => \$cfg,
                                "od:s"=>\$odir,
				"oneSam"=>\$oneSam,
                                ) or &USAGE;
&USAGE unless ($cfg  and $odir);


###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";


###############
my %para;
$cfg=&ABSOLUTE_DIR($cfg);
&para_load($cfg,\%para);

system "mkdir -p $odir" unless -d $odir;
$odir = &ABSOLUTE_DIR($odir);

#---unigene and pep.fa to Needed Data dir              
&runOrDie( "cp $para{Basic}/Cluster/Unigene/*.Unigene.fa  $odir/Unigene.fa" );
my %hash;
my $pep=(glob("$para{ANNO}/../Unigene_CDS_Predict/*.Unigene.pep.fa"))[0];
print "$pep";
open (PEP, "$pep") or die $!;
open (OUT,">$odir/Unigene.pep.fa") or die $!;
$/=">";
<PEP>;
while(<PEP>){
  chomp;
  s/^\s+//;s/\s+$//;s/\r+$//;
  next if (/^#/ || /^$/); 
  my ($head, $seq)=split(/\n+/,$_,2);
  my $id = (split(/\|/,$head))[0];
  $seq=~s/\s+//g;
  my $length=length($seq);
  if(exists $hash{$id}){my $len=length($hash{$id}); if($len<$length){$hash{$id}=$seq;}else{next;}}else{$hash{$id}=$seq;}
  print OUT ">$id\n$seq\n";
}
$/="\n";
close(PEP);
close(OUT);
#print Dumper(\%hash);

#---snp to Needed Data dir
unless(defined $oneSam){&runOrDie( "cp $para{SNP}/SNP/final.snp.list  $odir" );}

################Time
my $Time_End;
$Time_End = &sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
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

sub para_load {
	my ($file,$para)=@_;
	open IN,$file||die "$!";
	while (<IN>) {
		chomp;
		s/\r+//g;
		next if(/^$/||/^\#/);
		my ($key,$info)=(split/\s+/,$_)[0,1];
		if(!$key){print "$_\n";die;}
		$para->{$key}=$info;
	}
	close IN;
}
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "\nTotal elapsed time: ${t}s\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	#`rm -r $dir` if(-d $dir);
	`mkdir -p $dir` if(!-d $dir);
}

sub help{
	print << "	Usage End.";
	Description: Extract No_Ref Transcriptome Reaults for Html Process;
	version:$ver
	Usage:
		-cfg  <STR>   config file   force
		-od   <STR>   outdir        force
		-h            help
	Usage End.
		exit;
}
