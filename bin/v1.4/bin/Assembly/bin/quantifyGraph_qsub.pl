#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Simon Young <yangxh@biomarker.com.cn>
#Last modified  2015-05-14 
#Modifier        baij <baij@biomarker.com.cn>
##Last modified  2017-02-13 

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
GetOptions(\%opts,"id=s","key=s","cpu=s","worksh=s","interval=s","sh_num=s","h","queue=s");
if (!defined($opts{id}) || !defined($opts{worksh}) || defined($opts{h})) {
	&help;
	exit;
}

#######################

&timeLog("$Script Start.");

######################
my $notename=`hostname`;chomp $notename;
my $idir=&ABSOLUTE_DIR($opts{id});
$idir=~s/\/$//;

&MKDIR($opts{worksh});
my $worksh=&ABSOLUTE_DIR($opts{worksh});
my $qsub_dir="$worksh/Phase2Trinty_sh";
&MKDIR($qsub_dir);
my $cpu=$opts{cpu} || 60;

my $interval=$opts{interval} || 5;
my $sh_num=$opts{sh_num} || 1000;
my $key=$opts{key};

#------------------------------ by Simon Young 2015-05-14 --------------------------------
################################by songmm  2016-1-27
# open and close file handle only once when needed 
# to avoid endless write-in and optimize operation
&random("$idir/chrysalis/Phase2TrinityCommands") if(!-e "$idir/chrysalis/Phase2TrinityCommands.random");
&generate_sh("$idir/chrysalis/Phase2TrinityCommands.random",$qsub_dir,"Phase2Trinity");
#`rm $idir/chrysalis/quantifyGraph_commands.tmp`;

my @sh_file=glob "$qsub_dir/Phase2Trinity_*.sh";

foreach my $sh (@sh_file) {
	&qsubOrDie($sh,$opts{queue},$cpu,"20G");
	#system"sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $opts{queue} --resource vf=30.0G --maxproc $cpu --reqsub --independent --interva $interval $sh ";
}
my @sh_num=glob("$qsub_dir/*sh");
my @qsub_num=glob("$qsub_dir/*qsub");
my $cycle=0;
while($cycle<10){
	if($#sh_num != $#qsub_num){
		foreach my $sh (@sh_file) {
			&Shell_qsub($sh,$opts{queue},$cpu,"20G");
		}
	}
	$cycle++;
}
###############Time
&timeLog("$Script End.");
&Runtime($BEGIN);

###########subs
sub random
{
	my ($quantify_cmd)=@_;
	open (CMD1,"grep -v '#' $quantify_cmd |sed '/^\$/ d' |") or die $!;
	open OUT1,">$quantify_cmd.random";
	my @cmd1 = <CMD1>;
	my %tag;
	for ( 1 .. @cmd1 ) {
		my $num = 1 + int (rand (@cmd1));
		redo if $tag{$num}++;
		print OUT1 "$cmd1[$num-1]";
	}
	close OUT1;
	close CMD1;
}

sub generate_sh
{
	my ($cmd_tmp,$qsub_dir,$cmd_flag)=@_;
	my @commands;
	open CMD,"$cmd_tmp";
	while (<CMD>) {
		chomp;
		$commands[int(($.-1)/(10*$sh_num))][int((($.-1)%(10*$sh_num))/10)].= "$_ && ";
	}
	close CMD;
	my $decimals=int(log(@commands)/log(10))+1;
	for (my $i=0; $i<@commands; $i++) {
		my $index=sprintf "%0${decimals}d",($i+1);
		open (SPL, ">$qsub_dir/${cmd_flag}_$index.sh") or die;
		print SPL "hostname&&";
		for my $l (@{$commands[$i]}) {
			print SPL "$l\n";
		}
		close SPL;
	}
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

sub help
{
	print <<"	Usage End.";
	Description:
		Function : use Trinity Assembly Second Step quantifyGraph qsub proccess;
		Version  : $ver
				Update Extract The Failed quantifyGraph Commands To failed.sh, for Reqsub; 2011-12-13
		Writer   : mengf <mengf\@biomarker.com.cn>
		Usage    :
		-id
		    Trinity Denovo dir include the chrysalis dir;
		-worksh
		    shell line out dir;
		-key
		    prefix of OUT shell (Munro, Species Name);
		-cpu
		    The cpu for qsub;
		-interval
		    The Interval Time to qsub missions, default 5 second;
		-sh_num
		    per shell contains The line number of Commands, default 1000;
	Usage End.

	exit;
}
