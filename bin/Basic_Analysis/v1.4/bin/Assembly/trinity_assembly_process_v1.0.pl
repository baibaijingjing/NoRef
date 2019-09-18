#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $version="1.0.0";
my $BEGIN=time();
use threads;
use Thread::Semaphore;
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
my ($assembly_config,$data_config,$min_kmer_cov,$step,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"assembly_config:s"=>\$assembly_config,
				"data_config:s"=>\$data_config,
				"s:n"=>\$step,
                                "min_kmer_cov:i"=>\$min_kmer_cov,
				"od=s"=>\$od,
				) or &USAGE;
&USAGE unless ($assembly_config and $data_config and $od) ;

###############Time start
&timeLog("$programe_dir start.");

############### Mkdir
my $notename=`hostname`;chomp $notename;
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);
&MKDIR("$od/work_sh");

$step = $step || 1;
################ Load grid info
#my %middle_grid;
#my %great_grid;
#my %huge_grid;
#my %max_grid;
#&Grid_info (\%middle_grid,\%great_grid,\%huge_grid,\%max_grid);

############### Load Sample info
my %Samples;
&parse_sample_config ($data_config,\%Samples);

my %Assembly_Combination;
my %Assembly_Para;
&LOAD_PARA ($assembly_config,\%Assembly_Combination,\%Assembly_Para);

#######################################################
# Prepare Reads & Log file & Config files for Assembly
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
my %Data_log;


#########################################################
# Step1 : prepare cat reads shell & configs for Assembly
#         based on the assembly aombination
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
if ($step == 1) {
	&stepStart(1,"Prepare Config files for Assembly.");
	&MKDIR ("$od/Config");
	if ($Assembly_Combination{Separate}==1) {

		foreach my $sam (sort keys %Samples) {
			push @{$Data_log{$Assembly_Combination{Sep_Data}}}, $sam;
			my $config = "$od/Config/$sam.config";
			my $conf;
			&creat_trinity_config($config,$sam,$Samples{$sam}{FQ1},$Samples{$sam}{FQ2},\%Assembly_Para,$conf);
		}
	}
	elsif ($Assembly_Combination{SamG}==1) {
		
		foreach my $com (sort keys %{$Assembly_Combination{Samples}}) {
			my $data_size = "$com"."_Data";
			my @sams = split /,/, "$Assembly_Combination{Samples}{$com}";
			my ($fq1_list,$fq2_list);
			foreach my $sam (sort @sams) {
				$fq1_list .= "$Samples{$sam}{FQ1} ";
				$fq2_list .= "$Samples{$sam}{FQ2} ";
			}
			push @{$Data_log{$Assembly_Combination{$data_size}}}, $com;

			my $config = "$od/Config/$com.config";
			my $conf;
			&creat_trinity_config ($config,$com,$fq1_list,$fq2_list,\%Assembly_Para,$conf);
		}
	}

	&stepTime(1);
	$step = 2;
}

#######################################
# Step 2: Creat trinity assembly Shell
#         differ type grid shell creat
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
if ($step == 2) {
	&MKDIR ("$od/Trinity_assembly");
	&MKDIR ("$od/Trinity_sh");
	&stepStart(2,"creat Assembly Shell.");
	foreach my $data_s (sort {$a<=>$b} keys %Data_log) {
#		my $grid_type = 0;
#		$grid_type = 1 if ($data_s > 64);
#		$grid_type = 2 if ($data_s > 32 && $data_s <=64);
##		$grid_type = 3 if ($data_s > 4 && $data_s <= 32);
##		$grid_type = 4 if ($data_s <= 4);
## small grids (limit up to 31.5G) are no longer used.
## --------------------------- by Simon Young 2015-04-08 -------------------------------------
#        $grid_type = 3 if ($data_s <= 32);
#		$grid_type =4 if  ($data_s <= 12);
		foreach my $config_pre (@{$Data_log{$data_s}}) {
			open SH,">$od/Trinity_sh/$config_pre.trinity.$data_s.sh" || die $!;
			my $cmd="mkdir -p $od/Trinity_assembly/$config_pre && sh $Bin/export.sh && perl $Bin/bin/trinity_process.pl -cfg $od/Config/$config_pre.config -od $od/Trinity_assembly/$config_pre  ";
                       $cmd.="-min_kmer_cov $min_kmer_cov" if(defined $min_kmer_cov);
                       $cmd.="\n";
                       print SH $cmd;
		       close SH;
		}
	}
	&stepTime(2);
	$step = 3;
}

#################################################
# Step 3: run assembly step (run on differ grid)
#         search suitable grid
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
if ($step == 3) {
	&stepStart(3,"search compute grid for Trinity assembly.");

	my @shells = glob "$od/Trinity_sh/*trinity.*sh";
	my %cpu;
	&Grid_info(\%cpu);
	my $semaphore=Thread::Semaphore->new(5);
	foreach my $shell (@shells) {
		next if $shell =~/run.sh/;
		$shell =~ /.*\.(\d+)\.sh/;
		my $key = $1;

		my $status = "Null";
		my $num=0;
		while ($status eq "Null") {
			
			$status = &grid_run_cmd (\%cpu,$shell,$key);
			$num++;
			print "search compute num:$num, status:$status\n";
			if ($status eq "Null") {
				sleep (900);
			}
			else {
				
				threads->create(sub{
					$semaphore->down();
					print "$status start\n";
					system " sh $status ";
					$semaphore->up();
				});
				for my $thread (threads->list(threads::joinable)){
					$thread->join();
				}
			}
		}
	}
	for my $thread (threads->list(threads::all)){
		$thread->join();
	}
	$step = 4;
}


########################################################
# Step 4: Check assembly result (twice Interval 1 hour)
#          all assembly either runing or done!
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
if ($step == 4) {
	&stepStart(4,"Check Trinity assembly");
	my %unigene_list;
	my $done = 0;
	my $fasta_num = 0;
	while ($done == 0) {
		my $sam_num = 0;
		my $finish_num = 0;
		foreach my $num (sort keys %Data_log) {
			foreach my $combine (sort @{$Data_log{$num}}) {
				$sam_num++;
				my $fasta = "$od/Trinity_assembly/$combine/$combine"."_Result/Unigenes/$combine.Unigenes.fa";
				$unigene_list{$fasta} = 1;
				if (-f "$fasta") {
					$finish_num++;
				}
			}
		}
		$fasta_num = $sam_num;
		if ($finish_num == $fasta_num) {
			sleep (100);
			$done = 1;
			&stepTime(3);
		}
		else {
			sleep (900);
		}
	}

	open OUT,">$od/unigene.list" || die $!;
	foreach (keys %unigene_list) {
		print OUT "$_\n";
	}
	close OUT;

	&stepTime(4);
}

###############Time end
&totalTime();


####################subs
sub parse_sample_config
{ #load config file
	my $config_file= shift;
	my $SamData= shift;
	
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/\s+$//;s/\r$//;
		next if(/^$/ or /^\#/);
		if ( /^Sample/ ) {
			my $fq1_info = <IN> ; chomp $fq1_info;
			my $fq2_info = <IN> ; chomp $fq2_info;
			my $sample_name = (split /\s+/, $_)[1];
			my $fq1 = (split /\s+/, $fq1_info)[1];
			my $fq2 = (split /\s+/, $fq2_info)[1];

			if (! -e $fq1 || ! -e $fq2) {
				warn "Non-exist: Reads file dir for $sample_name isn't exists\n" ;
				exit;
			}
			$SamData->{$sample_name}->{FQ1} = $fq1;
			$SamData->{$sample_name}->{FQ2} = $fq2;
		}
	}
}

sub LOAD_PARA
{
	my $para_file= shift;
	my $combination = shift;
	my $para = shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	my $sam_group = 0;
	while (<IN>) {
		chomp;
		s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		if(/\w+_Data/){
			my $data_size = (split /\s+/,$_)[1];
			$data_size =~ s/G$//i;
			$combination->{$para_key} = $data_size;
		}
		elsif ($para_key =~ /^para/) {
			$para->{$para_key} = $para_value;
		}
		elsif ($para_key=~/^Group/){
			$combination->{Samples}{$para_key}=$para_value;
		}
		else{
			$combination->{$para_key} = $para_value;
		}
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
	close IN;
}

sub creat_trinity_config {
	my $config = shift;
	my $key = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $para = shift;
	my $conf = shift;

	$conf .= "Index\t$key\n";
	$conf .= "FQ1\t$fq1\n";
	$conf .= "FQ2\t$fq2\n\n";

	foreach my $para_meter (sort keys %{$para}) {
		$para_meter =~ /^para_(.*)/;
		$conf .= "$1\t$para->{$para_meter}\n";
	}
	open OUT,">$config" || die $!;
	print OUT "$conf";
	close OUT;
}

# --------------------------- by Simon Young 2014-12-05 -------------------------------------
sub Grid_info {
     my $compute=shift;
	my (%name,$host_cpu);
	for my $line(split/\n/,`qhost`){
		my @tmp=split/\s+/,$line;
		if($line=~/HOSTNAME/){
			for(my $i=0;$i<@tmp;$i++){
				$name{$tmp[$i]}=$i;
			}
			next;
		}
		next if ($line=~/------------/ or $tmp[$name{LOAD}] =~/-/ );
		$host_cpu.="$tmp[$name{HOSTNAME}]\t$tmp[$name{NCPU}]\t$tmp[$name{LOAD}]\t$tmp[$name{MEMTOT}]\t$tmp[$name{MEMUSE}]\n";
	}

    for my $line (split /\n/,$host_cpu) {
        my ($host, $cpu,$usedCpu,$ram,$usedRam) = split /\t/,$line;
        next if $ram=~/M$/;
        $ram =~s/G$//;
        $usedRam=~s/G$//;
        $usedRam=0 if $usedRam=~/M$/;
		my $load=$ram-$usedRam;
        $$compute{$host}=$load;
    }
}

sub grid_run_cmd {
	my %sys_cfg=%{&readconf("$Bin/../../../../../config/sys.cfg")};
	my $queue_g = shift;
	my $sh = shift;
	my $key = shift;
	my @run_grid = glob "$sys_cfg{trinity_dir}/*";
	#my @run_grid = glob "/home/yangxh/trinity_run_grid/*";

	my %use_grid;
	foreach my $grid (@run_grid) {
		my $grid_node = basename $grid;
		$use_grid{$grid_node} = 1;
	}

	my $stat = 0;
	foreach my $grid (sort {$$queue_g{$a}<=>$$queue_g{$b}} keys %{$queue_g}) {
		next if $$queue_g{$grid}<1.5*$key;
		if (!exists $use_grid{$grid}) {
			$stat = 1;
			system"touch $sys_cfg{trinity_dir}/$grid ";
			#system"touch /home/yangxh/trinity_run_grid/$grid ";
			open OUT,">$sh.run.sh" || die $!;
			my $headnode=$sys_cfg{headnode};
			print OUT "ssh $headnode 2>/dev/null <<TT\n";
			print OUT "ssh $grid 2>/dev/null <<EOF\n";
			print OUT " touch $sys_cfg{trinity_dir}/$grid\n";
			print OUT " sh $sh >$sh.run.sh.log 2>&1 \n";
			print OUT " rm $sys_cfg{trinity_dir}/$grid\n\n";
			print OUT "EOF\nTT\n";
			close OUT;
			return ("$sh.run.sh");
			last;
		}
	}
	return ("Null") if ($stat == 0);
}

sub Shell_qsub
{ # &Shell_qsub($sh,$qeue,$cpu);
	my $sh = shift;
	my $queue = shift;
	my $cpu = shift;
	#system"sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $qeue --reqsub --maxproc $cpu --independent $sh  ";
	&qsubOrDie($sh,$queue,$cpu,"5G");
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
		die "Warning just for file and dir $in\n";
	}
	
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: Pair End Reads Assembly With Trinity process
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:

Usage:

  -assembly_config                  <str>                  Index of Output                                 must be given;
  -data_config                      <str>                seqType for Assembly                            must be given;
  -od                               <str>                out dir                                         must be given;

  -s                        <int>                  step of process
                                                           1 load config & creat cat data shell          default;
                                                           2 cat reads prepare for assembly
                                                           3 creat assembly shell
                                                           4 search grid run trinity analysis
                                                           5 check final assembly result

USAGE
	print $usage;
	exit;
}
