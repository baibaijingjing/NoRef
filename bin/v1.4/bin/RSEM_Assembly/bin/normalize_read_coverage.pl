#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
use Cwd qw(abs_path);
use threads;


my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;

my %config = ();
#my $trinity_HOME = "/share/nas2/genome/biosoft/trinity/trinityrnaseq_r2013-02-25"; # 2014-12-25 ~ 
my %sys_cfg = %{&readconf("$Bin/../../../../../../config/sys.cfg")};
my $trinity_HOME  =$sys_cfg{trinity};
my $UTIL_DIR = "$trinity_HOME/util";

my $fIn;
my $fKey;
my $output_directory = "normalized_reads";

my $max_cov;
my $KMER_SIZE = 25;
my $SS_lib_type;
my $max_memory;
my $max_pct_stdev = 100;
my $data_s;
my $normalize_cpu;
GetOptions(
				"help|?" =>\&USAGE,
				
				"i:s"=>\$fIn,
				"k:s"=>\$fKey,
				"d:s"=>\$output_directory,
				"datas:s"=>\$data_s,
				# Jellyfish
				'JM=s'=> \$max_memory, # in GB

				# misc
				'max_pct_stdev=i' => \$max_pct_stdev,
				"max_cov=i" => \$max_cov,
				'kmer_size=i' => \$KMER_SIZE,
				"SS_lib_type=s" => \$SS_lib_type,
				"nlz_cpu=s"=>\$normalize_cpu,
				
				) or &USAGE;
&USAGE unless ($fIn and $fKey);

#
# check options 
#
if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(R|F|RF|FR)$/) {
        die "Error, unrecognized SS_lib_type value of $SS_lib_type. Should be: F, R, RF, or FR\n";
    }
}

unless ($max_cov && $max_cov >= 2) {
    die "Error, need to set --max_cov at least 2";
}

unless ($max_memory) {
    die "Error, must specify max memory for jellyfish to use, eg.  --JM 10G \n";
}
########cancel the set 2016-08-19
#else{
	## check memory, if too large, set to a small value 

#	$max_memory =~s/[MG]$//g;
#	$max_memory = 60 if ($max_memory > 60) ;
#
#	$max_memory .= "G";
#}

main: {

	$fIn = Cwd::abs_path($fIn);

	unless (-d $output_directory) {
		mkdir $output_directory or die "Error, cannot mkdir $output_directory";
	}

	$output_directory = Cwd::abs_path($output_directory);

#	chdir ($output_directory) or die "Error, cannot cd to $output_directory";


	#
	# step 1.parse original data config and generate list file
	#
	parse_data_config_file($fIn, \%config);

	#print Dumper %config;die;

	my ($fq1_list, $fq2_list) = make_list_file(\%config, "$output_directory/$fKey");


	#
	# statistic data size and choose the appropriate grid 
	#
#	my $data_s = 0;
#	foreach my $sample (keys %config) {
#		$data_s += (-s $config{$sample}{'fq1'}) + (-s $config{$sample}{'fq2'});
#	}
#	$data_s = int($data_s/(10**9)/3);

#	my $grid_type = 0;
#	$grid_type = 1 if ($data_s > 64);
#	$grid_type = 2 if ($data_s > 32 && $data_s <=64);
#	$grid_type = 3 if ($data_s > 4 && $data_s <= 32);
#	$grid_type = 3 if ($data_s <= 4);
#-------------------- by Simon Young 2015-01-05 --------------------------------
# #memory   base    state   time
# 188.9G    101G    DONE    
# 377.9G    213G    DONE    ~2d
# 757.4G    254G    DONE    
# 757.4G    622G    ABORT   
#	$grid_type = 1 if ($data_s > 250);
#	$grid_type = 2 if ($data_s > 60 && $data_s <= 250);
#	$grid_type = 3 if ($data_s <= 60);
#	$grid_type = 4 if  ($data_s <= 12);

#	#
#	# load grid info 
#	#
#	my %middle_grid;
#	my %great_grid;
#	my %huge_grid;
#	my %max_grid;
	my %cpu;
	&Grid_info (\%cpu);

	#
	# step 2: calling normalize_by_kmer_coverage.pl to downsampling reads 
	#
	if (-f "$output_directory/normalized.check") {
		#`rm -f $output_directory/normalized.check`;
		print "normalized has finished\n";
	}
	else{
	my $shell = "$output_directory/normalize_by_kmer_coverage.sh";
	open (SH,">$shell") or die $!;
	#my $cmd = "perl $UTIL_DIR/normalize_by_kmer_coverage.pl --seqType fq --left_list $fq1_list --right_list $fq2_list --JM $max_memory --max_cov $max_cov --KMER_SIZE $KMER_SIZE --max_pct_stdev $max_pct_stdev --output $output_directory --min_kmer_cov 2 --PARALLEL_STATS --pairs_together";
	my $cmd = "perl $UTIL_DIR/normalize_by_kmer_coverage.pl --seqType fq --left_list $fq1_list --right_list $fq2_list --JM $max_memory --max_cov $max_cov --KMER_SIZE $KMER_SIZE --max_pct_stdev $max_pct_stdev --output $output_directory  --PARALLEL_STATS --pairs_together";
	$cmd .= " --SS_lib_type $SS_lib_type " if ($SS_lib_type) ;

	$cmd .= " && touch $output_directory/normalized.check";

	print SH $cmd;
	close (SH) ;
	
	## distribute job to the grid
	my $status = "Null";
	$data_s||=$max_memory;
	$data_s=~s/G$//;
	$data_s=$data_s/2;
	while ($status eq "Null") {
		
		$status = &grid_run_cmd (\%cpu,$shell,$data_s);
	
		if ($status eq "Null") {
			sleep (1800);
		}
		else {
			my $cmd_res=system"sh $status";
			if($cmd_res){
				die"Error:normalize failed\n";
			}
		}
	}
	}
	# step 3: generate new data config file 
	&make_new_data_config_file("$output_directory/$fKey");
	print STDOUT "\nNormalize Done. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
}
#
# subs
#

# --------------------------- by Simon Young 2014-12-15 -------------------------------------
sub Grid_info {
    my $compute = shift;
	my (%name,$host_cpu);
	for my $line(split/\n/,`qhost`){
		my @tmp=split/\s+/,$line;
		if($line=~/HOSTNAME/){
			for(my $i=0;$i<@tmp;$i++){
				$name{$tmp[$i]}=$i;
			}
			next;
		}
		next if ($line=~/------------/ or $tmp[$name{LOAD}] =~/-/);
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
	my $queue_g = shift;
	my $sh = shift;
	my $data_s = shift;
	my @run_grid = glob "$sys_cfg{trinity_dir}/*";
	#my @run_grid = glob "/home/yangxh/trinity_run_grid/*";
    
    if($normalize_cpu){
        open OUT,">$sh.run.sh" || die $!;
        my $headnode=$sys_cfg{headnode};
		print OUT "ssh $headnode 2>/dev/null <<TT\n";
		print OUT "ssh $normalize_cpu 2>/dev/null <<EOF\n";
		print OUT " sh $sh >$sh.run.sh.log 2>&1\n";
		print OUT "EOF\nTT\n";
		close OUT;
		return ("$sh.run.sh");
		last;
    }
    
	my %use_grid;
	foreach my $grid (@run_grid) {
		my $grid_node = basename $grid;
		$use_grid{$grid_node} = 1;
	}

	my $stat = 0;
	foreach my $grid (sort {$$queue_g{$a}<=>$$queue_g{$b}} keys %{$queue_g}) {
		next if $$queue_g{$grid}<$data_s/2;          ####set normalize read compute by method of script normalize_by_kmer_coverage.pl
		if (!exists $use_grid{$grid}) {
			$stat = 1;
			#system"touch $sys_cfg{trinity_dir}/$grid ";
			#system"touch /home/yangxh/trinity_run_grid/$grid ";
			open OUT,">$sh.run.sh" || die $!;
			my $headnode=$sys_cfg{headnode};
			print OUT "ssh $headnode 2>/dev/null <<TT\n";
			print OUT "ssh $grid 2>/dev/null <<EOF\n";
			print OUT " touch $sys_cfg{trinity_dir}/$grid &&sh $sh >$sh.run.sh.log 2>&1&& rm $sys_cfg{trinity_dir}/$grid \n";
			#print OUT " sh $sh >$sh.run.sh.log 2>&1&& rm /home/yangxh/trinity_run_grid/$grid \n";
			print OUT "EOF\nTT\n";
			close OUT;
			return ("$sh.run.sh");
			last;
		}
	}
	return ("Null") if ($stat == 0);
}

sub make_new_data_config_file {#
	my ($file_prefix) = @_;

	my $normalized_file_suffix = "normalized_K${KMER_SIZE}_C${max_cov}_pctSD${max_pct_stdev}.fq";

	open (OUT,">$output_directory/$fKey.normalized.data.config") or die $!;
	
	print OUT "Sample\tAll\n";
	print OUT "fq1\t", "$file_prefix.left.fq.list.$normalized_file_suffix", "\n";
	print OUT "fq2\t", "$file_prefix.right.fq.list.$normalized_file_suffix", "\n";
	print OUT "\n";
	
	close (OUT) ;
}

sub make_list_file {#
	my ($config_aref, $file_prefix) = @_;

	my $left_file_list = "$file_prefix.left.fq.list";
	my $right_file_list = "$file_prefix.right.fq.list";

	open (LFQ,">$left_file_list") or die $!;
	open (RFQ,">$right_file_list") or die $!;

	foreach my $sample_id (sort {$config_aref->{$a}{'order'} <=> $config_aref->{$b}{'order'}} keys %{$config_aref}) {
		print LFQ "$config_aref->{$sample_id}{'fq1'}", "\n";
		print RFQ "$config_aref->{$sample_id}{'fq2'}", "\n";
#		print LFQ Cwd::abs_path("$config_aref->{$sample_id}{'fq1'}"), "\n";
#		print RFQ Cwd::abs_path("$config_aref->{$sample_id}{'fq2'}"), "\n";
	}
	
	close (RFQ) ;
	close (LFQ) ;

	return ($left_file_list, $right_file_list);
}

sub parse_data_config_file {#
	my ($infile, $config_aref) = @_;

	my $sample_id = '';
	my $sample_order = 0;

	open (IN,"$infile") or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ or /^\#/) ;
		
		if (/^Sample/) {
			($sample_id) = $_ =~/\S+\s+(\S+)/;
			$config_aref->{$sample_id}{'order'} = $sample_order++;
		}

		($config_aref->{$sample_id}{'fq1'}) = $_ =~/^fq1\s+(\S+)/ if (/^fq1/) ;
		($config_aref->{$sample_id}{'fq2'}) = $_ =~/^fq2\s+(\S+)/ if (/^fq2/) ;
		
	}
	close (IN) ;

	return;
}


sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";

    my $start_time = time();
    my $ret = system($cmd);
    my $end_time = time();

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    print "CMD finished (" . ($end_time - $start_time) . " seconds)\n";    

    return;
}

sub check_memory {#
	my %compute_node = ();

	my $qhost_str = `qhost |grep ^co `;chomp $qhost_str;
	my (%name,@line);
	for my $line(split/\n/,`qhost`){
		my @tmp=split/\s+/,$line;
		if($line=~/HOSTNAME/){
			for(my $i=0;$i<@tmp;$i++){
				$name{$tmp[$i]}=$i;
			}
			next;
		}
		next if ($tmp[$name{LOAD}] =~/-/ or $line=~/------------/);
		push @line,"$tmp[$name{HOSTNAME}]\t$tmp[$name{MEMTOT}]\t$tmp[$name{MEMUSE}]"
	}
	foreach my $line (@line) {
		my ($compute_node,$all_memory,$used_memory) = split /\s+/,$line;

		next if ($used_memory eq '-') ;

		$all_memory=~s/G//g;	
		next if ($all_memory < 70) ;

		$used_memory=~s/G//g;
		if ($used_memory =~/M$/) {
			$used_memory=~s/M//g;
			$used_memory /=1024;
		}

		$compute_node{$compute_node} = $used_memory / $all_memory;	
#		$compute_node{$compute_node} = $all_memory - $used_memory;	
	}

	my ($target_node) = sort {$compute_node{$a} <=> $compute_node{$b}} keys  %compute_node;
#	my ($target_node) = sort {$compute_node{$b} <=> $compute_node{$a}} keys  %compute_node;

	return $target_node;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription:
	
	format of data config is like this��
	
	Sample     Demo_T1
	fq1        /share/nas5/RNA_seq/Project/SOP_Test/New_Assembly_analysis/raw_Data/Demo_T1_1.fq
	fq2        /share/nas5/RNA_seq/Project/SOP_Test/New_Assembly_analysis/raw_Data/Demo_T1_2.fq

	Sample     Demo_T2
	fq1        /share/nas5/RNA_seq/Project/SOP_Test/New_Assembly_analysis/raw_Data/Demo_T2_1.fq
	fq2        /share/nas5/RNA_seq/Project/SOP_Test/New_Assembly_analysis/raw_Data/Demo_T2_2.fq

Usage:
  Options:
  -i              <file>    Fastq data config file, forced
  -k              <str>     Key of output file, forced
  -d              <str>     Directory where output file produced,optional,default [normalized_reads]
  -datas          <str>     the data size,(eg 10G)  include the "G" chr
  -JM             <str>     Jellyfish Memory, number of GB of system memory to use for 
                             k-mer counting by jellyfish  (eg. 10G) *include the 'G' char 
  -max_cov        <int>     Targeted maximum coverage for reads.
  -SS_lib_type    <str>     Strand-specific RNA-Seq read orientation.
                              if paired: RF or FR,
                              if single: F or R.   (dUTP method = RF)
  
  -kmer_size      <int>     Kmer size, default 25
  -max_pct_stdev  <int>     Maximum pct of mean for stdev of kmer coverage across read (default: 100)
  -nlz_cpu        <str>     the cpu name to run normalize step
  -h                        Help

USAGE
	print $usage;
	exit;
}

