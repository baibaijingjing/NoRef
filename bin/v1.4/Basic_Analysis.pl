#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.1";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($basic_config,$data_config,$i,$qphred,$od,$min_kmer_cov,$step);
GetOptions(
				"help|?" =>\&USAGE,
				"basic_config:s"=>\$basic_config,
				"data_config:s"=>\$data_config,
				"i:s"=>\$i,
                                "qphred:s"=>\$qphred,
				"od:s"=>\$od,
                                "min_kmer_cov:i"=>\$min_kmer_cov,
				"s:s"=>\$step,
				) or &USAGE;
&USAGE unless ($basic_config and $data_config and $i and $od);
my $notename=`hostname`;chomp $notename;

###############Time
&timeLog("$Script start.");

mkdir $od unless -d $od;
$basic_config=&ABSOLUTE_DIR($basic_config);
$data_config=&ABSOLUTE_DIR($data_config);
$od=&ABSOLUTE_DIR($od);
my $queue;
$step||=1;

my $combine;
open (IN,$basic_config) or die $!;
while (<IN>) {
	chomp;
	if (/^Combine/) {
		$combine=(split/\s+/,$_)[1];
	}
	elsif(/^para_queue/){
		$queue=(split/\s+/,$_)[1];
	}
}
close IN;

$combine ||= 0;

if ($combine==0) {
	if ($step==1) {
		&timeLog("Assembly Process:");
                if(defined $min_kmer_cov){
		  &timeLog("perl $Bin/bin/Assembly/trinity_assembly_process_v1.0.pl -assembly_config $basic_config -data_config $data_config -min_kmer_cov $min_kmer_cov -od $od/Assembly");
		  system "perl $Bin/bin/Assembly/trinity_assembly_process_v1.0.pl -assembly_config $basic_config -data_config $data_config -min_kmer_cov $min_kmer_cov -od $od/Assembly";
                  }else{
                  &timeLog("perl $Bin/bin/Assembly/trinity_assembly_process_v1.0.pl -assembly_config $basic_config -data_config $data_config -od $od/Assembly");
                  system "perl $Bin/bin/Assembly/trinity_assembly_process_v1.0.pl -assembly_config $basic_config -data_config $data_config -od $od/Assembly";
                  }
                &timeLog("Assembly is done.");
		$step=2;
	}

	if ($step==2) {
		my $list=join ",",(glob "$od/Assembly/Trinity_assembly/*/*_Result/Unigenes/*.Unigenes.fa");
		&timeLog("Cluster Process:");
		&timeLog("perl $Bin/bin/cluster2Unigene/cluster2unigene.pl -cfg $basic_config -i $i -l $list -od $od/Cluster");
		system"perl $Bin/bin/cluster2Unigene/cluster2unigene.pl -cfg $basic_config -i $i -l $list -od $od/Cluster";
#		`perl $Bin/bin/draw_length_distribution/draw_length_distribution.pl -i $od -od $od/Cluster/Unigene/";
        &timeLog("Cluster is done.");
		$step=3;
	}

	if ($step==3) {
		&timeLog("Remap Process:");
		&timeLog("perl $Bin/bin/Remap/remap_expression.pl -queue $queue -cfg $data_config -id $od/Cluster/Unigene -od $od/Remap");
		system"perl $Bin/bin/Remap/remap_expression.pl -queue $queue -cfg $data_config -id $od/Cluster/Unigene -od $od/Remap";
        &timeLog("Remap is done.");
	}
}

if ($combine==1) {
	if ($step==1) {
		&timeLog("Assembly Process:");
                if(defined $min_kmer_cov){
		&timeLog("perl $Bin/bin/RSEM_Assembly/trinity_assembly_process_v1.0.pl -assembly_config $basic_config -data_config $data_config -min_kmer_cov $min_kmer_cov -od $od/Assembly");
		system"perl $Bin/bin/RSEM_Assembly/trinity_assembly_process_v1.0.pl -assembly_config $basic_config -data_config $data_config -min_kmer_cov $min_kmer_cov -od $od/Assembly";
                }else{
                &timeLog("perl $Bin/bin/RSEM_Assembly/trinity_assembly_process_v1.0.pl -assembly_config $basic_config -data_config $data_config -od $od/Assembly");
                system"perl $Bin/bin/RSEM_Assembly/trinity_assembly_process_v1.0.pl -assembly_config $basic_config -data_config $data_config -od $od/Assembly";
                }
        &timeLog("Assembly is done.");
        $step=2;
	}

	if ($step==2) {
		my $list=join ",",(glob "$od/Assembly/Trinity_assembly/*/*_Result/Unigenes/*.Unigenes.fa");
		mkdir "$od/Cluster" unless -d "$od/Cluster";
		mkdir "$od/Cluster/Unigene" unless -d "$od/Cluster/Unigene";
		&timeLog("Cluster Process:");
        &timeLog("perl $Bin/bin/cluster2Unigene/bin/assembly_stat.pl -i $i -fa $od/Assembly/Trinity_assembly/All_Combination/All_Combination_Result/Unigenes/All_Combination.Unigenes.fa -s 300,500,1000,2000 -len 30000 -od $od/Cluster/Unigene ");
		system"perl $Bin/bin/cluster2Unigene/bin/assembly_stat.pl -i $i -fa $od/Assembly/Trinity_assembly/All_Combination/All_Combination_Result/Unigenes/All_Combination.Unigenes.fa -s 300,500,1000,2000 -len 30000 -od $od/Cluster/Unigene";
        &timeLog("perl $Bin/bin/cluster2Unigene/bin/Makegff.pl -i $od/Cluster/Unigene/$i.Unigene.fa -o $od/Cluster/Unigene/$i.Unigene.gff ");
		system"perl $Bin/bin/cluster2Unigene/bin/Makegff.pl -i $od/Cluster/Unigene/$i.Unigene.fa -o $od/Cluster/Unigene/$i.Unigene.gff";
        &timeLog("perl $Bin/bin/cluster2Unigene/bin/Getorf_Extract.pl -fa $od/Cluster/Unigene/$i.Unigene.fa -od $od/Cluster/Unigene/Unigene_Orf ");
		system"perl $Bin/bin/cluster2Unigene/bin/Getorf_Extract.pl -fa $od/Cluster/Unigene/$i.Unigene.fa -od $od/Cluster/Unigene/Unigene_Orf";
#		`perl $Bin/bin/cluster2Unigene/bin/assembly_stat.pl -i $i.orf -fa $od/Cluster/Unigene/Unigene_Orf/$i.Unigene.cds.fa -s 300,500,1000,2000 -od $od/Cluster/Unigene/Unigene_Orf/Graph`; #ORF len stat
        &timeLog("perl $Bin/bin/cluster2Unigene/bin/SSR_Analysis.pl -fa $od/Cluster/Unigene/$i.Unigene.fa -od $od/Cluster/Unigene/Unigene_SSR ");
        system"perl $Bin/bin/cluster2Unigene/bin/SSR_Analysis.pl -fa $od/Cluster/Unigene/$i.Unigene.fa -od $od/Cluster/Unigene/Unigene_SSR";
#		`perl $Bin/bin/draw_length_distribution/draw_length_distribution.pl -i $od -od $od/Cluster/Unigene/ -n $combine`;
        &timeLog("perl $Bin/bin/ppi_network/ppi_network.pl --qseq $od/Cluster/Unigene/$i.Unigene.fa --odir $od/Cluster/Unigene/Unigene_PPI ");
        system"perl $Bin/bin/ppi_network/ppi_network.pl --qseq $od/Cluster/Unigene/$i.Unigene.fa --odir $od/Cluster/Unigene/Unigene_PPI";
        &timeLog("Cluster and ppi analysis is done.");
		$step=3;
	}

	if ($step==3) {
		my %Sample;
		open (IN,$data_config) or die $!;
		open (OUT,">$od/remap_data.cfg") or die $!;
		while (<IN>) {
			if (/^Sample/) {
				print OUT $_;
				my $sam=(split/\s+/,$_)[1];
				<IN>;<IN>;
				my $line1=<IN>;
				my $fq1file=(split/\s+/,$line1)[1];
				$Sample{$sam}=(split/\s/,`wc -l $fq1file`)[0]/4;
				print OUT $line1;
				print OUT scalar <IN>;
			}
		}
		close IN;
		close OUT;
		&timeLog("Remap Process:");
		&timeLog("perl $Bin/bin/RSEM_Remap/v1.0/remap_expression.pl -queue $queue  -cfg $od/remap_data.cfg -qphred $qphred -i $od/Assembly/Trinity_assembly/All_Combination/All_Combination_Result/Transcripts/All_Combination.Transcripts.fa -od $od/Remap");
		system"perl $Bin/bin/RSEM_Remap/v1.0/remap_expression.pl -queue $queue  -cfg $od/remap_data.cfg -qphred $qphred -i $od/Assembly/Trinity_assembly/All_Combination/All_Combination_Result/Transcripts/All_Combination.Transcripts.fa -od $od/Remap";
		open SH,">$od/Remap/work_sh/detail.sh" || die "$!";
		foreach my $sam (keys %Sample) {
			print SH "perl $Bin/bin/bam_file_stat/v1.0/bam_file_stat.pl -t $Sample{$sam} -i $od/Remap/rmap/$sam/$sam.bowtie.bam -o $od/Remap/rmap/$sam/$sam.Mapped.stat.xls && ";
			print SH "perl $Bin/bin/RSEM_Remap/v1.0/bin/draw_insertlen_by_bam.pl -i $od/Remap/rmap/$sam/$sam.bowtie.bam -index $sam -od $od/Remap/rmap/$sam && ";
			print SH "perl $Bin/bin/RSEM_Remap/v1.0/bin/random_graph.pl -i $od/Remap/rmap/$sam/$sam.bowtie.bam -exp $od/Remap/geneExpression/$sam.geneExpression.xls -orf $od/Cluster/Unigene/Unigene_Orf/$i.Unigene.cds_pep.stat.xls -o $od/Remap/rmap/$sam/$sam && ";
			print SH "perl $Bin/bin/RSEM_Remap/v1.0/bin/gene_tagnum_graph.pl -i1 $Sample{$sam} -i2 $od/Remap/geneExpression/$sam.geneExpression.xls -value 0.1 -o $od/Remap/rmap/$sam/$sam \n";
		}
		#&runOrDie("$od/Remap/work_sh/detail.sh");
		&qsubOrDie("$od/Remap/work_sh/detail.sh",$queue,20,"10G");
		#my $cmd_res=system"sh $od/Remap/work_sh/detail.sh >$od/Remap/work_sh/detail.sh.log 2>&1";
		close SH;
		#&Shell_qsub ("$od/Remap/work_sh/detail.sh","general.q",20);
        &timeLog("Cluster is done.");
	}
	system"cp $od/Remap/rmap/*/*.png $od/Remap/geneExpression ";
	system"cp $od/Remap/rmap/*/*.isoforms.results $od/Remap/geneExpression ";
	system"cp $od/Remap/rmap/*/*.Mapped.stat.xls $od/Remap/geneExpression ";


}
&timeLog("perl $Bin/bin/draw_gene_tag/draw_gene_tag.pl -id $od/Remap/rmap -od $od/Remap/geneExpression");
system"perl $Bin/bin/draw_gene_tag/draw_gene_tag.pl -id $od/Remap/rmap -od $od/Remap/geneExpression";
&timeLog("perl $Bin/bin/draw_total_random/draw_total_random.pl -id $od/Remap/rmap -od $od/Remap/geneExpression");
system"perl $Bin/bin/draw_total_random/draw_total_random.pl -id $od/Remap/rmap -od $od/Remap/geneExpression";

#######################################################################################
###############Time
&timeLog("$Script end, elapsed time:".(time()-$BEGIN_TIME)."s.");
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		die "Warning just for file and dir $in\n";
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################

sub max{#&max(lists or arry);
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#锟斤拷取锟街凤拷锟斤拷锟斤拷锟叫的凤拷锟津互诧拷锟斤拷锟叫ｏ拷锟斤拷锟街凤拷锟斤拷锟斤拷式锟斤拷锟截★拷ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

sub Shell_qsub
{ # &Shell_qsub($sh,$qeue,$cpu);
	my $sh = shift;
	my $qeue = shift;
	my $cpu = shift;
	
	#system"sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --queue $qeue --reqsub --maxproc $cpu --independent $sh  ";
}

################################################################################################################
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Zhang XueChuan <zhangxc\@biomarker.com.cn> 
Usage:
  Options:
  -basic_config    <file>  basic analysis config file,forced 

  -data_config     <file>  data file,forced 

  -i               <str>   Index of clustered unigene

  -od              <dir>   output dir,forced 

  -s               <num>   default  1
                           1  Assembly
                           2  Cluster
                           3  Remap

  -h         Help

USAGE
	print $usage;
	exit;
}
