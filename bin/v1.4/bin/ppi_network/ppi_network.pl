#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# utils and data
my $STRING_db = ${&readconf("$Bin/../../../../../config/sys.cfg")}{STRING_db};
my $abstractFabyId = "$Bin/util/abstractFabyId.pl";
my $parallel_blast = "$Bin/util/parallel_blast_on_split_query.pl";
my $Y2H_in_silico = "$Bin/util/Y2H_in_silico.pl";
my $cmd;

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($refseq,  $clade_or_taxid, $index, $odir,$queue);

GetOptions(
        "help|?" =>\&USAGE,
        "qseq=s"=>\$refseq,
        "clade:s"=>\$clade_or_taxid,
        "index:s"=>\$index,
        "odir=s"=>\$odir,
        "queue=s"=>\$queue,
        ) or &USAGE;
&USAGE unless ( $refseq and $odir);

system "mkdir -p $odir" unless (-d $odir);
$odir = &ABSOLUTE_DIR($odir);
my @refseq = split /,/,$refseq;
for (my $i=0; $i<@refseq; $i++) {
    if (-f "$refseq[$i]") {
        $refseq[$i] = &ABSOLUTE_DIR($refseq[$i]);
    } else {
        print "ERROR: file $refseq[$i] not exists.\n";
        &USAGE;
    }
}
$refseq = join ",",@refseq;

$queue ||= "general.q";
$index ||= "ppi_qurey";
$clade_or_taxid ||= 'eukaryota';

&timeLog("$Script start.");
# ------------------------------------------------------------------
# load STRING database 
# ------------------------------------------------------------------
my ($protein_links, $protein_actions, $protein_seq);

if ($clade_or_taxid eq 'all') {
    $protein_links = "$STRING_db/protein.links.v10.txt";
    $protein_actions = "$STRING_db/protein.actions.v10.txt";
    $protein_seq = "$STRING_db/protein.sequences.v10.fa";
} elsif ($clade_or_taxid =~/^\d+$/) {
    $protein_links = "$STRING_db/split/protein_links/eukaryota.$clade_or_taxid.protein.links.txt";
    $protein_actions = "$STRING_db/split/protein_actions/eukaryota.$clade_or_taxid.protein.actions.txt";
    $protein_seq = "$STRING_db/split/protein_seq/eukaryota.$clade_or_taxid.protein.sequences.fa";
} else {
    $protein_links = "$STRING_db/split/protein_links/$clade_or_taxid.protein.links.txt";
    $protein_actions = "$STRING_db/split/protein_actions/$clade_or_taxid.protein.actions.txt";
    $protein_seq = "$STRING_db/split/protein_seq/$clade_or_taxid.protein.sequences.fa";
}

unless (-f $protein_actions) {
    print "ERROR: Illegal arguments of --clade is given! \n";
    print "ERROR: the taxid $clade_or_taxid is illegal, maybe it isn't an eukaryota, or not exist actually.\n\n" if ($clade_or_taxid =~/^\d+$/);
    &USAGE;
}

# ------------------------------------------------------------------
# load query lists and sequence 
# ------------------------------------------------------------------
&timeLog("load query lists and abstract sequences: ");

system "grep \">\" $refseq|sed  's/\>//' > $odir/$index.id.list ";

$cmd = "perl $abstractFabyId -l $odir/$index.id.list -f $refseq -o $odir/$index.seq.fa >/dev/null 2>&1 ";
&runOrDie($cmd) unless (-s "$odir/$index.seq.fa");

# ------------------------------------------------------------------
# query ids and sequences abstract & blast with STRING
# ------------------------------------------------------------------
&timeLog("blast with STRING: ");

$cmd = "perl $parallel_blast --query $odir/$index.seq.fa --database $protein_seq --odir $odir/ --index $index --cutnum 100 --evalue 1e-5 -queue $queue";
&runOrDie($cmd) unless (-s "$odir/$index.blast.tab.best");

system "rm -f $index.blast";    # this file is too large 
# ------------------------------------------------------------------
# abstract associations of query DEG, and output
# ------------------------------------------------------------------
&timeLog("abstract associations of query DEG, and output: ");

$cmd = "perl $Y2H_in_silico --bind $odir/$index.blast.tab.best --fusion $protein_actions --ppi $odir/PPI.txt ";
&runOrDie($cmd);

#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&timeLog("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
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
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

#############################################################################################################
sub abstract_interacts {
    my ($fId, $fIn, $fOut) = @_;

    my %ids;
    open (ID,$fId) or die $!;
    while (<ID>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        $ids{(split /\t/)[0]}=1;
    }
    close ID;

    open (IN,$fIn) or die $!;
    open (OUT,">$fOut") or die $!;
    while (<IN>) {
        chomp;
        next if (/^\s*$/);
        if (/^#/) {
            print OUT "$_\n";
        } else {
            my ($k1, $k2) =(split /\t/)[0,2];
            print OUT "$_\n" if ( exists $ids{$k1} && exists $ids{$k2} );
        }
    }
    close IN;
    close OUT;
}

#############################################################################################################
sub make_sif {
    my ($fIn, $fOut) = @_;

    my %ppi;
    open (PPI,$fIn) or die $!;
    while (<PPI>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        my ($a,$b,$c) = (split /\t/)[0,1,2];
        ($a,$c) = ($a lt $c) ? ($a,$c) : ($c,$a);
        $ppi{$a}{$c} = $b;
    }
    close PPI;

    open (OUT,">$fOut") or die $!;
    for my $i (sort keys %ppi) {
        for my $j (sort keys %{$ppi{$i}}) {
            print OUT "$i\tpp\t$j\n";
        }
    }
    close OUT;
}

#############################################################################################################
#&run_or_die($cmd);
sub run_or_die() {
    my ($cmd) = @_ ;
    &log_current_time($cmd);
    my $flag = system($cmd);
#    my $flag = 0;

    if ($flag){
        &log_current_time("Error: command fail: $cmd");
        exit(1);
    }
}

#############################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#############################################################################################################
sub USAGE {
	my $usage = <<"USAGE";
 ProgramName: $Script
     Version: $version
     Contact: Simon Young <yangxh\@biomarker.com.cn> 
Program Date: 2015-05-21
      Modify: 
 Description: This program is used to build protein-protein interactions for DEG, by interolog transfers with STRING database.
       Usage: 
        Options:
        --qseq  <FILE>  nucleotide sequence file that contain queries, FASTA format, required


        --clade <STR>   organism category or certain eukaryota\'s taxid to transfer interaction
                        'all' | 'archaea' | 'bacteria' | 'eukaryota' | TAXID,           ['eukaryota']
                        eg. 'eukaryota' or 9606 for Homo sapiens
                        [tip1: option 'all' is time-consuming! ]
                        [tip2: your can search TAXID of a species in file $STRING_db/species.taxid.txt!]

        --index <STR>   prefix of output files, optional,                               ['ppi_qurey']

        --odir  <DIR>   output directory, required

        Examples:
            perl $Script --qseq Maize.Known.longest_transcript.fa  --clade 4577 --odir DEG_Analysis/PPI_dir/
            perl $Script --qseq Basic_Analysis/Cluster/Unigene/Bee.Unigene.fa --clade eukaryota --odir DEG_Analysis/PPI_dir/

USAGE
	print $usage;
	exit;
}
