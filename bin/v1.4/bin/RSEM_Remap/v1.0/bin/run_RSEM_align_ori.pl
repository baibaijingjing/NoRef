#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<__EOUSAGE__;


#########################################################################
##
##  --transcripts <string>           transcript fasta file
##  --seqType <string>              fq|fa
##
## If Paired-end:
##  --left <string>
##  --right <string>
## or Single-end:
##  --single <string>
##
## Optional:
##
## --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
##
## --no_group_by_component          Trinity-mode, using 'components' as 'genes'
##
## --thread_count                   number of threads to use (default = 4)
##  
######################
##  Non-Trinity options:
##
##  --gene_trans_map <string>        file containing 'gene(tab)transcript' identifiers per line.
##
##
##########################################################################
##  
##  To pass additional parameters to rsem-prepare-reference, 
##    type ' -- ' followed by additional pass-through params
##
##########################################################################

__EOUSAGE__
;


my $help_flag;
my $transcripts;
my $SS_lib_type;
my $left;
my $right;
my $single;
my $no_group_by_component = 0;
my $thread_count = 4;
my $gene_trans_map_file;

&GetOptions ( 'h' => \$help_flag,
              'transcripts=s' => \$transcripts,
              'SS_lib_type=s' => \$SS_lib_type,
              'left=s' => \$left,
              'right=s' => \$right,
              'single=s' => \$single,
              'no_group_by_component' => \$no_group_by_component,
              'thread_count=i' => \$thread_count,
              'gene_trans_map=s' => \$gene_trans_map_file,    
              );


if ($help_flag) {
    die $usage;
}

unless ($transcripts) {
    die $usage;
}




if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(RF|FR|R|F)$/) {
        die "Error, do not recognize SS_lib_type: [$SS_lib_type]\n";
    }
    if ($left && $right && length($SS_lib_type) != 2 ) {
        die "Error, SS_lib_type [$SS_lib_type] is not compatible with paired reads";
    }
}

if ( $thread_count !~ /^\d+$/ ) {
    die "Error, --thread_count value must be an integer";
}

my $RSEM_dir = "/share/nas2/genome/biosoft/trinity/trinityrnaseq_r2013-02-25/trinity-plugins/rsem";

my $cmd = "$RSEM_dir/rsem-prepare-reference";

unless (-s "TRANS.1.ebwt" or -s "TRANS.1.ebwtl") { ## this step already run

        if ($gene_trans_map_file) {
            $cmd .= " --transcript-to-gene-map $gene_trans_map_file ";
        }
        elsif (! $no_group_by_component) {
            my $trans_to_gene_map_file = &write_gene_to_trans_map_file($transcripts);
            
            $cmd .= " --transcript-to-gene-map $trans_to_gene_map_file";
        
        }
        $cmd .= " $transcripts TRANS";
        
        &process_cmd($cmd);
    }


sub process_cmd {
    my ($cmd) = @_;

    print STDERR "**********************************************************************\n";
    print STDERR "**  Running RSEM_align Command:\n";
    print STDERR "**  $cmd\n";
    print STDERR "**********************************************************************\n\n\n";

    
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret: $ret";
    }
    
    return;
}


sub write_gene_to_trans_map_file {
    my ($transcripts_fasta_file) = @_;
    
        
    open (my $fh, $transcripts_fasta_file) or die "Error, cannot open file $transcripts_fasta_file";
    
    my $mapping_file = "$transcripts_fasta_file.component_to_trans_map";
    open (my $ofh, ">$mapping_file") or die "Error, cannot write to file: $mapping_file";
    
    while (<$fh>) {
        if (/>(comp\S+)/) {
            my $acc = $1;
            $acc =~ /^(comp\d+_c\d+)_seq\d+/ or die "Error, cannot parse the trinity component ID from $acc";
            my $comp_id = $1;
            print $ofh "$comp_id\t$acc\n";
        }
	elsif (/>(c\S+)/) {
	    my $acc = $1;
	    $acc =~ /^(c[^\_]+_c\d+)_seq\d+/ or die "Error, cannot parse the trinity component ID from $acc";
	    my $comp_id = $1;
	    print $ofh "$comp_id\t$acc\n";
	}
    }
    close $fh;
    close $ofh;

    return($mapping_file);
}

