#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<__EOUSAGE__;
#########################################################################
##
##  --seqType <string>              fq|fa
## 
##  If Paired-end:
##
##  --left <string>
##  --right <string>
##  
##    or Single-end:
##
##  --single <string>
##
##
##
## Optional:
## 
## --prefix <string>                prefix for RSEM output files (default: 'RSEM')
##
## --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
##
## --thread_count                   number of threads to use (default = 4)
##
## --debug                  retain intermediate files
##########################################################################
##  
##  To pass additional parameters to rsem-calculate-expression, 
##    type ' -- ' followed by additional pass-through params
##
##########################################################################


__EOUSAGE__

    ;

my $help_flag;
my $bam_file;
my $paired_flag;
my $DEBUG_flag = 0;
my $SS_lib_type;
my $thread_count = 4;
my $seqType;
my $left;
my $right;
my $single;
my $prefix = "RSEM";

&GetOptions ( 'h' => \$help_flag,
              'name_sorted_bam=s' => \$bam_file,
              'paired' => \$paired_flag,
              'debug' => \$DEBUG_flag,
              'SS_lib_type=s' => \$SS_lib_type,
              'thread_count=i' => \$thread_count,
              'seqType=s' => \$seqType,
              'left=s' => \$left,
              'right=s' => \$right,
              'single=s' => \$single,              
              'prefix=s' => \$prefix,
              
              );



if ($help_flag) {
    die $usage;
}

unless ($seqType && ($single || ($left && $right))) {
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

    my $keep_intermediate_files_opt = ($DEBUG_flag) ? "--keep-intermediate-files" : "";


    my $SS_opt = "";
    if ($SS_lib_type) {
        if ($SS_lib_type =~ /^F/) {
            $SS_opt = "--forward-prob 1.0";
        }
        else {
            $SS_opt = "--forward-prob 0";
        }
    }
    
 my  $cmd = "$RSEM_dir/rsem-calculate-expression @ARGV " ## allow for custom params
        . " -p $thread_count"
        . " $keep_intermediate_files_opt"
        . " $SS_opt";
    
    if ($seqType eq "fa")  {
        $cmd .= " --no-qualities";
    }

    if ($left && $right) {
        $cmd .= " --paired-end $left $right";
    }
    else {
        $cmd .= " $single";
    }
            
    $cmd .= " TRANS $prefix";
    
    &process_cmd($cmd);

sub process_cmd {
    my ($cmd) = @_;

    print STDERR "**********************************************************************\n";
    print STDERR "**  Running RSEM_estimate Command:\n";
    print STDERR "**  $cmd\n";
    print STDERR "**********************************************************************\n\n\n";

    
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret: $ret";
    }
    
    return;
}

