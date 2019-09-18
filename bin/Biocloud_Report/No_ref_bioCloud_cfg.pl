#!/usr/bin/perl -w
use strict;
use warnings;
use Config::General;
use Cwd qw(abs_path getcwd);
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#use Array::Diff;
my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($id,$data_cfg,$detail_cfg);

GetOptions(
				"help|h" =>\&USAGE,
				"id:s"=>\$id,
				"cfg1:s"=>\$data_cfg,
				"cfg2:s"=>\$detail_cfg,
				) or &USAGE;
&USAGE unless ($id and $data_cfg and $detail_cfg);
$id=abs_path($id);
$detail_cfg=abs_path($detail_cfg);
my $project_id;
die "$id does not exists,$!\n" unless -d $id ;

##################################################detail cfg and Create New Cfg
&load_data_cfg($data_cfg,"$id/data.cfg");
&load_detail_cfg($detail_cfg,"$id/detail.cfg");

sub load_data_cfg{
	my$data_cfg=shift;
	my$out=shift;
	my$Sample_name;
	my%Data;	
	open (IN,$data_cfg) or die $!;
	while (<IN>) {
	    if (/^Sample/) {
	        $Sample_name=(split/\s+/,$_)[1];
	    }
	    if ($_=~/^fq1/ || $_=~/^fq2/ || $_=~/^raw_fq1/ || $_=~/^raw_fq2/) {
	        my $file=(split/\s+/,$_)[1];
	        unless (-e $file) {
	            print "$file is not exist!\n";
	        }
	        $Data{$Sample_name}{fq1}=$file if $_=~/^fq1/;
	        $Data{$Sample_name}{fq2}=$file if $_=~/^fq2/;
	        $Data{$Sample_name}{raw_fq1}=$file if $_=~/^raw_fq1/;
	        $Data{$Sample_name}{raw_fq2}=$file if $_=~/^raw_fq2/;
	    }
	}
	close IN;
	&writeConfig(\%Data,$out);
	
	
}

sub load_detail_cfg{
	my$detail_cfg=shift;
	my$out=shift;
	my%detail;
	#`sed  -i 's/para_K-cov/para_K_cov/g'  $detail_cfg`;
	my@conItems=qw(projectInfo basicAnalysis geneAnn DEGAnalysis SNPAnalysis);
	$/="\n>>>>>";
	my $COM_limit=0;
	my $Check_num=0;
	my $key_word;
	open (IN,$detail_cfg) or die $!;
	open OUT,">$out" or die "$!";
	my@COM;
	while (<IN>) {
	    chomp;
	    $Check_num=$.;
		print OUT "<$conItems[$Check_num-1]>\n";
		print OUT "$_\n";
		print OUT "</$conItems[$Check_num-1]>\n";
		$project_id=$1 if(/Contract_NO\s+(\S+)/);
	}

	    
	$/="\n";
	close IN;
	if ($Check_num!=5) {
	    print "Check Your detail_cfg,it should have 4 '>>>>>' \n";die;
	}
	
	close OUT;
}

sub writeConfig{
    my($config,$outFile)=@_;
#    if(-e $outFile){
#        system(qq(mv $outFile  $outFile.bak));
#    }
    my $con=Config::General->new(-SplitPolicy=> 'custom' , -SplitDelimiter => '\s*=\s*',-StoreDelimiter=>'  ',-SaveSorted => 1);
     $con->save_file($outFile, \%$config);
}

sub USAGE {
        my $usage=<<"USAGE";
#-----------------------------------------------------------------------------------------
  Program: $Script
  Version: $version
    Usage:
      Options:
      --id             <dir> directory of Web report.
      --cfg1           <str> data cfg file
      --cfg2           <str> detail file 
#-----------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}
