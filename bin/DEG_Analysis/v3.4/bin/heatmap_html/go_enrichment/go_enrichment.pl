#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($od,$fin1,$dir);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fin1,
				"dir:s"=>\$dir,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($od and $fin1 and $dir);
unless (-d $od) {
	mkdir $od || die "can not mkdir $od!";
}
$dir = &ABSOLUTE_DIR($dir);
$od = &ABSOLUTE_DIR($od);
$fin1 = &ABSOLUTE_DIR($fin1);

my $line_count = (split " ",`wc -l $fin1`)[0];

if ($line_count <= 1) {
	print STDERR "There is no data(or only title) in $fin1.";
}
else {
	open IN,$fin1;
	open OUT,">$od/gene_ID.list";
	while (<IN>) {
		chomp;
		next if (/^#/ or /^\s*$/ or ($. == 1));
		my @line = split /\t/;
		print OUT $line[0],"\n";
	}
	close IN;
	close OUT;


	my $cmd = "perl $Bin/bin/KeggGo_enrich_map_webv1.pl -d $od/gene_ID.list -k go_enrichment ";
	$cmd .= " -i $dir -o $od -func go ";
	print $cmd,"\n";
	`$cmd`;
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact: 
Description: 
Usage:
  Options:
  -i   <file>        ID file, ID in first column, first line is title
  -dir <dir>         annotation directory
  -od  <dir>         output dir

  -h         Help

USAGE
	print $usage;
	exit;
}
