# Writer:         Mengf <mengf@biomarker.com.cn>
# Program Date:   2012.
# Modifier:       Mengf <mengf@biomarker.com.cn>
# Modifier:       lium <lium@biomarker.com.cn>
# Modifier:       baij <baij@biomarker.com.cn>
# Last Modified:  2016-01-12
my $ver="1.0.0";

use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $programe_dir=basename($0);
my $path=dirname($0);
######################


our %opts;
GetOptions(\%opts,"i=s","od=s","h" );


if(!defined($opts{i}) || !defined($opts{od}) || defined($opts{h})){
	&help();
	exit;
}


###############Time
&timeLog("$programe_dir start.");


###### set env variables
my $temp = `echo \$PATH`;
#print "PATH=$temp\n";
$ENV{"PATH"} = ${&readconf("$Bin/../../../config/sys.cfg")}{PATH}.":" . $ENV{"PATH"};
$temp = `echo \$PATH`;
#print "PATH=$temp\n";

################
&MKDIR($opts{od});
my $gene_exp_count =&ABSOLUTE_DIR($opts{i});
my $outdir=&ABSOLUTE_DIR($opts{od});
my $cmd;

##########################################################
##
##           fpkm correlation plot
##
##########################################################

##### fpkm free combination plot
&MKDIR("$outdir");
&MKDIR("$outdir/cor_plot");
$cmd = ${&readconf("$Bin/../../../config/sys.cfg")}{Rscript}." $Bin/bin/fpkm_plot/fpkm_cor_plot.r $gene_exp_count $outdir ";
$cmd .= ">$outdir/free_cor_plot.log 2>&1";
&runOrDie($cmd);

##########################################################
##
##           fpkm density plot
##
##########################################################
###### do fpkm density
&MKDIR("$outdir");
$cmd = ${&readconf("$Bin/../../../config/sys.cfg")}{Rscript}." $Bin/bin/fpkm_plot/fpkm_density_plot_func.r $opts{i} $outdir ";
$cmd .= ">$outdir/plot_density.log 2>&1";
&timeLog($cmd);
&runOrDie($cmd);

##########################################################
##
##           fpkm box plot
##
##########################################################
&MKDIR("$outdir");
$cmd = ${&readconf("$Bin/../../../config/sys.cfg")}{Rscript}." $Bin/bin/fpkm_plot/fpkm_box_plot_func.r $opts{i} $outdir ";
$cmd .= ">$outdir/plot_box.log 2>&1";
&timeLog($cmd);
&runOrDie($cmd);

###############Time
&timeLog("$programe_dir finish.");
&Runtime($BEGIN);

##############subs#################
sub LOAD_PARA {
	my $para_file= shift;
	my $para = shift;

	my $error_status = 0;
	open(IN, $para_file) || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		if ($para_key=~/_G/) {
			push @{$para->{$para_key}},$para_value;
			next;
		}
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub parse_config
{ # load config file
	my $config_file= shift;
	my $DataBase= shift;
	
	my $error_status = 0;
	open(IN,$config_file) || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if(/^$/ or /^\#/);
		my ($software_name,$software_address) = split(/\s+/,$_);
		$DataBase->{$software_name} = $software_address;
		if (! -e $software_address){
			warn "Non-exist:  $software_name  $software_address\n"; 
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
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

sub help{
print <<"Usage End.";
Version: $ver
Usage:
-i                All Samples Gene counts files                   must be given;
-od               Out dir                                         must be given;
-h                help document
Usage End.
exit;
}
