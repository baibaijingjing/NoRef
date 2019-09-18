# Writer:         Mengf <mengf@biomarker.com.cn>
# Program Date:   2012.
# Modifier:       Mengf <mengf@biomarker.com.cn>
# Modifier:       lium <lium@biomarker.com.cn>
# Last Modified:  2012-7-28
# Modifier:       baij  <baij@biomarker.com.cn>
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
######################get options


our %opts;
GetOptions(\%opts,"i=s","cfg=s","od=s","h" );


if(!defined($opts{i}) || !defined($opts{cfg}) || !defined($opts{od}) || defined($opts{h})){
	&help();
	exit;
}



###############Time
&timeLog("$programe_dir start.");


###### set env variables
my $temp = `echo \$PATH`;
print "PATH=$temp\n";
$ENV{"PATH"} = ${&readconf("$Bin/../../../../../config/sys.cfg")}{PATH}.":" . $ENV{"PATH"};
$temp = `echo \$PATH`;
print "PATH=$temp\n";

my $Rscript = ${&readconf("$Bin/../../../../../config/sys.cfg")}{Rscript};
################
my $in=&ABSOLUTE_DIR($opts{i});
my $cfg=&ABSOLUTE_DIR($opts{cfg});
&MKDIR($opts{od});
&MKDIR("$opts{od}/work_sh");
my $outdir=&ABSOLUTE_DIR($opts{od});
my $sh_dir=&ABSOLUTE_DIR("$opts{od}/work_sh");
my $notename=`hostname`;chomp $notename;

my %DEG_type;
&LOAD_PARA($cfg,\%DEG_type);



###### reading program from file
our %PROGRAM = ();
$PROGRAM{de_with_deseq}="$Bin/bin/deseq_analysis.r";
$PROGRAM{de_with_ebseq}="$Bin/bin/ebseq_analysis.r";
$PROGRAM{de_with_deseq2}="$Bin/bin/deseq2_analysis.r";
$PROGRAM{de_with_edger}="$Bin/bin/edger_analysis.r";
$PROGRAM{de_with_noRep_edger}="$Bin/bin/edger_withoutrep_analysis.r";
#$PROGRAM{de_with_noiseq}="$Bin/bin/noiseq_analysis.r";



###### read config file into hash
my %CONTROL = ();
read_config_into_hash($cfg, \%CONTROL);



open (SH,">$sh_dir/gene_counts_to_FC_FDR.sh") or die $!;
###### do DE analysis according to control
control_DE_analysis($in,$outdir, \%CONTROL);



###### do DE analysis according to control
sub control_DE_analysis{
	# get para
        my ($in,$outdir, $control) = @_;
        # get cutoff
        my $FDR_cutoff=0;
        my $Pvalue_cutoff=0;
        if (defined $CONTROL{FDR}) {
        $FDR_cutoff = $CONTROL{FDR};
        }
        my $FC_cutoff = $CONTROL{fold};

	# do combinate compare without replicates
	if( exists $$control{Com} ) {
                my $com_groups = $$control{Com};
                my $len = @$com_groups;
		for(my $i=0; $i<$len; $i++) {
			combinate_compare_two_condition_without_replicates($$com_groups[$i], $outdir, $in, $FDR_cutoff, $FC_cutoff,"$outdir/All_gene_fpkm.list");
		}
        }
	# do compare with replicates
	if( exists $$control{Sep} ) {
               my $sep_groups = $$control{Sep};
               my $len = @$sep_groups;
               for(my $i=0; $i<$len; $i++) {
                       compare_with_replicates($$sep_groups[$i], $outdir, $in, $FDR_cutoff, $FC_cutoff,"$outdir/All_gene_fpkm.list");
                }
        }

        close SH;
        &runOrDie("$sh_dir/gene_counts_to_FC_FDR.sh");
}



###### reading config file into hash
sub read_config_into_hash{
	# get config file name
	my ($in, $control) = @_;

	# open config file
	open(IN,"$in") || die "open file $in is failed: $!";

	# reading every line
	while( my $line = <IN> ) {
		chomp $line;
		next if($line!~/\w/); 	# ignore the empty line
		next if($line=~/^#/); 	# ignore the comment line

		# get parameters according the first word
		my @str = $line=~/(\S+)/g;

		# check the number of word in one line
		my $len = @str;
		if( $len != 2 ){ print "$in: $line: the number of world != 2\n"; exit; }

		# store the group info
		if( exists $$control{$str[0]} ){
			# check sample analysis
			if( ($str[0] ne "Com") && ($str[0] ne "Sep") ) {
				print "$in: repeat program word $str[0]\n"; exit;
			}
			# push
			my $tmp = $$control{$str[0]};
			push(@$tmp, $str[1]);
		}else{
			if( ($str[0] eq "Com") || ($str[0] eq "Sep") ) {
				my @tmp = ();
				$tmp[0] = $str[1];
				$$control{$str[0]} = \@tmp;
			} else {
				$$control{$str[0]} = $str[1];
			}
		}
	}

	# check key
	if( (! exists $$control{Com}) && (! exists $$control{Sep}) ) { 
		die "ERROR: key Com and Sep not exists in file $in"; 
	}

	# close file handle
	close IN;
}

=pod
###### reading program from file
sub read_program{
	# get config file name
	my ($in) = @_;

	# open config file
	open(IN,"$in") || die "open fine $in is failed: $!";

	# reading every line
	while( my $line = <IN> ) {
		chomp $line;
		next if($line!~/\w/); 	# ignore the empty line
		next if($line=~/^#/); 	# ignore the comment line

		# get parameters according the first word
		my @str = $line=~/(\S+)/g;

		# check the number of word in one line
		my $len = @str;
		if( $len != 2 ){ print "$in: $line: the number of world != 2\n"; exit; }

		# store the group info
		if( exists $PROGRAM{$str[0]} ){
			print "$in: repeat program word $str[0]\n"; exit;
		}else{
			$PROGRAM{$str[0]} = $str[1];
		}
	}

	# close file handle
	close IN;
}
=cut



# combinate all two condition, and do compare
# note: here the condition has no replicates
sub combinate_compare_two_condition_without_replicates {
	# get para
        my ($sample_id, $outdir, $read_count_file, $FDR_cutoff, $FC_cutoff,$fpkm_file) = @_;
	# abstract id names
	my @sample = split(',', $sample_id);
	# check sample number
	if($#sample < 1) {die "ERROR: the number of sample < 2\n"; }

	# combinate compare
	for (my $i=0;$i<=$#sample-1;$i++) {
		for (my $j=$i+1;$j<=$#sample;$j++) {
			# get output prefix
			my $vs_out = "$sample[$i]" . "_vs_" . "$sample[$j]";
                        my $groupfile="$sample[$i]" . ";" . "$sample[$j]";
			# create dir for output
			&MKDIR("$outdir/$vs_out");
                        &create_vs_config_noRep("$outdir/$vs_out/$vs_out", $sample[$i], $sample[$j]);
			# do compare
			my $cmd;
                        $CONTROL{Method_NR} ||= "EBSeq";
			if ($CONTROL{Method_NR} eq "EBSeq"){
                        $cmd = "$Rscript $PROGRAM{de_with_ebseq}  $read_count_file ";
                                $cmd .= "$sample[$i],$sample[$j] $FDR_cutoff $FC_cutoff $outdir/$vs_out/$vs_out $fpkm_file ";
                                $cmd .= ">$outdir/$vs_out/$vs_out.log 2>&1";
  #                              print SH "$cmd\n";
                          }
                        elsif($CONTROL{Method_NR} eq "edgeR"){
                                $cmd = "$Rscript $PROGRAM{de_with_noRep_edger} $read_count_file ";
                                $cmd .= "$outdir/$vs_out/$vs_out.de.config $FDR_cutoff $FC_cutoff $outdir/$vs_out/$vs_out $fpkm_file ";
                                $cmd .= ">$outdir/$vs_out/$vs_out.log 2>&1";
 #                               print SH "$cmd\n";
                           }
                   print SH  "cd $Bin/bin/ && $cmd \n";


		}
	}
}



# create compare config for two levels
sub create_vs_config {
	# get para
	my ($outprefix, $v1, $v2) = @_;
	# get number of replicate
	my $len_1 = @$v1;
	my $len_2 = @$v2;

	# print config
	# open file
	open(OUT,">$outprefix.de.config") || die "create file $outprefix.de.config is failed: $!";
	# print
	print OUT "sample\tprocess\n";
	# for level one
	for(my $i=0; $i<$len_1; $i++){
		print OUT "$$v1[$i]\tone\n";
	}
	# for level two
	for(my $i=0; $i<$len_2; $i++){
		print OUT "$$v2[$i]\ttwo\n";
	}
	close OUT;
}

sub create_vs_config_noRep {
        my ($outprefix, $v1, $v2) = @_;
        my $len_1 = 1;
        my $len_2 = 1;
        open(OUT,">$outprefix.de.config") || die "create file $outprefix.de.config is failed: $!";
        print OUT "sample\tprocess\n";
        # for level one
        for(my $i=0; $i<$len_1; $i++){
                print OUT "$v1\tone\n";
        }
        # for level two
        for(my $i=0; $i<$len_2; $i++){
                print OUT "$v2\ttwo\n";
        }
        close OUT;
}

# compare one factor which has two or more levels
# # for two levels, here only process one combinate compare
# # for multiple level, after process combinate compare, 
# # here will do one compare which contain all levels
sub compare_with_replicates {
        # get para
        my ($sample_id, $outdir, $read_count_file, $FDR_cutoff, $FC_cutoff,$fpkm_file) = @_;

        # remove empty character
        $sample_id =~ s/\s//g;
        # abstract id groups
        my @groups = split(';', $sample_id);
        # check group numbers
        if($#groups < 1) {die "ERROR: the number of group < 2\n"; }
        # abstract id names
        my @sample = ();
        for(my $i=0; $i<=$#groups; $i++) { 
                my @tmp = split(',', $groups[$i]);
                $sample[$i] = \@tmp; 
        }
       # combinate compare
        for (my $i=0;$i<=$#groups-1;$i++) {
                for (my $j=$i+1;$j<=$#groups;$j++) {
                        my $v1 = $sample[$i];
                        my $v2 = $sample[$j];
                        my $v1_name = join("_", @$v1);
                        my $v2_name = join("_", @$v2);
                        # get output prefix
                        my $vs_out = "$v1_name" . "_vs_" . "$v2_name";
                        my $groupfile=join(",", @$v1).";".join(",", @$v2);
                        # create dir for output
                        &MKDIR("$outdir/$vs_out");

                        # create vs config
                        &create_vs_config("$outdir/$vs_out/$vs_out", $v1, $v2);

                        # do compare
                        my $cmd;
                        $CONTROL{Method_RE} ||= "DESeq";
                        if ($CONTROL{Method_RE} eq "DESeq") {
                        $cmd = "$Rscript $PROGRAM{de_with_deseq} $read_count_file ";
                                $cmd .= "$outdir/$vs_out/$vs_out.de.config $FDR_cutoff $FC_cutoff $outdir/$vs_out/$vs_out $fpkm_file ";
                                $cmd .= ">$outdir/$vs_out/$vs_out.log 2>&1";
#                                print SH "$cmd\n";
                        } 
                        elsif($CONTROL{Method_RE} eq "edgeR") {
                                $cmd = "$Rscript $PROGRAM{de_with_edger} $read_count_file ";
                                $cmd .= "$outdir/$vs_out/$vs_out.de.config $FDR_cutoff $FC_cutoff $outdir/$vs_out/$vs_out $fpkm_file ";
                                $cmd .= ">$outdir/$vs_out/$vs_out.log 2>&1";
#                                print SH "$cmd\n";
                        }
                        elsif($CONTROL{Method_RE} eq "DESeq2") {
                                $cmd = "$Rscript $PROGRAM{de_with_deseq2} $read_count_file ";
                                $cmd .= "$outdir/$vs_out/$vs_out.de.config $FDR_cutoff $FC_cutoff $outdir/$vs_out/$vs_out $fpkm_file ";
                                $cmd .= ">$outdir/$vs_out/$vs_out.log 2>&1";
#                                print SH "$cmd\n";
                        }
                 print SH  "cd $Bin/bin/ &&  $cmd \n"; 
    }
  }
}

=pod
sub compare_with_replicates_with_NOISeq {
        my ($sample_id, $outdir, $read_count_file) = @_;
        $sample_id =~ s/\s//g;
        my @groups = split(';', $sample_id);
        if($#groups < 1) {die "ERROR: the number of group < 2\n"; }
        my @sample = ();
        for(my $i=0; $i<=$#groups; $i++) {
                my @tmp = split(',', $groups[$i]);
                $sample[$i] = \@tmp;
        }
        for (my $i=0;$i<=$#groups-1;$i++) {
                for (my $j=$i+1;$j<=$#groups;$j++) {
                        my $v1 = $sample[$i];
                        my $v2 = $sample[$j];
                        my $v1_name = join("_", @$v1);
                        my $v2_name = join("_", @$v2);
                        my $vs_out = "$v1_name" . "_vs_" . "$v2_name";
                        &MKDIR("$outdir/$vs_out");
                        create_vs_config("$outdir/$vs_out/$vs_out", $v1, $v2);
                        my $cmd = "cd $Bin/bin/ && $Rscript $PROGRAM{de_with_noiseq} $read_count_file ";
                        $cmd .= "$outdir/$vs_out/$vs_out.de.config $outdir/$vs_out/$vs_out";
                        print SH "$cmd\n";
                }
        }
}
=cut





###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);



sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Shell_run
{
	my $sh = shift;
	my $t=GetTime();
	print "\[$t\] $sh starts:\n";
	unless(-e "$sh.finish"){
		my $cmd_res=system"sh $sh >$sh.log 2>&1";
		if($cmd_res){
			die"Err:$sh run wrong,$!\n";
		}
		else{
			system"touch $sh.finish";
		}
	}
	my $e=GetTime();
	print "\[$e\] $sh ends.\n";
}

###########subs
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
		die"Warning just for file and dir $in in [sub ABSOLUTE_DIR]\n";
	}
	
	chdir $cur_dir;
	return $return;
}

sub help{
print <<"Usage End.";
Description: Differencial Express Analysis pipeline;
Version: $ver

Usage:
-i                All_gene_counts.list                            must be given;
-cfg              soft parameter to DE miRNA Analysis             must be given;
-od               Out dir                                         must be given;
-h                help document
Usage End.
exit;
}
