#!/usr/bin/perl -w
#
# Copyright (c) BMK 2012
# Writer:         mengf <mengf@biomarker.com.cn>
# Program Date:   2012.
# Modifier:       mengf <mengf@biomarker.com.cn>
# Last Modified:  2012.
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"cfg=s","od=s","h");
if (!defined($opts{cfg})||!defined($opts{od})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";


###############
my %para;
my $cfg=&ABSOLUTE_DIR($opts{cfg});
&para_load($cfg,\%para);

#my $group_num=0;
#foreach my $key (keys %para) {
#	if ($key=~/^Sample/) {
#		$group_num++;
#	}
#}

###############
&MKDIR($opts{od});
my $od=&ABSOLUTE_DIR($opts{od});
my $Unigene_od="$od/Unigene";
&MKDIR($Unigene_od);
my $DEG_od="$od/DEG_Analysis";
&MKDIR($DEG_od) if (defined $para{DEG});
my $cleandata_od="$od/cleandata";
&MKDIR($cleandata_od);
my $SNP_od="$od/SNP_Analysis";
&MKDIR($SNP_od) if (defined $para{SNP});

#open LOG,"$od/log.txt" || die;
############### Extract Assembly dir
if (defined $para{Rawdata}) {
#	system "cp $para{Rawdata}/AllSample.data.stat $cleandata_od";
	system "cp $para{Rawdata}/AllSample_GC_Q.stat $cleandata_od";
	system "cp -r $para{Rawdata}/PNG $cleandata_od";
}

if (defined $para{Basic}) {
	mkdir "$od/Assembly" unless -d "$od/Assembly";
	my @DIR=glob "$para{Basic}/Assembly/Trinity_assembly/*";
	foreach my $dir (@DIR) {
		my $dir_name=basename $dir;
		mkdir "$od/Assembly/$dir_name" unless -d "$od/Assembly/$dir_name";
		system "cp -r $dir/${dir_name}_Result/* $od/Assembly/$dir_name";
        system "rm $od/Assembly/$dir_name/*/*_normal_len.fa";
        system "rm $od/Assembly/$dir_name/Transcripts/*.component_to_trans_map" if (-f "$od/Assembly/$dir_name/Transcripts/*.component_to_trans_map");
	}
	system "cp -r $para{Basic}/Remap/geneExpression/ $od";
	system "cp -r $para{Basic}/Cluster/Unigene/Unigene_CDS $Unigene_od";
	system "cp $para{Basic}/Cluster/Unigene/*.Unigene.fa $Unigene_od";
	system "cp $para{Basic}/Cluster/All_trans.cd-hit_clu.fasta_cl_clusters $Unigene_od/$para{Project}.cluster.xls" if -f "$para{Basic}/Cluster/All_trans.cd-hit_clu.fasta_cl_clusters";
	system "cp $para{Basic}/Cluster/Unigene/*.xls $Unigene_od";
	system "cp $para{Basic}/Cluster/Unigene/*.png $Unigene_od";
	system "cp -r $para{Basic}/Cluster/Unigene/Unigene_SSR/ $Unigene_od";
	#system "rm $Unigene_od/Unigene_SSR/*.ini $Unigene_od/Unigene_SSR/*ssr.primer.xls $Unigene_od/Unigene_SSR/*detail.xls $Unigene_od/Unigene_SSR/*misa";
	system "rm $Unigene_od/Unigene_SSR/*.ini";
}




if (defined $para{ANNO}) {
	&MKDIR("$Unigene_od/Unigene_Anno");
	system "cp -r $para{ANNO}/*.stat $Unigene_od/Unigene_Anno";
	system "cp -r $para{ANNO}/*.txt $Unigene_od/Unigene_Anno";
	system "cp -r $para{ANNO}/*.Kegg.* $Unigene_od/Unigene_Anno";
	system "cp -r $para{ANNO}/Kegg_map $Unigene_od/Unigene_Anno";
	system "cp -r $para{ANNO}/*.xls $Unigene_od/Unigene_Anno";
	system "cp $para{ANNO}/*.png $Unigene_od/Unigene_Anno";
}

if (defined $para{DEG}) {
	my @DEG_dir=glob "$para{DEG}/*";
	my %Anno;
	my %Stat;
	system "mv $para{DEG}/All_gene_expression.list $DEG_od/All_gene_expression.list" unless (-f "$para{DEG}/All_gene_expression.list");
	foreach my $dir (@DEG_dir) {
		if (-f $dir) {
			next if $dir=~/\.svg$/;
			next if $dir=~/\.log$/;
			system "cp $dir $DEG_od";
		}
		if (-d $dir) {
			$dir=~m/.*\/(\S+)/;
			my $nam=$1;
			&MKDIR("$DEG_od/$nam") unless ($nam eq "Venn");
			if ($dir=~/\/All_DEG$/) {
				system "cp -r $dir/* $DEG_od/$nam/";
			}
			elsif ($dir=~/\/DEG_PPI$/) {
                system "cp -r $dir/* $DEG_od/$nam/";
                system "rm -r $DEG_od/DEG_PPI/work_sh/";
            }
			elsif ($dir=~/\/work_sh$/) {
				`rm -r $DEG_od/$nam`;
			}
			elsif ($dir=~/\/density$/) {
				system "cp -r $dir/*.png $od/geneExpression/";
				system "cp -r $dir/*.cor $od/geneExpression/";
				system "cp -r $dir/cor_plot $od/geneExpression/";
                system "rm -r $DEG_od/density/";
			}
			elsif ($dir=~/_vs_/){
				my $file=(glob "$dir/*.DEG_final.xls")[0];
                                if($file){
				open (IN,$file) or die $!;
				while (<IN>) {
					chomp;
					next if /^\#/;
					my $type=(split/\s+/,$_)[-1];
					$Stat{$nam}{up}++ if $type eq 'up';
					$Stat{$nam}{down}++ if $type eq 'down';
					$Stat{$nam}{total}++;
				}
				close IN;
				if(!defined $Stat{$nam}{total}){
					$Stat{$nam}{up}=0;$Stat{$nam}{down}=0;$Stat{$nam}{total}=0;
				}
				system "cp -r $dir/Anno_enrichment/* $DEG_od/$nam/";
				#system "cp -r $dir/DEG_Cluster $DEG_od/$nam";
				system "cp -r $dir/*final.* $DEG_od/$nam";
				system "cp $dir/*cluster.txt $DEG_od/$nam" if $Stat{$nam}{total}!=0;
				system "cp -r $dir/*.png $DEG_od/$nam" unless $dir=~/.*_vs_.*_vs_.*/;
                system "rm $DEG_od/$nam/go_enrichment/*.svg.list" if (-f "$DEG_od/$nam/go_enrichment/*.svg.list");

				if (-d "$dir/Graphic") {
					system "cp $dir/Graphic/*.png $DEG_od/$nam";
					system "cp $dir/Graphic/*.pdf $DEG_od/$nam";
				}
				my %Site;
				$file=(glob "$dir/Anno_enrichment/*.annotation.xls")[0];
				open (IN,$file) or die $!;
				while (<IN>) {
					chomp;
					if (/^\#/) {
						my @Anno=split/\s+/,$_;
						for (my $s=0;$s<@Anno ;$s++) {
							if ($Anno[$s] eq 'COG_class') {
								$Site{'COG'}=$s;
							}
							if ($Anno[$s] eq 'KOG_class') {
								$Site{'KOG'}=$s;
							}
                            if ($Anno[$s] eq 'Swissprot_annotation') {
								$Site{'Swiss-Prot'} = $s;
                            }
							elsif ($Anno[$s]=~/^([^_]+)_annotation/) {
								$Site{$1}=$s;
							}
						}
#                        $Anno{$nam}{'Annotated'} = 0;
					}
					else{
						my @Info=split /\t+/,$_;
						foreach my $key (keys %Site) {
							$Anno{$nam}{$key}||=0;
							$Anno{$nam}{$key}++ unless ($Info[$Site{$key}] eq '--');
						}
                        $Anno{$nam}{'Annotated'} ++;
					}
				}
				close IN;
			}else{
                                $Stat{$nam}{up}=0;$Stat{$nam}{down}=0;$Stat{$nam}{total}=0;            
                                system "cp -r $dir/*_vs_*.final.xls $DEG_od/$nam";
                                system "cp -r $dir/*.png $DEG_od/$nam" unless $dir=~/.*_vs_.*_vs_.*/;
                              }
                       }
		}
	}

      if(%Anno){
	open (OUT,">$DEG_od/DEG.anno.stat") or die $!;
	my $limit_anno=0;
	foreach my $key (sort keys %Anno) {#!
		if ($limit_anno==0) {
			print OUT "#DEG Set";
			foreach my $key1 (sort keys %{$Anno{$key}}) {
				print OUT "\t$key1";
			}
			print OUT "\n";
			$limit_anno++;
		}
		print OUT "$key";
		foreach my $key1 (sort keys %{$Anno{$key}}) {
			print OUT "\t$Anno{$key}{$key1}";
		}
		print OUT "\n";
	}
	close OUT;
       }

	open (OUT,">$DEG_od/DEG.stat") or die $!;
	print OUT "DEG Set\tAll DEG\tup-regulated\tdown-regulated\n";
	foreach my $key (sort keys %Stat) {
		$Stat{$key}{up}||=0;
		$Stat{$key}{down}||=0;
		print OUT "$key\t$Stat{$key}{total}\t$Stat{$key}{up}\t$Stat{$key}{down}\n";
	}
	close OUT;
}

#if (defined $para{SNP}) {
#	my %SNP;
#	my @SNP_dir=glob "$para{SNP}/*";
#	open STAT,">$SNP_od/SNP.stat.xls" || die;
#	print STAT "#Type\tS1.homo.S2.homo\tS1.hete.S2.homo\tS1.homo.S2.hete\tS1.hete.S2.hete\tTotal\n";
#	foreach my $dir (@SNP_dir) {
#		if (-d $dir) {
#			$dir=~m/.*\/(\S+)/;
#			my $nam=$1;
#			&MKDIR("$SNP_od/$nam");
#			system "cp $dir/*.xls $SNP_od/$nam";
#			system "cp $dir/*.snp.density.* $SNP_od/$nam";
#			system "rm $SNP_od/$nam/*.snp.stat.xls";
#			my $file;
#			foreach my $f (glob "$dir/*.homo.snp.xls") {
#				unless ($f=~/\.hete_/) {
#					$file=$f;
#					last;
#				}
#			}
#			my $mm=`wc -l $file`;
#			$SNP{$nam}{mm}=(split/\s+/,$mm)[0]-1;
#			my $tm=`wc -l $dir/*.hete_*.homo.snp.xls`;
#			$SNP{$nam}{tm}=(split/\s+/,$tm)[0]-1;
#			my $mt=`wc -l $dir/*.homo_*.hete.snp.xls`;
#			$SNP{$nam}{mt}=(split/\s+/,$mt)[0]-1;
#			foreach my $f (glob "$dir/*.hete.snp.xls") {
#				unless ($f=~/\.homo_/) {
#					$file=$f;
#					last;
#				}
#			}
#			my $tt=`wc -l $file`;
#			$SNP{$nam}{tt}=(split/\s+/,$tt)[0]-1;
#			$SNP{$nam}{total}=$SNP{$nam}{mm}+$SNP{$nam}{tm}+$SNP{$nam}{mt}+$SNP{$nam}{tt};
#		}
#	}
#	foreach my $name (keys %SNP) {
#		print STAT "$name\t$SNP{$name}{mm}\t$SNP{$name}{tm}\t$SNP{$name}{mt}\t$SNP{$name}{tt}\t$SNP{$name}{total}\n";
#	}
#	close STAT;
#}

if (defined $para{SNP}) {
    system "cp $para{SNP}/AllSample.SNP_density* $SNP_od/";
    system "cp $para{SNP}/SNP/final.snp.list $para{SNP}/SNP/final.*.vcf $SNP_od/";
    system "cp -r $para{SNP}/SNP/stat/* $SNP_od/";
    system "cp -r $para{SNP}/Pairwised_SNP/ $SNP_od/";
}

################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
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

sub para_load {
	my ($file,$para)=@_;
	open IN,$file||die "$!";
	while (<IN>) {
		chomp;
		s/\r+//g;
		next if(/^$/||/^\#/);
		my ($key,$info)=(split/\s+/,$_)[0,1];
		if(!$key){print "$_\n";die;}
		$para->{$key}=$info;
	}
	close IN;
}
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "\nTotal elapsed time: ${t}s\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub help{
	print << "	Usage End.";
	Description: Extract No_Ref Transcriptome Reaults for Html Process;
	version:$ver
	Usage:
		-cfg  <STR>   config file   force
		-od   <STR>   outdir        force
		-h            help
	Usage End.
		exit;
}
