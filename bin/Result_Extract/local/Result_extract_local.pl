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
use newPerlBase;
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
################# mkdir webReport Result dirctory
&MKDIR($opts{od});
my $od=&ABSOLUTE_DIR($opts{od});
my $rawData_od = "$od/BMK_1_rawData";
&MKDIR($rawData_od);
my $unigeneAssembly_od = "$od/BMK_2_Unigene_Assembly";
&MKDIR($unigeneAssembly_od);
my $geneExpression_od = "$od/BMK_3_geneExpression";
&MKDIR($geneExpression_od);
my $unigeneAnno_od = "$od/BMK_4_Unigene_Anno";
&MKDIR($unigeneAnno_od);
my $unigeneStructure_od = "$od/BMK_5_Unigene_Structure";
&MKDIR($unigeneStructure_od);
my $DEG_od="$od/BMK_6_DEG_Analysis";
&MKDIR($DEG_od) if (defined $para{DEG});
###########################################
############### Extract Assembly dir
&runOrDie("cp $Bin/readme.pdf $od");
if (defined $para{Rawdata}) {

	my $png_od = "$rawData_od/PNG";
	&MKDIR($png_od);
	&runOrDie( "cp $para{Rawdata}/AllSample_GC_Q.stat $Bin/BMK_1_rawData/readme.pdf $rawData_od");
	&runOrDie( "cp -r $para{Rawdata}/PNG/* $png_od");
}

if (defined $para{Basic}) {
	my $cluster_dir = "$para{Basic}/Cluster/Unigene";
	&MKDIR("$unigeneAssembly_od/Final_Unigene");
	&para_load("$para{Basic}/../Config/basicAnalysis.cfg",\%para);
	&runOrDie( "cp $cluster_dir/*.Unigene.distribution.png $cluster_dir/*.Unigene.fa $cluster_dir/*.Unigene.stat.xls $unigeneAssembly_od/Final_Unigene");
	&MKDIR("$unigeneAssembly_od/Map_stat");

	if($para{Combine}==1)
	{
		my $tran_od = "$para{Basic}/Assembly/Trinity_assembly/All_Combination/All_Combination_Result/Transcripts/";
		&MKDIR("$unigeneAssembly_od/Final_Transcript");
		&runOrDie( "cp $tran_od/*.distribution.png $tran_od/*.Transcripts.fa $tran_od/*.stat.xls $unigeneAssembly_od/Final_Transcript");
		&runOrDie( "cp $para{Basic}/Remap/rmap/*/*.Mapped.stat.xls $unigeneAssembly_od/Map_stat");
		&runOrDie("cp $Bin/BMK_2_Unigene_Assembly/readme_com.pdf $unigeneAssembly_od/readme.pdf");
	}
	if($para{Separate}==1)
	{
		my @Separate = glob("$para{Basic}/Assembly/Config/*.config");
		foreach my $sep(@Separate)
		{
			$sep=basename($sep);
			$sep=~s/\.config//;
			my $sep_dir = "$para{Basic}/Assembly/Trinity_assembly/$sep/${sep}_Result";
			&MKDIR("$unigeneAssembly_od/${sep}_Assembly/Unigene");
			&runOrDie("cp $sep_dir/Unigenes/*distribution.png $sep_dir/Unigenes/*.Unigenes.fa $sep_dir/Unigenes/*.Unigenes.stat.xls $unigeneAssembly_od/${sep}_Assembly/Unigene");
			&MKDIR("$unigeneAssembly_od/${sep}_Assembly/Transcript");
			&runOrDie("cp $sep_dir/Transcripts/*distribution.png $sep_dir/Transcripts/*.Transcripts.fa $sep_dir/Transcripts/*.stat.xls $unigeneAssembly_od/${sep}_Assembly/Transcript");
		}
		&runOrDie( "cp $para{Basic}/Remap/rmap/*/*.Mapped.stat.xls $unigeneAssembly_od/Map_stat");
		&runOrDie("cp $Bin/BMK_2_Unigene_Assembly/readme_sep.pdf $unigeneAssembly_od/readme.pdf");

	}elsif($para{SamG}==1)
	{
		my $num = 1;
		while(exists $para{"Group".$num})
		{
			my $group = $para{"Group".$num};
			my $group_dir = "$para{Basic}/Assembly/Trinity_assembly/Group$num";
			$group =~s/,/_/g;
			my $Unigene = "$unigeneAssembly_od/${group}_Assembly/Unigene";
			&MKDIR($Unigene);
			&runOrDie("cp $group_dir/*/Unigenes/*distribution.png $group_dir/*/Unigenes/*.Unigenes.fa $group_dir/*/Unigenes/*.Unigenes.stat.xls $Unigene");
			my $trans = "$unigeneAssembly_od/${group}_Assembly/Transcript";
			&MKDIR($trans);
			&runOrDie("cp $group_dir/*/Transcripts/*.distribution.png $group_dir/*/Transcripts/*.Transcripts.fa $group_dir/*/Transcripts/*.stat.xls $trans");
			$num++;
		}
		&runOrDie( "cp $para{Basic}/Remap/rmap/*/*.Mapped.stat.xls $unigeneAssembly_od/Map_stat");
		&runOrDie("cp $Bin/BMK_2_Unigene_Assembly/readme_sep.pdf $unigeneAssembly_od/readme.pdf");

	}

	####### BMK_3_geneExpression dir
	my $remap_dir = "$para{Basic}/Remap/geneExpression/";
	&MKDIR("$geneExpression_od/BMK_1_randcheck");
	&runOrDie( "cp $remap_dir/*.randcheck.png $geneExpression_od/BMK_1_randcheck");
	&MKDIR("$geneExpression_od/BMK_2_insertSize");
	&runOrDie( "cp $remap_dir/*.insertSize.png $geneExpression_od/BMK_2_insertSize");
	&MKDIR("$geneExpression_od/BMK_3_saturation");
	&runOrDie( "cp $remap_dir/*.express.gene_tag.png $geneExpression_od/BMK_3_saturation");

	###### BMK_5_Unigene_Structure/BMK_1_Unigene_CDS
        if(defined $para{ANNO}){
	&MKDIR("$unigeneStructure_od/BMK_1_Unigene_CDS");
	&runOrDie( "cp $para{ANNO}/../Unigene_CDS_Predict/* $unigeneStructure_od/BMK_1_Unigene_CDS && rm $unigeneStructure_od/BMK_1_Unigene_CDS/*.svg");
	&runOrDie( "cp $Bin/BMK_5_Unigene_Structure/BMK_1_Unigene_CDS/readme.pdf $unigeneStructure_od/BMK_1_Unigene_CDS");}
	###### BMK_5_Unigene_Structure/BMK_2_Unigene_SSR
	&MKDIR("$unigeneStructure_od/BMK_2_Unigene_SSR");
	&runOrDie( "cp $cluster_dir/Unigene_SSR/*.fa  $cluster_dir/Unigene_SSR/*.density.png $cluster_dir/Unigene_SSR/*.result.xls $cluster_dir/Unigene_SSR/*.stat.xls $unigeneStructure_od/BMK_2_Unigene_SSR");
	&runOrDie( "cp $Bin/BMK_5_Unigene_Structure/BMK_2_Unigene_SSR/readme.pdf $unigeneStructure_od/BMK_2_Unigene_SSR");
}
&runOrDie( "cp $Bin/BMK_3_geneExpression/readme.pdf $geneExpression_od");
############ BMK_5_Unigene_Structure/BMK_3_SNP_Analysis/
if (defined $para{SNP})
{
	&MKDIR("$unigeneStructure_od/BMK_3_SNP_Analysis/BMK_1_All_SNP");
	&MKDIR("$unigeneStructure_od/BMK_3_SNP_Analysis/BMK_2_Pairwised_SNP");
	&MKDIR("$unigeneStructure_od/BMK_3_SNP_Analysis/BMK_3_Single_Sample_SNP");
	&runOrDie( "cp $para{SNP}/*SNP_density.*  $para{SNP}/SNP/stat/*.snp.stat $unigeneStructure_od/BMK_3_SNP_Analysis/BMK_1_All_SNP");
	&runOrDie( "cp $para{SNP}/SNP/final.snp.list $unigeneStructure_od/BMK_3_SNP_Analysis/BMK_1_All_SNP/final.snp.xls");
	&runOrDie( "cp $para{SNP}/Pairwised_SNP/* $unigeneStructure_od/BMK_3_SNP_Analysis/BMK_2_Pairwised_SNP");
	&runOrDie( "cp $para{SNP}/SNP/stat/*.list $unigeneStructure_od/BMK_3_SNP_Analysis/BMK_3_Single_Sample_SNP");
	&runOrDie( "cp $Bin/BMK_5_Unigene_Structure/BMK_3_SNP_Analysis/readme.pdf $unigeneStructure_od/BMK_3_SNP_Analysis");

}




if (defined $para{ANNO}) {

	&runOrDie( "cp $Bin/BMK_4_Unigene_Anno/readme.pdf $unigeneAnno_od");
	&runOrDie( "cp $para{ANNO}/All_Database_annotation.xls $para{ANNO}/Function_Annotation.stat.xls $unigeneAnno_od");
	&MKDIR("$unigeneAnno_od/BMK_1_annotation");
	&runOrDie( "cp $para{ANNO}/*.anno.txt $para{ANNO}/*_class.txt $unigeneAnno_od/BMK_1_annotation");
	&runOrDie( "cp $Bin/BMK_4_Unigene_Anno/BMK_1_annotation/readme.pdf $unigeneAnno_od/BMK_1_annotation/");
	&MKDIR("$unigeneAnno_od/BMK_2_statistic");
	&runOrDie( "cp $para{ANNO}/*.png $para{ANNO}/*.pdf $para{ANNO}/*.stat $para{ANNO}/*.GO_enrichment.stat.xls $unigeneAnno_od/BMK_2_statistic/");
	&runOrDie( "cp $Bin/BMK_4_Unigene_Anno/BMK_2_statistic/readme.pdf $unigeneAnno_od/BMK_2_statistic/");
	&MKDIR("$unigeneAnno_od/BMK_3_KEGG_map");
	&runOrDie( "cp $para{ANNO}/Kegg_map/* $unigeneAnno_od/BMK_3_KEGG_map");
}

if (defined $para{DEG}) {
	my @DEG_dir=glob "$para{DEG}/*";
	my %Anno;
	my %Stat;
	my $count = 2;
	&runOrDie( "cp $para{DEG}/All_gene_expression.list $geneExpression_od/All_gene_expression.xls");
	&runOrDie( "cp $Bin/BMK_6_DEG_Analysis/readme.pdf $DEG_od");
	foreach my $dir (@DEG_dir) {
		if (-d $dir) {
			$dir=~m/.*\/(\S+)/;
			my $nam=$1;
			if ($dir=~/\/All_DEG$/) {
				&MKDIR("$DEG_od/BMK_1_All_DEG");
				&runOrDie( "cp $Bin/BMK_6_DEG_Analysis/BMK_1_All_DEG/readme.pdf $DEG_od/BMK_1_All_DEG");
				&runOrDie( "cp $dir/All* $DEG_od/BMK_1_All_DEG && mv $DEG_od/BMK_1_All_DEG/All.DEG_final.xls $DEG_od/BMK_1_All_DEG/All.DEG.Expression.xls");
			}
			elsif ($dir=~/\/DEG_PPI$/) {
                		&runOrDie( "cp  $dir/ppi_qurey.ppi.cytoscapeInput.sif $DEG_od/BMK_1_All_DEG");

			}
			elsif ($dir=~/\/density$/) {
				&MKDIR("$geneExpression_od/BMK_4_density");
				&runOrDie( "cp $dir/*fpkm_density.png $dir/all.fpkm_box.png $geneExpression_od/BMK_4_density/");
				&MKDIR("$geneExpression_od/BMK_5_correlation");
				&runOrDie( "cp $dir/sample_cluster.png $geneExpression_od/BMK_5_correlation/cor.cluster.png && cp $dir/cor_plot/*.png $geneExpression_od/BMK_5_correlation/");
				&runOrDie("awk '{for(i=1;i<=NF/2+1;i++) {if(NR==1){if(i==1) printf(\"\%s\\t\%s\\t\",\"sample\",\$i);if(i>1 && i<NF/2+1) printf(\"\%s\\t\",\$i) } else {printf(\"\%s\\t\",\$i)}} printf \"\\n\"}' $dir/free_com.cor|sed 's/decor.//g' > $geneExpression_od/BMK_5_correlation/correlation.txt");
			}
			elsif ($dir=~/_vs_/){

				my $BMK_1 = "$DEG_od/BMK_$count"."_$nam/BMK_1_Statistics_Visualization";
				my $BMK_2 = "$DEG_od/BMK_$count"."_$nam/BMK_2_DEG_Annotation";
				my $BMK_3 = "$DEG_od/BMK_$count"."_$nam/BMK_3_GO_Enrichment";
				my $BMK_4 = "$DEG_od/BMK_$count"."_$nam/BMK_4_Pathway_Enrichment";
				my $BMK_5 = "$DEG_od/BMK_$count"."_$nam/BMK_5_PPI";
				my $DEG_file = (glob "$dir/*.DEG_final.xls")[0];
                                if(!$DEG_file){
                                  $Stat{$nam}{up}=0;$Stat{$nam}{down}=0;$Stat{$nam}{total}=0;
                                  &MKDIR($BMK_1);
                                  &runOrDie( "cp  $dir/*.final.xls   $dir/*.FC_count.png   $dir/*.FC_FDR.png $dir/Anno_enrichment/*.annotation.xls  $Bin/BMK_6_DEG_Analysis/BMK_1_Statistics_Visualization/readme.pdf $BMK_1  && mv $BMK_1/$nam.final.xls $BMK_1/$nam.all.gene.xls ");
                                  $Stat{$nam}{up}=0;$Stat{$nam}{down}=0;$Stat{$nam}{total}=0;
}else{
                                &MKDIR($BMK_1);
				&runOrDie( "cp $dir/*.cluster.* $dir/*.final.xls $dir/*.DEG_final.xls $dir/*.FC_count.png   $dir/*.FC_FDR.png $dir/*.FC_FDR.png $dir/Anno_enrichment/*.annotation.xls $Bin/BMK_6_DEG_Analysis/BMK_1_Statistics_Visualization/readme.pdf $BMK_1 && mv $BMK_1/$nam.DEG_final.xls $BMK_1/$nam.DEG.xls && mv $BMK_1/$nam.final.xls $BMK_1/$nam.all.gene.xls ");

				&MKDIR($BMK_2);
				&runOrDie( "cp $dir/Anno_enrichment/Kog_Anno/*.classfy.* $dir/Anno_enrichment/Cog_Anno/*.classfy.* $Bin/BMK_6_DEG_Analysis/BMK_2_DEG_Annotation/readme.pdf $BMK_2");

				&MKDIR($BMK_3);
				&runOrDie( "cp $dir/Anno_enrichment/go_enrichment/*.GO.png  $dir/Anno_enrichment/Graph/topGO.map $Bin/BMK_6_DEG_Analysis/BMK_3_GO_Enrichment/readme.pdf $BMK_3 && mv $BMK_3/topGO.map $BMK_3/GO.map");

				&MKDIR("$BMK_3/$nam.BP");
				&runOrDie( "cp $dir/Anno_enrichment/Graph/*.topGO_BP* $BMK_3/$nam.BP");

				&MKDIR("$BMK_3/$nam.CC");
				&runOrDie( "cp $dir/Anno_enrichment/Graph/*.topGO_CC* $BMK_3/$nam.CC");

				&MKDIR("$BMK_3/$nam.MF");
				&runOrDie( "cp $dir/Anno_enrichment/Graph/*.topGO_MF* $BMK_3/$nam.MF");

                                my @mapfile = glob "$dir/Anno_enrichment/pathway/kegg_map/*";
				if(@mapfile){
                                &MKDIR("$BMK_4/KEGG_map");
				&runOrDie( "cp $Bin/BMK_6_DEG_Analysis/BMK_4_Pathway_Enrichment/readme.pdf $BMK_4");
				&runOrDie( "cp $dir/Anno_enrichment/pathway/kegg_map/* $BMK_4/KEGG_map");}

                                my @enrichfile = glob "$dir/Anno_enrichment/pathway/kegg_enrichment/*";
				if(@enrichfile){
                                &MKDIR("$BMK_4/KEGG_enrichment");
				&runOrDie( "cp $dir/Anno_enrichment/Graph/*.KEGG.list $dir/Anno_enrichment/Graph/*.KEGG.Phase.png $BMK_4/KEGG_enrichment && mv $BMK_4/KEGG_enrichment/$nam.KEGG.Phase.png $BMK_4/KEGG_enrichment/$nam.KEGG.enrichment_factor.png && cp $dir/Anno_enrichment/pathway/kegg_enrichment/*.Kegg.ko $dir/Anno_enrichment/pathway/kegg_enrichment/*.KEGG.png $dir/Anno_enrichment/pathway/kegg_enrichment/*.KEGG.stat $BMK_4/KEGG_enrichment && cut -f1-4,7,8 $dir/Anno_enrichment/pathway/kegg_enrichment/$nam.KEGG.xls > $BMK_4/KEGG_enrichment/$nam.KEGG.xls");
				open I,"$BMK_4/KEGG_enrichment/$nam.KEGG.stat";
				my %pathway_enrich;
				while(<I>){
					next if /^\#/;
					chomp;
					my @tmp=split/\t/;
					$pathway_enrich{$tmp[1]}="$tmp[-2]\t$tmp[-1]";
				}
				close I;
				open I,"$BMK_4/KEGG_enrichment/$nam.KEGG.xls";
				open O,">$BMK_4/KEGG_enrichment/$nam.KEGG.xls.final";
				while(<I>){
					chomp;
					my @tmp=split/\t/;
					print O "$_\n" and next if($.==1);
					if($.==2){
						print O join("\t",@tmp[0..3])."\tP-value\tCorrected_P-value\t$tmp[4]\t$tmp[5]\n";
						next;
					}
					next unless exists $pathway_enrich{$tmp[1]};
					print O join("\t",@tmp[0..3])."\t$pathway_enrich{$tmp[1]}\t$tmp[4]\t$tmp[5]\n";
				}
				close I;close O;
				#&runOrDie(("grep -r \"#\" $BMK_4/KEGG_enrichment/$nam.KEGG.xls|sed 's/gene\$/gene\\tP-value\\tCorrected_P-value/' > $BMK_4/KEGG_enrichment/$nam.KEGG.xls.final"));

				#&runOrDie(("gawk -F \"\\t\" '{if(NR==FNR) A[NR]=\$1\"\\t\"\$5\"\\t\"\$6;else {for(i in A) {if(A[i]~\$1) print \$0\"\\t\"A[i]}}}' $BMK_4/KEGG_enrichment/$nam.KEGG.stat $BMK_4/KEGG_enrichment/$nam.KEGG.xls |gawk -F \"\\t\" '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$7\"\\t\"\$8\"\\t\"\$5\"\\t\"\$6}' >> $BMK_4/KEGG_enrichment/$nam.KEGG.xls.final"));

				&runOrDie(("rm $BMK_4/KEGG_enrichment/$nam.KEGG.xls $BMK_4/KEGG_enrichment/$nam.KEGG.stat && mv $BMK_4/KEGG_enrichment/$nam.KEGG.xls.final $BMK_4/KEGG_enrichment/$nam.KEGG.xls"));}

				&MKDIR($BMK_5);
				&runOrDie( "cp $dir/../DEG_PPI/$nam.ppi.cytoscapeInput.sif $BMK_5");
				my $file=(glob "$dir/*.DEG_final.xls")[0];
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
							if ($Anno[$s] eq 'eggNOG_class') {
                                                                $Site{'eggNOG'} = $s;
                                                        }
							if ($Anno[$s] eq 'Swissprot_annotation') {
								$Site{'Swiss-Prot'} = $s;
							}
							elsif ($Anno[$s]=~/^([^_]+)_annotation/) {
								$Site{$1}=$s;
							}
						}
					}
					else{
						my @Info=split /\t/,$_;
						foreach my $key (keys %Site) {
							$Anno{$nam}{$key}||=0;
							$Anno{$nam}{$key}++ unless ($Info[$Site{$key}] eq '--');
						}
                        			$Anno{$nam}{'Annotated'} ++;
					}
				}
				close IN;
                          }
				$count++;
		}
	}
        if(%Anno){
	open (OUT,">$DEG_od/BMK_1_All_DEG/All.DEG.anno.stat.xls") or die $!;
	my $limit_anno=0;
	foreach my $key (sort keys %Anno) {#!
		if ($limit_anno==0) {
			print OUT "#DEG_Set";
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
        mkdir "$DEG_od/BMK_1_All_DEG" unless -d "$DEG_od/BMK_1_All_DEG";
	open (OUT,">$DEG_od/BMK_1_All_DEG/All.DEG.stat.xls") or die $!;
	print OUT "DEG_Set\tAll_DEG\tup-regulated\tdown-regulated\n";
	foreach my $key (sort keys %Stat) {
		$Stat{$key}{up}||=0;
		$Stat{$key}{down}||=0;
		print OUT "$key\t$Stat{$key}{total}\t$Stat{$key}{up}\t$Stat{$key}{down}\n";
	}
	close OUT;
}
}
##### rename
&runOrDie("rename GO_enrichment.stat.xls GO.stat $unigeneAnno_od/BMK_2_statistic/*");
&rename_recursion($od,"Kegg","KEGG");
&rename_recursion($od,"Kog","KOG");
&rename_recursion($od,"Cog","COG");
&rename_recursion($od,".nr.",".NR.");
&rename_recursion("$od/BMK_6_DEG_Analysis/","classfy","classification");
################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
sub rename_recursion
{
	my($dir,$key,$replace)=@_;
	my %options;
	#目录路径
  	my @cases;
  	if (-d $dir) {#判断目录是否存在
		my @files;
    		my $dh;
    		push(@files, $dir);
    		while (@files) {
      			if (-d $files[0])
			{#若是目录执行以下操作
        			opendir $dh, $files[0] or die $!;#打开目录句柄,若失败打印错误信息
        			@_ = grep { /^[^\.]/ } readdir $dh;#过滤掉以"."和".."的文件,即UNIX下的隐藏文件
        			foreach (@_) {
        				if ($_ !~ /\./) {
          					push @files,"$files[0]/".$_."/";
          					push @cases,"$files[0]/".$_."/"; ##
          				}
				}
        			closedir $dh;
			}
			shift @files;
		}
	}else {
		@cases = ($dir);
	}


&runOrDie("rename $key $replace $_*") foreach @cases;#打印文件列表
}

#####################################
sub add_col
{
	my($in,$out)=@_;

}
####################################
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
	#`rm -r $dir` if(-d $dir);
	`mkdir -p $dir` if(!-d $dir);
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
