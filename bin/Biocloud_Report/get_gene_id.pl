#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;
use newPerlBase;
my $USAGE = qq{
Name:
    $0
Function:
    1.get all deg 
Usage:
    perl $0 -inputdir /share/bioCloud/cloud/project_report/2014-09-22/BMK140701-E16_Litchi/BMK140701-E16_Litchi -out /share/bioCloud/cloud/project_report/2014-09-22/BMK140701-E16_Litchi
Options:
    -inputdir    <string>    input webreport  dir
    -out    <string>    output path;
    -oneSam
Author:

    huangls
Version:
    v1.0;2014-10-27
Modifier: baij;2017-02-04
};




my ($deg_dir,$input_dir,$cfg,$fout,$all_deg,@deg,$oneSam);

GetOptions(
				"cfg:s" => \$cfg,
				"out:s"=>\$fout,
				"inputdir:s"=>\$input_dir,
				"oneSam"=>\$oneSam,
				
				) ;
die $USAGE unless ( $cfg and $fout);
if(-d $fout){}else{mkdir "$fout";}
$fout="$fout/no_ref_trans_all_table.xls";
my %config=%{readconf("$Bin/../../config/sys.cfg")};

my %para;
$cfg=&ABSOLUTE_DIR($cfg);
&para_load($cfg,\%para);

my %DEG_para;
&para_load("$para{DEG}/../Config/DEGAnalysis.cfg",\%DEG_para);
my $FDR=$DEG_para{FDR};
my $FC=$DEG_para{fold};
my $Method_RE=$DEG_para{Method_RE};
my $Method_NR=$DEG_para{Method_NR};


	$deg_dir="$para{DEG}";
	$all_deg="$para{DEG}/All_gene_expression.list";#All_gene_fpkm.list
	if(-e $all_deg){
		open INALL,"$all_deg" or die "$!:$all_deg";
		#print STDOUT "aaa:$all_deg\n";
	}elsif(-e "$deg_dir/All_gene_expression.list"){
		$all_deg="$deg_dir/All_gene_expression.list";
		#print STDOUT "bbb:$all_deg\n";
		open INALL,"$all_deg" or die "$!:$all_deg";
	}else{
		$all_deg="$deg_dir/All_geneExpression.list";
		#print STDOUT "bbb:$all_deg\n";
		open INALL,"$all_deg" or die "$!:$all_deg";
		
	}
	my %id_all;
	my ($head_line1,@head_line2);
	while (my $line=<INALL>){
		chomp $line;
		if ($line=~/^#/){
			my @colNum = split(/\t/,$line);
			my @col;
			for (my $i=1;$i< @colNum-1 ;$i++){push @col,"\t";}
			$head_line2[0]="$line";
			$head_line1="#All_gene_expression\t".join(" ",@col)."\t";
			next;
		}
		my @tmp=split(/\t/,$line);
		my $tmp="@tmp[1..$#tmp]";
		$id_all{$tmp[0]}=$line;
		#if (($tmp=~s/0//g)==$#tmp){
		#	next;
		#}else{
		#	$id_all{$tmp[0]}=$line;
		#}
	}
	my %hash_files;
	my $ann = "$para{ANNO}/Integrated_Function.annotation.xls";
unless($oneSam){	
	foreach my $dir (glob "$deg_dir/*_vs_*") {
	    my $deg = basename($dir);
            my $software;
	    my ($group1, $group2) = split(/_vs_/,$deg);
            if($group1=~/_/ || $group2=~/_/){$software=$Method_RE;}else{$software=$Method_NR;}
            my $degname="${deg}_${software}";
####Modify the venn licq 20160616
            my $final = "$dir/$deg.final.xls";
            my $deg_final = (-f "$dir/$deg.DEG.final.xls") ? "$dir/$deg.DEG.final.xls" :
                                        (-f "$dir/$deg.DEG_final.xls") ? "$dir/$deg.DEG_final.xls" :
                                        'unknown';
            
            if (-f $final && -f $deg_final) {
                push @deg,$deg;
                # get deg ids
                my %deg_ids;
                open DEG,$deg_final or die $!;
                while (<DEG>) {
                    chomp;
                    next if (/^#/ or /^\s*$/);
                    my $id = (split /\t/)[0];
                    $deg_ids{$id} = 1;
                }
                close DEG;

                #switch "up" or "down" to "normal" if needed
                open IN1,$final or die $!;
                while (my$line=<IN1>) {
                    chomp $line;
                    my @tmp=split(/\t/,$line);
                    my$c=pop @tmp;
                    my$b=pop @tmp;
                    my$a=pop @tmp;
                    if ($line=~/^#/ ){
                            #$head_line2=$head_line2. "${deg}_$a\t${deg}_$b\t${deg}_$c";
                            push @head_line2,"${degname}_$a\t${degname}_$b\t${degname}_(FDR_${FDR}_FC_${FC})_$c";
                            #$head_line1=$head_line1."$tmp[1]_vs_$tmp[2]\t \t \t";
                            next;
                    }
                    $c = 'normal' unless (exists $deg_ids{$tmp[0]});
                    $hash_files{$deg}{$tmp[0]}="$a\t$b\t$c";
                }
            close IN1;

            } elsif (!-f $final && -f $deg_final) {
                push @deg,$deg;
                open IN1,$deg_final or die $!;
                while (my$line=<IN1>) {
                    chomp $line;
                    my @tmp=split(/\t/,$line);
                    my$c=pop @tmp;
                    my$b=pop @tmp;
                    my$a=pop @tmp;
                    if ($line=~/^#/ ){
                            #$head_line2=$head_line2. "${deg}_$a\t${deg}_$b\t${deg}_$c";
                            push @head_line2,"${degname}_$a\t${degname}_$b\t${degname}_(FDR_${FDR}_FC_${FC})_$c";
                            #$head_line1=$head_line1."$tmp[1]_vs_$tmp[2]\t \t \t";
                            next;
                    }
                    $hash_files{$deg}{$tmp[0]}="$a\t$b\t$c";
                }
            close IN1;
            } else {
                print STDOUT "file $dir/$deg.DEG_final.xls or $dir/$deg.DEG.final.xls not exsist.\n";
            }
	
	}
 }
	ANN:
	open INANN,"$ann" or die "$!";
	my $anncolNum;
        my $go_anno_index;
        my %go_link_ref;
        my %go_link;
	while (my$line=<INANN>){
		 chomp$line;
	        my @tmp=split(/\t/,$line);
	         if ($line=~/^#/ || $line=~/^\s*$/){
				$anncolNum=@tmp;
                                for my $col_index(0..$#tmp){
 				    if ($tmp[$col_index] eq "GO_annotation") {
                                                   $go_anno_index=$col_index;
                                                   last;
                                    }
				}
	            push @head_line2,join("\t",@tmp[1..$#tmp]);
	            #$head_line1=$head_line1."Ann";
	            next;
	         };
	        $hash_files{basename($ann)}{$tmp[0]}=join("\t",@tmp[1..$#tmp]);
                my $go_anno=$tmp[$go_anno_index];
                if ($go_anno eq "--") {
                        $go_link_ref{$tmp[0]}="--";
                        next;
                }
                my @go_ids=($go_anno=~/(GO\:\d+)/g);
                for my $go_id(@go_ids) {
                        push @{$go_link{$tmp[0]}},$go_id;
                }
	}
        close (INANN);
	my@annNull;
	for (2..$anncolNum) {
		push @annNull,"--";
	}

       #by baij time 2017-2-7 Add GO second level annotation
        my $component_file =$config{component_file};
        my $function_file =$config{function_file};
        my $process_file =$config{process_file};
        my @go_tree_file = ($component_file, $function_file, $process_file);
        my %go_sub_id;
        my %go_tree_hash;
        foreach my $file (@go_tree_file) {
                open IN,$file or die "cannot open file $file, $!\n";
                my @array = ();
                my $name;
                my $key;
                while (my $line = <IN>) {
                        chomp($line);
                        next if ($line =~ /^$|^\#/);
                        next if ($line !~ /(^\s+>)|(^\s+%)|(^\s+<)/);
                        if ($line =~ /^\s>(.*)\s;/) {
                                $name = $1;
                                $name =~ s/[_\-]/ /;
                                next;
                        }
                        if ($line=~/^\s\s%/) {
                                my @key_in=split(/\s+;\s+/,$line);
                                $key_in[0]=~s/^\s\s%//;
                                $key = $key_in[0];
                                my ($go_id) = $key_in[1] =~ /(GO:\d+)/;
                                push @{$go_tree_hash{$name}{$go_id}},$key;
                                if (exists $go_sub_id{$name}{$key} and $go_sub_id{$name}{$key} ne $go_id) {
                                        print "two go id for $name:$key\n";
                                }
                                $go_sub_id{$name}{$key}=$go_id;
                        } else {
                                my ($go_id) = ($line =~ /(GO:\d+)/);
                                next unless (defined $go_id);
                                push @{$go_tree_hash{$name}{$go_id}},$key;
                        }
                }
                close IN;
        }

        foreach my $gene(keys %go_link){
                $go_link_ref{$gene}="";
                my %record;
                my @go_ids=@{$go_link{$gene}};
                my $mark=0;
                foreach my $go_id(@go_ids){
                        my $mark=0;
                        foreach my $main_term(keys %go_tree_hash){
                                if (exists $go_tree_hash{$main_term}{$go_id}) {
                                        $mark=1;
                                        my @sub_items=@{$go_tree_hash{$main_term}{$go_id}};
                                        foreach my $sub_item (@sub_items){
                                                next if (exists $record{$main_term}{$sub_item});
                                                $record{$main_term}{$sub_item}=1;
                                                $go_link_ref{$gene}.="$main_term: $sub_item ($go_sub_id{$main_term}{$sub_item});; ";
                                        }
                                }
                        }
                        if ($mark==0) {
                                print "$go_id has no secondary categary\n";
                        }
                }
                if($go_link_ref{$gene}) {
                        $go_link_ref{$gene}=~s/;;\s+$//;
                } else {
                        $go_link_ref{$gene}="--";
                }
        }
        push @head_line2,"GO_second_level_annotation" ;

	###by xugl time:2016-10-09 添加KEGG_pathway注释
	my %pathway_ko;
	my @pathwayFile=glob("$para{ANNO}/*Kegg.pathway");
	for my $file (@pathwayFile){
		open P,"$file" or die "$!";
		while(<P>){
			chomp;
			next if (/^\s*$/ || /^\#/);
			my ($pathway,$ko,$gene_id) = (split /\t+/)[0,1,3];
			$gene_id=~s/;$//;
			my @gene_ids = split /;/,$gene_id;
			my $current_ko_pathway = $pathway." ($ko)";
			foreach my $id (@gene_ids)
			{
				if (exists $pathway_ko{$id}) {
					my $ko_pathway = $pathway_ko{$id};
					if ($current_ko_pathway !~ /$ko_pathway/) {
						$pathway_ko{$id}.=";; ".$current_ko_pathway;
						}
					}
				else
				{
					$pathway_ko{$id} = $current_ko_pathway;
				}
			}
		}
		close P;
	}
        push @head_line2, "KEGG_pathway_annotation";
        push @deg,basename($ann);



        #by baij time: 2017-2-7 Add SSR Annotation
        my %SSR;
        my $SSR_file=(glob("$para{Basic}/Cluster/Unigene/Unigene_SSR/*SSR.result.xls"))[0];
	open S,"$SSR_file" or die "$!\n";
	my $SSR_col_num;
	my @SSR_header_blank;
        while(<S>){
             next if /^\s*$/;
             chomp;
             my @tmp=split /\t+/;
             if (/^\#/){
                    for(my $i=0;$i<$#tmp;$i++){
                           $SSR_header_blank[$i]="--";
                           $SSR_col_num=$#tmp;
                            }
                       next;
                     }
             for my $i (1..$SSR_col_num){
                      $tmp[$i]="--" if $i>$#tmp;
                     }
             my @t=@tmp[1..$SSR_col_num];
             $SSR{$tmp[0]}=\@t;
          }
        close S;
        push @head_line2,qw(SSR_SSR_nr SSR_SSR_type SSR_SSR SSR_Size SSR_SSR_Start SSR_SSR_End SSR_FPr1(5'-3') SSR_Tm_1F SSR_Size_1F SSR_RPr1(5'-3') SSR_Tm_1R SSR_Size_1R SSR_PSize1 SSR_PStart1 SSR_PEnd1F SSR_Pr2(5'-3') SSR_Tm_2F SSR_FSize_2F SSR_RPr2(5'-3') SSR_Tm_2R SSR_Size_2R SSR_PSize2 SSR_PStart2 SSR_PEnd2 SSR_FPr3(5'-3') SSR_Tm_3 SSR_FSize_3F SSR_RPr3(5'-3') SSR_Tm_3R SSR_Size_3R SSR_PSize3 SSR_PStart3 SSR_PEnd3) ;    
	#print OUT "$head_line1\n$head_line2\n";
        open OUT,">$fout" or die "$!";
	print OUT join("\t",@head_line2),"\n";
#	push @deg,basename($ann);
	foreach my $key  (keys %id_all){
		my @aa=();
		push @aa,$id_all{$key};
	    foreach my $file (@deg){
	    	if (exists ($hash_files{$file}{$key})){
	    		push @aa, $hash_files{$file}{$key};
	    	}else{
	    		if ($file eq basename($ann)){
	    			push @aa,join("\t",@annNull);
	    		}else{
	    			push @aa,"--\t--\t--";
	    		}
	    		
	    	}
	    	
	    }
#	    my $out=join("\t",@aa);
            if (exists $go_link_ref{$key}) {
                        push @aa,$go_link_ref{$key};
                } else {
                        push @aa,"--";
                }
	   
	    if (exists $pathway_ko{$key}){
			push @aa, $pathway_ko{$key};
		}
		else{
			push @aa, "--";
		}
            if(exists $SSR{$key}){
                        push @aa, @{$SSR{$key}};
               }
               else{
     			push @aa,@SSR_header_blank;
               }
            my $out = join("\t",@aa);
	    print OUT "$out\n";
	}
	close OUT;



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
		warn "Warning just for file and dir $in\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
