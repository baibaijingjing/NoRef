#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.1";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($Index,$fGFF,$orf,$l_len,$l_FPKM,$BLOCK_SIZE,$SEARCH_DEPTH);
GetOptions(
				"help|?" =>\&USAGE,
				"g:s"=>\$fGFF,
				"i:s"=>\$Index,
				"orf:s"=>\$orf,
				"l_len:s"=>\$l_len,
				"l_FPKM:s"=>\$l_FPKM,
				"BLOCK_SIZE:i"=>\$BLOCK_SIZE,
				"SEARCH_DEPTH:i"=>\$SEARCH_DEPTH,
				) or &USAGE;
&USAGE unless ($Index and $fGFF and $orf) ;

$l_len||=500;
$l_FPKM||=10;
$fGFF=&AbsolutePath("file",$fGFF);
$Index=&AbsolutePath("file",$Index);
my $TMR=&GetTMR("$Index.Mapped.stat.xls");
print STDERR "limit of length: $l_len\n";
print STDERR "limit of FPKM: $l_FPKM\n";
print STDERR "Gff file: $fGFF\n";
print STDERR "orf file: $orf\n";
print STDERR "Rmap: $Index.rmap\n";
print STDERR "Total Mapped Reads(TMR): $TMR\n";
print STDERR "Out file: $Index.geneExpression.xls\n\n";
#######################################################################################


$BLOCK_SIZE=1000;
$SEARCH_DEPTH||=30;


# Get gene info form GFF file.
my %Gff;
&getGff($fGFF,\%Gff);

# Cluster Genes in order to save time later.
my %ClusterGff;
&clusterGene(\%Gff,$BLOCK_SIZE,\%ClusterGff);


# Get the reads and calculate the coverage of genes.
&getCover($Index,$BLOCK_SIZE,\%ClusterGff);


# Calculate gene coverage.
&geneCoverage(\%ClusterGff,$TMR);


# Output the geneExpression
&Output(\%ClusterGff,$Index);

#randCheck
my %ST;
$/=">";
open (IN,$orf) or die $!;
<IN>;
while (<IN>) {
	chomp;
	my ($genename,$strand)=(split/\s+/,$_)[0,2];
	if ($strand=~/-/) {
		$ST{$genename}=1;
	}
}
$/="\n";
close IN;
my %randCheck;
&randCheck($Index,\%ClusterGff,\%ST);







#######################################################################################
&timeLog("$Script Done. Total elapsed time : ".time()-$BEGIN_TIME."s");
#######################################################################################

# ------------------------------------------------------------------
# Get gene info form GFF file.
# ------------------------------------------------------------------
sub getGff{
	my ($fGFF,$refGff)=@_;
	open (IN,"<",$fGFF) or die $!;
	while (<IN>) {
		chomp;
		next if (/^\#/);
		next if (/^$/);
		my ($chro,$class,$start,$end,$strand,$id)=(split /\t/,$_)[0,2,3,4,6,8];
		($start,$end)=sort {$a <=> $b} ($start,$end);
		if ($class=~/mRNA/) {
			$id=~/ID=(\S+)/;$id=$1;
			$refGff->{$chro}{$id}{"start"}=$start;
			$refGff->{$chro}{$id}{"end"}=$end;
			$refGff->{$chro}{$id}{"strand"}=$strand;
		}elsif($class=~/CDS/){
			$id=~/Parent=(\S+)/;$id=$1;
			push @{$refGff->{$chro}{$id}{"arrExon"}},$start,$end;
		}
	}
	close (IN);
}

# ------------------------------------------------------------------
# Cluster Genes in order to save time later.
# ------------------------------------------------------------------
sub clusterGene{
	#We cluster the genes based by the start position of the genes;
	my ($refGff,$SIZE,$refCluster)=@_;
	foreach my $chro (keys %{$refGff}) {
		foreach my $gene (keys %{$refGff->{$chro}}) {
			my $blockNum=int($refGff->{$chro}{$gene}{"start"}/$SIZE);
			$refCluster->{$chro}{$blockNum}{$gene}=$refGff->{$chro}{$gene};
		}
	}
}

# ------------------------------------------------------------------
# Get the reads and calculate the coverage of genes.
# ------------------------------------------------------------------
sub getCover{
	my ($Key,$SIZE,$refCluster)=@_;

	open (IN,"<","$Key.rmap") or die $!;
	while (<IN>) {
		chomp;
		my ($Id,$Mode,$Chro,$SO12,$PlusStart,$PlusEnd,$MinusStart,$MinusEnd,$Cut1,$Cut2,$OnePlus,$MisMatch,$InDel,$TotalHits,$PEhits,$MPhits,$SOhits,$SThits)=split /\t/,$_;

		if ($Mode eq "PE") {
			&belongGen($PlusStart,$MinusEnd,$Chro,$TotalHits,$SIZE,$refCluster);
		}elsif($Mode eq "MP"){
			&belongGen($PlusStart,$PlusEnd,$Chro,$TotalHits,$SIZE,$refCluster);
			&belongGen($MinusStart,$MinusEnd,$Chro,$TotalHits,$SIZE,$refCluster);
		}elsif($Mode eq "SO" or $Mode eq "ST"){
			&belongGen($PlusStart,$Cut1,$Chro,$TotalHits,$SIZE,$refCluster);
			&belongGen($Cut2,$MinusEnd,$Chro,$TotalHits,$SIZE,$refCluster);
		}
	}
	close (IN);
}

# ------------------------------------------------------------------
# Find which gene the reads belong to.
# ------------------------------------------------------------------
sub belongGen{
	my ($start,$end,$chro,$totalReadsNum,$SIZE,$refCluster)=@_;
	my $blockNum=int($start/$SIZE);

	#This reads may belong to this block or the prior N block
	for (my $i=0;$i<$SEARCH_DEPTH ;$i++) {
		if (exists $refCluster->{$chro}{$blockNum}) {
			foreach my $gene (keys %{$refCluster->{$chro}{$blockNum}}) {
				if ($refCluster->{$chro}{$blockNum}{$gene}{"start"}<=$start and $refCluster->{$chro}{$blockNum}{$gene}{"end"}>=$end or
					$refCluster->{$chro}{$blockNum}{$gene}{"start"}>$start and $refCluster->{$chro}{$blockNum}{$gene}{"start"}<$end or
					$refCluster->{$chro}{$blockNum}{$gene}{"end"}>$start and $refCluster->{$chro}{$blockNum}{$gene}{"end"}<$end
				){

					#gene count add one;
					if ($totalReadsNum==1) {
						$refCluster->{$chro}{$blockNum}{$gene}{"UniqReads"}++;
					}else{
						$refCluster->{$chro}{$blockNum}{$gene}{"MultiReads"}++;
					}

					#gene coverage adjust
					my $Lmax=$refCluster->{$chro}{$blockNum}{$gene}{"start"}>=$start?$refCluster->{$chro}{$blockNum}{$gene}{"start"}:$start;
					my $Rmax=$refCluster->{$chro}{$blockNum}{$gene}{"end"}<=$end?$refCluster->{$chro}{$blockNum}{$gene}{"end"}:$end-1;
					push @{$refCluster->{$chro}{$blockNum}{$gene}{"arrReads"}},$Lmax,$Rmax;
					return 0;
				}
			}
		}
		$blockNum--;
	}
}


# ------------------------------------------------------------------
# Calculate gene coverage.
# ------------------------------------------------------------------
sub geneCoverage{
	my ($refCluster,$TMR)=@_;
	foreach my $chro (keys %{$refCluster}) {
		foreach my $block (keys %{$refCluster->{$chro}}) {
			foreach my $gene (keys %{$refCluster->{$chro}{$block}}) {

				#calculate the sum len of exons
				my @arrExon=&cat(0,@{$refCluster->{$chro}{$block}{$gene}{"arrExon"}});
				my $sumLenExons=0;
				while (@arrExon) {
					my $l=shift @arrExon;
					my $r=shift @arrExon;
					$sumLenExons+=$r-$l;
				}
				$refCluster->{$chro}{$block}{$gene}{"sumLenExons"}=$sumLenExons;

				my @arrReads=&cat(0,@{$refCluster->{$chro}{$block}{$gene}{"arrReads"}});
				my $sum=0;
				while (@arrReads) {
					my $l=shift @arrReads;
					my $r=shift @arrReads;
					$sum+=$r-$l+1;
				}
				$refCluster->{$chro}{$block}{$gene}{"cov"}=$sum/($refCluster->{$chro}{$block}{$gene}{"end"}-$refCluster->{$chro}{$block}{$gene}{"start"}+1);

				@arrReads=@{$refCluster->{$chro}{$block}{$gene}{"arrReads"}};
				my $DepthSum=0;
				while (@arrReads) {
					my $s=shift @arrReads;
					my $e=shift @arrReads;
					next unless ($s<$e) ;
					$DepthSum+=$e-$s+1;
				}
				if ($DepthSum<0) {
					print STDERR "Warning:$gene 's depth is below 0, the progarm will change it to 0!\n";
					$DepthSum=0;
				}
				$refCluster->{$chro}{$block}{$gene}{"depth"}=$DepthSum/$sumLenExons;

				#calculate the FPKM of gene
				my $uniqPosiReads=exists $refCluster->{$chro}{$block}{$gene}{"UniqReads"}?$refCluster->{$chro}{$block}{$gene}{"UniqReads"}:0;
				my $mulitPosiReads=exists $refCluster->{$chro}{$block}{$gene}{"MultiReads"}?$refCluster->{$chro}{$block}{$gene}{"MultiReads"}:0;
				$refCluster->{$chro}{$block}{$gene}{"FPKM"}=1e9*($uniqPosiReads+$mulitPosiReads)/$sumLenExons/$TMR;
			}
		}
	}
}


# ------------------------------------------------------------------
# Output the geneExpression
# ------------------------------------------------------------------
sub Output{
	my ($refCluster,$Key)=@_;
	open (OUT,">","$Key.geneExpression.xls") or die $!;
	print OUT  "#Gene_ID\tLength\tDepth\tCoverage\tFPKM\tTotalReads\tUniqReads\tMultiReads\n";
	foreach my $chro (sort {$a cmp $b} keys %{$refCluster}) {
		foreach my $block (sort {$a <=> $b} keys %{$refCluster->{$chro}}) {
			foreach my $gene (sort {$a cmp $b} keys %{$refCluster->{$chro}{$block}}) {
				my ($strand,$start,$end,$length,$exonLen,$introLen,$depth,$coverage,$FPKM,$totalReads,$uniqPosiReads,$mulitPosiReads);
				$strand=$refCluster->{$chro}{$block}{$gene}{"strand"};
				$start=$refCluster->{$chro}{$block}{$gene}{"start"};
				$end=$refCluster->{$chro}{$block}{$gene}{"end"};
				$length=$end-$start;
				$strand=$refCluster->{$chro}{$block}{$gene}{"strand"};
				$uniqPosiReads=exists $refCluster->{$chro}{$block}{$gene}{"UniqReads"}?$refCluster->{$chro}{$block}{$gene}{"UniqReads"}:0;
				$mulitPosiReads=exists $refCluster->{$chro}{$block}{$gene}{"MultiReads"}?$refCluster->{$chro}{$block}{$gene}{"MultiReads"}:0;
				$totalReads=$uniqPosiReads+$mulitPosiReads;
				$coverage=$refCluster->{$chro}{$block}{$gene}{"cov"};
				$depth=$refCluster->{$chro}{$block}{$gene}{"depth"};
				$FPKM=$refCluster->{$chro}{$block}{$gene}{"FPKM"};
				$exonLen=$refCluster->{$chro}{$block}{$gene}{"sumLenExons"};
				$introLen=$length-$exonLen;
				print OUT $gene,"\t",join("\t",$length,$depth,$coverage,$FPKM,$totalReads,$uniqPosiReads,$mulitPosiReads),"\n";
			}
		}
	}
}


# ------------------------------------------------------------------
# RandCheck
# ------------------------------------------------------------------
sub randCheck {#
	my ($Key,$refCluster,$st)=@_;
	my %refRandCheck;
	my %Extreme;
	my %Low;
	open (OUT,">","$Key.Extreme.list") or die $!;
	my $Gene_num=0;
	foreach my $chro (keys %{$refCluster}) {
		foreach my $block (keys %{$refCluster->{$chro}}) {
			LINE:foreach my $gene (keys %{$refCluster->{$chro}{$block}}) {
				my %GeneCheck;
				my $Gene_total=0;
				my $len=$refCluster->{$chro}{$block}{$gene}{"end"}-$refCluster->{$chro}{$block}{$gene}{"start"}+1;
				for (my $i=0;$i<@{$refCluster->{$chro}{$block}{$gene}{"arrReads"}} ;$i+=2) {
#					$refRandCheck{int(100*(($refCluster->{$chro}{$block}{$gene}{"arrReads"}[$i]-$refCluster->{$chro}{$block}{$gene}{"start"})/$len))}++;
					my $j=$i+1;
					my $ZXC_MID=($refCluster->{$chro}{$block}{$gene}{"arrReads"}[$i]+$refCluster->{$chro}{$block}{$gene}{"arrReads"}[$j])/2;
					$GeneCheck{int(100*(($ZXC_MID-$refCluster->{$chro}{$block}{$gene}{"start"})/$len))}++;
					$Gene_total++;
				}
				if($ClusterGff{$chro}{$block}{$gene}{"FPKM"}<$l_FPKM){
					$Low{$gene}=1;
					next LINE;
				}
				if($len<$l_len){
					next LINE;
				}
				for (my $zxc_i=int(10000/$len);$zxc_i<=int(($len-100)*80/$len) ;$zxc_i++) {
					my $limit=0;
					for (my $zxc_j=0;$zxc_j<=int (($len-100)*20/$len) ;$zxc_j++) {
						my $zxc_now=$zxc_i+$zxc_j;
						$GeneCheck{$zxc_now}||=0;
						$limit+=$GeneCheck{$zxc_now};
					}
					if ($limit>$Gene_total*0.8) {
						$Extreme{$gene}=1;
						print OUT "$gene\n";
						next LINE;
					}
				}
				$Gene_num++;
				foreach my $key (keys %GeneCheck) {
					my $st_ne=100-$key;
					$refRandCheck{$key}+=$GeneCheck{$key}/$Gene_total unless exists $st->{$gene};
					$refRandCheck{$st_ne}+=$GeneCheck{$key}/$Gene_total if exists $st->{$gene};
				}
			}
		}
	}
	close OUT;
	
	my ($max)=sort {$b <=> $a} values %refRandCheck;
	$max=sprintf "%.2f",100*$max/$Gene_num;
	my ($y,$ystep)=&define_axis ($max);
	open (OUT,">","$Key.randcheck.list") or die $!;

	print OUT "Type:Line","\n";
	print OUT "Width:600","\n";
	print OUT "Height:400","\n";
	print OUT "WholeScale:0.9","\n";
	print OUT "FontSize:25","\n";
	print OUT "X:Relative Position in Genes(5'-3')","\n";
	print OUT "Y:Percent of Reads","\n";
	print OUT "XStep:10","\n";
	print OUT "YStep:",$ystep,"\n";
	print OUT "XStart:0","\n";
	print OUT "YStart:0","\n";
	print OUT "XEnd:100","\n";
	print OUT "YEnd:",$y,"\n\n";
	print OUT "Color:red","\n";

	foreach my $posi (sort {$a <=> $b }keys %refRandCheck) {
		my $per=sprintf "%.2f",100*$refRandCheck{$posi}/$Gene_num;
		print OUT $posi,":",$per,"\n";
	}
	close (OUT);
	my ${distributing_svg}=${&readconf("$Bin/../../../../../../config/sys.cfg")}{"distributing_svg"};
	system "$distributing_svg $Key.randcheck.list $Key.randcheck.svg";
	my $bn=basename $Index;
	my $svg2xxx=${&readconf("$Bin/../../../../../../config/sys.cfg")}{"svg2xxx"};
	`$svg2xxx $bn.randcheck.svg`;

	open (IN,"<","$Key.rmap") or die $!;
	open (OUT,">","$Key.Extreme.rmap") or die $!;
	while (<IN>) {
		my $G_now=(split/\t/,$_)[2];
		print OUT $_ if exists $Extreme{$G_now};
	}
	close IN;
	close OUT;


}

sub define_axis () {
        my $i=0;
        my ($max)=@_;
        my $time=1;
        my @ret=();
        my @limit=(1,2,3,4,5,6,8,10,12,15,16,20,24,25,30,40,50,60,80,100,120);
        my @unlim=(0,1,2,3,4,5,6,8 ,10,11,14,15,18.5,21,23,29,37,47,56,76 ,92 ,110);

        while ($max >$unlim[21]) {
                 $max=$max/10;
                 $time=$time*10;
        }
        for ($i=0;$i<=20 ;$i++) {
                 if ($max>$unlim[$i] && $max<=$unlim[$i+1]) {
                         $ret[0]=$limit[$i]*$time;

                         if ($i==2 || $i==5 || $i==9 || $i==14) {
                                 $ret[1]=$ret[0]/3;
                         }
                         elsif ($i==4 || $i==7 || $i==13 || $i==16){
                                 $ret[1]=$ret[0]/5;
                         }
                         else {
                                 $ret[1]=$ret[0]/4;
                         }

                 }
        }
        return @ret;
}

sub cat
#function:quit redundance
#input:($para,@array), para is the merge length
#output:(@array),
#for example (0,1,3,4,7,5,8)->(1,3,4,8) (1,1,3,4,7,5,8)->(1,8)
{
	my($merge,@input) = @_;
	my $i = 0;
	my @output = ();
	my %hash = ();
	my $each = 0;
	my $begin = "";
	my $end = 0;


	for ($i=0;$i<@input;$i+=2)
	{
		my $Qb = $input[$i];
		my $Qe = $input[$i+1];

		if($Qb > $Qe) { next; }
		if(defined($hash{$Qb}))	{ if($hash{$Qb} < $Qe) { $hash{$Qb} = $Qe; } }
		else { $hash{$Qb} = $Qe; }
		$Qb = 0;
	}

	foreach $each (sort {$a <=> $b} keys %hash)
	{
		if($begin eq "")
		{
			$begin = $each;
			$end = $hash{$each};
		}
		else
		{
			if($hash{$each} > $end)
			{
				if($each > $end + $merge)
				{
					push(@output,$begin);
					push(@output,$end);
					$begin = $each;
					$end = $hash{$each};
				}
				else { $end = $hash{$each}; }
			}
		}
	}
	if(keys %hash > 0)
	{
		push(@output,$begin);
		push(@output,$end);
	}

	%hash = ();

	return(@output);
}

sub AbsolutePath
{		#获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/^Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: GeneExpression.pl (Calculate the gene expression profile using pe, mp, so, st, gff files)
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>

Description:
	This program is a memeber of Transcript Analysis Kit, it can Calculate the gene expression profile
	using pe, mp, so, st reads and gff file. the Total Mapped Reads(TMR) which can be found in
	.Mapped.stat Info file generated by Rmap.pl also should be give in commond line.

	The program will calculate the FPKM, Depth, Coverage, Length ... of each genes.

Usage:
  -g               <file>   Gff File, forced;
  -i               <str>    Index of inFiles and outFiles. Index.rmap, Index.geneExpression.xls, Index.randCheck
  -orf             <file>   orf result,*.Unigene.cds_pep.stat.xls,forced
  -l_len           <num>    limit of length,default 500
  -l_FPKM          <num>    limit of FPKM,default 10

  -BLOCK_SIZE      <int>    Block Size, used to save time, default 1000;
  -SEARCH_DEPTH    <int>    Search Depth, used to save time with BLOCK_SIZE, default 30;


USAGE
	print $usage;
	exit;
}
