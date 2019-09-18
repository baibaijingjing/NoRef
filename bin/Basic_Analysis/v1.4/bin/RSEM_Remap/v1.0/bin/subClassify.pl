#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.2.0";
#######################################################################################

my ($in1,$in2,$out1,$out2,$MAX_MISMATCH,$MAX_INDEL,$MAX_LEN_OF_INTRO,$MIN_LEN_OF_INTRO,$MAX_INSERT_SIZE,$NOISE)=@ARGV;

$MAX_MISMATCH||=3;
$MAX_INDEL||=2;

$MAX_LEN_OF_INTRO||=8_000;
$MIN_LEN_OF_INTRO||=30;

$MAX_INSERT_SIZE||=250;

$NOISE||=10;

my %ID;
my %uniq;

my $total_mapped_reads_1;
my $satisfied_Max_Len_Intro_1;
my $satisfied_Maxmismatch_1;
my $satisfied_QueryIndel_1;
my $satisfied_MaxIndel_1;
my $satisfied_Intro_1;


open (IN,"<",$in1) or die $!;
while (<IN>) {
	chomp;
	my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts)=split(/\t/,$_);

	$total_mapped_reads_1++;
	next if ($tEnd-$tStart>$MAX_LEN_OF_INTRO) ;                                         # modified by mengf 2011-8-18;ֱ�ӹ��˵�reads�ȶԵ�Target������ $MAX_LEN_OF_INTRO �����У�
	$satisfied_Max_Len_Intro_1++;
	#�ж�û�бȶ��ϵļ�������������$MAX_MISMATCH,�͹��˵�
	my $meMisMatches=$qSize-$matches-$misMatches;
	next if ($meMisMatches>$MAX_MISMATCH) ;
	$satisfied_Maxmismatch_1++;

	#������Query Reads�����������Gap(��������query��������1�����ϵļ����INDEL:������Gap�Ĵ���,����������������������ϵ�Gap)
	next unless ($qBaseInsert==0 or $qBaseInsert==$qNumInsert) ;
	$satisfied_QueryIndel_1++;

	my @BlockSizes=split ",",$blockSizes;
	my @QStarts=split ",",$qStarts;
	my @TStarts=split ",",$tStarts;

	#�޸ı߽�,���ܶ˵�ȶԲ��ϵ�Reads��(query�ıȶ���ʼλ�ò��ǵ�һ���������qStart��tStart����������ĩβҲ����ͬ������������Ϊ��Read
	#�Ǵӿ�ʼ�ȵ�β����)
	if ($qStart>0) {
		$misMatches+=$qStart;
		$BlockSizes[0]++;
		$QStarts[0]-=$qStart;
		$TStarts[0]-=$qStart;
		$tStart-=$qStart;
		$qStart-=$qStart;
	}
	if ($qEnd<$qSize) {
		$misMatches+=$qSize-$qEnd;
		$BlockSizes[-1]++;
		$qEnd+=$qSize-$qEnd;
		$tEnd+=$qSize-$qEnd;
	}

	#�ж�InDel���
	my $InDel=0;
	my @Intro=();
	for (my $i=0;$i<$#BlockSizes ;$i++) {
		my $QInDel=$QStarts[$i+1]-($QStarts[$i]+$BlockSizes[$i]);
		my $TInDel=$TStarts[$i+1]-($TStarts[$i]+$BlockSizes[$i]);
		next if ($QInDel>1 ) ;
		next unless ($TInDel<=1 or $TInDel>=$MIN_LEN_OF_INTRO and  $TInDel<=$MAX_LEN_OF_INTRO ) ;
		$InDel++ if ($QInDel==1) ;
		$InDel++ if ($TInDel==1) ;
		if ($TInDel>=$MIN_LEN_OF_INTRO and  $TInDel<=$MAX_LEN_OF_INTRO ) {
			#�Ѵ�Intro�������������ڷ������á�
			push @Intro,$TStarts[$i]+$BlockSizes[$i],$TStarts[$i+1];
		}
	}
	next if ($InDel>$MAX_INDEL) ;
	$satisfied_MaxIndel_1++;
	
	#���һ��READ���ж��Intro���͹��˵���READ.
	next if (@Intro>2) ;
	$satisfied_Intro_1++;

	#���ڣ�Alignment�Ѿ����ˣ�MisMatch��InDel�Ѿ�����Ҫ��READ������Gap��READ����Ӧ1��Intro�������˱߽�ȶԲ��϶�������Read���֡�

	#2012-10-8 Adjust for Hiseq sequencing ID;
	my $seq_pre;
	if ($qName=~/\//) {
		$qName=~/^(\S+)\/(\d)$/;
		$seq_pre=$1;
	}
	else {
		$seq_pre=$qName;
	}
	my %ref;
	$ref{"qSize"}=$qSize;
	$ref{"strand"}=$strand;
	$ref{"tName"}=$tName;

	$ref{"tStart"}=$tStart;
	$ref{"tEnd"}=$tEnd;

	@{$ref{"Intro"}}=();
	push @{$ref{"Intro"}},@Intro;

	$ref{"meMisMatches"}=$meMisMatches;
	$ref{"meIndel"}=$InDel;
	push @{$ID{$seq_pre}{1}},\%ref;
}
close (IN);

my $total_mapped_reads_2;
my $satisfied_Max_Len_Intro_2;
my $satisfied_Maxmismatch_2;
my $satisfied_QueryIndel_2;
my $satisfied_MaxIndel_2;
my $satisfied_Intro_2;

open (IN,"<",$in2) or die $!;
while (<IN>) {
	chomp;
	my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts)=split(/\t/,$_);
	$total_mapped_reads_2++;

	next if ($tEnd-$tStart>$MAX_LEN_OF_INTRO) ;                                         # modified by mengf 2011-8-18;ֱ�ӹ��˵�reads�ȶԵ�Target������ $MAX_LEN_OF_INTRO �����У�
	$satisfied_Max_Len_Intro_2++;

	#�ж�MisMatch�����������$MAX_MISMATCH,�͹��˵�
	my $meMisMatches=$qSize-$matches-$misMatches;
	next if ($meMisMatches>$MAX_MISMATCH) ;
	$satisfied_Maxmismatch_2++;

	#������Reads����Gap(������1�����ϵļ����INDEL)
	next unless ($qBaseInsert==0 or $qBaseInsert==$qNumInsert) ;
	$satisfied_QueryIndel_2++;
	
	my @BlockSizes=split ",",$blockSizes;
	my @QStarts=split ",",$qStarts;
	my @TStarts=split ",",$tStarts;

	#�޸ı߽�,���ܶ˵�ȶԲ��ϵ�Reads��
	if ($qStart>0) {
		$misMatches+=$qStart;
		$BlockSizes[0]++;
		$QStarts[0]-=$qStart;
		$TStarts[0]-=$qStart;
		$tStart-=$qStart;
		$qStart-=$qStart;
	}
	if ($qEnd<$qSize) {
		$misMatches+=$qSize-$qEnd;
		$BlockSizes[-1]++;
		$qEnd+=$qSize-$qEnd;
		$tEnd+=$qSize-$qEnd;
	}

	#�ж�InDel���
	my $InDel=0;
	my @Intro=();
	for (my $i=0;$i<$#BlockSizes ;$i++) {
		my $QInDel=$QStarts[$i+1]-($QStarts[$i]+$BlockSizes[$i]);
		my $TInDel=$TStarts[$i+1]-($TStarts[$i]+$BlockSizes[$i]);
		next if ($QInDel>1 ) ;
		next unless ($TInDel<=1 or $TInDel>=$MIN_LEN_OF_INTRO and  $TInDel<=$MAX_LEN_OF_INTRO ) ;
		$InDel++ if ($QInDel==1) ;
		$InDel++ if ($TInDel==1) ;
		if ($TInDel>=$MIN_LEN_OF_INTRO and  $TInDel<=$MAX_LEN_OF_INTRO ) {
			#�Ѵ�Intro�������������ڷ������á�
			push @Intro,$TStarts[$i]+$BlockSizes[$i],$TStarts[$i+1];
		}
	}
	next if ($InDel>$MAX_INDEL) ;
	$satisfied_MaxIndel_2++;
	
	#���һ��READ���ж��Intro���͹��˵���READ.
	next if (@Intro>2) ;
	$satisfied_Intro_2++;

	#���ڣ�Alignment�Ѿ����ˣ�MisMatch��InDel�Ѿ�����Ҫ��READ������Gap��READ����Ӧ1��Intro�������˱߽�ȶԲ��϶�������Read���֡�

	#2012-10-8 Adjust for Hiseq sequencing ID;
	my $seq_pre;
	if ($qName=~/\//) {
		$qName=~/^(\S+)\/(\d)$/;
		$seq_pre=$1;
	}
	else {
		$seq_pre=$qName;
	}
	my %ref;
	$ref{"qSize"}=$qSize;
	$ref{"strand"}=$strand;
	$ref{"tName"}=$tName;

	$ref{"tStart"}=$tStart;
	$ref{"tEnd"}=$tEnd;

	@{$ref{"Intro"}}=();
	push @{$ref{"Intro"}},@Intro;

	$ref{"meMisMatches"}=$meMisMatches;
	$ref{"meIndel"}=$InDel;
	push @{$ID{$seq_pre}{2}},\%ref;
}
close (IN);

###########################################################################################
my $satisfied_filt=0;
my $satisfied_pair=0;
my $satisfied_Rmap=0;
my $skipped_strand=0;
my $skipped_Nosie=0;
my $skipped_ST=0;
my $skipped_MP=0;

#Classify the reads
open (RMAP,">",$out1) or die $!;
foreach my $id (keys %ID) {
	my @line;
	my $line;
	$satisfied_filt++;

	if ($ID{$id}{1} and $ID{$id}{2}) {
		$satisfied_pair++;
		foreach my $ref_read1 (@{$ID{$id}{1}}) {
			foreach my $ref_read2 (@{$ID{$id}{2}}) {
				if ($ref_read1->{"strand"} ne $ref_read2->{"strand"} and $ref_read1->{"tName"} eq $ref_read2->{"tName"}) {

					my $plus=$ref_read1->{"strand"} eq "+"?$ref_read1:$ref_read2;
					my $minus=$ref_read1->{"strand"} eq "-"?$ref_read1:$ref_read2;
					my $oneStrand=$ref_read1->{"strand"} eq "+"?"1":"0";

					my $PlusStart=$plus->{"tStart"};
					my $MinusStart=$minus->{"tStart"};
					my $PlusEnd=$plus->{"tEnd"};
					my $MinusEnd=$minus->{"tEnd"};

					#���˵������Ŀ�ʼλ��С��������ʼλ�õ�READS
					if ($MinusStart<$PlusStart) {
						$skipped_strand++;
						next;
					}
					my $Distance=$MinusStart-$PlusEnd;

					#���˵�1��READ��ȱ�Ĳ��ֱ���һ��READ���ǵ�READS��������������ڣ��϶���READS������Ե�READS
					#�˴����������READS���ಿ����

					if ($Distance <$MAX_INSERT_SIZE) {
						if (scalar @{$plus->{"Intro"}}==0 and scalar @{$minus->{"Intro"}}==0) {
							$uniq{$id}{"PE"}++;
							$line=$id."\tPE\t".$plus->{"tName"}."\t-\t".$PlusStart."\t".$PlusEnd."\t".$MinusStart."\t".$MinusEnd."\t-\t-\t".$oneStrand."\t".($plus->{"meMisMatches"}+$minus->{"meMisMatches"})."\t".($plus->{"meIndel"}+$minus->{"meIndel"});
							push @line,$id,$line;
							$satisfied_Rmap++;
						}elsif(scalar @{$plus->{"Intro"}}==2 and scalar @{$minus->{"Intro"}}==0){
							#Output this SO record which the plus strand was splited.
							$uniq{$id}{"SO"}++;
							my ($S,$E)=@{$plus->{"Intro"}};

							#���˵�1��READ��ȱ�Ĳ��ֱ���һ��READ���ǵ�READS��������������ڣ��϶���READS������Ե�READS,����Ϊ10bp.
							if (&N1($S+1,$E-1,$MinusStart,$MinusEnd)>$NOISE) {
#								print $id,"\t",&N1($S+1,$E-1,$MinusStart,$MinusEnd),"\n";
								$skipped_Nosie++;
								next;
							}

							$line=$id."\tSO\t".$plus->{"tName"}."\t1\t".$PlusStart."\t".$PlusEnd."\t".$MinusStart."\t".$MinusEnd."\t".$S."\t".$E."\t".$oneStrand."\t".($plus->{"meMisMatches"}+$minus->{"meMisMatches"})."\t".($plus->{"meIndel"}+$minus->{"meIndel"});
							push @line,$id,$line;
							$satisfied_Rmap++;
						}elsif(scalar @{$plus->{"Intro"}}==0 and scalar @{$minus->{"Intro"}}==2){
							#Output this SO record which the minus strand was splited.
							$uniq{$id}{"SO"}++;
							my ($S,$E)=@{$minus->{"Intro"}};

							#���˵�1��READ��ȱ�Ĳ��ֱ���һ��READ���ǵ�READS��������������ڣ��϶���READS������Ե�READS
							if (&N1($S+1,$E-1,$PlusStart,$PlusEnd)>$NOISE) {
#								print $id,"\t",&N1($S+1,$E-1,$PlusStart,$PlusEnd),"\n";
								$skipped_Nosie++;
								next;
							}

							$line=$id."\tSO\t".$plus->{"tName"}."\t2\t".$PlusStart."\t".$PlusEnd."\t".$MinusStart."\t".$MinusEnd."\t".$S."\t".$E."\t".$oneStrand."\t".($plus->{"meMisMatches"}+$minus->{"meMisMatches"})."\t".($plus->{"meIndel"}+$minus->{"meIndel"});
							push @line,$id,$line;
							$satisfied_Rmap++;
						}elsif(scalar @{$plus->{"Intro"}}==2 and scalar @{$minus->{"Intro"}}==2){
							#Output this SO record which both the plus and the minus strand was splited.
							$uniq{$id}{"ST"}++;
							my ($S,$E)=@{$plus->{"Intro"}};
							my ($Minus_S,$Minus_E)=@{$minus->{"Intro"}};
							if ($S!=$Minus_S or $E!=$Minus_E) {
								$skipped_ST++;
								next;
							}
							$line=$id."\tST\t".$plus->{"tName"}."\t-\t".$PlusStart."\t".$PlusEnd."\t".$MinusStart."\t".$MinusEnd."\t".$S."\t".$E."\t".$oneStrand."\t".($plus->{"meMisMatches"}+$minus->{"meMisMatches"})."\t".($plus->{"meIndel"}+$minus->{"meIndel"});
							push @line,$id,$line;
							$satisfied_Rmap++;
						}
					}else{
						#Output this MP record.
						if (scalar @{$plus->{"Intro"}}==2 or scalar @{$minus->{"Intro"}}==2) {
							$skipped_MP++;
							next;
						}
						if ($MinusStart-$PlusEnd<$MIN_LEN_OF_INTRO or $MinusStart-$PlusEnd>$MAX_LEN_OF_INTRO) {
							$skipped_MP++;
							next;
						}
						$uniq{$id}{"MP"}++;
						$line=$id."\tMP\t".$plus->{"tName"}."\t-\t".$PlusStart."\t".$PlusEnd."\t".$MinusStart."\t".$MinusEnd."\t-\t-\t".$oneStrand."\t".($plus->{"meMisMatches"}+$minus->{"meMisMatches"})."\t".($plus->{"meIndel"}+$minus->{"meIndel"});
						push @line,$id,$line;
						$satisfied_Rmap++;
					}
				}
			}
		}
	}

	while (@line) {
		my $id=shift @line;
		my $line=shift @line;
		my $cpe=exists $uniq{$id}{"PE"}?$uniq{$id}{"PE"}:0;
		my $cmp=exists $uniq{$id}{"MP"}?$uniq{$id}{"MP"}:0;
		my $cso=exists $uniq{$id}{"SO"}?$uniq{$id}{"SO"}:0;
		my $cst=exists $uniq{$id}{"ST"}?$uniq{$id}{"ST"}:0;
		print RMAP $line,"\t",$cpe+$cmp+$cso+$cst,"\t",$cpe,"\t",$cmp,"\t",$cso,"\t",$cst,"\n";
	}

}

close (RMAP);

open STAT,">$out2" || die;
print STAT "fq1:\nTotal_mapped_Reads\t$total_mapped_reads_1\nSatisfied_Max_Len_Intro\t$satisfied_Max_Len_Intro_1\nSatisfied_Maxmismatch\t$satisfied_Maxmismatch_1\n";
print STAT "Satisfied_QueryIndel\t$satisfied_QueryIndel_1\nSatisfied_MaxIndel\t$satisfied_MaxIndel_1\nSatisfied_Intro_num\t$satisfied_Intro_1\n";
print STAT "\n\n";
print STAT "fq2:\nTotal_mapped_Reads\t$total_mapped_reads_2\nSatisfied_Max_Len_Intro\t$satisfied_Max_Len_Intro_2\nSatisfied_Maxmismatch\t$satisfied_Maxmismatch_2\n";
print STAT "Satisfied_QueryIndel\t$satisfied_QueryIndel_2\nSatisfied_MaxIndel\t$satisfied_MaxIndel_2\nSatisfied_Intro_num\t$satisfied_Intro_2\n";
print STAT "\n\n";
print STAT "Rmap:\nSatisfied_filt\t$satisfied_filt\nSatisfied_pair\t$satisfied_pair\nSatisfied_Rmap\t$satisfied_Rmap\n";
print STAT "Skipped_strand\t$skipped_strand\nSkipped_Nosie\t$skipped_Nosie\nSkipped_ST\t$skipped_ST\nSkipped_MP\t$skipped_MP\n";

close STAT;

#############################
sub N1{
	#�󹫹����ϵĴ�С��N(1,10)=10;N(1,10,5,10)=6;N(1,10,5,10,8,10)=3;�����ǡ���ʼλ�ã�����λ�á� ���ظ���
	my $sum=0;
	my $L=shift @_;
	my $R=shift @_;
	while (@_) {

		my $s=shift @_;
		my $e=shift @_;

		$L=$L>$s?$L:$s;
		$R=$R<$e?$R:$e;
	}

	$sum=$R>=$L?($R-$L+1):0;
	return $sum;
}
