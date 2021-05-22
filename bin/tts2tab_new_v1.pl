#!/usr/bin/env perl 
# First created by SPMane on 11/15/2020
# parses CDS.tab and transtermHP output and creates regulartory tab file with TTS.
# terminator locus_tag is matched to the upstream gene it terminates.
# Only prints terminators within 500 bp downstream of a gene.

use strict;
use warnings;
use Getopt::Long;
#my $ttsfile='/lrlhps/users/c240616/P_lentus/final/assembly/1240/PLEN1240_g_11232016/PLEN1240_g_11232016.tt';
#my $cdsfile='/lrlhps/users/c240616/P_lentus/final/assembly/1240/PLEN1240_g_11232016/gbk/cds.tab';
#my $ttsfile='/lrlhps/users/c240616/P_lentus/final/assembly/161/PLEN161_g_11232016/PLEN161_g_11232016.tt';
#my $cdsfile='/lrlhps/users/c240616/P_lentus/final/assembly/161/PLEN161_g_11232016/gbk/cds.tab';

#my $ttsfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU170331_09182017/LREU170331_09182017.tt';
#my $cdsfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU170331_09182017/gbk/cds.tab';

#my $ttsfile='/lrlhps/users/c240616/projects/Streptomyces/Assembly/SA73515/SA73515.tt';
#my $cdsfile='/lrlhps/users/c240616/projects/Streptomyces/Assembly/SA73515/gbk/cds.tab';

#my $cdsfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU3632_09182017/gbk/new/cds.tab';
#my $ttsfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU3632_09182017/LREU3632_09182017.tt';
GetOptions( 'cdsfile=s' => \my $cdsfile,
            'ttsfile=s' => \my $ttsfile,
          , 'tag=s' => \my $tag  # same for --user or -u
          );
#print "$cdsfile\n";
#print "$ttsfile\n";
#print "$tag\n";
#my $cdsfile='gbk/cds.tab';
#my $ttsfile='PTA-126787.tt';

my %contigs;
my %postiveEnd2LT;
my %negativeEnd2LT;
open my $cfh, '<', $cdsfile or die $!;
my $cnt=0;
while(<$cfh>){
	chomp;
	next if /^\#/;
	next if /^chromosome/;
	$cnt++;
	my @col= split /\t/;
	my $cont=$col[0];my $s=$col[2];my $e=$col[3];
	$contigs{$cont}++;
	#($e,$s)=($s,$e) if($col[4] == -1);
	#print "CDS: $cont\t$s\t$e\t$col[7]\n";
	#$coord2locustag{"$cont\t$s\t$e"}="$col[7]($col[4])";
	if($col[4] == 1){
		$postiveEnd2LT{$cont}{$e}=$col[7];
	}else{
		$negativeEnd2LT{$cont}{$s}=$col[7];
		#print "$cont\t$s = $col[7]\n";
	}
	#print "$cont\t$e\t$col[4] = $col[7]\n";
	#last if $cnt >=10;
}
close $cfh;
#exit;

my %codes;
$codes{'G'} = "Terminator in the interior of a gene (at least 50bp from an end)";
$codes{'F'} = "Terminator between two +strand genes";
$codes{'R'} = "Terminator between two -strand genes";
$codes{'T'} = "Terminator between a +strand gene and a -strand gene";
# $codes{'N'} = "none of the above";
my $det=" 5' tail, 5' stem, loop, 3' stem, 3' tail: ";

my @geneStart;
my $lastGene="";
my $gene='-';
my $operon='-';



print join("\t", qw/chromosome FeatureType start end strand regulatory_class locus_tag operon note/),"\n";



my $featureCnt=0;
my $chrom='-';

open my $tfh, '<', $ttsfile or die $!;
while(<$tfh>){
	chomp;
	if(/$tag\S+/){ #******************IMPORTANT******************************
		s/ - /\t/;
		s/ +/\t/g;
		my @col= split /\t/;
		$chrom=$col[0];
		my $gen='-';
		#$gen=$start2gene{$col[1]} if exists $start2gene{$col[1]};
		#print "****$col[0]\t$col[1]\t**$col[2]\n";
		#print scalar(@col)," *** ", join(" ",@col),"\n";
		push @geneStart,$col[1];
		$gene='-';
		$operon='-';
		next;
	}
	if (/ TERM/){
		$featureCnt++;
		s/TERM\s+/TERM/;
		s/ - /\t/;
		s/^ +//;
		s/ +/\t/g;
		#next;
		my @col= split /\t/;
		my $line=<$tfh>;
		chomp($line);
		#print join("\t",@col),"\n";
		my $start=$col[1];
		my $end=$col[2];
		my $strand=$col[3];
		my $loc=$col[4];
		($start,$end) = ($end,$start) if ($start>$end);
		my @letter=split("",$loc);
		my $noteDesc=$line;
		foreach my $let(@letter){
			next if ! defined $let;
			$let=uc $let;
			my $desc='-';
			$desc=$codes{$let} if exists $codes{$let};
			#print "$let => $desc\n" if $desc ne '-';
			$noteDesc.="; $desc" if $desc ne '-';
		}

		#print "$l_id\t$start\t$end\t$strand\t$loc\t$gene\t$operon\t$line\n";
		#print "$chrom\tterminator\t$start\t$end\t$strand\tterminator\t$gene\t$operon\t$l_id\t$det$noteDesc\n";
		#print "$chrom\tterminator\t$start\t$end\t$strand\tterminator\t-\t$operon\t$det$noteDesc\n"; #\t$det$noteDesc

		if($strand eq "+"){
			my $tts2locus=getSmall($chrom,$start);
			print "$chrom\tterminator\t$start\t$end\t$strand\tterminator\t$postiveEnd2LT{$chrom}{$tts2locus}\t$operon\t$det$noteDesc\n" if $tts2locus >0;
		}
		if($strand eq "-"){
			my $tts2locus=getLarge($chrom,$end);
			print "$chrom\tterminator\t$start\t$end\t$strand\tterminator\t$negativeEnd2LT{$chrom}{$tts2locus}\t$operon\t$det$noteDesc\n" if $tts2locus >0;
		}
		
		
		@geneStart=();
	}
	

}
close $tfh;

sub getSmall{
	my $chrom=shift;
	my $num=shift;
	require Number::Closest;
	use List::MoreUtils qw(first_index);
	my @array= sort keys %{$postiveEnd2LT{$chrom}};
	my $finder = Number::Closest->new(number => $num, numbers => \@array) ;
	my $closest5 = $finder->find(5) ;
	push @{$closest5}, $num;
	my @sorted=sort { $a <=> $b }@{$closest5};
	my $idx = first_index { $_ eq $num }  @sorted;
	my $newnum=$sorted[$idx-1];
	if(abs($newnum-$num)>500){
		$newnum=0;
	}
	return $newnum;
}


sub getLarge{
	my $chrom=shift;
	my $num=shift;
	require Number::Closest;
	use List::MoreUtils qw(first_index);
	my @array= sort keys %{$negativeEnd2LT{$chrom}};
	my $finder = Number::Closest->new(number => $num, numbers => \@array) ;
	my $closest5 = $finder->find(8) ;
	push @{$closest5}, $num;
	my @sorted=sort { $a <=> $b }@{$closest5};
	my $idx = first_index { $_ eq $num }  @sorted;
	my $newnum=$sorted[$idx+1];
	if(abs($newnum-$num)>500){
		$newnum=0;
	}
	return $newnum;

}
