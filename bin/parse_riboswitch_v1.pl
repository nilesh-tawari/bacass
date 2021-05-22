#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
#use Data::Dumper;
#my $file='/lrlhps/users/c240616/P_lentus/final/assembly/1240/PLEN1240_g_11232016/PLEN1240_g_11232016.cmscan';
#my $file='/lrlhps/users/c240616/P_lentus/final/assembly/161/PLEN161_g_11232016/PLEN161_g_11232016.cmscan';
#my $file='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU170331_09182017/LREU170331_09182017.cmscan';
#my $file='/lrlhps/users/c240616/projects/Streptomyces/Assembly/SA73515/SA73515.cmscan';
#my $file='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU2091_09182017/LREU2091_09182017.cmscan';
#my $file='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU3632_09182017/LREU3632_09182017.cmscan';

GetOptions( 'file=s' => \my $file
          , 'tag=s' => \my $tag  # same for --user or -u
          );
#print "$file\n";
#print "$tag\n";

open my $fh, '<', $file or die $!;


my @types=qw/SMK_box_riboswitch c-di-GMP-I Cobalamin crcB FMN L10_leader Lysine pan pfl PreQ1 Purine SAM T-box TPP ydaO-yuaA ykkC-yxkD ylbH/;
my %rnatype;
$rnatype{$_}++ for (@types);

my %rnadesc=("c-di-GMP-I" => "Cyclic di-GMP-I riboswitch",
"Cobalamin" => "Cobalamin riboswitch",
"crcB" => "crcB RNA motif/fluoride riboswitch",
"FMN" => "FMN riboswitch (RFN element)",
"L10_leader" => "Ribosomal protein L10 leader",
"Lysine" => "Lysine riboswitch",
"pan" => "pan motif",
"pfl" => "ZMP/ZTP riboswitch",
"PreQ1" => "PreQ1 riboswitch",
"Purine" => "Purine riboswitch",
"SAM" => "SAM riboswitch (S box leader)",
"T-box" => "T-box leader",
"TPP" => "TPP riboswitch (THI element)",
"ydaO-yuaA" => "ydaO/yuaA leader",
"ykkC-yxkD" => "ykkC-yxkD leader",
"ylbH" => "ylbH leader",
"SMK_box_riboswitch" => 'SMK box riboswitch'
);

# NEED TO ADD locus tag
my $genePrefix=$tag; #'LREU3632_ribo';
$genePrefix="$genePrefix";
my $digits=5;
my $d="d";
my $featureCnt=0;





print join("\t", qw/chromosome FeatureType start end strand regulatory_class note locus_tag/),"\n";

# print join("\t", qw/FeatureType start end strand ncRNA_class note locus_tag/),"\n";
my %sorted;
my $chrom='-';
while(<$fh>){
	chomp;
	s/(ID\=)*gnl\|Elanco\|//g;
	$chrom=$1 if /Query:\s+(\S+)/;
	last if /^Hit alignments/;
	next if ! /\(\d+\) \!/;
	s/\s+/\t/g;
	my @col=split (/\t/);
	#print join("\t",@col),"\n";
	#print "$col[7]\n";
	#@print Dumper @col,"\n";
	my $start=$col[7]; my $end=$col[8];
	my $strand=1;
	if ($start>$end){
		($start,$end) = ($end,$start);
		$strand=-1;
	}
	if (exists $rnatype{$col[6]}){
		#print "ncRNA\t$start\t$end\t$strand\t$rnatype{$col[6]}\t$col[6]\t$l_id\n" ;
		$sorted{$start}="$chrom\triboswitch\t$start\t$end\t$strand\triboswitch\t$rnadesc{$col[6]}";
	}
	

}
close $fh;
foreach my $pos(sort { $a <=> $b } keys %sorted){
	$featureCnt++;
	my $l_id=sprintf("%.$digits$d",$featureCnt);
	$l_id="$genePrefix$l_id";
	print "$sorted{$pos}\t$l_id\n";
}
