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

my %rnatype=(
	"6S"	=> 	"other",
"Bacteria_large_SRP"	=> 	"SRP_RNA",
"Bacteria_small_SRP"	=> 	"SRP_RNA",
"BsrC"	=> 	"other",
"glmS"	=> 	"ribozyme",
"group-II-D1D4-4"	=> 	"ribozyme",
"Intron_gpII"	=> 	"ribozyme",
"PK-G12rRNA"	=> 	"other",
"PyrR"	=> 	"other",
"sau-50" => "other",
"RNaseP_bact_a"	=> 	"ribozyme",
"ASpks" => "other",
"ASdes" => "other",
"msiK" => "other",
"Actino-pnp" => "other",
"5_ureB_sRNA" => "other",
"6C" => "other"

);

# NEED TO ADD locus tag
my $genePrefix=$tag; #'LREU3632_ncRNA';
$genePrefix="$genePrefix";
my $digits=5;
my $d="d";
my $featureCnt=0;
my $chrom='-';
print join("\t", qw/chromosome FeatureType start end strand ncRNA_class note locus_tag/),"\n";
my %sorted;
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
		$sorted{$start}="$chrom\tncRNA\t$start\t$end\t$strand\t$rnatype{$col[6]}\t$col[6]";
	}
	

}
close $fh;
foreach my $pos(sort { $a <=> $b } keys %sorted){
	$featureCnt++;
	my $l_id=sprintf("%.$digits$d",$featureCnt);
	$l_id="$genePrefix$l_id";
	print "$sorted{$pos}\t$l_id\n";
}
