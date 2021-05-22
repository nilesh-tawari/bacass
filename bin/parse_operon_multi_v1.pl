#!/usr/bin/env perl 
###########################################
# Parse Rockhopper operon output and print it in tab format. Also, assign new IDs to operons
# NOTE: The code needs to be changed to generate different outputs and operon prefixes
# Created : 1/18/2017
# Shrinivasrao P. Mane
###########################################
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

GetOptions('tag=s' => \my $tag  
                    );

my @files=glob("*_Rockhopper/_operons.txt");

#print join("\n",@files),"\n";

my $chrom="";
my $genePrefix=$tag;
$genePrefix="$genePrefix";
my $digits=5;
my $d="d";
print join("\t", qw/chromosome FeatureType start end strand locus_tag note/),"\n";
my $featureCnt=0;
foreach my $file(@files){
	my $f_pre=$file;
	$f_pre=~s/\_Rockhopper\/\_operons\.txt// ;
	#print $f_pre,"\n";
	$chrom=$f_pre;
	open my $fh , '<' , $file or die $!;
	
	while(<$fh>){
		chomp;
		next if /Start/;
		$featureCnt++;
		#print "$_\n";
		my @col=split /\t/;
		my $l_id=sprintf("%.$digits$d",$featureCnt);
		$l_id="$genePrefix$l_id";
		my $start=$col[0];
		my $end=$col[1];
		my $glist=$col[4];
		$glist=~s/ +//g;
		 print "$chrom\toperon\t$start\t$end\t$col[2]\t$l_id\tNum. CDS:$col[3]; $glist\n";
		my @genes=split /,/, $glist;
		foreach my $gene(@genes){
		
			#print "$gene\t$l_id\n";
		}
	}
	close $fh;
}

exit;

#my $file =shift or die $!;


