# parse ISEscan.py prediction
# Initial script written by 
# Shrinivasrao P. Mane on 11/16/2020
# Important: does not parse or print "terminal_inverted_repeat" from ISEscan output;
# the reason for that is that gbk2tbl.pl cannot parse >1 feature from the same feature type sharing the same chromosome location
# nned to fix that at a later time 
# ###################################
use strict;
use warnings;
use Getopt::Long;

#my $file="ISEscan/prediction/PTA-126787.fna.gff"; 
#/home/NILESH_RAMESH.TAWARI/Projects/18_Lreuteri/Publication_data/01_annotation/results/PTA-126787/PTA-126787_annotation/predictions/ISEscan/prediction/PTA-126787.fna.gff
GetOptions( 'file=s' => \my $file
#          , 'tag=s' => \my $tag  # same for --user or -u
          );
#print "$file\n";
#print "$tag\n";


my @header=qw/chromosome FeatureType start end strand inference mobile_element_type rpt_type note locus_tag/;
my $inference="COORDINATES:profile:ISEScan:1.7.2.1";

print join("\t", @header),"\n";
open my $fh, '<', $file or die $!;
while(<$fh>){
	chomp;
	next if /^\#/;
	#print "$_\n";
	s/(ID\=)*gnl\|Elanco\|//g;
	my @col=split /\t/;
	#print $col[-1],"\n";
	my @desc=split/;/, $col[-1];
	my $locus_tag=$desc[0];
	my $chrom= $locus_tag;
	$chrom=~s/\_IS\S+//;
	if ($col[-1]=~/family/){
		
		print "$chrom\tmobile_element\t$col[3]\t$col[4]\t$col[6]\t$inference\tinsertion sequence\t-\t$desc[1];$desc[2]\t$locus_tag\n";
	}
	# "terminal_inverted_repeat" parsing code
	#if ($col[-1]=~/parent/){
	#	my $locus_tag=$desc[0];
	#	my $chrom= $locus_tag;
	#	$chrom=~s/\_IS\S+//;
	#	print "$chrom\tmobile_element\t$col[3]\t$col[4]\t$col[6]\t$inference\tinsertion sequence\tinverted\t$desc[1];$col[2]\t$locus_tag\n";
	#}
}

# ##gff-version 3
# gnl|Elanco|IU404_1	ISEScan	insertion_sequence	93	1350	.	+	.	ID=gnl|Elanco|IU404_1_IS_1;family=IS30;cluster=IS30_223
# gnl|Elanco|IU404_1	ISEScan	terminal_inverted_repeat	93	109	.	+	.	ID=gnl|Elanco|IU404_1_IS_1_TIR;parent=gnl|Elanco|IU404_1_IS_1
# gnl|Elanco|IU404_1	ISEScan	terminal_inverted_repeat	1335	1350	.	+	.	ID=gnl|Elanco|IU404_1_IS_1_TIR;parent=gnl|Elanco|IU404_1_IS_1
# gnl|Elanco|IU404_1	ISEScan	insertion_sequence	23144	24357	.	-	.	ID=gnl|Elanco|IU404_1_IS_2;family=IS110;cluster=IS110_139
# gnl|Elanco|IU404_1	ISEScan	terminal_inverted_repeat	23144	23181	.	-	.	ID=gnl|Elanco|IU404_1_IS_2_TIR;parent=gnl|Elanco|IU404_1_IS_2
# gnl|Elanco|IU404_1	ISEScan	terminal_inverted_repeat	24321	24357	.	-	.	ID=gnl|Elanco|IU404_1_IS_2_TIR;parent=gnl|Elanco|IU404_1_IS_2
