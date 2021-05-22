#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
# my $operonfile='/lrlhps/users/c240616/P_lentus/final/assembly/1240/PLEN1240_g_11232016/gbk/operon.tab';
# my $cdsfile='/lrlhps/users/c240616/P_lentus/final/assembly/1240/PLEN1240_g_11232016/gbk/cds.tab';
# my $rbsfile='/lrlhps/users/c240616/P_lentus/final/assembly/1240/PLEN1240_g_11232016/PLEN1240_g_11232016.rbs';
# 

# my $operonfile='/lrlhps/users/c240616/P_lentus/final/assembly/161/PLEN161_g_11232016/gbk/operon.tab';
# my $cdsfile='/lrlhps/users/c240616/P_lentus/final/assembly/161/PLEN161_g_11232016/gbk/cds.tab';
# my $rbsfile='/lrlhps/users/c240616/P_lentus/final/assembly/161/PLEN161_g_11232016/PLEN161_g_11232016.rbs';

# my $operonfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU170331_09182017/gbk/operons.tab';
# my $cdsfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU170331_09182017/gbk/cds.tab';
# my $rbsfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU170331_09182017/LREU170331_09182017.rbs';

#my $operonfile='/lrlhps/users/c240616/projects/Streptomyces/Assembly/SA73515/gbk/operons.tab';
#my $cdsfile='/lrlhps/users/c240616/projects/Streptomyces/Assembly/SA73515/gbk/cds.tab';
#my $rbsfile='/lrlhps/users/c240616/projects/Streptomyces/Assembly/SA73515/SA73515.rbs';

#my $operonfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU3632_09182017/gbk/new/operon.tab';
#my $cdsfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU3632_09182017/gbk/new/cds.tab';
#my $rbsfile='/lrlhps/users/c240616/projects/Lactobacillus/FINAL/prokka/LREU3632_09182017/LREU3632_09182017.rbs';
GetOptions( 'cdsfile=s' => \my $cdsfile,
            'rbsfile=s' => \my $rbsfile,
          , 'tag=s' => \my $tag  # same for --user or -u
          );
#print "$cdsfile\n";
#print "$rbsfile\n";
#print "$tag\n";

#my $cdsfile='gbk/cds.tab';
#my $rbsfile='PTA-126787.rbs';


my %start2gene;
my %gene2operon;
#my %gene2symbol;
# my %start2chrom;
# my %chroms;

open my $cfh, '<', $cdsfile or die $!;
while(<$cfh>){
	chomp;
	next if /^chromosome/;
	my @col= split /\t/;
	my $chrom=shift @col;
	#$chroms{$chrom}++;
	#$start2chrom{$col[6]}=$chrom;
	$start2gene{$col[1]}=$col[6];
	$gene2operon{$col[6]}=$col[-1];
	#$gene2symbol{$col[6]}=$col[4];

}
close $cfh;





print join("\t", qw/chromosome FeatureType start end strand regulatory_class locus_tag operon note/),"\n";

my $genePrefix=$tag; #'IU404_rbs';
$genePrefix="$genePrefix";
my $digits=5;
my $d="d";
my $featureCnt=0;

open my $fh, '<', $rbsfile or die $!;
while(<$fh>){
	chomp;
	next if !/^\s*$genePrefix/; #******************IMPORTANT******************************
	s/^\s+//g;
	my @col= split /\s+/;
	#print join("\t", @col),"\n";
	my $note='-';
	if ($col[4] >0){
		$featureCnt++;
		my $l_id=sprintf("%.$digits$d",$featureCnt);
		$l_id="$genePrefix$l_id";
		my $start=$col[-1]; my $end=$col[2];
		my $gene='-';
		my $operon='-';
		#my @out;

		print "$col[0]\tRBS\t";
		#push @out, 'RBS';
		if($col[1]<$col[2]){
			print "$col[4]\t", $col[4]+4, "\t+\t";
			# push @out,$col[4];
			# push @out,$col[4]+4;
			# push @out,'+';

		}else{
			print $col[4]-4, "\t" ,$col[4],"\t-\t";
			# push @out,$col[4]-4;
			# push @out,$col[4];
			# push @out, '-';
			($start,$end)=($end,$start);
			 
		}
		$gene=$start2gene{$start} if exists $start2gene{$start};
		$operon=$gene2operon{$gene} if exists $gene2operon{$gene};
		#$gene=$gene2symbol{$gene} if (exists $gene2symbol{$gene} and $gene2symbol{$gene} ne '-');
		$note="New CDS start predicted based on RBSFinder:$col[1]" if ($col[1]!=$col[-1]);
		print "ribosome_binding_site\t$gene\t$operon\t$note\n";
		#print "$col[0]\tRBS\t", join("\t",@out),"\tribosome_binding_site\t$gene\t$operon\t$l_id\t$note\n";
	}
	
}
close $fh;

#perl -F"\s+" -lane ' if ($F[4] >0){if($F[1]<$F[2]){print "perl getSubSeq -i ../PLEN1240_g_11232016.fna -t \"",$F[3],":",$F[1],"-",$F[2],"\" -s $F[4] -e ", $F[4]+4 }else{print "perl getSubSeq -i ../PLEN1240_g_11232016.fna -t \"",$F[3],":",$F[1],"-",$F[2],"\"  -o - -s ",$F[4]-4," -e ", $F[4] }}' ../PLEN1240_g_11232016.rbs   >test ;bash test |perl -nle 's/\:/\t/g;s/\-/\t/g;print' > test.res
