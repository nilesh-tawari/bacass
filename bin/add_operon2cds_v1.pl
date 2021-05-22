use strict;
use warnings;
use Getopt::Long;
use File::Copy 'move';

GetOptions( 'cdsfile=s' => \my $cdsfile,
            'opfile=s' => \my $opfile
          );
#my $cdsfile='gbk/cds.tab';
#my $opfile='Rockhopper/operon.tab.txt';
# operon.tab header format: chromosome	FeatureType	start	end	strand	locus_tag	note
my $backup=$cdsfile.".bak";

move $cdsfile, $backup or die "can't move $cdsfile to $backup!\n";

my %gene2op;
open my $opfh ,'<', $opfile or die $!;

while(<$opfh>){
	chomp;
	my @col=split /\t/;
	#print "$col[-2]\t$col[-1]\n";
	$col[-1]=~s/Num. CDS:\d+;\s*//;
	#print "$col[-1]\n";
	my @loc_tags=split(",",$col[-1]);
	foreach my $loc (@loc_tags){
		$gene2op{$loc}=$col[-2];
		#print "$loc\t$col[-2]\n";
	}

}
close $opfh;
open my $cdsfh ,'<', $backup or die $!;
#open my $cdsfh ,'<', $cdsfile or die $!;
open my $cdsout ,'>', $cdsfile or die $!;
# chromosome	FeatureType	start	end	strand	gene	inference	locus_tag	product	transl_table	EC_number	note	pseudo
my @header;
my $startIdx;
while(<$cdsfh>){
	chomp;
	my %data;
	my @col=split /\t/;
	if (/^chromosome/){
		@header=@col;
		print $cdsout join("\t", @header),"\t","operon\n";
	}else{
		for(my $i=0;$i<@col;$i++){
			#print "$header[$i]\t=>\t$col[$i]\n";
			$data{$header[$i]}=$col[$i];
			
		}
		print $cdsout join("\t",@col);
		if (exists $gene2op{$data{'locus_tag'}}){
			#print $gene2op{$data{'locus_tag'}},"\toperon\n"
			print  $cdsout "\t", $gene2op{$data{'locus_tag'}},"\n";
		}else{
			print $cdsout "\t-\n";
		}
	}

}
close $cdsfh;
close $cdsout;
