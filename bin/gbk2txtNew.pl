use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Std;
getopts('i:f:h');
our($opt_i,$opt_f,$opt_h);

HelpMsg() if defined $opt_h;

my $ftype=$opt_f or die Usage();

my $file = $opt_i or die Usage() ;

#my @pTags;
#my @allTags;

my %allTags;
push @{$allTags{'CDS'}}, qw/chromosome FeatureType start end strand gene inference locus_tag product transl_table EC_number note/;
push @{$allTags{'tRNA'}}, qw/chromosome FeatureType start end strand inference locus_tag note product/;
push @{$allTags{'tmRNA'}}, qw/chromosome FeatureType start end strand gene inference locus_tag note product/;
push @{$allTags{'repeat_regions'}}, qw/chromosome FeatureType start end strand rpt_family/;# CRISPR features from prokka

print join("\t",@{$allTags{$ftype}}),"\n";
#print join("\t", qw/chromosome FeatureType start end strand gene inference locus_tag product transl_table EC_number note/),"\n"; # for CDS features
#print join("\t", qw/chromosome FeatureType start end strand inference locus_tag note product/),"\n"; # for tRNA features
#print join("\t", qw/chromosome FeatureType start end strand gene inference locus_tag note product/),"\n"; # for tmRNA features
#print join("\t", qw/chromosome FeatureType start end strand rpt_family/),"\n"; # CRISPR features from prokka

my %allFeatureTags;
push @{$allFeatureTags{'CDS'}}, qw/gene inference locus_tag product transl_table EC_number note/;
push @{$allFeatureTags{'tRNA'}}, qw/inference locus_tag note product/;
push @{$allFeatureTags{'tmRNA'}}, qw/gene inference locus_tag note product/;
push @{$allFeatureTags{'repeat_region'}}, 'rpt_family';

my $gbk = Bio::SeqIO->new(-file => $file, -format=>'genbank');
while(my $seq = $gbk->next_seq){
	my $chromosome=$seq->id();
	my $description=$seq->desc;
	my $sequence=$seq->seq();
	#print "$chromosome\t$description\n";
	my @features  = $seq->get_SeqFeatures;
	foreach my $feature (@features){
		if ($feature->primary_tag eq $ftype){ #CDS
			my @tags=$feature->get_all_tags;
			print $chromosome,"\t",$feature->primary_tag,"\t",$feature->start,"\t",$feature->end,"\t",$feature->strand;
			foreach my $t (@{$allFeatureTags{$feature->primary_tag}}){
				print "\t",tag($feature,$t); 
			}
			print "\n";
		}elsif($feature->primary_tag eq $ftype){ #tRNA
			print $chromosome,"\t",$feature->primary_tag,"\t",$feature->start,"\t",$feature->end,"\t",$feature->strand;
			foreach my $t (@{$allFeatureTags{$feature->primary_tag}}){
				print "\t",tag($feature,$t); 
			}
			print "\n";
		}elsif($feature->primary_tag eq $ftype){ #tmRNA
			print $chromosome,"\t",$feature->primary_tag,"\t",$feature->start,"\t",$feature->end,"\t",$feature->strand;
			foreach my $t (@{$allFeatureTags{$feature->primary_tag}}){
				print "\t",tag($feature,$t); 
			}
			print "\n";
		}elsif($feature->primary_tag eq $ftype){ #tmRNA
			print $chromosome,"\t",$feature->primary_tag,"\t",$feature->start,"\t",$feature->end,"\t",$feature->strand;
			foreach my $t (@{$allFeatureTags{$feature->primary_tag}}){
				print "\t",tag($feature,$t);  # tag($f, 'locus_tag');
			}
			print "\n";
		}

	}


}



sub tag {
   my($f, $tag) = @_;
   return '-' unless $f->has_tag($tag);
   return join(' ', $f->get_tag_values($tag));
}


sub Usage{
        return "Usage: $0 (-f) -i infile.gbk > outfile.tab\n";
}

sub HelpMsg{
        print STDERR Usage();
        print STDERR <<EOM;
        -i      input genbank file
        -f      feature type to parse (CDS,tRNA,tmRNA, CRISPR)
        NOTE: Type perldoc $0 for documentation.
EOM
        exit(0);

}
