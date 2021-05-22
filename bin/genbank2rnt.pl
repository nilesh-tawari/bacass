#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

# This script takes a GenBank file as input, and produces a
# NCBI PTT file (protein table) as output. A PTT file is
# a line based, tab separated format with fixed column types.
#
# Written by Torsten Seemann
# 18 September 2006
#genbank_to_ptt.pl < infile.gbk > outfile.ptt
# Modified my Shrinivasrao P. Mane to extract tRNAs on 11/23/2016

my $gbk = Bio::SeqIO->new(-fh=>\*STDIN, -format=>'genbank');
my $seq = $gbk->next_seq;
my @cds = grep { $_->primary_tag eq 'tRNA' } $seq->get_SeqFeatures;

print $seq->description, " - 0..",$seq->length,"\n";
print scalar(@cds)," RNAs\n";
print join("\t", qw(Location Strand Length PID Gene Synonym Code COG 
Product)),"\n";

for my $f (@cds) {
   my $gi = '-';
   $gi = $1 if tag($f, 'db_xref') =~ m/\bGI:(\d+)\b/;
   my $cog = '-';
   # $cog = $1 if tag($f, 'product') =~ m/^(COG\S+)/;
   my @col = (
     $f->start.'..'.$f->end,
     $f->strand >= 0 ? '+' : '-',
     ($f->length)-1,
#     $gi,
#   I need to change what would be the GID to the Synonym Code for a
#   file which lacks official GB GIDs
     tag($f, 'locus_tag'),
     tag($f, 'gene'),
     tag($f, 'locus_tag'),
     $cog,
     '-',
     tag($f, 'product'),
   );
   print join("\t", @col), "\n";
}

sub tag {
   my($f, $tag) = @_;
   return '-' unless $f->has_tag($tag);
   return join(' ', $f->get_tag_values($tag));
}
