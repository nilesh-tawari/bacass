#!/usr/bin/perl
use strict;
use warnings;


my $infile=shift or die "Usage: $0 infile.gbk\n";
use Bio::SeqIO;
my $seq_in = Bio::SeqIO->new(
	 -file   => "<$infile",
	-format => "genbank",
);
my $seq_out = Bio::SeqIO->new(
	-fh   => \*STDOUT, -format => "fasta" );
while (my $inseq = $seq_in->next_seq) {
    $seq_out->write_seq($inseq);
}
