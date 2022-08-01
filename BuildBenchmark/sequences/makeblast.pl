#!/usr/bin/perl

use strict;
use warnings;

die "USE $0 sequences_fasta_file blast_output" if (scalar @ARGV != 2);

(my $fasta, my $output) = @ARGV;

# 1.- if there is no BLAST database build, we make it...
if (! -e $fasta . ".phr")
{
    print "Making BLAST database...\n";
    `makeblastdb -dbtype prot -in $fasta`;
}

# 2.- we run BLAST and generate the "blast_output" file. Note that
# Smith-Waterman alignment are generated, and therefore the running
# time could be a bit tedious

print "Running Smith-Waterman alignment. Please, be patient\n";

`blastp -query $fasta -db $fasta -use_sw_tback -out $output -outfmt 6`;

print "BLAST output file ", $output, " generated\n";


