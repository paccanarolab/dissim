#!/usr/bin/perl

use strict;
use warnings;
use LWP::UserAgent;

######################################################
# Gets a mimtoprot file as an input, and writes a 
# Fasta file and a file with the unique identifiers
# as an output
######################################################

die "USE $0 mimtoprot_file fasta_file unique_ids_file" if (scalar @ARGV != 3);

(my $mimtoprot, my $fasta, my $unique_ids) = @ARGV;

# 1.- read of the mimtosp file to gather the uniprot ids

my %uniprot_ids = ();

open (my $in1, "<", $mimtoprot) or die "The input file \"$mimtoprot\" does not exist";
while (<$in1>)
{
	chomp;
	(my $omim_id, my $uniprot_id) = split;
	$uniprot_ids{$uniprot_id}++;
}

close ($in1);

# 2.- print of the unique ids file
print "Printing unique ids...\n";

open (my $out2, ">", $unique_ids) or die "The output file \"$unique_ids\" does not exist";

for my $id (sort keys %uniprot_ids)
{
	print $out2 $id, "\n";
}

close ($out2);

# 3.- we write the UniProt fasta file

print "Fetching sequences...\n";
open (my $out, ">", $fasta) or die "The output file \"$fasta\" does not exist";

my $list = $unique_ids; # File containg list of UniProt identifiers.
my $base = 'http://www.uniprot.org';
my $tool = 'batch';

my $contact = 'aeromero@cs.rhul.ac.uk'; # Please set your email address here to help us debug in case of problems.
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
push @{$agent->requests_redirectable}, 'POST';
my $response = $agent->post("$base/$tool/",
                            [ 'file' => [$list],
                              'format' => 'fasta',
                            ],
                            'Content_Type' => 'form-data');

while (my $wait = $response->header('Retry-After')) {
  sleep $wait;
  $response = $agent->get($response->base);
}

$response->is_success ?
  print $out $response->content :
  die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";

close($out);

print "Sequences successfully retrieved from UniProtKB\n";








