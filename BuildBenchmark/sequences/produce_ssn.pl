#!/usr/bin/perl

####################################################
#  Produces the SSN
#
#  It takes two files as an input:
#       1.- a table relating diseases an protein ids
#           disease_id    protein_acc_id
#           ...
#       2.- the output from a BLAST/SW result in 
#       tabular format of the table with itself
#
#  The third parameter will be the output file for
#  the binary graph. Two diseases will be related if
#  at least one of their corresponding proteins are
#  similar in sequence with evalue < 1e-6. This value
#  can be changed inputting an optional fourth parameter
#  to the script
####################################################

use strict;
use warnings;

die "USE $0 table_diseases blast_table output_file [evalue]" if (scalar @ARGV < 3 || scalar @ARGV > 4);

my $evalue = 1e-6;
(my $table_diseases, my $blast_table, my $output_file);
if (scalar @ARGV == 3)
{
    ($table_diseases, $blast_table, $output_file) = @ARGV;    
} else {
    ($table_diseases, $blast_table, $output_file, $evalue) = @ARGV;
    if ($evalue < 0.0)
    {
        print "'evalue' should be a positive number\n";
        exit;    
    }
}


# we read the disease table
open (my $in_disease, "<", $table_diseases);

my %disease = ();

while (<$in_disease>)
{
    chomp;
    (my $dis, my $prot) = split;
    if (not exists $disease{$dis})
    {
        $disease{$dis} = $prot;
    } else {
        $disease{$dis} = $disease{$dis} . "_" . $prot;    
    }
}
close ($in_disease);

my %pair_prots = ();
open (my $prots, "<", $blast_table);

while (<$prots>)
{
    chomp;
    my @fields = split;
    (my $prot1, my $prot2, my $_evalue) = @fields[0,1,10];
    next if ($_evalue >= $evalue);
    #remove the identical proteins.
    next if($prot1 eq $prot2);

    my $pr1 = (split(/\|/, $prot1))[1];
    my $pr2 = (split(/\|/, $prot2))[1];

    $pair_prots{ $pr1 . "_" . $pr2} = 1;
    $pair_prots{ $pr2 . "_" . $pr1} = 1;
}

close ($prots);

my @diseases = sort keys %disease;
my $n = scalar @diseases;

open (my $out, ">", $output_file);

for (my $i=0; $i<$n; ++$i)
{
    my $dis_i = $diseases[$i];
    my @prot_i = split(/_/, $disease{$dis_i});
    for (my $j=$i+1; $j<$n; ++$j)        
    {
        my $dis_j = $diseases[$j];
        my @prot_j = split(/_/, $disease{$dis_j});

        my $found = 0;

        for my $pi (@prot_i)
        {
            for my $pj (@prot_j)
            {
                if (exists $pair_prots{$pi . "_" . $pj} or exists $pair_prots{$pj . "_" . $pi})
                {
                    $found = 1;
                    last;    
                }
            }
            last if ($found == 1);
        }


        if ($found)
        {
            print $out $dis_i, "\t", $dis_j, "\t", 1.0, "\n";
        }
    }
}

close ($out);




