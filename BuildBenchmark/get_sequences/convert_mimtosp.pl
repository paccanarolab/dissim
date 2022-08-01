#!/usr/bin/perl

#######################################################
#
# This script requires the file mimtosp,
# available here:
#   http://www.uniprot.org/docs/mimtosp.txt
# The file is downloaded to the current directory
#
#######################################################

use strict;
use warnings;
use File::Temp qw/ tempfile tempdir /;

# file download
`wget -N http://www.uniprot.org/docs/mimtosp.txt`;

open (my $in, "<", "mimtosp.txt");

my $flag = 0;
($_, my $filename) = tempfile();

open (my $tmp, ">", $filename);

while (<$in>)
{
    if (m/^_/)
    {
       $flag = 1;   
       next;
     }

    last if (m/^\s*$/ and $flag == 1);

    next if ($flag == 0);

    s/^\s+,/,/g;
    s/,\s+$/, /;

    print $tmp $_;
}
close ($tmp);
close ($in); 

open ($tmp, "<", $filename);
open (my $out, ">", "mimtoprot.txt");
while (<$tmp>)
{
    chomp;
    my @fields = split(/:/, $_);
    
    my $omim = $fields[0];

    my @protz = split(/,/, $fields[1]);
    for my $prot (@protz)
    {
        if ($prot =~ m/(\()(\w+)(\))/)
        {
            print $out $omim, "\t", $2, "\n";
        }
        
    }
    
}

close ($tmp);
close ($out);

print $filename;

