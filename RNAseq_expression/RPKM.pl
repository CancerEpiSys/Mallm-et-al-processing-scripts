#!/usr/bin/perl

# Barbara Hutter (b.hutter@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Benedikt Brors, Barbara Hutter

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################


# calculate RPKM

use strict;
use warnings;

if (@ARGV < 2)
{
	die "$0 (1)BED file with read numbers and coordinates to calculate RPKM - (2)total number of mapped reads\n";
}

my $file = shift;
my $total = shift;

# in million:
$total = $total/1_000_000;

my $ctr = 0;
my @help = ();
my $fpkm = 0;
my $rncol = 0;
open (FH, $file) or die "could not open $file: $!\n";
while (<FH>)
{
	$ctr++;
	#if ($ctr > 3){last;}
	chomp $_;
	@help = split ("\t", $_);
	# the column with the read count is the 3rd-to-last one
	$rncol = $#help-3; # original ist -3
	# fragments per kilobase of feature (exon/repeat/...) per million mapped reads
	$fpkm = sprintf ("%.4f", ($help[$rncol]/(($help[2]-$help[1])/1000)/$total));
	print $_, "\t", $fpkm, "\n";
}
close FH;
print STDERR "$ctr lines\n";
