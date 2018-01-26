#!/usr/bin/perl

# Author: Barbara Hutter (b.hutter@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Benedikt Brors, Barbara Hutter

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

use strict;
use warnings;

if (@ARGV < 2)
{
	die "Usage: $0: 1:readcount-coverage of BedTools on RefSeq exons for summing counts per gene - 2:Refseq bed file - optional 3:total read number (default:sum of reads from first file\n";
}
my $file = shift;
my $refseq = shift;
my $total = shift;

my @help;
my $lines = 0;
my %gene_reads = ();	# gene name => reads 
my %gene_sums = ();	# gene name => sum exon lengths for fpkm (has to be divided by 1000 to yield "per kb exon"!)
my $readsum = 0;	# to calculate FPKM per gene
my $fpkm = 0;
# there is no obvious order of genes/coordinates so have to collect per gene name
# format: chr start end name strand exon_nr exon_length reads sum_covered_bases exon_length coverage
#chr17   41267742        41267796        BRCA1   -       22      54      0       0       54      0.0000000
#chr17   41276033        41276132        BRCA1   -       23      99      1       51      99      0.5151515

open (FH, $file) or die "cannot open $file: $!\n";

while (<FH>)
{
	if ($_ =~ /#/){next;}
	$lines++;
	@help = split ("\t", $_);
	$gene_reads{$help[3]} += $help[7];
	$gene_sums{$help[3]} += $help[6];
	$readsum+=$help[7];
}
close FH;
# total reads mapped to exons (can occur multiple times if overlap > 1 exon!) in million:
if (defined $total && $total > 1)
{
	$readsum = $total;
}
print STDERR "$readsum reads\n";
$readsum = $readsum/1_000_000;

# now I need to get the gene coordinates!
my %refgenes;
open (RG, $refseq) or die "cannot open $refseq: $!\n";
my $header = <RG>;
chomp $header;
while (<RG>)
{
	chomp;
	@help = split ("\t", $_);
	$refgenes{$help[3]} = $_;
}
close RG;


print "$header\treads\texonsum\tRPKM\n";
foreach my $key (sort keys %gene_reads)
{
	print $refgenes{$key}, "\t", $gene_reads{$key}, "\t", $gene_sums{$key}, "\t";
	#$fpkm = sprintf ("%.4f", ($gene_reads{$key}/$gene_sums{$key}/$readsum));
	$fpkm = sprintf ("%.2f", ($gene_reads{$key}/($gene_sums{$key}/1000)/$readsum));
	print $fpkm, "\n";
}
