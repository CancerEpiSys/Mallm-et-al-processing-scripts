#

# Author: Naveed Ishaque (n.ishaqaue@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Naveed Ishaque

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

use Getopt::Long;
use strict;

my $rose_directory="/path/to/ROSE";
my $usage = "\t$0 -histone_bam=H3K27ac.bam -control_bam=H3.bam -bed=H3K27ac_peak_calls.bed -gff=H3K27ac_peak_calls.gff -output_dir=super_enhancer_calling -pathToRose=\"/path/to/Rose\"\n\n\t\tNOTE: please provide either BED or GFF files\n\n";

my ($r, $c, $bed, $gff, $o);

GetOptions ("histone_bam=s" => \$r,
            "control_bam=s" => \$c,
            "bed=s"         => \$bed,
            "gff=s"         => \$gff,
            "output_dir=s"  => \$o,
            "pathToRose=s"  => \$rose_directory)
or die("ERROR: error in command line arguments\n\n\t$usage");

if (!($r||$c||($bed||$gff)||$o)){
  print "ERROR: please define all inputs\n\n\t$usage";
  exit;
}

my $rose_call = "cd $rose_directory && python ROSE_main.py -s 10000 -t 3000  -g HG19_GENCODE17";

if (-e "$gff"){
    print "\n$rose_call -i $gff -r $r -c $c -o $o\n\n";
}

else {
  chomp $bed;
  print "\nawk -F'\\t' '{print \$1\"\\t\"\$4\"\\t\\t\"\$2\"\\t\"\$3\"\\t\\t.\\t\\t\"\$4}' $bed > $bed.gff\n";
  print "\n$rose_call -i $bed.gff -r $r -c $c -o $o\n\n";
}
