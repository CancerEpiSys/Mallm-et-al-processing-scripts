#

# Author: Naveed Ishaque (n.ishaqaue@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Naveed Ishaque

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

#

use strict;
use List::Util qw(max min);

sub largest_value_mem (\%) {
    my $hash   = shift;
    my ($key, @keys) = keys   %$hash;
    my ($big, @vals) = values %$hash;

    for (0 .. $#keys) {
        if ($vals[$_] > $big) {
            $big = $vals[$_];
            $key = $keys[$_];
        }
    }
    $key
}

#######################################
###### GLOBAL PARAMS ##################

my $min_recurrence = 2;
$min_recurrence = $min_recurrence-1;

my $min_stable = 2;
$min_stable = $min_stable-1;

my $afile = shift or die;
my $bfile = shift or die;

my $acount=`wc -l $afile |  cut -f 1 -d " "`;
my $bcount=`wc -l $bfile |  cut -f 1 -d " "`;

my $observation = 1 / ($acount * $bcount);

######################################

# Translation of states ....

my %states;

$states{1}=1;
$states{2}=2;
$states{3}=3;
$states{4}=4;
$states{5}=5;
$states{6}=6;
$states{7}=7;
$states{8}=8;
$states{9}=9;
$states{10}=10;
$states{11}=11;
$states{12}=12;

#######################################

# init count matrix

my %transition_matrix_c;

foreach  my $i (0..15){  
  foreach  my $j (0..15){
    $transition_matrix_c{$i}{$j}=0; 
    $transition_matrix_c{$i}{$j}=0; 
  }
}

my $total_sum = 0;

#######################################

#open files 

open(X1, "cat $afile | xargs paste |") or die "Cannot open files in A file: $afile\n\n";
open(Y1, "cat $bfile | xargs paste |") or die "Cannot open files in B file: $bfile\n\n";

my $line_counter = 0;
while (<X1>){
  # READ LINES, remove new line; translate to old states
  
  $line_counter++;
  my $x1_line = $_;
  my $y1_line = <Y1>;

  chomp($x1_line);
  chomp($y1_line);

  my @a_array = split /\t/, $x1_line;
  my @b_array = split /\t/, $y1_line;

  # convert asthma states to smoking states

  foreach my $a_index (1 .. $acount){
    $a_array[$a_index-1]=$states{$a_array[$a_index-1]};
  }

  foreach my $b_index (1 .. $bcount){
    $b_array[$b_index-1]=$states{$b_array[$b_index-1]};
  }

  # count the number of each state at position

  my %a_array_counts;
  foreach my $a_index (1 .. $acount){
    $a_array_counts{$a_array[$a_index-1]}++;
  }

  my %b_array_counts;
  foreach my $b_index (1 .. $bcount){
    $b_array_counts{$b_array[$b_index-1]}++;
  } 

  # IF max in A is not max in B
  # AND
  # IF max in A is stable (above 0.5) and IF max in B is stable (above 0.5)

  my ($max_a_count, $max_b_count, $max_a_value, $max_b_value);

  $max_a_value = largest_value_mem %a_array_counts;
  $max_b_value = largest_value_mem %b_array_counts;

  $max_a_count = $a_array_counts{$max_a_value};
  $max_b_count = $b_array_counts{$max_b_value}; 

  #print "\#SANITY\tA : @a_array\t$max_a_value x $max_a_count\;\tB : @b_array\t$max_b_value x $max_b_count\n";

  if (!($max_a_value eq $max_b_value) && ($max_a_count > ($acount/2)) && ($max_b_count > ($bcount/2)) ){
    #print "\#PASS\tA : @a_array\t$max_a_value x $max_a_count\;\tB : @b_array\t$max_b_value x $max_b_count\n";

    foreach my $a_index (1 .. $acount){  
      foreach my $b_index (1 .. $bcount){
         $transition_matrix_c{$a_array[$a_index-1]}{$b_array[$b_index-1]} = $transition_matrix_c{$a_array[$a_index-1]}{$b_array[$b_index-1]} + $observation;
         print "$line_counter\t$a_index\t$a_array[$a_index-1]\t$b_index\t$b_array[$b_index-1]\n" if ($a_array[$a_index-1] eq 12 && ($max_b_value eq 1 || $max_b_value eq 9 || $max_b_value eq 8 || $max_b_value eq 11));
      }
    }  
    $total_sum = $total_sum + 1;
  }
}

# Print out matrix
 
foreach  my $i (1..12){
  foreach  my $j (1..12){
#    print "".($transition_matrix_c{$i}{$j}*200)."\t";
  }
  print "\n";
}

