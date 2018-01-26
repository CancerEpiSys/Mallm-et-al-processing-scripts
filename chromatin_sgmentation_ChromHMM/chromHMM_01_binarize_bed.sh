#!/bin/bash

# Author: Naveed Ishaque (n.ishaqaue@dkfz.de)

# Citation: Linking aberrant chromatin features in chronic lymphocytic leukemia to deregulated transcription factor networks. Jan-Philipp Mallm, Murat Iskar, Naveed Ishaque et al.

############################################

# Copyright 2018 DKFZ, Karsten Rippe, Daniel Mertens, Roland Eils, Naveed Ishaque

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################################

clear

## SOFTWARE STACK - these binaries were used for the analysis

CHROMHMM_PATH="/path/to/chomrhmm/folder/"
JAVA_BINARY="java8 -mx2000M"

## PATHS and FILES - this has to change with your files structure

DATABASE_PATH="/path/to/references/and/indexes"
THREADS=8
SAMPLE="SAMPLE_XYZ"
OUTPUT_DIR="/path/to/results/$SAMPLE/chromHMM"
INPUT_DIR="/path/to/results/$SAMPLE/beds"

chromosomelengthfile inputbeddir cellmarkfiletable outputbinarydir

## ChromHMM binarize BED files

echo
echo "#############################"
echo "### ChromHMM binarize bed ###"
echo "#############################"
echo

read -p "#INPUT REQUIRED: please define path to ChromHMM folder (default $CHROMHMM_PATH): " CHROMHMM_PATH_temp
if [ -d $CHROMHMM_PATH_temp ]
then
  CHROMHMM_PATH=$CHROMHMM_PATH_temp
else
  echo "#ERROR: could not find directory '$CHROMHMM_PATH_temp'... exiting"
  echo
  exit
fi

read -p "#INPUT REQUIRED: please define peak file list (no default): " PEAK_FILE_LIST_temp
if [ -e $PEAK_FILE_LIST_temp ]
then                                                                                                                                                                                                                                               
  PEAK_FILE_LIST=$PEAK_FILE_LIST_temp
else
  echo "#ERROR: peakfile list '$PEAK_FILE_LIST_temp' could not be found... exitting"
  echo
  exit
fi

read -p "#INPUT REQUIRED: please define peak file folder (no default): " INPUT_DIR_temp
if [ -d $INPUT_DIR_temp ]
then
  INPUT_DIR=$INPUT_DIR_temp
else
  echo "#ERROR: peakfile folder '$INPUT_DIR_temp' could not be found... exitting"
  echo
  exit
fi                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                   
read -p "#INPUT REQUIRED: please define path to output folder for alignments (default $OUTPUT_DIR): " OUTPUT_DIR_temp                                                                                                                              
if [ -d $OUTPUT_DIR_temp ]                                                                                                                                                                                                                         
then                                                                                                                                                                                                                                               
  OUTPUT_DIR=$OUTPUT_DIR_temp                                                                                                                                                                                                                      
else                                                                                                                                                                                                                                               
  echo "#ERROR: cannot find '$OUTPUT_DIR_temp'... using default STAR index folder $OUTPUT_DIR"                                                                                                                                                     
  echo                                                                                                                                                                                                                                             
fi                                                                      

echo
echo "$JAVA_BINARY -jar $CHROMHMM_PATH/ChromHMM.jar BinarizeBed -peaks $CHROMHMM_PATH/CHROMSIZES/hg19.txt $INPUT_DIR $PEAK_FILE_LIST $OUTPUT_DIR" 
echo
