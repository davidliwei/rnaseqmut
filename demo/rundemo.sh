#!/bin/bash

# This demo illustrates the basic usage of rnaseqmut, including calling de-novo mutations, merging mutations from different samples, calling mutations based on a given list of mutations, and filtering mutations.

##############
# Data files

BAMFILELIST=`ls data/*.bam`

# user defined labels
LABELS="NORMAL1,NORMAL2,TUMOR1,TUMOR2"

##############
# setting up environments, including paths
# the path to the DEMO dir
CURRENTPATH=`pwd`
# the path to the upper level
cd ../
BASEPATH=`pwd`
cd $CURRENTPATH

# add bin directory to the PATH
export PATH=$BASEPATH/bin:$PATH
# add script directory to the PATH
export PATH=$BASEPATH/script:$PATH

####################
# step 0, cleaning
echo "####### cleaning ##########"
if [ ! -d results ]; then
  mkdir results
fi
rm -rf results/*.txt results/*.vcf
# step 1, de-novo mutation calling

echo ""
echo "####### Step 1, de-novo mutation calling ##########"
for file in $BAMFILELIST; do
  filebase=`basename $file`
  CMD="rnaseqmut $file > results/$filebase.1st.txt"
  echo "#### COMMAND LINE: $CMD"
  eval $CMD
done

# step 2, merge the 1st pass of vcf files

echo ""
echo "####### Step 2, merging mutations in Step 1 into a candidate mutation list  ##########"
CMD="merge1stfile results/*.1st.txt > results/ALLMUTLIST.txt"
echo "#### COMMAND LINE: $CMD"
eval $CMD


# step 3, mutation calling from the merged lists

echo ""
echo "####### Step 3, calling mutations again using the given list in Step 2  ##########"
for file in $BAMFILELIST; do
  filebase=`basename $file`
  CMD="rnaseqmut -l results/ALLMUTLIST.txt $file > results/$filebase.2nd.txt"
  echo "#### COMMAND LINE: $CMD"
  eval $CMD
done

# step 4, merge the second pass of mutations into a big table
# only mutations with at least 4 alternative reads support is kept (by default of merge2ndvcf.py)
echo ""
echo "####### Step 4, merging mutations in Step 3 into a big table  ##########"
CMD="python3 $BASEPATH/script/merge2ndvcf.py -l $LABELS results/*.2nd.txt > results/ALLMUT.txt" 
echo "#### COMMAND LINE: $CMD"
eval $CMD


# step 5, filter mutations based on user-defined parameters
echo ""
echo "####### Step 5, custom filtering based on mutations in Step 4  ##########"
# defining the two normal samples as control groups
CONTROLGROUP="0,1"

# the following command keep mutations that occur in at least 1 non-control sample with at least 10 alternative read support. 

# By default, filtermut.py will only keep mutations that occur in at least 1 non-control sample (-t option) with 20% frequency (-f) and 10 alternative read support (-d) , excluding those that also occur in control samples (-a) or does not have enough read coverage in control samples (-b)

CMD="python3 $BASEPATH/script/filtermut.py -d 10 -f 0.0 -b 0 -c $CONTROLGROUP  -l $LABELS < results/ALLMUT.txt > results/ALLMUT_FILTERED.vcf"
echo "#### COMMAND LINE: $CMD"
eval $CMD


 
echo ""
echo "####### DEMO completed succesfully.  Check results/ALLMUT_FILTERED.vcf for detected mutations.  ##########"

