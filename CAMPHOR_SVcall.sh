#!/bin/bash

BAM=$1
SORTED_BAM=$2
OUTPUT=$3

CONFIG=./pram.config

. $CONFIG

SRC=./src
DATA=./data

if [ ! -d $OUTPUT ]; then mkdir $OUTPUT; fi
#echo "$SAMTOOLS view $BAM|python $SRC/read_format.py /dev/stdin $MIN_INDEL_LENGTH 0.3|python $SRC/SV_format.py /dev/stdin 0.05 $MIN_INDEL_LENGTH 0 10 0.3 1000 20 > $OUTPUT/SV_read.txt"
$SAMTOOLS view $BAM|python3 $SRC/read_format.py /dev/stdin $MIN_INDEL_LENGTH 0.3|python3 $SRC/SV_format.py /dev/stdin 0.05 $MIN_INDEL_LENGTH 0 10 0.3 1000 20 > $OUTPUT/SV_read.txt

#echo "python $SRC/SV_candidate.py $OUTPUT/SV_read.txt|sort -k2,2g -k3,3 -k4,4g|perl $SRC/SV_format_change.5.pl /dev/stdin > $OUTPUT/SV_read.txt.2"
python3 $SRC/SV_candidate.py $OUTPUT/SV_read.txt|sort -k2,2g -k3,3 -k4,4g|perl $SRC/SV_format_change.5.pl /dev/stdin > $OUTPUT/SV_read.txt.2

#echo "python $SRC/SV_candidate_INS.py $OUTPUT/SV_read.txt 100000|sort -k2,2g -k3,3 -k4,4g|perl $SRC/SV_format_change.5.pl /dev/stdin|python $SRC/SV_support_read.py /dev/stdin 50 0.2 20 50 1 30 0.3|sort -k1,1 -k2,2g|python $SRC/merge_BP_INS.py /dev/stdin 300 0.8|sort -k1,1 -k2,2g > $OUTPUT/INS_candidate.txt"
python3 $SRC/SV_candidate_INS.py $OUTPUT/SV_read.txt 100000|sort -k2,2g -k3,3 -k4,4g|perl $SRC/SV_format_change.5.pl /dev/stdin|python3 $SRC/SV_support_read.py /dev/stdin 50 0.2 20 50 1 30 0.3|sort -k1,1 -k2,2g|python3 $SRC/merge_BP_INS.py /dev/stdin 300 0.8|sort -k1,1 -k2,2g > $OUTPUT/INS_candidate.txt

#echo "python $SRC/SV_support_read.2.py $OUTPUT/SV_read.txt.2 50 DEL|sort -k1,1 -k2,2g -k3,3 -k4,4g|python $SRC/merge_BP.py /dev/stdin 1000 0.8|sort -k1,1 -k2,2g -k3,3 -k4,4g > $OUTPUT/DEL_candidate.txt"
python3 $SRC/SV_support_read.2.py $OUTPUT/SV_read.txt.2 50 DEL|sort -k1,1 -k2,2g -k3,3 -k4,4g|python3 $SRC/merge_BP.py /dev/stdin 1000 0.8|sort -k1,1 -k2,2g -k3,3 -k4,4g > $OUTPUT/DEL_candidate.txt

#echo "python $SRC/SV_support_read.2.py $OUTPUT/SV_read.txt.2 1500 CHR > $OUTPUT/CHR_candidate.txt"
python3 $SRC/SV_support_read.2.py $OUTPUT/SV_read.txt.2 1500 CHR > $OUTPUT/CHR_candidate.txt

#echo "python $SRC/SV_support_read.2.py $OUTPUT/SV_read.txt.2 100 DIF|sort -k1,1 -k2,2g -k3,3 -k4,4g|python $SRC/merge_BP.py /dev/stdin 1000 0.8|sort -k1,1 -k2,2g -k3,3 -k4,4g > $OUTPUT/INV_candidate.txt"
python3 $SRC/SV_support_read.2.py $OUTPUT/SV_read.txt.2 100 DIF|sort -k1,1 -k2,2g -k3,3 -k4,4g|python3 $SRC/merge_BP.py /dev/stdin 1000 0.8|sort -k1,1 -k2,2g -k3,3 -k4,4g > $OUTPUT/INV_candidate.txt

#echo "python $SRC/SV_support_read.2.py $OUTPUT/SV_read.txt.2 100 TRS|sort -k1,1 -k2,2g -k3,3 -k4,4g|python $SRC/merge_BP.py /dev/stdin 1000 0.8|sort -k1,1 -k2,2g -k3,3 -k4,4g > $OUTPUT/TRS_candidate.txt"
python3 $SRC/SV_support_read.2.py $OUTPUT/SV_read.txt.2 100 TRS|sort -k1,1 -k2,2g -k3,3 -k4,4g|python3 $SRC/merge_BP.py /dev/stdin 1000 0.8|sort -k1,1 -k2,2g -k3,3 -k4,4g > $OUTPUT/TRS_candidate.txt

rm $OUTPUT/SV_read.txt $OUTPUT/SV_read.txt.2

python3 ./src/SV_selection.py $OUTPUT/INS_candidate.txt 2 > $OUTPUT/INS_candidate.txt2
python3 ./src/SV_selection.py $OUTPUT/DEL_candidate.txt 2 > $OUTPUT/DEL_candidate.txt2
python3 ./src/SV_selection.py $OUTPUT/CHR_candidate.txt 2 > $OUTPUT/CHR_candidate.txt2
python3 ./src/SV_selection.py $OUTPUT/INV_candidate.txt 2 > $OUTPUT/INV_candidate.txt2
python3 ./src/SV_selection.py $OUTPUT/TRS_candidate.txt 2 > $OUTPUT/TRS_candidate.txt2

python3 ./src/SV_filter_for_comparison.py $OUTPUT/DEL_candidate.txt 0.05 20 800 1 > $OUTPUT/DEL_candidate.txt0

python3 ./src/SV_filter_for_comparison.py $OUTPUT/INS_candidate.txt 0.05 20 800 1 > $OUTPUT/INS_candidate.txt0

python3 ./src/SV_filter_for_comparison.py $OUTPUT/INV_candidate.txt 0.05 20 800 1 > $OUTPUT/INV_candidate.txt0

python3 ./src/SV_filter_for_comparison.py $OUTPUT/TRS_candidate.txt 0.05 20 800 1 > $OUTPUT/TRS_candidate.txt0

python3 ./src/SV_filter_for_comparison.py $OUTPUT/CHR_candidate.txt 0.05 20 800 1 > $OUTPUT/CHR_candidate.txt0

