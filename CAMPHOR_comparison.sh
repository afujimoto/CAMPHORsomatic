#!/bin.bash

CANCER_OUTPUT=$1
NORMAL_OUTPUT=$2
CANCER_BAM=$3
NORMAL_BAM=$4
FASTQ=$5
OUTPUT=$6

if [ ! -d $OUTPUT ]; then mkdir $OUTPUT; fi

CONFIG=./pram.config

. $CONFIG

#sh /share/amed_snt/WORK/fujimoto/src/nanopore/CAMPHOR_somatic/CAMPHOR_comparison2.sh /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_C_WGS.All/out/All.50.DEL  /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_B_WGS.All/out/All.50.no_read_num_flt.DEL.filt /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_C_WGS.All/out/All.100.DIF  /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_B_WGS.All/out/All.100.no_read_num_flt.DIF.filt /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_C_WGS.All/out/All.1500.CHR  /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_B_WGS.All/out/All.1500.no_read_num_flt.CHR.filt /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_C_WGS.All/out/All.100.TRS  /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_B_WGS.All/out/All.100.no_read_num_flt.TRS.filt /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_C_WGS.All/out/All.50.INS  /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_B_WGS.All/out/All.50.no_read_num_flt.INS.filt /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_C_WGS.All/RK019_C_WGS.All.fastq.minimap2.merge.sort.bam /archive/data/amed_snt/WORK/fujimoto/nanopore/180615/RK019_B_WGS.All/RK019_B_WGS.All.fastq.minimap2.merge.sort.bam /archive/data/amed_snt/WORK/fujimoto/nanopore/All_fastq/RK019_C_WGS.All.fastq /share/amed_snt/WORK/fujimoto/nanopore/201228/RK019_C_SV_test_new_other_SV_CHR1500/

#[afujrcb:amed_snt@gc019 ~]$ ls /share/amed_snt/WORK/fujimoto/nanopore/201228/RK019_C_WGS_test
#CHR_candidate.txt   CHR_candidate.txt2  DEL_candidate.txt0  INS_candidate.txt   INS_candidate.txt2  INV_candidate.txt0  SV_read.txt    TRS_candidate.txt   TRS_candidate.txt2
#CHR_candidate.txt0  DEL_candidate.txt   DEL_candidate.txt2  INS_candidate.txt0  INV_candidate.txt   INV_candidate.txt2  SV_read.txt.2  TRS_candidate.txt0
#[afujrcb:amed_snt@gc019 ~]$ ls /share/amed_snt/WORK/fujimoto/nanopore/201228/RK019_B_WGS_test
#CHR_candidate.txt   CHR_candidate.txt2  DEL_candidate.txt0  INS_candidate.txt   INS_candidate.txt2  INV_candidate.txt0  SV_read.txt    TRS_candidate.txt   TRS_candidate.txt2
#CHR_candidate.txt0  DEL_candidate.txt   DEL_candidate.txt2  INS_candidate.txt0  INV_candidate.txt   INV_candidate.txt2  SV_read.txt.2  TRS_candidate.txt0


CANCER_DEL=$CANCER_OUTPUT/DEL_candidate.txt2
NORMAL_DEL=$NORMAL_OUTPUT/DEL_candidate.txt0
CANCER_INV=$CANCER_OUTPUT/INV_candidate.txt2
NORMAL_INV=$NORMAL_OUTPUT/INV_candidate.txt0
CANCER_CHR=$CANCER_OUTPUT/CHR_candidate.txt2
NORMAL_CHR=$NORMAL_OUTPUT/CHR_candidate.txt0
CANCER_TRS=$CANCER_OUTPUT/TRS_candidate.txt2
NORMAL_TRS=$NORMAL_OUTPUT/TRS_candidate.txt0
CANCER_INS=$CANCER_OUTPUT/INS_candidate.txt2
NORMAL_INS=$NORMAL_OUTPUT/INS_candidate.txt0
#CANCER_BAM
#NORMAL_BAM=${12}
#FASTQ=${13}
#OUTPUT=${14}

#CANCER_DEL=$1
#NORMAL_DEL=$2
#CANCER_INV=$3
#NORMAL_INV=$4
#CANCER_CHR=$5
#NORMAL_CHR=$6
#CANCER_TRS=$7
#NORMAL_TRS=$8
#CANCER_INS=$9
#NORMAL_INS=${10}
#CANCER_BAM=${11}
#NORMAL_BAM=${12}
#FASTQ=${13}
#OUTPUT=${14}

#echo CANCER_DEL $CANCER_DEL
#echo NORMAL_DEL $NORMAL_DEL
#echo CANCER_INV $CANCER_INV
#echo NORMAL_INV $NORMAL_INV
#echo CANCER_CHR	$CANCER_CHR
#echo NORMAL_CHR $NORMAL_CHR
#echo CANCER_TRS $CANCER_TRS
#echo NORMAL_TRS	$NORMAL_TRS
#echo CANCER_INS $CANCER_INS
#echo NORMAL_INS $NORMAL_INS
#echo CANCER_BAM $CANCER_BAM
#echo NORMAL_BAM $NORMAL_BAM
#echo FASTQ $FASTQ
#echo OUTPUT $OUTPUT

echo "Deletion"
#echo "python ./src/compare_with_normal.py $NORMAL_DEL $CANCER_DEL 300 $NORMAL_BAM 7|python ./src/normal_filter.py /dev/stdin 7 8 10 > $OUTPUT/somatic_DEL_candidate1"
python3 ./src/compare_with_normal.py $NORMAL_DEL $CANCER_DEL 300 $NORMAL_BAM 7|python3 ./src/normal_filter.py /dev/stdin 7 8 10 > $OUTPUT/somatic_DEL_candidate1

#echo "python ./src/somatic_DEL_filt.py $OUTPUT/somatic_DEL_candidate1 30 5000 10|python ./src/somatic_DEL_filter2.py /dev/stdin 20 0.1 50 0.3 0.3,1000,0.1 0.1 800 10 10 2 0 1 > $OUTPUT/somatic_DEL_candidate2"
python3 ./src/somatic_DEL_filt.py $OUTPUT/somatic_DEL_candidate1 30 5000 10|python3 ./src/somatic_DEL_filter2.py /dev/stdin 20 0.1 50 0.3 0.3,1000,0.1 0.1 800 10 10 2 0 1 > $OUTPUT/somatic_DEL_candidate2

#echo "python ./src/VAF.py $OUTPUT/somatic_DEL_candidate2 $CANCER_BAM 0 30 1 4 11 > $OUTPUT/somatic_DEL_candidate3"
python3 ./src/VAF.py $OUTPUT/somatic_DEL_candidate2 $CANCER_BAM 0 30 1 4 11 > $OUTPUT/somatic_DEL_candidate3

#echo "python ./src/add_rep3.py $SR $OUTPUT/somatic_DEL_candidate3 2000|python ./src/add_rep3.py $RMSK /dev/stdin 2000 > $OUTPUT/somatic_DEL_candidate4"
python3 ./src/add_rep3.py $SR $OUTPUT/somatic_DEL_candidate3 2000|python3 ./src/add_rep3.py $RMSK /dev/stdin 2000 > $OUTPUT/somatic_DEL_candidate4

#echo "python ./src/somatic_DEL_filter3.py $OUTPUT/somatic_DEL_candidate4 0.8 > $OUTPUT/somatic_DEL_candidate5"
python3 ./src/somatic_DEL_filter3.py $OUTPUT/somatic_DEL_candidate4 0.8 > $OUTPUT/somatic_DEL_candidate5

#echo "python ./src/add_rep.py $SEGDUP $OUTPUT/somatic_DEL_candidate5 > $OUTPUT/somatic_DEL_candidate6"
python3 ./src/add_rep.py $SEGDUP $OUTPUT/somatic_DEL_candidate5 > $OUTPUT/somatic_DEL_candidate6

#echo "python ./src/add_rep_somatic2.py $SELF_CHAIN $OUTPUT/somatic_DEL_candidate6|python ./src/filter_same_repeat.py /dev/stdin -2 -1 > $OUTPUT/somatic_DEL_candidate7"
python3 ./src/add_rep_somatic2.py $SELF_CHAIN $OUTPUT/somatic_DEL_candidate6|python3 ./src/filter_same_repeat.py /dev/stdin -2 -1 > $OUTPUT/somatic_DEL_candidate7

#echo "python ./src/somatic_SB_CH.py $OUTPUT/somatic_DEL_candidate7 300 $CANCER_BAM $SR $RMSK $FASTQ 14 14|python ./src/somatic_DEL_filter4.py /dev/stdin 0,500,1000 4,3,2 14 15 50 0.03 12 18 13 0.3 > $OUTPUT/somatic_DEL_candidate8"
python3 ./src/somatic_SB_CH.py $OUTPUT/somatic_DEL_candidate7 300 $CANCER_BAM $SR $RMSK $FASTQ 14 14|python3 ./src/somatic_DEL_filter4.py /dev/stdin 0,500,1000 4,3,2 14 15 50 0.03 12 18 13 0.3 > $OUTPUT/somatic_DEL_candidate8

#echo "python ./src/compare_with_normal.2.py $NORMAL_DEL $OUTPUT/somatic_DEL_candidate8 100000 0.8 1000 7 > $OUTPUT/somatic_DEL_candidate_filtered.txt"
python3 ./src/compare_with_normal.2.py $NORMAL_DEL $OUTPUT/somatic_DEL_candidate8 100000 0.8 1000 7 > $OUTPUT/somatic_DEL_candidate_filtered.txt

rm $OUTPUT/somatic_DEL_candidate1 $OUTPUT/somatic_DEL_candidate2 $OUTPUT/somatic_DEL_candidate3 $OUTPUT/somatic_DEL_candidate4 $OUTPUT/somatic_DEL_candidate5 $OUTPUT/somatic_DEL_candidate6 $OUTPUT/somatic_DEL_candidate7 $OUTPUT/somatic_DEL_candidate8

echo "Inversion"
#echo "python ./src/INV_filter.py $CANCER_INV|python ./src/INV_filter2.py /dev/stdin 20 500 0.1 0.1 800 1 7 2 0 1 > $OUTPUT/somatic_INV_candidate"
python3 ./src/INV_filter.py $CANCER_INV|python3 ./src/INV_filter2.py /dev/stdin 20 500 0.1 0.1 800 1 7 2 0 1 > $OUTPUT/somatic_INV_candidate

#echo "python ./src/compare_with_normal.py $NORMAL_INV $OUTPUT/somatic_INV_candidate 300 $NORMAL_BAM 8 > $OUTPUT/somatic_INV_candidate1"
python3 ./src/compare_with_normal.py $NORMAL_INV $OUTPUT/somatic_INV_candidate 300 $NORMAL_BAM 8 > $OUTPUT/somatic_INV_candidate1

#echo "python ./src/normal_filter.py $OUTPUT/somatic_INV_candidate1 8 9 10 > $OUTPUT/somatic_INV_candidate2"
python3 ./src/normal_filter.py $OUTPUT/somatic_INV_candidate1 8 9 10 > $OUTPUT/somatic_INV_candidate2

#echo "python ./src/add_rep_somatic3.py $SEGDUP $OUTPUT/somatic_INV_candidate2 > $OUTPUT/somatic_INV_candidate3"
python3 ./src/add_rep_somatic3.py $SEGDUP $OUTPUT/somatic_INV_candidate2 > $OUTPUT/somatic_INV_candidate3

#echo "python ./src/add_rep_somatic4.py $SELF_CHAIN $OUTPUT/somatic_INV_candidate3|python ./src/filter_same_repeat2.py /dev/stdin -4 -3|python ./src/filter_same_repeat2.py /dev/stdin -2 -1 > $OUTPUT/somatic_INV_candidate4"
python3 ./src/add_rep_somatic4.py $SELF_CHAIN $OUTPUT/somatic_INV_candidate3|python3 ./src/filter_same_repeat2.py /dev/stdin -4 -3|python3 ./src/filter_same_repeat2.py /dev/stdin -2 -1 > $OUTPUT/somatic_INV_candidate4

#echo "python ./src/VAF.py $OUTPUT/somatic_INV_candidate4 $CANCER_BAM 0.03 30 0.5 4 10 > $OUTPUT/somatic_INV_candidate5"
python3 ./src/VAF.py $OUTPUT/somatic_INV_candidate4 $CANCER_BAM 0.03 30 0.5 4 10 > $OUTPUT/somatic_INV_candidate5

#echo "python ./src/SB_CH2.py $OUTPUT/somatic_INV_candidate5 $FASTQ 2 13 13 > $OUTPUT/somatic_INV_candidate6"
python3 ./src/SB_CH2.py $OUTPUT/somatic_INV_candidate5 $FASTQ 2 13 13 > $OUTPUT/somatic_INV_candidate6

#echo "python ./src/add_rep3.py $SR $OUTPUT/somatic_INV_candidate6 2000|python ./src/add_rep3.py $RMSK /dev/stdin 2000|python ./src/TRS_filt3.py /dev/stdin 0.8 50|python ./src/INV_filt.py /dev/stdin 0.8 14|python ./src/TRS_filt4.py /dev/stdin 12 0.3 > $OUTPUT/somatic_INV_candidate7"
python3 ./src/add_rep3.py $SR $OUTPUT/somatic_INV_candidate6 2000|python3 ./src/add_rep3.py $RMSK /dev/stdin 2000|python3 ./src/TRS_filt3.py /dev/stdin 0.8 50|python3 ./src/INV_filt.py /dev/stdin 0.8 14|python3 ./src/TRS_filt4.py /dev/stdin 12 0.3 > $OUTPUT/somatic_INV_candidate7

#echo "python ./src/INV_filter3.py $OUTPUT/somatic_INV_candidate7 14|python ./src/compare_with_normal.2.py $NORMAL_INV /dev/stdin 0 0.8 1000 7 > $OUTPUT/somatic_INV_candidate_filtered.txt"
python3 ./src/INV_filter3.py $OUTPUT/somatic_INV_candidate7 14|python3 ./src/compare_with_normal.2.py $NORMAL_INV /dev/stdin 0 0.8 1000 7 > $OUTPUT/somatic_INV_candidate_filtered.txt

rm $OUTPUT/somatic_INV_candidate $OUTPUT/somatic_INV_candidate1 $OUTPUT/somatic_INV_candidate2 $OUTPUT/somatic_INV_candidate3 $OUTPUT/somatic_INV_candidate4 $OUTPUT/somatic_INV_candidate5 $OUTPUT/somatic_INV_candidate6 $OUTPUT/somatic_INV_candidate7 

echo "Chromosoal translocation"
#echo "python ./src/INV_filter.py $CANCER_CHR|python ./src/CHR_filter.py /dev/stdin 20 1500 0.1 0.1 800 1 7 2 0 1 > $OUTPUT/somatic_CHR_candidate"
python3 ./src/INV_filter.py $CANCER_CHR|python3 ./src/CHR_filter.py /dev/stdin 20 1500 0.1 0.1 800 1 7 2 0 1 > $OUTPUT/somatic_CHR_candidate

#echo "python ./src/compare_with_normal.py $NORMAL_CHR $OUTPUT/somatic_CHR_candidate 1500 $NORMAL_BAM 8|python ./src/normal_filter.py /dev/stdin 8 9 10 > $OUTPUT/somatic_CHR_candidate1"
python3 ./src/compare_with_normal.py $NORMAL_CHR $OUTPUT/somatic_CHR_candidate 1500 $NORMAL_BAM 8|python3 ./src/normal_filter.py /dev/stdin 8 9 10 > $OUTPUT/somatic_CHR_candidate1

#echo "python ./src/add_rep_somatic3.py $SEGDUP $OUTPUT/somatic_CHR_candidate1 > $OUTPUT/somatic_CHR_candidate2"
python3 ./src/add_rep_somatic3.py $SEGDUP $OUTPUT/somatic_CHR_candidate1 > $OUTPUT/somatic_CHR_candidate2

#echo "python ./src/add_rep_somatic4.py $SELF_CHAIN $OUTPUT/somatic_CHR_candidate2 > $OUTPUT/somatic_CHR_candidate3"
python3 ./src/add_rep_somatic4.py $SELF_CHAIN $OUTPUT/somatic_CHR_candidate2 > $OUTPUT/somatic_CHR_candidate3

#echo "python ./src/filter_same_repeat2.py $OUTPUT/somatic_CHR_candidate3 -4 -3|python ./src/filter_same_repeat2.py /dev/stdin -2 -1 > $OUTPUT/somatic_CHR_candidate4"
python3 ./src/filter_same_repeat2.py $OUTPUT/somatic_CHR_candidate3 -4 -3|python3 ./src/filter_same_repeat2.py /dev/stdin -2 -1 > $OUTPUT/somatic_CHR_candidate4

#echo "python ./src/VAF.py $OUTPUT/somatic_CHR_candidate4 $CANCER_BAM 0.03 30 0.5 4 10 > $OUTPUT/somatic_CHR_candidate5"
python3 ./src/VAF.py $OUTPUT/somatic_CHR_candidate4 $CANCER_BAM 0.03 30 0.5 4 10 > $OUTPUT/somatic_CHR_candidate5

#echo "python ./src/SB_CH2.py $OUTPUT/somatic_CHR_candidate5 $FASTQ 2 13 13|python ./src/TRS_filt4.py /dev/stdin 12 0.3 > $OUTPUT/somatic_CHR_candidate_filtered.txt"
python3 ./src/SB_CH2.py $OUTPUT/somatic_CHR_candidate5 $FASTQ 2 13 13|python3 ./src/TRS_filt4.py /dev/stdin 12 0.3 > $OUTPUT/somatic_CHR_candidate_filtered.txt

rm $OUTPUT/somatic_CHR_candidate $OUTPUT/somatic_CHR_candidate1 $OUTPUT/somatic_CHR_candidate2 $OUTPUT/somatic_CHR_candidate3 $OUTPUT/somatic_CHR_candidate4 $OUTPUT/somatic_CHR_candidate5

echo "Translocation"
#echo "python ./src/INV_filter.py $CANCER_TRS|python ./src/INV_filter2.py /dev/stdin 20 500 0.1 0.1 800 1 7 2 0 1 > $OUTPUT/somatic_TRS_candidate"
python3 ./src/INV_filter.py $CANCER_TRS|python3 ./src/INV_filter2.py /dev/stdin 20 500 0.1 0.1 800 1 7 2 0 1 > $OUTPUT/somatic_TRS_candidate

#echo "python ./src/compare_with_normal.py $NORMAL_TRS $OUTPUT/somatic_TRS_candidate 300 $NORMAL_BAM 8|python ./src/normal_filter.py /dev/stdin 8 9 10 > $OUTPUT/somatic_TRS_candidate1"
python3 ./src/compare_with_normal.py $NORMAL_TRS $OUTPUT/somatic_TRS_candidate 300 $NORMAL_BAM 8|python3 ./src/normal_filter.py /dev/stdin 8 9 10 > $OUTPUT/somatic_TRS_candidate1

#echo "python ./src/add_rep_somatic3.py $SEGDUP $OUTPUT/somatic_TRS_candidate1 > $OUTPUT/somatic_TRS_candidate2"
python3 ./src/add_rep_somatic3.py $SEGDUP $OUTPUT/somatic_TRS_candidate1 > $OUTPUT/somatic_TRS_candidate2

#echo "python ./src/add_rep_somatic4.py $SELF_CHAIN $OUTPUT/somatic_TRS_candidate2 > $OUTPUT/somatic_TRS_candidate3"
python3 ./src/add_rep_somatic4.py $SELF_CHAIN $OUTPUT/somatic_TRS_candidate2 > $OUTPUT/somatic_TRS_candidate3

#echo "python ./src/filter_same_repeat2.py $OUTPUT/somatic_TRS_candidate3 -4 -3|python ./src/filter_same_repeat2.py /dev/stdin -2 -1 > $OUTPUT/somatic_TRS_candidate4"
python3 ./src/filter_same_repeat2.py $OUTPUT/somatic_TRS_candidate3 -4 -3|python3 ./src/filter_same_repeat2.py /dev/stdin -2 -1 > $OUTPUT/somatic_TRS_candidate4

#echo "python ./src/VAF.py $OUTPUT/somatic_TRS_candidate4 $CANCER_BAM 0.03 30 0.5 4 10 > $OUTPUT/somatic_TRS_candidate5"
python3 ./src/VAF.py $OUTPUT/somatic_TRS_candidate4 $CANCER_BAM 0.03 30 0.5 4 10 > $OUTPUT/somatic_TRS_candidate5

#echo "python ./src/SB_CH2.py $OUTPUT/somatic_TRS_candidate5 $FASTQ 2 13 13 > $OUTPUT/somatic_TRS_candidate6"
python3 ./src/SB_CH2.py $OUTPUT/somatic_TRS_candidate5 $FASTQ 2 13 13 > $OUTPUT/somatic_TRS_candidate6

#echo "python ./src/add_rep3.py $SR $OUTPUT/somatic_TRS_candidate6 2000|python ./src/add_rep3.py $RMSK /dev/stdin 2000 > $OUTPUT/somatic_TRS_candidate7"
python3 ./src/add_rep3.py $SR $OUTPUT/somatic_TRS_candidate6 2000|python3 ./src/add_rep3.py $RMSK /dev/stdin 2000 > $OUTPUT/somatic_TRS_candidate7

#echo "python ./src/TRS_filt3.py $OUTPUT/somatic_TRS_candidate7 0.8 50 > $OUTPUT/somatic_TRS_candidate8"
python3 ./src/TRS_filt3.py $OUTPUT/somatic_TRS_candidate7 0.8 50 > $OUTPUT/somatic_TRS_candidate8

#echo "python ./src/TRS_filt4.py $OUTPUT/somatic_TRS_candidate8 12 0.3 > $OUTPUT/somatic_TRS_candidate9"
python3 ./src/TRS_filt4.py $OUTPUT/somatic_TRS_candidate8 12 0.3 > $OUTPUT/somatic_TRS_candidate9

#echo "python ./src/compare_with_normal.2.py $NORMAL_TRS $OUTPUT/somatic_TRS_candidate9 0 0.8 1000 7 > $OUTPUT/somatic_TRS_candidate_filtered.txt"
python3 ./src/compare_with_normal.2.py $NORMAL_TRS $OUTPUT/somatic_TRS_candidate9 0 0.8 1000 7 > $OUTPUT/somatic_TRS_candidate_filtered.txt

rm $OUTPUT/somatic_TRS_candidate $OUTPUT/somatic_TRS_candidate1 $OUTPUT/somatic_TRS_candidate2 $OUTPUT/somatic_TRS_candidate3 $OUTPUT/somatic_TRS_candidate4 $OUTPUT/somatic_TRS_candidate5 $OUTPUT/somatic_TRS_candidate6 $OUTPUT/somatic_TRS_candidate7 $OUTPUT/somatic_TRS_candidate8 $OUTPUT/somatic_TRS_candidate9

echo "Insertion"
#echo "python ./src/INS_filt.py $CANCER_INS 100|python ./src/INS_filt1.py /dev/stdin 20 50 0.1 800 1 8 2 1 1 > $OUTPUT/somatic_INS_candidate"
python3 ./src/INS_filt.py $CANCER_INS 100|python3 ./src/INS_filt1.py /dev/stdin 20 50 0.1 800 1 8 2 1 1 > $OUTPUT/somatic_INS_candidate 

#echo "python ./src/VAF.INS.py $OUTPUT/somatic_INS_candidate $CANCER_BAM 0.03 30 2 4 9 > $OUTPUT/somatic_INS_candidate1"
python3 ./src/VAF.INS.py $OUTPUT/somatic_INS_candidate $CANCER_BAM 0.03 30 2 4 9 > $OUTPUT/somatic_INS_candidate1

#echo "python ./src/INS_filt2.py $OUTPUT/somatic_INS_candidate1 1000 3 2,2 0.03 0.3 > $OUTPUT/somatic_INS_candidate2"
python3 ./src/INS_filt2.py $OUTPUT/somatic_INS_candidate1 1000 3 2,2 0.03 0.3 > $OUTPUT/somatic_INS_candidate2

#echo "python ./src/compare_with_normal.INS.py $NORMAL_INS $OUTPUT/somatic_INS_candidate1 500 0.5 3 1000 $NORMAL_BAM 11 $SAMTOOLS > $OUTPUT/somatic_INS_candidate2"
python3 ./src/compare_with_normal.INS.py $NORMAL_INS $OUTPUT/somatic_INS_candidate1 500 0.5 3 1000 $NORMAL_BAM 11 $SAMTOOLS > $OUTPUT/somatic_INS_candidate2

#echo "python ./src/INS_filt3.py $OUTPUT/somatic_INS_candidate2 10 0.3 4 > $OUTPUT/somatic_INS_candidate3"
python3 ./src/INS_filt3.py $OUTPUT/somatic_INS_candidate2 10 0.3 4 > $OUTPUT/somatic_INS_candidate3

#echo "python ./src/add_rep.INS.py $SR $OUTPUT/somatic_INS_candidate3 500 > $OUTPUT/somatic_INS_candidate4"
python3 ./src/add_rep.INS.py $SR $OUTPUT/somatic_INS_candidate3 500 > $OUTPUT/somatic_INS_candidate4

#echo "python ./src/add_rep.INS2.py $RMSK $OUTPUT/somatic_INS_candidate4 2000 > $OUTPUT/somatic_INS_candidate5"
python3 ./src/add_rep.INS2.py $RMSK $OUTPUT/somatic_INS_candidate4 2000 > $OUTPUT/somatic_INS_candidate5

#echo "python ./src/INS_filt4.py $OUTPUT/somatic_INS_candidate5 2 -2|python ./src/INS_filt4.py /dev/stdin 1 -1 > $OUTPUT/somatic_INS_candidate_filtered.txt"
python3 ./src/INS_filt4.py $OUTPUT/somatic_INS_candidate5 2 -2|python3 ./src/INS_filt4.py /dev/stdin 1 -1 > $OUTPUT/somatic_INS_candidate_filtered.txt

rm $OUTPUT/somatic_INS_candidate $OUTPUT/somatic_INS_candidate1 $OUTPUT/somatic_INS_candidate2 $OUTPUT/somatic_INS_candidate3 $OUTPUT/somatic_INS_candidate4 $OUTPUT/somatic_INS_candidate5

#echo "python ./src/vcf.py $OUTPUT/somatic_DEL_candidate_filtered.txt $OUTPUT/somatic_INS_candidate_filtered.txt $OUTPUT/somatic_INV_candidate_filtered.txt $OUTPUT/somatic_TRS_candidate_filtered.txt $OUTPUT/somatic_CHR_candidate_filtered.txt > $OUTPUT/somatic_SV.vcf"
python3 ./src/vcf.py $OUTPUT/somatic_DEL_candidate_filtered.txt $OUTPUT/somatic_INS_candidate_filtered.txt $OUTPUT/somatic_INV_candidate_filtered.txt $OUTPUT/somatic_TRS_candidate_filtered.txt $OUTPUT/somatic_CHR_candidate_filtered.txt > $OUTPUT/somatic_SV.vcf

echo "Output file;" $OUTPUT/somatic_SV.vcf
