#chr10   83897886        chr10   83897882        2       TRS     -3      25,25   20,19   0.103   1/2     0/19,0/19       0       -       -       2       NS;2/DIS;2/MQ;2/CN;2/RL;2
#chr10   115808358       chr10   115808267       2       TRS     -90     29,28   26,27   0.075   1/2     0/26,0/24       0       -       -       2       NS;2/DIS;2/MQ;2/CN;2/RL;2
#chr10   131176224       chr10   131177158       2       TRS     935     23,24   16,19   0.114   2/2     0/17,0/20       0       -       -       2       NS;2/DIS;2/MQ;2/CN;2/RL;2
import sys

#chr8,76672371,76672439,trf,34,2,GCAACATCTTCCTTTATCCTTTTGGAACTCCAGA,0.16,within|chr8,76672483,76672529,trf,2,23,AC,0.11,within chr8,76672483,76672529,(AC)n,1169,54,0,0,0,chr8,76672483,76672529,-68466107,+,(AC)n,Simple_repeat,Simple_repeat,1,46,0,3,0.11,within

f = open(sys.argv[1])
min_prop_rep = float(sys.argv[2])
min_length = int(sys.argv[3])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if int(line_l[6]) <= min_length:
		continue

	rep1_filt = 0
	if line_l[-2] != "-":
		rep1 = line_l[-2].split('\|')
		for rep1_tmp in rep1:
			rep1_tmp_l = rep1_tmp.split(",")
#			print(rep1_tmp_l)
			if float(rep1_tmp_l[-2]) >= min_prop_rep:
				rep1_filt = 1
			
	rep2_filt = 0
	if line_l[-1] != "-":
		rep2 = line_l[-1].split('\|')
		for rep2_tmp in rep2:
			rep2_tmp_l = rep2_tmp.split(",")
#			print(rep2_tmp_l)
			if float(rep2_tmp_l[-2]) >= min_prop_rep:
				rep2_filt = 1

	
	if rep1_filt == 0 and rep2_filt == 0:
		print(line)
#	else:
#		print("FILTER!!", line)
