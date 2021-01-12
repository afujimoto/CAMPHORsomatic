import sys

#chr10   930674  chr10   934257  2       DEL     3583    20,26   -       0       2,2,2   14,14   0.143   0,14/0,14       2       0       0/2     |||     54460fbf-e16d-432e-9709-819efedf7d06;within;chr10;930602;chr10;934162;8298;46;8278;926799;938174;15;60;-;-;-;560|7c74f7b2-ce59-42c0-8983-575eb1e8eb3e;within;chr10;930674;chr10;934257;9153;32;9114;926386;938475;17.3;60;-;-;-;652     chr10;930602;chr10;934162;GATCTCACGTTTGAAAGAATG;3560;54460fbf-e16d-432e-9709-819efedf7d06;46,8278;chr10;8298;926799,938174;-;15;60;20S13M5D2M1D23M1D2M2D1M1I26M1I7M1D10M2I17M3I2M1I10M1I5M3D11M1I6M1I17M1I25M1D5M2D9M1I11M2I15M3D5M7D17M1I15M2D15M1D6M1I23M4

f = open(sys.argv[1])
len_cutoff_l = sys.argv[2].split(",")
min_dis_cutoff_l = sys.argv[3].split(",")
support_read_col = int(sys.argv[4])
strand_bias_col = int(sys.argv[5])
distance_btw_BP = int(sys.argv[6])
VAF_cutoff = float(sys.argv[7])
VAF_col = int(sys.argv[8])
segment_col = int(sys.argv[9])
low_mq_read_col = int(sys.argv[10])
low_mq_read_prop = float(sys.argv[11])

#print("len_cutoff_l", len_cutoff_l)

len_cutoff_l = [int(i) for i in len_cutoff_l]
min_dis_cutoff_l = [int(i) for i in min_dis_cutoff_l]

for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

#	print("########################")
#	print(line)
#	print(line_l[VAF_col], VAF_cutoff)
	if float(line_l[VAF_col]) < VAF_cutoff:
		continue

	low_mq_prop_l = line_l[low_mq_read_col].split("/")
	low_mq_prop_l_1 = low_mq_prop_l[0].split(",")
	low_mq_prop_l_2 = low_mq_prop_l[1].split(",")

	total_read_num1, total_read_num2 = low_mq_prop_l_1[1], low_mq_prop_l_2[1]
	if float(total_read_num1) == 0:
		total_read_num1 = 1
	if float(total_read_num2) == 0:
		total_read_num2 = 1

#	print(line)
#	print("low_mq_prop_l", low_mq_prop_l)
#	print("total_read_num1", total_read_num1, "total_read_num2", total_read_num2)

#	if int(low_mq_prop_l_1[0]) >= 3 and int(low_mq_prop_l_2[0]) >= 3 and float(low_mq_prop_l_1[0])/float(low_mq_prop_l_1[1]) >= low_mq_read_prop or float(low_mq_prop_l_2[0])/float(low_mq_prop_l_2[1]) >= low_mq_read_prop:
	if int(low_mq_prop_l_1[0]) >= 3 and int(low_mq_prop_l_2[0]) >= 3 and float(low_mq_prop_l_1[0])/float(total_read_num1) >= low_mq_read_prop or float(low_mq_prop_l_2[0])/float(total_read_num2) >= low_mq_read_prop:
		continue

	support_read_filt = 0
	for i in range(len(len_cutoff_l) - 1, -1, -1):
		if int(line_l[6]) >= len_cutoff_l[i] and int(line_l[support_read_col]) >= min_dis_cutoff_l[i]:
			support_read_filt += 1
			break

#	print("support_read_filt", support_read_filt, line_l[strand_bias_col])

	if support_read_filt == 0:
		continue
	
	if int(line_l[strand_bias_col]) == 1:
		continue

#54460fbf-e16d-432e-9709-819efedf7d06;within;chr10;930602;chr10;934162;8298;46;8278;926799;938174;15;60;-;-;-;560
	overlap_read_pair = 0
	segment_l = line_l[segment_col].split('|')
	for i in range(len(segment_l)):
		segment_tmp_l1 = segment_l[i].split(";")
		for j in range(i + 1, len(segment_l)):
			segment_tmp_l2 = segment_l[j].split(";")
#			print("1", segment_tmp_l1)
#			print("2", segment_tmp_l2)
			if abs(int(segment_tmp_l1[3]) - int(segment_tmp_l2[3])) <= distance_btw_BP and abs(int(segment_tmp_l1[5]) - int(segment_tmp_l2[5])) <= distance_btw_BP:
				overlap_read_pair += 1

#	print("overlap_read_pair", overlap_read_pair)

	if overlap_read_pair >= 1:
		print(line)
