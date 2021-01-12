#chr10   196947  chr10   196947  2       INS     159.0   2,0,0   2,2,2   21,23   0.091   0,21/0,23       aec8eb37-c388-4650-a17e-303cd67f3f04;within;chr10;196947;chr10;196
import sys
import re

f = open(sys.argv[1])
length_border = int(sys.argv[2])
filt1 = sys.argv[3]
filt2 = sys.argv[4]
VAF_cutoff = float(sys.argv[5])
low_mq_read_ratio_cutoff = float(sys.argv[6])


filt1_l = filt1.split(",")
filt2_l = filt2.split(",")

for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	read_num_l = line_l[7].split(",")
	
	length = line_l[6]
	length = length.replace(">", "")

	VAF = float(line_l[10])
	if VAF < VAF_cutoff:
#		print("VAF", line)
		continue

	low_mq_read_ratio_l = line_l[11].split("/")
	low_mq_read_ratio_l_1 = low_mq_read_ratio_l[0].split(",")
	low_mq_read_ratio_l_2 = low_mq_read_ratio_l[1].split(",")
	ratio1, ratio2 = 0, 0
	if float(low_mq_read_ratio_l_1[1]) > 0:
		ratio1 = float(low_mq_read_ratio_l_1[0])/float(low_mq_read_ratio_l_1[1])

	if float(low_mq_read_ratio_l_2[1]) > 0:
		ratio2 = float(low_mq_read_ratio_l_2[0])/float(low_mq_read_ratio_l_2[1])

	if ratio1 >= low_mq_read_ratio_cutoff or ratio2 >= low_mq_read_ratio_cutoff:
#		print("low_mq_read_ratio_cutoff", line)
		continue

	if float(length) <= length_border:
		if int(read_num_l[0]) >= int(filt1_l[0]):
			print(line)
#		else:
#			print("smaller", line)
	else:
		if int(read_num_l[0]) >= int(filt2_l[0]) or (int(read_num_l[1]) >= int(filt2_l[1]) and int(read_num_l[2]) >= int(filt2_l[1])):
			print(line)
#		else:
#			print("larger", line)
