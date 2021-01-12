import sys
import re

def compare_with_repeat(line_l, rep_col):
	SR_l = re.split('\|', line_l[rep_col])

	max_prop = 0
	max_rep = "NA"
	max_rep_length = 0
	if line_l[rep_col] != "-":
		for sr_tmp in SR_l:
			sr_tmp_l = re.split(",", sr_tmp)
			if float(sr_tmp_l[7]) > float(max_prop):
				max_prop = float(sr_tmp_l[7])
				max_rep = len(sr_tmp_l[6])
				max_rep_length = int(sr_tmp_l[2]) - int(sr_tmp_l[1])

	start = int(line_l[1])
	end = int(line_l[3])
	if int(start) > int(end):
		start, end = end, start

	covered_length = 0
	if max_prop == 1:
		max_rep_length = int(end) - int(start)
		covered_length = int(end) - int(start)
	else:
		if line_l[rep_col] != "-":
			boundary = []
			rep_start = []
			rep_end = []
			for sr_tmp in SR_l:
				sr_tmp_l = re.split(",", sr_tmp)
				if int(sr_tmp_l[1]) < int(start):
					rep_start.append(int(line_l[1]))
				else:
					rep_start.append(int(sr_tmp_l[1]))
				if int(end) < int(sr_tmp_l[2]):
					rep_end.append(int(line_l[3]))
				else:
					rep_end.append(int(sr_tmp_l[2]))
                        
			boundary.extend(rep_start)
			boundary.extend(rep_end)
        
			boundary.sort()
			for i in range(len(boundary) - 1):
				covered = 0
				for j in range(len(rep_start)):
					if rep_start[j] <= boundary[i] <= rep_end[j] and rep_start[j] <= boundary[i + 1] <= rep_end[j]:
						covered = 1
						break
				if covered == 1:
					covered_length += boundary[i + 1] - boundary[i]

	rep_prop = 0
	if abs(int(end) - int(start)):
		rep_prop = round(float(covered_length)/float(abs(int(end) - int(start))), 2)
	
	return rep_prop, max_rep, max_rep_length

def get_max_length(rep):
	if rep == "-":
		return 0

	rep1_l = rep.split('|')
	max_length = 0
	for tmp in rep1_l:
		tmp_l = tmp.split(",")
		length = int(tmp_l[2]) - int(tmp_l[1]) + 1
		if length >= max_length:
			max_length = length

	return max_length

f = open(sys.argv[1])
rep_len_ratio = float(sys.argv[2])
rep_col = int(sys.argv[3])
#min_rep_prop = float(sys.argv[2])

for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)

#	rep_prop_trf, max_rep_trf, max_rep_length_trf = compare_with_repeat(line_l, -2)
#	rep_prop_rmsk, max_rep_rmsk, max_rep_length_rmsk = compare_with_repeat(line_l, -1)

#chr5	157682248	chr5	157682248	19	INS	>8358
#chr12,27240377,27241846,trf,39,41.4,TATATATATTACATATATGACTATGTATAGTTGTATAGT,1.0,complete|chr12,27240642,27241828,trf,31,36.5,TGTATAGTTATATATATTACATATATGACTC,1	chr5,157682245,157682272,(T)n,1788,28,0,0,0,chr5,157682245,157682272,-23855987,+,(T)n,Simple_repeat,Simple_repeat,1,27,0,3,1.0,complete

#	rep1_l = line_l[-2].split('|')
#	rep2_l = line_l[-1].split('|')

	max_rep_len1 = get_max_length(line_l[rep_col])
#	max_rep_len2 = get_max_length(line_l[-1])

#	print("max_rep_len1", max_rep_len1)
#	print("max_rep_len2", max_rep_len2)

	INS_length = line_l[6]
	INS_length = INS_length.replace(">", "")

#	if float(INS_length) >= max_rep_len1*rep_len_ratio and float(INS_length) >= max_rep_len2*rep_len_ratio:
	if float(INS_length) >= max_rep_len1*rep_len_ratio:
		print(line)
