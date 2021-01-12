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

f = open(sys.argv[1])
min_rep_prop = float(sys.argv[2])

for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)

#	print(line)

	rep_prop_trf, max_rep_trf, max_rep_length_trf = compare_with_repeat(line_l, -2)
	rep_prop_rmsk, max_rep_rmsk, max_rep_length_rmsk = compare_with_repeat(line_l, -1)

	if rep_prop_trf < min_rep_prop and rep_prop_rmsk < min_rep_prop:
		print(line)
#	else:
#		print(line_l[-2])
#		print("rep_prop_trf", rep_prop_trf)

#		print(line_l[-1])
#		print("rep_prop_rmsk", rep_prop_rmsk)
