import re
import sys
from operator import itemgetter
import statistics

#chr10   37451   chr10   37564   2       DEL     113     0       chr10;37409;chr10;37564;TTGAGAGCCTGTCCATTGATCTATATCCTGTTTTAAACAGTACCAGTTGTTTTAATATCACAGCCTTGCATAAATAGTCAGGCAGTATTCCAGCTGTTCTGGTTTAGACTGACTTATACGCGGCCTTTTTGTCCATATAATTTCAGCGTTTTCCACCCTGCATAACCATTATAGTTTGACGAGACAAATTGACCTAAATTACCTTGGCAAAATATTTTCAATATTGATTCTTCCTAACGCATGAGCATAATGTTCTTCATTTGTTTGTATCCTATTCTGGAGCAGGTTTGTCTCCTGAAGAGTATCCCACATCCCTGTTATCTTATGTATTTACCTTTGCGCAATCAGAATGGGAGTTCATCTGACGGCCTGCTTCATCAGTGGCGGCAGCATCGTTTACT;155;653b180d-754b-4b4a-89b7-af1b4885fd65_Basecall_1D_template;3445,5339/340,1949/11,339/1974,2900;chr10/chr10/chr10/chr10;5393;109799256,109801664/35682,37761/109795057,109795430/109797362,109798491;-/+/+/+;10.3/10.3/10.3/10.3;11/0/60/5;3444S13M4D12M10D14M1D13M1D36M2D17M1D21M1D9M4D10M1I3M2D14M2D2M2D4M154D2M1D8M2D3M1D1M1D2M1D4M2I3M3D20M2D5M1D3M2I2M5D3M1D2M4D7M2D5M1D3M3D5M1D7M5D3M1D9M1D7M3D6M2D7M1D5M5D5M3D2M3D14M3D8M4D11M1D4M1D13M1I5M4D6M1D12M8D5M1D24M5D7M3D6M1D10M1I1M2I1M3D1M4D6M1D8M2D1M1D9M10D1M1D6M4I7M2D1M1D3M3D2M1D19M4D1M1D8M4D23M3D4M3D7M9D6M2D22M1D5M1D11M4D5M2I15M2D28M2D9M1D3M3D7M4D7M3D4M3D2M2D6M1D26M1D6M2D16M1D24M1D4M1D3M2I15M1D7M2I7M2D9M2D4M1D3M2D17M2D2M2D4M6D9M3D9M1I2M3D3M2D6M1I1M1I22M1D9M3D2M2D2M1D15M2D5M1D12M1D9M3D17M2D2M2D2M1D2M1D9M1D8M2D3M2D2M2D3M4D11M3D6M1D13M1D5M2D2M1I23M2D7M1D4M2D1M1D3M4D9M1D3M4D12M1D2M1D19M1D22M5D2M2D7M1D8M3I8M1D6M1D3M1D6M1D6M1D31M1I6M1D5M1I9M1D8M2I1M1I3M2D6M2D38M1D1M2D3M2D19M1D11M2I28M1D9M2D6M2I7M1D9M1D17
def compare_with_repeat(line_l, rep_col):
	SR_l = re.split('\|', line_l[rep_col])

#['chr10', '64273', '64367', 'trf', '49', '1.9', 'CCTAGCCCTAGCCCTAACCCTAACCCTAACCCCAAACCCCAACCCTAGT', '0.24', 'within']
	max_prop = 0
	max_rep = "NA"
	max_rep_length = 0
	if line_l[rep_col] != "-":
		for sr_tmp in SR_l:
			sr_tmp_l = re.split(",", sr_tmp)
#			print("sr_tmp_l", sr_tmp_l)
#			print(sr_tmp_l[7])
			if float(sr_tmp_l[7]) > float(max_prop):
				max_prop = float(sr_tmp_l[7])
				max_rep = len(sr_tmp_l[6])
				max_rep_length = int(sr_tmp_l[2]) - int(sr_tmp_l[1])

#	print("max_prop>>>", max_prop)

	start = int(line_l[1])
	end = int(line_l[3])

	if int(start) > int(end):
		start, end = end, start

	covered_length = 0
	if max_prop == 1:
#		print("max_prop>>> 1 !!!!")
		max_rep_length = int(end) - int(start)
		covered_length = int(end) - int(start)
	else:
#		print("max_prop>>>i>>>>>>>>>>>>", max_prop)
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
					if rep_start[j]	<= boundary[i] <= rep_end[j] and rep_start[j] <= boundary[i + 1] <= rep_end[j]:
						covered = 1
						break
				if covered == 1:
					covered_length += boundary[i + 1] - boundary[i]

#	print("covered_length => ", covered_length)

	rep_prop = 0
#	if float(line_l[6]) > 0:
	if abs(int(end) - int(start)):
#		rep_prop = round(float(covered_length)/float(line_l[6]), 2)
		rep_prop = round(float(covered_length)/float(abs(int(end) - int(start))), 2)

#	print(rep_prop, max_rep, max_rep_length)
	return rep_prop, max_rep, max_rep_length

#def find_close_ins(cigar, mapping_start, del_pos, del_size, close_ins_distance, close_ins_len):
def find_close_ins(cigar, mapping_start, del_pos, del_size, close_ins_distance):
	cigar_type1 = re.split(r'\d+', cigar)
	cigar_len1 = re.split(r'[a-zA-Z]+', cigar)
	del cigar_type1[0]
	del cigar_len1[-1]

	ins_pos = {}
	target_del_found = 0
	map_end = mapping_start - 1
	for i in range(0, len(cigar_type1)):
		if cigar_type1[i] == "I" or cigar_type1[i] == "S" or cigar_type1[i] == "H":
			if cigar_type1[i] == "I":
				ins_pos[map_end] = cigar_len1[i]
			continue

		if map_end + 1 == int(del_pos) and int(cigar_len1[i]) == del_size:
			target_del_found = 1

		map_end += int(cigar_len1[i])

	close_ins_l = []
	for ins_pos_tmp in ins_pos:
		if 0 <= del_pos - ins_pos_tmp <= close_ins_distance or 0 <= ins_pos_tmp - (del_pos + del_size) <= close_ins_distance:
			close_ins_l.append(int(ins_pos[ins_pos_tmp]))
#			total_ins_len += int(ins_pos[ins_pos_tmp])

	max_close_ins = 0
	if len(close_ins_l) > 0:
		max_close_ins = max(close_ins_l)

	if max_close_ins:
		return max_close_ins
	elif target_del_found == 1:
		return 0
	else:
		return -1

def calculate_edit_distance(cigar, max_indel_len):
#	print("cigar => ", cigar)
	cigar_type1 = re.split(r'\d+', cigar)
	cigar_len1 = re.split(r'[a-zA-Z]+', cigar)
	del cigar_type1[0]
	del cigar_len1[-1]

	read_length = 0
	indel_dis = 0
	for i in range(0, len(cigar_type1)):
		if cigar_type1[i] == "S" or cigar_type1[i] == "H":
			continue
		
		if cigar_type1[i] == "M" or cigar_type1[i] == "I":
			read_length += int(cigar_len1[i])

		if cigar_type1[i] == "I" or cigar_type1[i] == "D":
			if int(cigar_len1[i]) <= max_indel_len:
				indel_dis += int(cigar_len1[i])

	return read_length, indel_dis

f = open(sys.argv[1])
#min_distance_prop = float(sys.argv[2])
#min_distance = float(sys.argv[3])
#min_mq = float(sys.argv[4])
#r_dis_g_dis_ratio = float(sys.argv[5])
close_ins_distance = int(sys.argv[2]) #30
#close_ins_len = float(sys.argv[3]) #20
max_del_len_for_close_ins_analysis = int(sys.argv[3]) #5000
read_col_num = int(sys.argv[4])

for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)
	read_l = re.split(':', line_l[read_col_num])

#	print(read_l)

	strange_read = []
	breakpoint_info = []
	for j in range(len(read_l)):
		read_tmp_l = re.split(";", read_l[j])
		mq_bq_found = 0
		if "/" in read_tmp_l[8]:
			read_pos = re.split("/", read_tmp_l[7])
			chr_l = re.split("/", read_tmp_l[8])
			genome_pos = re.split("/", read_tmp_l[10])
			bq_l = re.split("/", read_tmp_l[12])
			mq_l = re.split("/", read_tmp_l[13])
			cigar_l = re.split("/", read_tmp_l[14])
			strand_l = re.split("/", read_tmp_l[11])
			mapper_l = re.split("/", read_tmp_l[15])

			seq_align_pos = []
			for i in range(len(read_pos)):
				read_pos_tmp_l = re.split(",", read_pos[i])
				genome_pos_tmp_l = re.split(",", genome_pos[i])
				seq_align_pos.append([int(read_pos_tmp_l[0]), int(read_pos_tmp_l[1]), int(genome_pos_tmp_l[0]), int(genome_pos_tmp_l[1]), chr_l[i], bq_l[i], mq_l[i], cigar_l[i], strand_l[i], mapper_l[i]])

			seq_align_pos.sort(key=itemgetter(0))

			for i in range(len(seq_align_pos)):
#				print(i, seq_align_pos[i])
				r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand, mapper = seq_align_pos[i]
				if chr == read_tmp_l[0] and chr == read_tmp_l[2] and int(g_start) <= int(read_tmp_l[1]) <= int(g_end) and int(g_start) <= int(read_tmp_l[3]) <= int(g_end):
				
					total_ins_len = -1 #180319
					if int(read_tmp_l[5]) <= max_del_len_for_close_ins_analysis:
						total_ins_len = find_close_ins(cigar, g_start, int(read_tmp_l[1]), int(read_tmp_l[3]) - int(read_tmp_l[1]), close_ins_distance)
					if total_ins_len == -1:
						continue

					mq_bq_found = 1
#					breakpoint_info.append([read_tmp_l[6], "within", chr, read_tmp_l[1], chr, read_tmp_l[3], r_start, r_end, g_start, g_end, bq, mq, cigar, strand, mapper, total_ins_len, cigar])
					breakpoint_info.append([read_tmp_l[6], "within", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], r_start, r_end, g_start, g_end, bq, mq, strand, mapper, "-", total_ins_len])
					break

				if i < len(seq_align_pos) - 1:
					r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar2, strand2, mapper2 = seq_align_pos[i + 1]
					if chr2 == read_tmp_l[0]  and chr2 == read_tmp_l[2] and int(g_end) == int(read_tmp_l[1]) and int(g_start2) == int(read_tmp_l[3]):
						g_dis = int(g_start2) - int(g_end)
						r_dis = int(r_start2) - int(r_end)

						r_dis_g_dis_ratio = abs(float(r_dis)/float(g_dis))
#						print(g_dis, r_dis, r_dis_g_dis_ratio)
#						breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr, read_tmp_l[3], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper]), r_dis_g_dis_ratio, "-", ",".join([cigar, cigar2])])
						breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), r_dis_g_dis_ratio, "-"])

						mq_bq_found = 1
						break

			if mq_bq_found == 0:
				for i in range(len(seq_align_pos)):
					r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand, mapper = seq_align_pos[i]
					for k in range(len(seq_align_pos)):
						if i == k:
							continue
						r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar2, strand2, mapper2 = seq_align_pos[k]
						if chr2 == read_tmp_l[0]  and chr2 == read_tmp_l[2] and int(g_end) == int(read_tmp_l[1]) and int(g_start2) == int(read_tmp_l[3]):
							g_dis = int(g_start2) - int(g_end)

							r_dis = int(r_start2) - int(r_end)
							if strand2 == "-":
								r_dis = int(r_start) - int(r_end2)
							
							r_dis_g_dis_ratio = abs(float(r_dis)/float(g_dis))
#							print(r_start, r_end, r_start2, r_end2)
#							print(g_dis, r_dis, r_dis_g_dis_ratio)
#							breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr, read_tmp_l[3], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([str(bq), str(bq2)]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper]), r_dis_g_dis_ratio, "-", ",".join([cigar, cigar2])])
							breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), r_dis_g_dis_ratio, "-"])
#							breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_end)]), ",".join([str(r_start2), str(r_end2)]), ",".join([str(g_start), str(g_end)]), ",".join([str(g_start2), str(g_end2)]), ",".join([str(bq), str(bq2)]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper]), r_dis_g_dis_ratio, "-"])
							mq_bq_found = 1
							break
					if mq_bq_found == 1:
						break
		else:
#			print(read_tmp_l[10])
			read_start, read_end = read_tmp_l[10].split(",")
			total_ins_len = find_close_ins(read_tmp_l[14], int(read_start), int(read_tmp_l[1]), int(read_tmp_l[3]) - int(read_tmp_l[1]), close_ins_distance)
#			breakpoint_info.append([read_tmp_l[6], "within", chr, read_tmp_l[1], chr, read_tmp_l[3], r_start, r_end, g_start, g_end, bq, mq, cigar, strand, mapper, "-", total_ins_len, cigar])
#			breakpoint_info.append([read_tmp_l[6], "within", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], r_start, r_end, g_start, g_end, bq, mq, strand, mapper, "-", total_ins_len])
			r_start, r_end = re.split(",", read_tmp_l[7])
			g_start, g_end = re.split(",", read_tmp_l[10])
			bq = read_tmp_l[12]
			mq = read_tmp_l[13]
			strand = read_tmp_l[11]
			mapper = read_tmp_l[15]

			breakpoint_info.append([read_tmp_l[6], "within", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], r_start, r_end, g_start, g_end, bq, mq, strand, mapper, "-", total_ins_len])

		if mq_bq_found == 0:
			strange_read.append(j)

	breakpoint_start_l = []
	breakpoint_end_l = []
	for breakpoint_info_tmp in breakpoint_info:
		breakpoint_start_l.append(int(breakpoint_info_tmp[3]))
		breakpoint_end_l.append(int(breakpoint_info_tmp[5]))

	start_median = int(statistics.median(breakpoint_start_l))
	end_median = int(statistics.median(breakpoint_end_l))

	best_start, best_end, difference_from_median = 0, 0, 0
	for breakpoint_info_tmp in breakpoint_info:
		difference_from_median_tmp = abs(start_median - int(breakpoint_info_tmp[3])) + abs(end_median - int(breakpoint_info_tmp[5]))
		if breakpoint_info_tmp[1] == "within":
			if difference_from_median == 0:
				best_start, best_end, difference_from_median = int(breakpoint_info_tmp[3]), int(breakpoint_info_tmp[5]), difference_from_median_tmp
			elif difference_from_median_tmp <= difference_from_median:
				best_start, best_end, difference_from_median = int(breakpoint_info_tmp[3]), int(breakpoint_info_tmp[5]), difference_from_median_tmp

	difference_from_median = 0
	if difference_from_median == 0:
		for breakpoint_info_tmp in breakpoint_info:
			difference_from_median_tmp = abs(start_median - int(breakpoint_info_tmp[3])) + abs(end_median - int(breakpoint_info_tmp[5]))
			if difference_from_median == 0: 
				best_start, best_end, difference_from_median = int(breakpoint_info_tmp[3]), int(breakpoint_info_tmp[5]), difference_from_median_tmp
			elif difference_from_median_tmp < difference_from_median:
				best_start, best_end, difference_from_median = int(breakpoint_info_tmp[3]), int(breakpoint_info_tmp[5]), difference_from_median_tmp

	line_l[1] = str(best_start)
	line_l[3] = str(best_end)

	breakpoint_info_out = []
	for breakpoint_info_tmp in breakpoint_info:
		breakpoint_info_out.append(";".join(map(str, breakpoint_info_tmp)))

	print("\t".join(line_l[0:read_col_num]), "|".join(breakpoint_info_out), "\t".join(line_l[read_col_num:len(line_l)]), sep="\t") 

