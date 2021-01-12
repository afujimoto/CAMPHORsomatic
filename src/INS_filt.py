import re
import sys
from operator import itemgetter

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
			if float(sr_tmp_l[7]) > float(max_prop):
				max_prop = sr_tmp_l[7]
				max_rep = len(sr_tmp_l[6])
				max_rep_length = int(sr_tmp_l[2]) - int(sr_tmp_l[1])

	covered_length = 0
	if max_prop == 1:
		max_rep_length = line_l[6]
	else:
		if line_l[rep_col] != "-":
			boundary = []
			rep_start = []
			rep_end = []
			for sr_tmp in SR_l:
				sr_tmp_l = re.split(",", sr_tmp)
				if int(sr_tmp_l[1]) < int(line_l[1]):
					rep_start.append(int(line_l[1]))
				else:
					rep_start.append(int(sr_tmp_l[1]))
				if int(line_l[3]) < int(sr_tmp_l[2]):
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

	rep_prop = 0
	if float(line_l[6]) > 0:
		rep_prop = round(float(covered_length)/float(line_l[6]), 2)

	return rep_prop, max_rep, max_rep_length

f = open(sys.argv[1])
ins_range = int(sys.argv[2])

for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)

	read_l = re.split(':', line_l[8])

#	print("read_l", read_l)

	sc_max_prop = 0
	max_rep = "NA"
	max_rep_length = 0
	number_of_non_sc_read = 0

	output_read_BP = {}
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
				r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand, mapper = seq_align_pos[i]
				if (line_l[5] == "DEL" or line_l[5] == "INS") and chr == read_tmp_l[0] and chr == read_tmp_l[2] and int(g_start) < int(read_tmp_l[1]) < int(g_end) and int(g_start) < int(read_tmp_l[3]) < int(g_end):
					mq_bq_found = 1
					breakpoint_info.append([read_tmp_l[6], "within", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], r_start, r_end, g_start, g_end, bq, mq, strand, mapper, "-", "-"])
					break

				if i < len(seq_align_pos) - 1:
					r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar2, strand2, mapper2 = seq_align_pos[i + 1]
#					if chr == line_l[0] and chr2 == line_l[2] and abs(int(g_end) - int(line_l[1])) < ins_range and abs(int(g_start2) - int(line_l[3])) < ins_range:
					if chr == line_l[0] and chr2 == line_l[2] and abs(int(g_end) - int(line_l[1])) < ins_range and abs(int(g_start2) - int(line_l[3])) < ins_range and strand == strand2:
						g_start_out, g_end_out = g_end, g_start2
						if g_start_out > g_end_out:
							g_start_out, g_end_out = g_start2, g_end
						BP_name = "/".join([read_tmp_l[6], chr, str(g_start_out), chr, str(g_end_out)])
						if not BP_name in output_read_BP:
							breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr2, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), "NA", "-"])
#							print("2", [read_tmp_l[6], "split", chr, read_tmp_l[1], chr2, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), "NA", "-"])
							output_read_BP[BP_name] = 1
						mq_bq_found = 1	
						break
#					elif chr2 == line_l[0] and chr == line_l[2] and abs(int(g_start2) - int(line_l[1])) < ins_range and abs(int(g_end) - int(line_l[3])) < ins_range:
					elif chr2 == line_l[0] and chr == line_l[2] and abs(int(g_start2) - int(line_l[1])) < ins_range and abs(int(g_end) - int(line_l[3])) < ins_range and strand == strand2:
						g_start_out, g_end_out = g_end, g_start2
						if g_start_out > g_end_out:
							g_start_out, g_end_out = g_start2, g_end
						BP_name = "/".join([read_tmp_l[6], chr, str(g_start_out), chr, str(g_end_out)])
						if not BP_name in output_read_BP:
							breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr2, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), "NA", "-"])
#							print("2", [read_tmp_l[6], "split", chr, read_tmp_l[1], chr2, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), "NA", "-"])
							output_read_BP[BP_name] = 1
						mq_bq_found = 1
						break

			if mq_bq_found == 0:
				for i in range(len(seq_align_pos)):
					r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand, mapper = seq_align_pos[i]
					for k in range(len(seq_align_pos)):
						if i == k:
							continue
						r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar2, strand2, mapper2 = seq_align_pos[k]
#						if chr == line_l[0] and chr2 == line_l[2] and abs(int(g_end) - int(line_l[1])) < ins_range and abs(int(g_start2) - int(line_l[3])) < ins_range:
						if chr == line_l[0] and chr2 == line_l[2] and abs(int(g_end) - int(line_l[1])) < ins_range and abs(int(g_start2) - int(line_l[3])) < ins_range and strand == strand2:
							g_start_out, g_end_out = g_end, g_start2
							if g_start_out > g_end_out:
								g_start_out, g_end_out = g_start2, g_end
							BP_name = "/".join([read_tmp_l[6], chr, str(g_start_out), chr, str(g_end_out)])
							if not BP_name in output_read_BP:
								breakpoint_info.append([read_tmp_l[6], "split", chr, g_start_out, chr, g_end_out, read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), "NA", "-"])
#								print([read_tmp_l[6], "split", chr, g_start_out, chr, g_end_out, read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), "NA", "-"])
								output_read_BP[BP_name] = 1
							mq_bq_found = 1
							if abs(i - k) >= 2:
								strange_read.append(j)
							break
#						elif chr2 == line_l[0] and chr == line_l[2] and abs(int(g_start2) - int(line_l[1])) < ins_range and abs(int(g_end) - int(line_l[3])) < ins_range:
						elif chr2 == line_l[0] and chr == line_l[2] and abs(int(g_start2) - int(line_l[1])) < ins_range and abs(int(g_end) - int(line_l[3])) < ins_range and strand == strand2:
							g_start_out, g_end_out = g_start2, g_end
							if g_start_out > g_end_out:
								g_start_out, g_end_out = g_start2, g_end
							BP_name = "/".join([read_tmp_l[6], chr, str(g_start_out), chr, str(g_end_out)])
							if not BP_name in output_read_BP:
								breakpoint_info.append([read_tmp_l[6], "split", chr, g_start_out, chr, g_end_out, read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), "NA", "-"])
#								print([read_tmp_l[6], "split", chr, g_start_out, chr, g_end_out, read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq,bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper2]), "NA", "-"])
								output_read_BP[BP_name] = 1
							mq_bq_found = 1
							if abs(i - k) >= 2:
								strange_read.append(j)
							break
						if mq_bq_found == 1:
							break
			if mq_bq_found == 0:
				breakpoint_chr, breakpoint = read_tmp_l[0], int(read_tmp_l[1])
				for i in range(len(seq_align_pos)):
					r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand, mapper = seq_align_pos[i]

					if line_l[0] == read_tmp_l[0] and int(read_tmp_l[1]) - int(line_l[1]) <= ins_range:
						breakpoint_chr, breakpoint = read_tmp_l[0], int(read_tmp_l[1])

					if chr == breakpoint_chr and (int(g_start) == breakpoint or int(g_end) == breakpoint):
						type = "split"
						if int(g_start) == breakpoint:
							type = "lsplit"
						if int(g_end) == breakpoint:
							type = "rsplit"
						breakpoint_info.append([read_tmp_l[6], type, breakpoint_chr, breakpoint, "NA", "NA", read_tmp_l[9], ",".join([str(r_start), "NA"]), ",".join([str(r_end), "NA"]), ",".join([str(g_start), "NA"]), ",".join([str(g_end), "NA"]), ",".join([bq, "NA"]), ",".join([mq, "NA"]), ",".join([strand, "NA"]), ",".join([mapper, "NA"]), "-", "-"])
		else:
			chr1, chr2 = read_tmp_l[8], read_tmp_l[8]
			read_start, read_end = read_tmp_l[10].split(",")
			r_start, r_end = re.split(",", read_tmp_l[7])
			g_start, g_end = re.split(",", read_tmp_l[10])
			bq = read_tmp_l[12]
			mq = read_tmp_l[13]
			strand = read_tmp_l[11]
			mapper = read_tmp_l[15]

#			if chr1 == line_l[0] and chr2 == line_l[2] and int(g_start) < int(line_l[1]) < int(g_end):
			if chr1 == line_l[0] and chr2 == line_l[2] and int(g_start) < int(read_tmp_l[1]) < int(g_end):
				mq_bq_found = 1
				breakpoint_info.append([read_tmp_l[6], "within", chr1, read_tmp_l[1], chr2, read_tmp_l[3], read_tmp_l[9], r_start, r_end, g_start, g_end, bq, mq, strand, mapper, "-", "-"])
			if mq_bq_found == 0:
#				print(">>>>>>>", read_tmp_l)
#				print("read_tmp_l[0], int(read_tmp_l[1])", read_tmp_l[0], int(read_tmp_l[1]))
#				print("line_l[0]", line_l[0])
				breakpoint_chr, breakpoint  = "", ""
				if line_l[0] == read_tmp_l[0] and int(read_tmp_l[1]) - int(line_l[1]) <= ins_range:
					breakpoint_chr, breakpoint = read_tmp_l[0], int(read_tmp_l[1])

				if chr1 == breakpoint_chr and (int(g_start) == breakpoint or int(g_end) == breakpoint):
					type = "split"
					if int(g_start) == breakpoint:
						type = "lsplit"
					if int(g_end) == breakpoint:
						type = "rsplit"
					breakpoint_info.append([read_tmp_l[6], type, chr1, str(read_tmp_l[1]), "NA", "NA", read_tmp_l[9], ",".join([str(r_start), "NA"]), ",".join([str(r_end), "NA"]), ",".join([str(g_start), "NA"]), ",".join([str(g_end), "NA"]), ",".join([bq, "NA"]), ",".join([mq, "NA"]), ",".join([strand, "NA"]), ",".join([mapper, "NA"]), "-", "-"])
			strange_read.append(j)

	breakpoint_info_out = []
	for breakpoint_info_tmp in breakpoint_info:
#		print("breakpoint_info_tmp", breakpoint_info_tmp)
		breakpoint_info_out.append(";".join(map(str, breakpoint_info_tmp)))

	if len(breakpoint_info_out):
		print("\t".join(line_l[0:8]), "|".join(breakpoint_info_out), "\t".join(line_l[8:len(line_l)]), sep="\t") 
