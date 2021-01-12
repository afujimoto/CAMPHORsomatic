import re
import sys
from operator import itemgetter
import subprocess

#chr10   37451   chr10   37564   2       DEL     113     0       chr10;37409;chr10;37564;TTGAGAGCCTGTCCATTGATCTATATCCTGTTTTAAACAGTACCAGTTGTTTTAATATCACAGCCTTGCATAAATAGTCAGGCAGTATTCCAGCTGTTCTGGTTTAGACTGACTTATACGCGGCCTTTTTGTCCATATAATTTCAGCGTTTTCCACCCTGCATAACCATTATAGTTTGACGAGACAAATTGACCTAAATTACCTTGGCAAAATATTTTCAATATTGATTCTTCCTAACGCATGAGCATAATGTTCTTCATTTGTTTGTATCCTATTCTGGAGCAGGTTTGTCTCCTGAAGAGTATCCCACATCCCTGTTATCTTATGTATTTACCTTTGCGCAATCAGAATGGGAGTTCATCTGACGGCCTGCTTCATCAGTGGCGGCAGCATCGTTTACT;155;653b180d-754b-4b4a-89b7-af1b4885fd65_Basecall_1D_template;3445,5339/340,1949/11,339/1974,2900;chr10/chr10/chr10/chr10;5393;109799256,109801664/35682,37761/109795057,109795430/109797362,109798491;-/+/+/+;10.3/10.3/10.3/10.3;11/0/60/5;3444S13M4D12M10D14M1D13M1D36M2D17M1D21M1D9M4D10M1I3M2D14M2D2M2D4M154D2M1D8M2D3M1D1M1D2M1D4M2I3M3D20M2D5M1D3M2I2M5D3M1D2M4D7M2D5M1D3M3D5M1D7M5D3M1D9M1D7M3D6M2D7M1D5M5D5M3D2M3D14M3D8M4D11M1D4M1D13M1I5M4D6M1D12M8D5M1D24M5D7M3D6M1D10M1I1M2I1M3D1M4D6M1D8M2D1M1D9M10D1M1D6M4I7M2D1M1D3M3D2M1D19M4D1M1D8M4D23M3D4M3D7M9D6M2D22M1D5M1D11M4D5M2I15M2D28M2D9M1D3M3D7M4D7M3D4M3D2M2D6M1D26M1D6M2D16M1D24M1D4M1D3M2I15M1D7M2I7M2D9M2D4M1D3M2D17M2D2M2D4M6D9M3D9M1I2M3D3M2D6M1I1M1I22M1D9M3D2M2D2M1D15M2D5M1D12M1D9M3D17M2D2M2D2M1D2M1D9M1D8M2D3M2D2M2D3M4D11M3D6M1D13M1D5M2D2M1I23M2D7M1D4M2D1M1D3M4D9M1D3M4D12M1D2M1D19M1D22M5D2M2D7M1D8M3I8M1D6M1D3M1D6M1D6M1D31M1I6M1D5M1I9M1D8M2I1M1I3M2D6M2D38M1D1M2D3M2D19M1D11M2I28M1D9M2D6M2I7M1D9M1D17

def get_depth(chr, pos, bam, min_mq):
	cmd = "samtools view " + bam + " " + chr + ":" + pos + "-" + pos
	contents = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
	contents = str(contents)
	contents = contents.replace("b\'", "")
	contents_l = re.split("\\\\n", contents)
	total_read_num = 0
	low_mq_read_num = 0
	for read in contents_l:
		if read == "\'":
			continue
		read_l = read.split("\\t")
		if int(read_l[4]) < min_mq:
			low_mq_read_num += 1
		total_read_num += 1
	return total_read_num, low_mq_read_num

f = open(sys.argv[1])
minimum_length_prop = float(sys.argv[2])
min_mq = int(sys.argv[3])
min_segment_length = int(sys.argv[4])
min_num_read = int(sys.argv[5])

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)
	
	print(line)
	if line[0] == "#":
		continue

	read_l = re.split(':', line_l[8])

	sc_max_prop = 0
	max_rep = "NA"
	max_rep_length = 0
	number_of_non_sc_read = 0

#chr10   37451   chr10   37564   2       DEL     113     0       chr10;37409;chr10;37564;TTGAGAGCCTGTCCATTGATCTATATCCTGTTTTAAACAGTACCAGTTGTTTTAATATCACAGCCTTGCATAAATAGTCAGGCAGTATTCCAGCTGTTCTGGTTTAGACTGACTTATACGCGGCCTTTTTGTCCATATAATTTCAGCGTTTTCCACCCTGCATAACCATTATAGTTTGACGAGACAAATTGACCTAAATTACCTTGGCAAAATATTTTCAATATTGATTCTTCCTAACGCATGAGCATAATGTTCTTCATTTGTTTGTATCCTATTCTGGAGCAGGTTTGTCTCCTGAAGAGTATCCCACATCCCTGTTATCTTATGTATTTACCTTTGCGCAATCAGAATGGGAGTTCATCTGACGGCCTGCTTCATCAGTGGCGGCAGCATCGTTTACT;155;653b180d-754b-4b4a-89b7-af1b4885fd65_Basecall_1D_template;3445,5339/340,1949/11,339/1974,2900;chr10/chr10/chr10/chr10;5393;109799256,109801664/35682,37761/109795057,109795430/109797362,109798491;-/+/+/+;10.3/10.3/10.3/10.3;11/0/60/5;3444S13M4D12M10D14M1D13M1D36M2D17M1D21M1D9M4D10M1I3M2D14M2D2M2D4M154D2M1D8M2D3M1D1M1D2M1D4M2I3M3D20M2D5M1D3M2I2M5D3M1D2M4D7M2D5M1D3M3D5M1D7M5D3M1D9M1D7M3D6M2D7M1D5M5D5M3D2M3D14M3D8M4D11M1D4M1D13M1I5M4D6M1D12M8D5M1D24M5D7M3D6M1D10M1I1M2I1M3D1M4D6M1D8M2D1M1D9M10D1M1D6M4I7M2D1M1D3M3D2M1D19M4D1M1D8M4D23M3D4M3D7M9D6M2D22M1D5M1D11M4D5M2I15M2D28M2D9M1D3M3D7M4D7M3D4M3D2M2D6M1D26M1D6M2D16M1D24M1D4M1D3M2I15M1D7M2I7M2D9M2D4M1D3M2D17M2D2M2D4M6D9M3D9M1I2M3D3M2D6M1I1M1I22M1D9M3D2M2D2M1D15M2D5M1D12M1D9M3D17M2D2M2D2M1D2M1D9M1D8M2D3M2D2M2D3M4D11M3D6M1D13M1D5M2D2M1I23M2D7M1D4M2D1M1D3M4D9M1D3M4D12M1D2M1D19M1D22M5D2M2D7M1D8M3I8M1D6M1D3M1D6M1D6M1D31M1I6M1D5M1I9M1D8M2I1M1I3M2D6M2D38M1D1M2D3M2D19M1D11M2I28M1D9M2D6M2I7M1D9M1D17
	all_bq_l = []
	all_mq_l = []
	strange_read = []
	number_of_chr = []
	all_strand = []
	breakpoint_info = []
	read_length_l = []
	for j in range(len(read_l)):
		read_tmp_l = re.split(";", read_l[j])
		mq_bq_found = 0
		if "/" in read_tmp_l[8]:
			read_pos = re.split("/", read_tmp_l[7])
			chr_l = re.split("/", read_tmp_l[8])
			genome_pos = re.split("/", read_tmp_l[10])
			bq_l = re.split("/", read_tmp_l[12])
			mq_l = re.split("/", read_tmp_l[13])
			number_of_chr.append(len(set(chr_l)))
			cigar_l = re.split("/", read_tmp_l[14])
			strand_l = re.split("/", read_tmp_l[11])

			seq_align_pos = []
			for i in range(len(read_pos)):
				read_pos_tmp_l = re.split(",", read_pos[i])
				genome_pos_tmp_l = re.split(",", genome_pos[i])
				if int(read_pos_tmp_l[1]) - int(read_pos_tmp_l[0]) < minimum_length_prop*int(read_tmp_l[9]):
					continue
				seq_align_pos.append([int(read_pos_tmp_l[0]), int(read_pos_tmp_l[1]), int(genome_pos_tmp_l[0]), int(genome_pos_tmp_l[1]), chr_l[i], bq_l[i], mq_l[i], cigar_l[i], strand_l[i]])

			seq_align_pos.sort(key=itemgetter(0))

			for i in range(len(seq_align_pos)):
				r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand = seq_align_pos[i]
				if strand == "-":
					g_start, g_end = g_end, g_start
				if (line_l[5] == "DEL" or line_l[5] == "INS") and chr == read_tmp_l[0] and chr == read_tmp_l[2] and int(g_start) <= int(read_tmp_l[1]) <= int(g_end) and int(g_start) <= int(read_tmp_l[3]) <= int(g_end):
					all_bq_l.append(read_tmp_l[12])
					all_bq_l.append(bq)
					all_mq_l.append(mq)
					mq_bq_found = 1
					breakpoint_info.append([read_tmp_l[6], "within", read_tmp_l[1], read_tmp_l[3], read_tmp_l[9], mq])
					break

				if i < len(seq_align_pos) - 1:
					r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar, strand2 = seq_align_pos[i + 1]
					if strand2 == "-":
						g_start2, g_end2 = g_end2, g_start2
					if (chr == read_tmp_l[0] and chr2 == read_tmp_l[2] and int(g_end) == int(read_tmp_l[1]) and int(g_start2) == int(read_tmp_l[3])) or (chr2 == read_tmp_l[0] and chr == read_tmp_l[2] and int(g_start2) == int(read_tmp_l[1]) and int(g_end) == int(read_tmp_l[3])):
						r_len = float((int(r_end) - int(r_start) + int(r_end2) - int(r_start2)))/2
						r_dis = abs(int(r_end) - int(r_start2))
#						all_bq_l.append((float(bq)+ float(bq2))/2)
						all_bq_l.append("-")
						all_mq_l.append(mq + "," + mq2)
						mq_bq_found = 1
						segment_length1 = abs(int(g_end) - int(g_start))
						segment_length2 = abs(int(g_end2) - int(g_start2))

						breakpoint_info.append([read_tmp_l[6], "split", read_tmp_l[1], read_tmp_l[3], read_tmp_l[9], read_tmp_l[0], read_tmp_l[2], r_len, r_dis, segment_length1, segment_length2, mq, mq2])
						break

			if mq_bq_found == 0:
				for i in range(len(seq_align_pos)):
					r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand = seq_align_pos[i]
					if strand == "-":
						g_start, g_end = g_end, g_start
					for k in range(len(seq_align_pos)):
						if i == k:
							continue
						r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar, strand2 = seq_align_pos[k]
						if strand2 == "-":
							g_start2, g_end2 = g_end2, g_start2
						if (chr == read_tmp_l[0] and chr2 == read_tmp_l[2] and int(g_end) == int(read_tmp_l[1]) and int(g_start2) == int(read_tmp_l[3])) or (chr2 == read_tmp_l[0] and chr == read_tmp_l[2] and int(g_start2) == int(read_tmp_l[1]) and int(g_end) == int(read_tmp_l[3])):
							r_len = float((int(r_end) - int(r_start) + int(r_end2) - int(r_start2)))/2
							r_dis = abs(int(r_end) - int(r_start2))
							segment_length1 = abs(int(g_end) - int(g_start))
							segment_length2 = abs(int(g_end2) - int(g_start2))

							breakpoint_info.append([read_tmp_l[6], "split", read_tmp_l[1], read_tmp_l[3], read_tmp_l[9], read_tmp_l[0], read_tmp_l[2] ,r_len, r_dis, segment_length1, segment_length2, mq, mq2])
#							all_bq_l.append((float(bq)+ float(bq2))/2)
							all_bq_l.append("-")
							all_mq_l.append(mq + "," + mq2)
							mq_bq_found = 1
							if abs(i - k) >= 2:
								strange_read.append(j)
							break
						if mq_bq_found == 1:
							break
		else:
			read_start, read_end = read_tmp_l[10].split(",")
			all_bq_l.append(read_tmp_l[12])
			all_mq_l.append(read_tmp_l[13])
			mq_bq_found = 1
			number_of_chr.append(1)
			breakpoint_info.append([read_tmp_l[6], "within", read_tmp_l[1], read_tmp_l[3], read_tmp_l[9], read_tmp_l[13]])

		if mq_bq_found == 0:
			all_bq_l.append("NA")
			all_mq_l.append("NA")
			strange_read.append(j)

#['1f59ac08-d4a0-4161-b0e7-8c72962a7e67', 'split', '124329', '127265', '8941', 'chr10', 'chr10', 4427.5, 27, 4143, 5101, '60', '60']
#['0a0f3af7-fb67-414c-b967-b9cebfc9accb', 'within', '35245', '35516', '3184', '23']
	num_read_with_long_segment = 0
	num_split_read = 0
	name_of_read_l = []
	for i in range(len(breakpoint_info)):
#		print(breakpoint_info[i])
		if breakpoint_info[i][1] == "within":
			if int(breakpoint_info[i][4]) >= min_segment_length and int(breakpoint_info[i][5]) >= min_mq:
				num_read_with_long_segment += 1

		else:
			if int(breakpoint_info[i][9]) >= min_segment_length and int(breakpoint_info[i][10]) >= min_segment_length and int(breakpoint_info[i][11]) >= min_mq and int(breakpoint_info[i][12]) >= min_mq:
				num_read_with_long_segment += 1

#	print(num_read_with_long_segment)
	if num_read_with_long_segment >= min_num_read:
		print(line)

#	print("\t".join(line_l[0:insert_col_num]), str(num_read_with_long_segment) + "/" + str(num_split_read), num_read, number_of_read_from_same_ch, "\t".join(line_l[insert_col_num:len(line_l)]), name_of_read, sep="\t")
