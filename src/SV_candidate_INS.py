import re
import sys
from operator import itemgetter

#chr4    54593956        chr4    54593921        TRS     35      TTCGTAACTACGCATGCGCTAAACTATATAAAGTGTTTTTTTCACCAAAATGCGGAGATTCGCAATTTAATTTGCCTAATAGACGGATAAAATCTAACGCCTGTAGGATTGCTCTTCTCTCTGGAGTGAGAGAAATTTCTGAACAAAAGAATCCATAAAATGGTGCAGGGGTAAAGTATGGTGTGATCTCTGTTCTTAAAGAATGCTCCTCTTCAGCACTCCCAAGGGCTGGCTATGCTCAAGAAAACATAACTAACAGAGTGTGCCCAGCAGAGATGTATGGAAAAAATCTCTTCATTAATTAGCAGCTTACATGTGGCA       00025171-ed0a-4db2-bead-907ec5e53024    7182    663,4831;4865,7063      chr4;chr4       54589109,54593921;54593956,54596243     +;+     8.4;8.4 60;18   662H10M4D12M2I8M2D10M3D7M3D5M1I25M1I18M1I12M1D4M1I9M3D8M2D5M1I12M1I2M2I9M3I7M1I7M1I2M2I2M2I15M1D12M2D21M1I2M1I11M1I3M2I9M1I4M1I2M1I6M1I2M1I7M2I4M1D9M7D2M1D4M1I9M1D7M2D1M1D4M1D5M3D1M3D4M2D2M2D13M1D2M2D23M3D10M1I5M1I2M1D7M1D6M1D8M1D20M8D2M1D4M2D1M4D1M1D5M1D7M10D1M1D13M1D3M5D6M2D6M2I7M1I3M1I4M1D2M1I4M1I18M1I5M2D10M1D9M1D7M3D26M1D10M1D8M8D11M2D20M2D9M1I6M2I1M1I5M2I3M1I4M1D19M2D7M1D10M5I10M1I4M1D1M3D9M3D7M1D3M4D4M1D5M1I5M1D3M1D15M2D7M2D12M1D2M3D23M3D11M5I6M2D2M1D21M5I10M4D17M1I9M1I1M1I15M3D4M2I13M1D10M1I13M1I5M1I4M1I2M2I7M2I3M1D3M2D15M1D13M2D16M4D2M1D2M3D6M5D8M2D13M3D9M2D7M1D4M1D3M3D5M6D4M2D8M1D11M4D3M1I12M1D4M2D7M1D4M1I15M1I12M1D5M2D31M3D1M5D4M1D3M1D4M1D15M1I8M1D22M4D6M1D6M2D4M1I4M1D11M1D1M1D31M1D1M2D9M3D11M1D25M1D6M1I23M1D3M2D9M1D8M4D5M1I20M3D4M2D11M1D2

def get_strand(line):
	read_tmp_l = re.split("\t", line)
	bp_chr1, bp_pos1, bp_chr2, bp_pos2 = read_tmp_l[0], read_tmp_l[1], read_tmp_l[2], read_tmp_l[3]
	dir1, dir2 = "", ""
	if ";" in read_tmp_l[9]:
		chr_l = re.split(";", read_tmp_l[10])
		bq_l = re.split(";", read_tmp_l[13])
		mq_l = re.split(";", read_tmp_l[14])
		cigar_l = re.split(";", read_tmp_l[15])
		strand_l = re.split(";", read_tmp_l[12])
		read_pos = re.split(";", read_tmp_l[9])
		genome_pos = re.split(";", read_tmp_l[11])

#		print("chr_l", chr_l)
#		print("bq_l", bq_l)
#		print("mq_l", mq_l)
#		print("cigar_l", cigar_l)
#		print("strand_l", strand_l)
#		print("read_pos", read_pos)
#		print("genome_pos", genome_pos)

		seq_align_pos = []
		for i in range(len(read_pos)):
			read_pos_tmp_l = re.split(",", read_pos[i])
			genome_pos_tmp_l = re.split(",", genome_pos[i])
#			print("read_pos_tmp_l => ", read_pos_tmp_l)
#			print("genome_pos_tmp_l => ", genome_pos_tmp_l)
#			print("read_tmp_l[0]", read_tmp_l[0])
#			print("read_pos_tmp_l[0]", read_pos_tmp_l[1])
#			print("genome_pos_tmp_l[0]", genome_pos_tmp_l[0])
#			print("genome_pos_tmp_l[1]", genome_pos_tmp_l[1])
			seq_align_pos.append([int(read_pos_tmp_l[0]), int(read_pos_tmp_l[1]), int(genome_pos_tmp_l[0]), int(genome_pos_tmp_l[1]), chr_l[i], bq_l[i], mq_l[i], cigar_l[i], strand_l[i]])

		seq_align_pos.sort(key=itemgetter(0))
		
		for i in range(len(seq_align_pos)):
			r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand = seq_align_pos[i]
			if strand == "-":
				g_start, g_end = g_end, g_start
			if (type == "DEL" or type == "INS") and chr == bp_chr1 and chr == bp_chr2 and int(g_start) <= int(bp_pos1) <= int(g_end) and int(g_start) <= int(bp_pos2) <= int(g_end):
				return ["within", strand_l[i]]

			if i < len(seq_align_pos) - 1:
				r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar, strand2 = seq_align_pos[i + 1]
				if strand2 == "-":
					g_start2, g_end2 = g_end2, g_start2
				if (chr == bp_chr1 and chr2 == bp_chr2 and int(g_end) == int(bp_pos1) and int(g_start2) == int(bp_pos2)) or (chr2 == bp_chr1 and chr == bp_chr2 and int(g_start2) == int(bp_pos1) and int(g_end) == int(bp_pos2)):
					if int(g_end) == int(bp_pos1) and strand == "+":
						dir1 = "r"
					elif int(g_end) == int(bp_pos1) and strand == "-":
						dir1 = "l"
					elif int(g_start2) == int(bp_pos1) and strand2 == "+":
						dir1 = "l"
					elif int(g_start2) == int(bp_pos1) and strand2 == "-":
						dir1 = "r"

					if int(g_end) == int(bp_pos2) and strand == "+":
						dir2 = "r"
					elif int(g_end) == int(bp_pos2) and strand == "-":
						dir2 = "l"
					elif int(g_start2) == int(bp_pos2) and strand2 == "+":
						dir2 = "l"
					elif int(g_start2) == int(bp_pos2) and strand2 == "-":
						dir2 = "r"
			
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
				if (chr == bp_chr1 and chr2 == bp_chr2 and int(g_end) == int(bp_pos1) and int(g_start2) == int(bp_pos2)) or (chr2 == bp_chr1 and chr == bp_chr2 and int(g_start2) == int(bp_pos1) and int(g_end) == int(bp_pos2)):
					if int(g_end) == int(bp_pos1) and strand == "+":
						dir1 = "r"
					elif int(g_end) == int(bp_pos1) and strand == "-":
						dir1 = "l"
					elif int(g_start2) == int(bp_pos1) and strand2 == "+":
						dir1 = "l"
					elif int(g_start2) == int(bp_pos1) and strand2 == "-":
						dir1 = "r"

					if int(g_end) == int(bp_pos2) and strand == "+":
						dir2 = "r"
					elif int(g_end) == int(bp_pos2) and strand == "-":
						dir2 = "l"
					elif int(g_start2) == int(bp_pos2) and strand2 == "+":
						dir2 = "l"
					elif int(g_start2) == int(bp_pos2) and strand2 == "-":
						dir2 = "r"

	return ["split", dir1, dir2]

#target_pos = re.split(":|-", sys.argv[2])

#min_deletion_length = int(sys.argv[3])
min_deletion_length = int(sys.argv[2])

#chr = target_pos[0]
#start = int(target_pos[1])
#end = int(target_pos[2])

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)

#	print(line)

	if line_l[4] == "DEL" and int(line_l[5]) <= min_deletion_length:
		continue

	if line_l[4] == "BR" and line_l[0] == "NA":
		line_l[4] = "lBR"
		line_l[0], line_l[1] = line_l[2], line_l[3]
		print("\t".join(line_l))
	elif line_l[4] == "BR" and line_l[2] == "NA":
		line_l[4] = "rBR"
		line_l[2], line_l[3] = line_l[0], line_l[1]
		print("\t".join(line_l))
	elif line_l[4] == "INS":
		print(line)
	elif line_l[0] == line_l[2]:
		strand_l = get_strand(line)
#		print("strand_l1", strand_l)
		line_l[4] = strand_l[1] + line_l[4]
		print("\t".join(line_l))
	elif line_l[0] != line_l[2]:
		strand_l = get_strand(line)
		print("\t".join(line_l))
#		print("strand_l2", strand_l)
		chr1, pos1 = line_l[0], line_l[1]
		chr2, pos2 = line_l[2], line_l[3]
		line_l[0], line_l[1] = chr2, pos2
		line_l[2], line_l[3] = chr1, pos1
		line_l[4] = strand_l[2] + line_l[4] #strand_l[1] => strand_l[2] 180707
		print("\t".join(line_l))

"""
	if line_l[4] == "BR" and line_l[0] == "NA" and line_l[2] == chr and start <= int(line_l[3]) <= end:
		line_l[4] = "lBR"
		line_l[0], line_l[1] = line_l[2], line_l[3]
		print("\t".join(line_l))
	elif line_l[4] == "BR" and line_l[2] == "NA" and line_l[0] == chr and start <= int(line_l[1]) <= end:
		line_l[4] = "rBR"
		line_l[2], line_l[3] = line_l[0], line_l[1]
		print("\t".join(line_l))
	elif line_l[4] == "INS" and line_l[0] == chr and start <= int(line_l[1]) <= end:
		print(line)
	elif line_l[0] == chr and start <= int(line_l[1]) <= end:
		strand_l = get_strand(line)
#		print("strand_l1", strand_l)
		line_l[4] = strand_l[1] + line_l[4]
		print("\t".join(line_l))
	elif line_l[2] == chr and start <= int(line_l[3]) <= end:
		strand_l = get_strand(line)
#		print("strand_l2", strand_l)
		chr1, pos1 = line_l[0], line_l[1]
		chr2, pos2 = line_l[2], line_l[3]
		line_l[0], line_l[1] = chr2, pos2
		line_l[2], line_l[3] = chr1, pos1
		line_l[4] = strand_l[2] + line_l[4] #strand_l[1] => strand_l[2] 180707
		print("\t".join(line_l))
"""
