import sys
import subprocess
import re
import pysam

#chr10   84226   chr10   84406   2       DEL     180     NA      chr10;84180;chr10;84432;CTGAGACCTGAGCACACTCAC;252;3a13e1f5-aa6
#chr10   264615  chr10   264728  2       DEL     128     NA      chr10;264600;chr10;264728;AACTGTAGAATTTTAGAAGAT;128;d306c1ea-1
#chr10   264911  chr10   265000  17      DEL     122     NA      chr10;264879;chr10;265001;AGATTTTGACTGCTGTAGGGT;122;2b6ad657-d


def get_depth(chr, pos, bam):
	bamfile = pysam.AlignmentFile(bam, "rb")
	depth = 0
	for read in bamfile.fetch(chr, int(pos), int(pos) + 1):
		depth += 1
	return depth

f = open(sys.argv[1])
distance = int(sys.argv[2])
bam = sys.argv[3]
slice_col = int(sys.argv[4])

SV_list_N, SV_list_C = [], []
pre_chr, pre_pos = 0, ""
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line_l[1] != pre_chr or pre_pos > int(line_l[2]) + distance:
		for tmp_C in SV_list_C:
			normal_SV = []
			for tmp_B in SV_list_N:
				if abs(int(tmp_C[1]) - int(tmp_B[1])) <= distance and abs(int(tmp_C[3]) - int(tmp_B[3])) <= distance:
					normal_SV.append(",".join(tmp_B[0:7]))
		
			if len(normal_SV) == 0:
				normal_SV = ["-"]

			depth1 = 100
			depth2 = 100
			if "NC_" in line_l[0]:
				depth1 = 1
			else:
				depth1 = get_depth(tmp_C[0], str(int(tmp_C[1]) - 100), bam)
			if "NC_" in line_l[2]:
				depth2 = 1
			else:
				depth2 = get_depth(tmp_C[2], str(int(tmp_C[3]) + 100), bam)

			print("\t".join(tmp_C[0:slice_col]), str(depth1) + "," + str(depth2), ";".join(normal_SV), "\t".join(tmp_C[slice_col:]), sep="\t")	

		SV_list_N, SV_list_C = [], []

	if line_l[0] == "N":
		SV_list_N.append(line_l[1:])

	if line_l[0] == "C":
		SV_list_C.append(line_l[1:])
	pre_chr, pre_pos = line_l[1], int(line_l[2])

"""
	cancer_SV_out = []
	if line_l[0] == line_l[2]:
		if line_l[0] in cancer_SV_list:
			for cancer_SV in cancer_SV_list[line_l[0]]:
				if abs(int(cancer_SV[1]) - int(line_l[1])) <= distance and abs(int(cancer_SV[3]) - int(line_l[3])) <= distance:
#					print("#################")
					SV_tmp = ",".join((line_l[0:7]))
					cancer_SV_comparison["\t".join(cancer_SV[0:6])].append(SV_tmp)
#					cancer_SV_out.append(SV_tmp)
					if int(line_l[1]) + distance < int(cancer_SV[1]):
						break

	if line_l[0] == "N":
		SV_list_N.append(line_l[1:])

	if line_l[0] == "C":
		SV_list_C.append(line_l[1:])
	pre_chr, pre_pos = line_l[1], int(line_l[2])

cancer_f = open(sys.argv[1])
for line in cancer_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	
	depth1 = 100
	depth2 = 100
	if "NC_" in line_l[0]:
		depth1 = 1
	else:
		depth1 = get_depth(line_l[0], str(int(line_l[1]) - 100), bam)
	if "NC_" in line_l[2]:
		depth2 = 1
	else:
		depth2 = get_depth(line_l[2], str(int(line_l[3]) + 100), bam)

#	print(line)
#	print(cancer_SV_comparison["\t".join(line_l[0:6])])

	repeat = []
	if len(cancer_SV_comparison["\t".join(line_l[0:6])]) == 0:
		repeat = ["-"]
	else:
		repeat = cancer_SV_comparison["\t".join(line_l[0:6])]

#normal_SV_out.append(SV_tmp)
#	print("\t".join(line_l[0:slice_col]), str(depth1) + "," + str(depth2), ";".join(normal_SV_out), "\t".join(line_l[slice_col:len(line_l)]), sep="\t")
	print("\t".join(line_l[0:slice_col]), str(depth1) + "," + str(depth2), ";".join(repeat), "\t".join(line_l[slice_col:len(line_l)]), sep="\t")
"""
