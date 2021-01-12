import sys
import subprocess
import re
import pysam

#chr10   84226   chr10   84406   2       DEL     180     NA      chr10;84180;chr10;84432;CTGAGACCTGAGCACACTCAC;252;3a13e1f5-aa6
#chr10   264615  chr10   264728  2       DEL     128     NA      chr10;264600;chr10;264728;AACTGTAGAATTTTAGAAGAT;128;d306c1ea-1
#chr10   264911  chr10   265000  17      DEL     122     NA      chr10;264879;chr10;265001;AGATTTTGACTGCTGTAGGGT;122;2b6ad657-d


def get_depth(chr, pos, bam):
#	cmd = "/share/amed_snt/WORK/fujimoto/src/tools/samtools-0.1.6/samtools view " + bam + " " + chr + ":" + pos + "-" + pos
##	print(cmd)
#	contents = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
#	contents = str(contents)
#	contents = contents.replace("b\'", "")
#	contents_l = re.split("\\\\n", contents)
#	num = len(contents_l)
#	return num
	bamfile = pysam.AlignmentFile(bam, "rb")
	depth = 0

	if int(pos) <= 0:
		pos = 0

	#print(chr, int(pos))
	for read in bamfile.fetch(chr, int(pos), int(pos) + 1):
	#	print(read)
		depth += 1

	#print("depth", depth)

	return depth

cancer_f = open(sys.argv[2])
cancer_SV_list = {}
cancer_SV_comparison = {}
for line in cancer_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line[0] == "#":
		continue

	if line_l[0] == line_l[2]:
		cancer_SV_list.setdefault(line_l[0], []).append(line_l)
	else:
		cancer_SV_list.setdefault(line_l[0], []).append(line_l)
		cancer_SV_list.setdefault(line_l[2], []).append(line_l)

	cancer_SV_comparison["\t".join(line_l[0:6])] = []

normal_f = open(sys.argv[1])
distance = int(sys.argv[3])
bam = sys.argv[4]
slice_col = int(sys.argv[5])

for line in normal_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line[0] == "#":
		continue

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
"""
	else:
		if line_l[0] in cancer_SV_list:
			for cancer_SV in cancer_SV_list[line_l[0]]:
				if cancer_SV[0] == line_l[0] and cancer_SV[2] == line_l[2] and abs(int(cancer_SV[1]) - int(line_l[1])) <= distance and abs(int(cancer_SV[3]) - int(line_l[3])) <= distance:
					
					SV_tmp = ",".join((cancer_SV))
					cancer_SV_out.append(SV_tmp)
#					if int(line_l[1]) + distance < int(normal_SV[1]):
#						break

		if line_l[2] in cancer_SV_list:
			for cancer_SV in cancer_SV_list[line_l[2]]:
				if cancer_SV[2] == line_l[0] and cancer_SV[0] == line_l[2] and abs(int(cancer_SV[3]) - int(line_l[1])) <= distance and abs(int(cancer_SV[1]) - int(line_l[3])) <= distace:
					SV_tmp = ",".join((cancer_SV))
					cancer_SV_out.append(SV_tmp)
#					if int(line_l[3]) + distance < int(normal_SV[1]):
#						break
"""
cancer_f = open(sys.argv[2])
for line in cancer_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line[0] == "#":
		continue
	
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
