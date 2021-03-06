import sys
import subprocess
import re
import pysam
#chr1    47280   chr1    49577   1       DEL     2297
#chr1    50682   chr1    51052   1       DEL     370

def get_depth(chr, pos, bam):
#	cmd = "samtools view " + bam + " " + chr + ":" + pos + "-" + pos
#	contents = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
#	contents = str(contents)
#	contents = contents.replace("b\'", "")
#	contents_l = re.split("\\\\n", contents)
#	num = len(contents_l)
#	return num

	bamfile = pysam.AlignmentFile(bam, "rb")
	depth = 0
	for read in bamfile.fetch(chr, int(pos), int(pos) + 1):
		depth += 1

	return depth

normal_f = open(sys.argv[1])
normal_SV_list = {}
for line in normal_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line[0] == "#":
		continue

	if line_l[0] == line_l[2]:
		normal_SV_list.setdefault(line_l[0], []).append(line_l[0:7])
	else:
		normal_SV_list.setdefault(line_l[0], []).append(line_l[0:7])
		normal_SV_list.setdefault(line_l[2], []).append(line_l[0:7])

cancer_f = open(sys.argv[2])
max_length_for_prop_filt = int(sys.argv[3])
overlap_prop = float(sys.argv[4])
max_distance_btw_BPs = int(sys.argv[5])
slice_col = int(sys.argv[6])

for line in cancer_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if max_length_for_prop_filt and int(line_l[3]) - int(line_l[1]) >= max_length_for_prop_filt:
		print("\t".join(line_l[0:slice_col]), "-", "\t".join(line_l[slice_col:len(line_l)]), sep="\t")
		continue
#del_OL_length = min(int(read_l[i][3]), int(read_l[j][3])) - max(int(read_l[i][1]), int(read_l[j][1])) + 1       
#if del_OL_length > 0 and float(del_OL_length)/float(read_l[i][6]) >= overlap_prop and float(del_OL_length)/float(read_l[j][6]) >= overlap_prop:
	
#	print(">>>>>", line_l[0:5])

	normal_SV_out = []
	if line_l[0] == line_l[2]:
		if line_l[0] in normal_SV_list:
			for normal_SV in normal_SV_list[line_l[0]]:
#				print("#############", normal_SV[0:6])
				if abs(int(normal_SV[1]) - int(line_l[1])) > max_distance_btw_BPs and abs(int(normal_SV[3]) - int(line_l[3])) > max_distance_btw_BPs:
					continue
				del_OL_length = min(int(normal_SV[3]), int(line_l[3])) - max(int(normal_SV[1]), int(line_l[1])) + 1
				if del_OL_length > 0 and float(del_OL_length)/float(line_l[6]) >= overlap_prop and float(del_OL_length)/float(normal_SV[6]) >= overlap_prop:
					SV_tmp = ",".join((normal_SV))
					normal_SV_out.append(SV_tmp)
#					if int(line_l[1]) + int(line_l[6])*overlap_prop < int(normal_SV[1]):
#						break
	else:
		if line_l[0] in normal_SV_list:
			for normal_SV in normal_SV_list[line_l[0]]:
				if normal_SV[0] == line_l[0] and normal_SV[2] == line_l[2] and abs(int(normal_SV[1]) - int(line_l[1])) <= distance and abs(int(normal_SV[3]) - int(line_l[3])) <= distance:
					
					SV_tmp = ",".join((normal_SV))
					normal_SV_out.append(SV_tmp)
#					if int(line_l[1]) + distance < int(normal_SV[1]):
#						break

		if line_l[2] in normal_SV_list:
			for normal_SV in normal_SV_list[line_l[2]]:
				if normal_SV[2] == line_l[0] and normal_SV[0] == line_l[2] and abs(int(normal_SV[3]) - int(line_l[1])) <= distance and abs(int(normal_SV[1]) - int(line_l[3])) <= distace:
					SV_tmp = ",".join((normal_SV))
					normal_SV_out.append(SV_tmp)
#					if int(line_l[3]) + distance < int(normal_SV[1]):
#						break

#	depth1 = get_depth(line_l[0], str(int(line_l[1]) - 100), bam)
#	depth2 = get_depth(line_l[2], str(int(line_l[3]) + 100), bam)

	if len(normal_SV_out) == 0:
		normal_SV_out = ["-"]
		print("\t".join(line_l[0:slice_col]), ";".join(normal_SV_out), "\t".join(line_l[slice_col:len(line_l)]), sep="\t")

#	print("\t".join(line_l[0:slice_col]), ";".join(normal_SV_out), "\t".join(line_l[slice_col:len(line_l)]), sep="\t")
