import sys
import subprocess
import re

#chr1    47280   chr1    49577   1       DEL     2297
#chr1    50682   chr1    51052   1       DEL     370

def get_depth(chr, pos, bam):
	cmd = "/share/amed_snt/WORK/fujimoto/src/tools/samtools-0.1.6/samtools view " + bam + " " + chr + ":" + pos + "-" + pos
#	print(cmd)
	contents = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
	contents = str(contents)
	contents = contents.replace("b\'", "")
	contents_l = re.split("\\\\n", contents)
	num = len(contents_l)
	return num

normal_f = open(sys.argv[1])
normal_SV_list = {}
for line in normal_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line[0] == "#":
		continue

	if line_l[0] == line_l[2]:
		normal_SV_list.setdefault(line_l[0], []).append(line_l)
	else:
		normal_SV_list.setdefault(line_l[0], []).append(line_l)
		normal_SV_list.setdefault(line_l[2], []).append(line_l)

cancer_f = open(sys.argv[2])
distance = int(sys.argv[3])
range_ratio = float(sys.argv[4])
length_ratio = float(sys.argv[5])
length_filter_min = int(sys.argv[6])
bam = sys.argv[7]
slice_col = int(sys.argv[8])

for line in cancer_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line[0] == "#":
		continue

	INS_length = float(line_l[6].replace(">", ""))

	if INS_length*range_ratio < distance:
		distance2 = distance
	else:
		distance2 = INS_length*range_ratio

	normal_SV_out = []
	if line_l[0] == line_l[2]:
		if line_l[0] in normal_SV_list:
			for normal_SV in normal_SV_list[line_l[0]]:
				normal_INS_length = float(normal_SV[6].replace(">", ""))

				distance3 = 0
				if normal_INS_length*range_ratio < distance:
					distance3 = distance
				else:
					distance3 = normal_INS_length*range_ratio

				distance4 = distance2
				if distance2 < distance3:
					distance4 = distance3

				if abs(int(normal_SV[1]) - int(line_l[1])) <= distance4 and abs(int(normal_SV[3]) - int(line_l[1])) <= distance4:
					normal_INS_length = float(normal_SV[6].replace(">", ""))
#					print(INS_length, normal_INS_length, max(INS_length/normal_INS_length, normal_INS_length/INS_length))
					if max(INS_length/normal_INS_length, normal_INS_length/INS_length) <= length_ratio:
						SV_tmp = ",".join((normal_SV))
						normal_SV_out.append(SV_tmp)
					elif INS_length >= length_filter_min and normal_INS_length >= length_filter_min:
						SV_tmp = ",".join((normal_SV))
						normal_SV_out.append(SV_tmp)

					if int(line_l[1]) + distance4 < int(normal_SV[1]):
						break

				elif abs(int(normal_SV[1]) - int(line_l[3])) <= distance4 and abs(int(normal_SV[3]) - int(line_l[3])) <= distance4:
					normal_INS_length = float(normal_SV[6].replace(">", ""))
#					print(INS_length, normal_INS_length, max(INS_length/normal_INS_length, normal_INS_length/INS_length))
					if max(INS_length/normal_INS_length, normal_INS_length/INS_length) <= length_ratio:
						SV_tmp = ",".join((normal_SV))
						normal_SV_out.append(SV_tmp)
					elif INS_length >= length_filter_min and normal_INS_length >= length_filter_min:
						SV_tmp = ",".join((normal_SV))
						normal_SV_out.append(SV_tmp)
					if int(line_l[1]) + distance4 < int(normal_SV[1]):
						break
	else:
		if line_l[0] in normal_SV_list:
			for normal_SV in normal_SV_list[line_l[0]]:
				if normal_SV[0] == line_l[0] and normal_SV[2] == line_l[2] and abs(int(normal_SV[1]) - int(line_l[1])) <= distance4 and abs(int(normal_SV[3]) - int(line_l[3])) <= distance4:
					
					SV_tmp = ",".join((normal_SV))
					normal_SV_out.append(SV_tmp)
#					if int(line_l[1]) + distance < int(normal_SV[1]):
#						break

		if line_l[2] in normal_SV_list:
			for normal_SV in normal_SV_list[line_l[2]]:
				if normal_SV[2] == line_l[0] and normal_SV[0] == line_l[2] and abs(int(normal_SV[3]) - int(line_l[1])) <= distance4 and abs(int(normal_SV[1]) - int(line_l[3])) <= distance4:
					SV_tmp = ",".join((normal_SV))
					normal_SV_out.append(SV_tmp)
#					if int(line_l[3]) + distance < int(normal_SV[1]):
#						break

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

	if len(normal_SV_out) == 0:
		normal_SV_out = ["-"]

	print("\t".join(line_l[0:slice_col]), str(depth1) + "," + str(depth2), ";".join(normal_SV_out), "\t".join(line_l[slice_col:len(line_l)]), sep="\t")
