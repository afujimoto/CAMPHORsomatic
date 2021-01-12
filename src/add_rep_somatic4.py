#chr17   60000   100819  chr7:56364573-56412260
#chr17   60000   104243  chr5:208519-249668
#chr10   10869   chr15   1281062 2       CHR     1270194 0       0       NS;2/DIS;2/MQ;0/CN;1/RL;2/RD;2  NA/NA/NA/NA     NA/NA/NA/NA     2,7     NA      53-69,0-0       11.4,10.2
#       5,0,18,0        3543.0,2861.0   7,3     chr10   10869   chr15   1280993 2       CHR     1270124 0       chr10;10816;chr7;1280993;GGTGTGAGCTTGTATGTGCATGAGGCTTGGCATGGATGGCCAA
import re
import sys

selfchain_f = open(sys.argv[1])

analysis_set = 100000

num_lines = 0
for line in selfchain_f:
	num_lines += 1

candidate2chain_1 = {}
candidate2chain_2 = {}
file_pos = 0
num_lines2 = 0
for i in range(0, int(num_lines/analysis_set) + 1):
#	print("#i = ", i)

	start = i*analysis_set
	end = (i + 1)*analysis_set - 1

#	print("#start = ", start, "end = ", end)
#	print("#file_pos = ", file_pos)
	selfchain_f = open(sys.argv[1])
	selfchain_f.seek(int(file_pos))

	chain_list = {}
	for line in selfchain_f:
		num_lines2 += 1
		file_pos += len(line)
#		print(line.replace("\n", ""))

#		if num_lines2 < start or num_lines2 > end:
#			continue

		line2 = line.replace('\n', '')
		chain_l = re.split("\t", line2)

		chain_list.setdefault(str(chain_l[0]), []).append(chain_l)

		if num_lines2 == end:
#			print("#break file_pos", file_pos)
			break

	selfchain_f.close()

	deletion_candidate_f = open(sys.argv[2])
	for row in deletion_candidate_f:
		row = row.replace('\n', '')
		row_l = re.split('\t', row)

		cover = ""
		coverage_prop = ""
		output = ""
		out = []
		out2 = []
		if row_l[0] == chr:
			for line in chain_list[row_l[0]]:
				if int(line[2]) > int(row_l[1]) >= int(line[1]):
					coverage = "within"
					coverage_prop = str(1.00)
					out_tmp = [",".join(line), coverage_prop, coverage]
					out_tmp_str = ",".join(out_tmp)
					out.append(out_tmp_str)

				if int(line[1]) > int(row_l[1]):
					break

		if row_l[2] == chr:
			for line in chain_list[row_l[2]]:
				if int(line[2]) > int(row_l[3]) >= int(line[1]):
					coverage = "within"
					coverage_prop = str(1.00)
					out_tmp = [",".join(line), coverage_prop, coverage]
					out_tmp_str = ",".join(out_tmp)
					out2.append(out_tmp_str)

				if int(line[1]) > int(row_l[3]):
					break

		if len(out) == 0:
			out = ["-"]
		if len(out2) == 0:
			out2 = ["-"]

		if row_l[0] == chr:
			candidate2chain_1[row] = out
		if row_l[2] == chr:
			candidate2chain_2[row] = out2

	deletion_candidate_f.close()

deletion_candidate_f = open(sys.argv[2])
for row in deletion_candidate_f:
	row = row.replace('\n', '')
	row_l = re.split('\t', row)

	if not row in candidate2chain_1:
		candidate2chain_1[row] = ["-"]
	if not row in candidate2chain_2:
		candidate2chain_2[row] = ["-"]
	print(row, "|".join(candidate2chain_1[row]), "|".join(candidate2chain_2[row]), sep="\t")
