#chr1    10000   10468   trf,6,77.2,TAACCC
#chr1    10627   10800   trf,29,6,AGGCGCGCCGCGCCGGCGCAGGCGCAGAG

#chr10   10014858        chr10   10015285        2       DEL     427     chr10;10014824;chr10;10015288;d2ad20f5-4ed8-4d06-892a-e6c66dc7129b_Basecall_1D_template;464;4877;c
#chr10   10169954        chr10   10170424        2       DEL     470     chr10;10169906;chr10;10170472;e74c59d0-a94b-4208-8673-78da9cc41f38_Basecall_1D_template;566;12364;
import re
import sys

selfchain_f = open(sys.argv[1])
deletion_candidate_f = open(sys.argv[2])
max_del_len = int(sys.argv[3])

chain_list = {}

for line in selfchain_f:
	line = line.replace('\n', '')
	chain_l = re.split("\t", line)
	chain_list.setdefault(str(chain_l[0]), []).append(chain_l)

for row in deletion_candidate_f:
	row = row.replace('\n', '')
	row_l = re.split('\t', row)

	if row_l[0] == "#chr1":
		continue
	
	if int(row_l[3]) - int(row_l[1]) >= max_del_len:
#		print(">>>>>>>>>>>", row_l)
		out = ["-"]
		print(row, "|".join(out), sep="\t")
		continue

	cover = ""
	coverage_prop = ""
	output = ""
	out = []

	if row_l[0] in chain_list:
		for line in chain_list[row_l[0]]:
#			print(line)
			if int(line[2]) > int(row_l[1]) >= int(line[1]) and int(line[1]) < int(row_l[3]) <= int(line[2]):
				coverage = "complete"
		#		coverage_prop = str(round((abs(int(line[2]) - int(line[1])))/(abs(int(row_l[3]) - int(row_l[1]))), 2))
				coverage_prop = str(1.00)
				out_tmp = [",".join(line), coverage_prop, coverage]
				out_tmp_str = ",".join(out_tmp)
				out.append(out_tmp_str)

			elif int(line[1]) < int(row_l[1]) < int(line[2]) and int(row_l[3]) > int(line[2]):
				coverage = "partial"
				coverage_prop = str(round((abs(int(line[2]) - int(row_l[1])))/(abs(int(row_l[3]) - int(row_l[1]))), 2))
				out_tmp = [",".join(line), coverage_prop, coverage]
				out_tmp_str = ",".join(out_tmp)
				out.append(out_tmp_str)

			elif int(row_l[1]) < int(line[1]) and int(line[1]) < int(row_l[3]) < int(line[2]):
				coverage = "partial"
				coverage_prop = str(round((abs(int(row_l[3]) - int(line[1])))/(abs(int(row_l[3]) - int(row_l[1]))), 2))
				out_tmp = [",".join(line), coverage_prop, coverage]
				out_tmp_str = ",".join(out_tmp)
				out.append(out_tmp_str)

			elif int(row_l[3]) > int(line[1]) > int(row_l[1]) and int(row_l[1]) < int(line[2]) < int(row_l[3]):
				coverage = "within"
				coverage_prop = str(round((abs(int(line[2]) - int(line[1])))/(abs(int(row_l[3]) - int(row_l[1]))), 2))
				out_tmp = [",".join(line), coverage_prop, coverage]
				out_tmp_str = ",".join(out_tmp)
				out.append(out_tmp_str)
		
			else:
				coverage = "none"	

			if int(line[1]) > int(row_l[3]):
				break

	if len(out) == 0:
		out = ["-"]

	print(row, "|".join(out), sep="\t")

