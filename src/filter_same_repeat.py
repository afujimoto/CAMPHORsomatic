import sys
import re

#chr9,6727603,6728097,chr9:30690509-30690990,1.0,within
#chr10   1293159 chr10   1293297 2       DEL     139     NA      2       NS;2/DIS;2/MQ;2/CN;2/CI;2       0.8/12/176/0.69 0.0/NA/0/NA     1,1     NA

f = open(sys.argv[1])
rep1_col = int(sys.argv[2])
rep2_col = int(sys.argv[3])

for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	
#	print(line)

	chr1, pos1, chr2, pos2 = line_l[0], int(line_l[1]), line_l[2], int(line_l[3])

	rep1 = line_l[rep1_col].split("\|")
	rep2 = line_l[rep2_col].split("\|")

	ol_rep = []

	found_in_rep1 = 0
	if line_l[rep1_col] != "-":
		for rep_tmp in rep1:
			rep_tmp_l = rep_tmp.split(",")
			partner_l = re.split(r':|-', rep_tmp_l[3])

#			if rep_tmp_l[0] == partner_l[0] and rep_tmp_l[1] == partner_l[1] and rep_tmp_l[2] == partner_l[2]:
#				continue
#			elif rep_tmp_l[0] == partner_l[0] and int(partner_l[1]) <= int(rep_tmp_l[1]) <= int(partner_l[2]):
#				continue
#			elif rep_tmp_l[0] == partner_l[0] and int(partner_l[1]) <= int(rep_tmp_l[2]) <= int(partner_l[2]):
#				continue
#			elif rep_tmp_l[0] == partner_l[0] and int(rep_tmp_l[1]) <= int(partner_l[1]) <= int(rep_tmp_l[2]):
#				continue

			if (chr1 == rep_tmp_l[0] and int(rep_tmp_l[1]) <= pos1 <= int(rep_tmp_l[2])) and (chr2 == partner_l[0] and int(partner_l[1]) <= pos2 <= int(partner_l[2])):
				found_in_rep1 += 1
				ol_rep.append(rep_tmp)
			if (chr2 == rep_tmp_l[0] and int(rep_tmp_l[1]) <= pos2 <= int(rep_tmp_l[2])) and (chr1 == partner_l[0] and int(partner_l[1]) <= pos1 <= int(partner_l[2])):
				found_in_rep1 += 1
				ol_rep.append(rep_tmp)

			if (chr1 == rep_tmp_l[0] and int(rep_tmp_l[1]) <= pos1 <= int(rep_tmp_l[2])) and (chr2 == rep_tmp_l[0] and int(rep_tmp_l[1]) <= pos2 <= int(rep_tmp_l[2])):
				found_in_rep1 += 1
				ol_rep.append(rep_tmp)

	found_in_rep2 = 0
	if line_l[rep2_col] != "-":
		for rep_tmp in rep2:
			rep_tmp_l = rep_tmp.split(",")
			partner_l = re.split(r'-|:', rep_tmp_l[3])

#			if rep_tmp_l[0] == partner_l[0] and rep_tmp_l[1] == partner_l[1] and rep_tmp_l[2] == partner_l[2]:
#				continue
#			elif rep_tmp_l[0] == partner_l[0] and int(partner_l[1]) <= int(rep_tmp_l[1]) <= int(partner_l[2]):
#				continue
#			elif rep_tmp_l[0] == partner_l[0] and int(partner_l[1]) <= int(rep_tmp_l[2]) <= int(partner_l[2]):
#				continue
#			elif rep_tmp_l[0] == partner_l[0] and int(rep_tmp_l[1]) <= int(partner_l[1]) <= int(rep_tmp_l[2]):
#				continue

#			if rep_tmp_l[0] == partner_l[0] and rep_tmp_l[1] == partner_l[1] and rep_tmp_l[2] == partner_l[2]:
#				continue

			if (chr1 == rep_tmp_l[0] and int(rep_tmp_l[1]) <= pos1 <= int(rep_tmp_l[2])) and (chr2 == partner_l[0] and int(partner_l[1]) <= pos2 <= int(partner_l[2])):
				found_in_rep1 += 1
				ol_rep.append(rep_tmp)
			if (chr2 == rep_tmp_l[0] and int(rep_tmp_l[1]) <= pos2 <= int(rep_tmp_l[2])) and (chr1 == partner_l[0] and int(partner_l[1]) <= pos1 <= int(partner_l[2])):
				found_in_rep1 += 1
				ol_rep.append(rep_tmp)

			if (chr1 == rep_tmp_l[0] and int(rep_tmp_l[1]) <= pos1 <= int(rep_tmp_l[2])) and (chr2 == rep_tmp_l[0] and int(rep_tmp_l[1]) <= pos2 <= int(rep_tmp_l[2])):
				found_in_rep1 += 1
				ol_rep.append(rep_tmp)
	uniq = {}
	for rep in ol_rep:
		uniq[rep] = 1
	ol_rep2 = uniq.keys()

	if len(ol_rep2) == 0:
		print (line, sep="\t")
