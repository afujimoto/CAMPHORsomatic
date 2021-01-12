import re
import sys

#chr7    100958268       chr7    100958513       DEL     245     45a434d9-692a-4367-a8ad-3dd95b1f8af3_Basecall_1D_template/12.2  452     137,322;27,45;370,389;85,103;347,3
#chr7    100958519       chr7    100958871       DEL     352     45a434d9-692a-4367-a8ad-3dd95b1f8af3_Basecall_1D_template/12.2  452     137,322;27,45;370,389;85,103;347,3

#target_pos = re.split(":|-", sys.argv[2])

#chr = target_pos[0]
#start = int(target_pos[1])
#end = int(target_pos[2])

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)

	if line_l[4] == "BR" and line_l[0] == "NA":
		line_l[4] = "lBR"
		line_l[0], line_l[1] = line_l[2], line_l[3]
		print("\t".join(line_l))
	elif line_l[4] == "BR" and line_l[2] == "NA":
		line_l[4] = "rBR"
		line_l[2], line_l[3] = line_l[0], line_l[1]
		print("\t".join(line_l))
	else:
		print(line)
#	elif line_l[0] == chr and start <= int(line_l[1]) <= end:
#		print(line)

"""
	if line_l[4] == "BR" and line_l[0] == "NA" and line_l[2] == chr and start <= int(line_l[3]) <= end:
		line_l[4] = "lBR"
		line_l[0], line_l[1] = line_l[2], line_l[3]
		print("\t".join(line_l))
	elif line_l[4] == "BR" and line_l[2] == "NA" and line_l[0] == chr and start <= int(line_l[1]) <= end:
		line_l[4] = "rBR"
		line_l[2], line_l[3] = line_l[0], line_l[1]
		print("\t".join(line_l))
	elif line_l[0] == chr and start <= int(line_l[1]) <= end:
		print(line)
"""
