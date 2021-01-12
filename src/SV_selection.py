import sys

##chr1   BP1     chr2    BP2     support_read    SV_type read_name
#chr1    10154   chr1    10385   1       DEL     231     0       chr1;10154;chr1;10385;;231;ERR3219857.2117920;3381,14702;chr1;1;10040,21666;+;-;0;3380S12M
#chr1    10193   chr1    10470   1       DEL     277     0       chr1;10193;chr1;10470;;277;ERR3219855.1240223;1877,3605;chr1;1;10054,11448;+;-;0;1876S4M1D

f = open(sys.argv[1])
read_number = int(sys.argv[2])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line[0] == "#":
#		print(line)
		continue

	if int(line_l[4]) >= read_number:
		print(line)
