import sys

#585	chr1	10000	10468	trf	6	77.2	6	95	3	789	33	51	0	15	1.43	TAACCC
#585	chr1	10627	10800	trf	29	6	29	100	0	346	13	38	47	0	1.43	AGGCGCGCCGCGCCGGCGCAGGCGCAGAG
#585	chr1	10757	10997	trf	76	3.2	76	95	2	434	17	30	45	6	1.73	GGCGCAGGCGCAGAGAGGCGCGCCGCGCCG

#chr1    10000   10468   trf,6,77.2,TAACCC
#chr1    10627   10800   trf,29,6,AGGCGCGCCGCGCCGGCGCAGGCGCAGAG

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	print(line_l[1], line_l[2], line_l[3], ",".join(line_l[4:7]) + "," + line_l[-1], sep="\t")
