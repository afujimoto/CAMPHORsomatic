import sys

#chr10   930674  chr10   934257  2       DEL     3583    20,26   -       0       2,2,2   14,14   0.143   0,14/0,14       2       0       0/2     |||     54460fbf-e16d-432e-9709-819efedf7d06;within;chr10;930602;chr10;934162;8298;46;8278;926799;938174;15;60;-;-;-;560|7c74f7b2-ce59-42c0-8983-575eb1e8eb3e;within;chr10;930674;chr10;934257;9153;32;9114;926386;938475;17.3;60;-;-;-;652     chr10;930602;chr10;934162;GATCTCACGTTTGAAAGAATG;3560;54460fbf-e16d-432e-9709-819efedf7d06;46,8278;chr10;8298;926799,938174;-;15;60;20S13M5D2M1D23M1D2M2D1M1I26M1I7M1D10M2I17M3I2M1I10M1I5M3D11M1I6M1I17M1I25M1D5M2D9M1I11M2I15M3D5M7D17M1I15M2D15M1D6M1I23M4

f = open(sys.argv[1])
low_mq_read_col = int(sys.argv[2])
low_mq_read_prop = float(sys.argv[3])

for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	low_mq_prop_l = line_l[low_mq_read_col].split("/")
	low_mq_prop_l_1 = low_mq_prop_l[0].split(",")
	low_mq_prop_l_2 = low_mq_prop_l[1].split(",")
	if (float(low_mq_prop_l_1[1]) and float(low_mq_prop_l_1[0])/float(low_mq_prop_l_1[1]) >= low_mq_read_prop) or (float(low_mq_prop_l_2[1]) and float(low_mq_prop_l_2[0])/float(low_mq_prop_l_2[1]) >= low_mq_read_prop):
		continue

	print(line)
