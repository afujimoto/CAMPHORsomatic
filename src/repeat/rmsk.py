import sys

#585     463     13      6       17      chr1    10000   10468   -248945954      +       (TAACCC)n       Simple_repeat   Simple_repeat   1       471     0
#585     18      232     0       19      chr1    15797   15849   -248940573      +       (TGCTCC)n       Simple_repeat   Simple_repeat   1       52      0

#chr1	10000	10468	(TAACCC)n	585,463,13,6,17,chr1,10000,10468,-248945954,+,(TAACCC)n,Simple_repeat,Simple_repeat,1,471,0,1
#chr1	15797	15849	(TGCTCC)n	585,18,232,0,19,chr1,15797,15849,-248940573,+,(TGCTCC)n,Simple_repeat,Simple_repeat,1,52,0,6

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	print("\t".join(line_l[5:8]), line_l[10], ",".join(line_l), sep="\t")
