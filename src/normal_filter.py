#chr10   264662  chr10   264834  3       DEL     172     22,22   -       NA      chr10;264594;chr10;264834;GTTGACTCATAACTGCCGTTGCACAGGTATGGGAGTGATGCAAAAGGTGTCTGATTACGCTGTTACTAG
#chr10   264928  chr10   265001  14      DEL     73      22,23   -       NA      chr10;264826;chr10;265005;CAATGGTATGTCAATAGGCTGACATAACTGTAGGGTACAGGCATGGAGTACGCAAGACCGCGTGGCTGT
import sys

f = open(sys.argv[1])
normal_depth_col = int(sys.argv[2])
normal_SV_col = int(sys.argv[3])
min_normal_depth = int(sys.argv[4])

for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	depth_l = line_l[normal_depth_col].split(",")
	
	if line_l[normal_SV_col] == "-" and (int(depth_l[0]) + int(depth_l[1]))/2 >= min_normal_depth:
		print(line)
