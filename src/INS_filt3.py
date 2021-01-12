import sys

#perl -ne '@l=split("\t"); @depth=split(",", $l[11]); @rep=split("/", $l[13]); @rep1=split(",", $rep[0]); @rep2=split(",", $rep[1]); if(($depth[0] + $depth[1])/2 < 10){next;} if($l[12] eq "-" and ($rep1[1] and $rep1[0]/$rep1[1] < 0.3) and ($rep2[1] and $rep2[0]/$rep2[1] < 0.3) and $l[4] >= 4){print}'

f = open(sys.argv[1])
min_depth = int(sys.argv[2]) #10
min_prop_lowq_reads = float(sys.argv[3]) #0.3
min_support_read = int(sys.argv[4]) #4
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	depth_l = line_l[11].split(",")

	rep_l = line_l[13].split("/")

	rep_l1 = rep_l[0].split(",")
	rep_l2 = rep_l[1].split(",")

	if (int(depth_l[0]) + int(depth_l[1]))/2 < min_depth:
#		print("depth", (int(depth_l[0]) + int(depth_l[1]))/2, line)
		continue

	if line_l[12] == "-" and (rep_l1[1] and int(rep_l1[0])/int(rep_l1[1]) < min_prop_lowq_reads) and (rep_l2[1] and int(rep_l2[0])/int(rep_l2[1]) < min_prop_lowq_reads) and int(line_l[4]) >= min_support_read:
#		print(line_l[12], int(rep_l1[0])/int(rep_l1[1]), int(rep_l2[0])/int(rep_l2[1]), line_l[4])
		print(line)
#	else:
#		print(">>>", line_l[12], int(rep_l1[0])/int(rep_l1[1]), int(rep_l2[0])/int(rep_l2[1]), line_l[4])
#		print("!!!", line)
