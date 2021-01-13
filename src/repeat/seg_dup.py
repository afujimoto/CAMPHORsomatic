import sys

##bin    chrom   chromStart      chromEnd        name    score   strand  otherChrom      otherStart      otherEnd        otherSize       uid     posBasesHit     testResult
#1       chr1    16761934        16799163        chr1:234783386  0       -       chr1    234783386       234820600       37214   289     1000    N/A     N/A     N/A     N/A
#23      chr1    120533238       120705592       chr1:148804200  0       -       chr1    148804200       148977664       173464  661     1000    N/A     N/A     N/A     N/A
#23      chr1    120575596       120705592       chr1:149194882  0       +       chr1    149194882       149326044       131162  664     1000    N/A     N/A     N/A     N/A

#chr1    10000   19844   chrX:156020216-156030574
#chr1    10000   20818   chr12:10043-20853

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if line[0] == "#":
		continue

	print(line_l[1], line_l[2], line_l[3], line_l[7] + ":" + line_l[8] + "-" + line_l[9], sep="\t") 
