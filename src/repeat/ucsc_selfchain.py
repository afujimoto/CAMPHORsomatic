#dL_f
#bin    score   tName   tSize   tStart  tEnd    qName   qSize   qStrand qStart  qEnd    id      normScore
#9      27375   chr1    248956422       1048486 1056733 chr22   50818468        +       39616480        39618495        2933611 17.9
#11     5215    chr1    248956422       17825763        17826504        chr5    181538259       +       77785246        77785973        52571646        17.7
#12     74963   chr1    248956422       28310607        28312400        chr5    181538259       +       74881392        74882758        503299  66.2
#12     30071   chr1    248956422       28310619        28311678        chr15   101991189       -       70062229        70062669        2475342 70.6

#585     897863  chr1    248956422       10000   19844   chrY    57227415        -       10435   20679   14932   92
#585     897863  chr1    248956422       10000   19844   chrX    156040895       -       10435   20679   14931   92
#585     2691807 chr1    248956422       10000   39238   chr12   133275309       +       10043   41283   3727    92.6
import sys
import re

query_region_dict = {}

ucsc_self_chain_f = open(sys.argv[1])
#target_chr = sys.argv[2]

for line in ucsc_self_chain_f:
	line = line.replace('\n', '')
	line_l = re.split('\t', line)
	
	if line_l[0] == "#bin":
		continue

#	if line_l[2] != target_chr:
#		continue

	region = str(line_l[2]) + "\t" + str(line_l[4]) + "\t" + str(line_l[5])
	query_pos_tmp = str(line_l[6]) + ":" + str(line_l[9]) + "-" + str(line_l[10])

	query_region_dict.setdefault(region,[]).append(query_pos_tmp)
	
#	query_region_dict.setdefault("q_pos", []).append(query_pos_tmp)
#	query_region_dict.setdefault("q_size", []).append(str(line_l[7]))
#	query_region_dict.setdefault("q_strand",[]).append(line_l[8])
#	query_region_dict.setdefault("q_id", []).append(str(line_l[11]))
#	query_region_dict.setdefault("normscore", []).append(str(line_l[12]))


#out = (line_l[0], line_l[1], line_l[2], line_l[3], line_l[4], line_l[5], ";".join(query_region_dict["q_pos"]), ",".join(query_region_dict["q_strand"]), ",".join(query_region_dict["normscore"]))
#out = map(str, out)
#print("\t".join(out), sep="\t")

for region in query_region_dict.keys():
	query_out = ",".join(query_region_dict[region])
	print(region, query_out, sep="\t")
