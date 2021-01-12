import re 
import sys
import subprocess
#from Bio.Align.Applications import MafftCommandline

#chr1   63980559        chr1    63982331        FF
#chr1   63980566        chr1    63982331        FF

def mafft_align(read_name, fasta_file, match_cutoff):
#	print(read_name)
	read_name_l = re.split(":", read_name)
	fasta_f = open(fasta_file, "w")
	seq_len = 0
	for read_name_tmp in read_name_l:
		read_name_tmp_l = re.split(";", read_name_tmp)
		seq_len = len(read_name_tmp_l[4])
		print(">", read_name_tmp_l[6], "\n", read_name_tmp_l[4], sep="", file=fasta_f)
	fasta_f.close()

	cmd = ""
	if len(read_name_l) < 100:
		cmd = "/usr/local/bin/mafft --ep 0.0 --op 1 --maxiterate 1000 --globalpair " + fasta_file
	else:
		cmd = "/usr/local/bin/mafft --retree 2 --maxiterate 2 " + fasta_file

	contents = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
	contents = str(contents)
	contents = contents.replace("b\'", "")
#	contents = contents.replace('\'', '')
	contents_l = re.split("\\\\n", contents)
#	print(contents_l)

	name = ""
	seq_dic = {}
	for tmp in contents_l:
		tmp = tmp.replace("\n", "")
		if tmp == "\'":
			continue
		if tmp[0] == ">":
			name = tmp
			seq_dic[name] = ""
		else:
			seq_dic[name] += tmp

	seq_len = 0
	for name in seq_dic:
		if len(seq_dic[name]) > seq_len:
			seq_len = len(seq_dic[name])

	consensus = ""
	for i in range(seq_len):
		seq_tmp = {}
		for name in seq_dic:
			seq_tmp[seq_dic[name][i]] = seq_tmp.get(seq_dic[name][i], 0) + 1

		consensus_tmp = ""
		for tmp in seq_tmp:
			if seq_tmp[tmp] > len(seq_dic)/2:
				consensus_tmp += tmp
				break

		if consensus_tmp == "":
			consensus_tmp = "-"

		consensus += consensus_tmp

	print(consensus)	

	diff_rate = {}
	for name in seq_dic:
		same = 0
		diff = 0
		
		for i in range(len(consensus)):
			if consensus[i] == seq_dic[name][i]:
				same += 1
			else:
				diff += 1

		diff_rate_tmp = float(diff)/float(same + diff)

		diff_rate[name] = diff_rate_tmp

#	print("\n")
#	print(diff_rate)
#	print("\n")

	num_match_read = 0
	for name in diff_rate:
		if diff_rate[name] <= match_cutoff:
			num_match_read += 1

	for name in seq_dic:
		num_base = 0
		for s in seq_dic[name]:
			if s != "-":
				num_base += 1
		seq2 = ""
		num_base2 = 0
		for s in seq_dic[name]:
			if s != "-":
				num_base2 += 1
			seq2 += s
			if num_base2 == 200:
				seq2 += "<>"

		print(name, num_base, seq2)

	return num_match_read

	
def define_BP(BP_grp2, max_or_min, num_element):
        candidate_l = re.split('\t', BP_grp2[0])

        if max_or_min == "max":
                max_BP = candidate_l[num_element]
                for candidate in BP_grp2:
                        candidate_l = re.split('\t', candidate)
                        if max_BP <= candidate_l[num_element]:
                                max_BP = candidate_l[num_element]
                return max_BP
        elif max_or_min == "min":
                min_BP = candidate_l[num_element]
                for candidate in BP_grp2:
                        candidate_l = re.split('\t', candidate)
                        if min_BP >= candidate_l[num_element]:
                                min_BP = candidate_l[num_element]
                return min_BP

def get_BP2(BP_grp):
        BP_grp_dic = {}
        for BP_tmp in BP_grp:
                BP_tmp_l = re.split("\t", BP_tmp)
                BP_grp_dic.setdefault(BP_tmp, {})["chr1"] = BP_tmp_l[1]
                BP_grp_dic.setdefault(BP_tmp, {})["bp1"] = int(BP_tmp_l[2])
                BP_grp_dic.setdefault(BP_tmp, {})["chr2"] = BP_tmp_l[3]
                BP_grp_dic.setdefault(BP_tmp, {})["bp2"] = int(BP_tmp_l[4])

        pre_BP_tmp = []
        BP_grp2 = []
        breakpoints = []
        for BP_tmp_t in sorted(BP_grp_dic.items(), key=lambda x:x[1]['bp1']):
                BP_tmp = re.split("\t", BP_tmp_t[0])
#                print("###########", BP_tmp)
                if len(pre_BP_tmp) > 0 and abs(int(BP_tmp[2]) - int(pre_BP_tmp[2])) > SV_range:
#                        for tmp in BP_grp2:
#                                print("<><><><>", tmp)
                        defined_BP1 = define_BP(BP_grp2, "max", 2)
                        defined_BP2 = define_BP(BP_grp2, "min", 4)
                        support_read = len(BP_grp2)
                        support_read_name = []
                        for read_tmp in BP_grp2:
                                read_tmp = re.split("\t", read_tmp)
                                support_read_name.append(read_tmp[0])
                        breakpoint_tmp = str(defined_BP1) + "\t" + str(defined_BP2) + "\t" + str(support_read) + "\t" + ":".join(support_read_name)
                        breakpoints.append(breakpoint_tmp)
                        BP_grp2_dic = {}
                        BP_grp2_2 = []
                        BP_grp2 = []

                BP_grp2.append("\t".join(BP_tmp))
                pre_BP_tmp = BP_tmp

        defined_BP1 = define_BP(BP_grp2, "max", 2)
        defined_BP2 = define_BP(BP_grp2, "min", 4)
        support_read = len(BP_grp2)
        support_read_name = []
        for read_tmp in BP_grp2:
                read_tmp = re.split("\t", read_tmp)
 #               print("read_tmp",read_tmp)
                support_read_name.append(read_tmp[0])
        breakpoint_tmp = str(defined_BP1) + "\t" + str(defined_BP2) + "\t" + str(support_read) + "\t" + ":".join(support_read_name)
        breakpoints.append(breakpoint_tmp)
#        for tmp in breakpoints:
#                print("!!!!!!", tmp)
        return(breakpoints)
		
def get_BP(BP_grp):
        BP_grp_dic = {}
        for BP_tmp in BP_grp:
                BP_tmp_l = re.split("\t", BP_tmp)
#                print("BP_tmp_l>>>>", BP_tmp_l[1], BP_tmp_l[2], BP_tmp_l[3], BP_tmp_l[4])
                BP_grp_dic.setdefault(BP_tmp, {})["chr1"] = BP_tmp_l[1]
                BP_grp_dic.setdefault(BP_tmp, {})["bp1"] = int(BP_tmp_l[2])
                BP_grp_dic.setdefault(BP_tmp, {})["chr2"] = BP_tmp_l[3]
                BP_grp_dic.setdefault(BP_tmp, {})["bp2"] = int(BP_tmp_l[4])
#        print("BP_grp_dic", BP_grp_dic)
        pre_BP_tmp = []
        BP_grp2 = []

        breakpoints_all = []
        for BP_tmp_t in sorted(BP_grp_dic.items(), key=lambda x:x[1]['bp2']):
                BP_tmp = re.split("\t", BP_tmp_t[0])
#                print("BP_tmp========>", BP_tmp[1:6])
                if len(pre_BP_tmp) > 0 and abs(int(BP_tmp[4]) - int(pre_BP_tmp[4])) > SV_range:
                        breakpoints = get_BP2(BP_grp2)
#                        for tmp in breakpoints:
#                                 print(">>>>>>", tmp)
                        breakpoints_all.extend(breakpoints)
                        BP_grp2 = []

                BP_grp2.append("\t".join(BP_tmp))
                pre_BP_tmp = BP_tmp

        breakpoints = get_BP2(BP_grp2)
        breakpoints_all.extend(breakpoints)
        return (breakpoints_all)

def get_distance(read):
	read_l = read.split(":")
	max_len = 0
	distance = 0
	for read_tmp in read_l:
		read_tmp_l = read_tmp.split(";")
		if max_len < int(read_tmp_l[9]):
			max_len = int(read_tmp_l[9])
			distance = read_tmp_l[5]
	return distance

SV_range = int(sys.argv[2])
SV_target_type = sys.argv[3]
#mafft_tmp_dir = sys.argv[4]
#match_cutoff = float(sys.argv[5])
SV_target_type_dic = {}
for SV_target_tmp in re.split(",", SV_target_type):
        SV_target_type_dic[SV_target_tmp] = 1

#print (SV_target_type_dic)

"""
chr20;23114936;chr20;23115121;ACTATATAGCAACTACTAAGCGAAACGTTACAAATATGATGAACAAGCCTAGAACACCGTGGTTCAATTGCTGCAGTGATGTAACACAACAACAAAAATGCAGCGCCTTACCAATCTCAGCACTAAAACAGCATGGGCTCTATGCCGTTCTGGTCAGTGCTGACTCCATTTGCCAGCAGGCACTGGGCTATCACGGCCTTCGCCGCAATGTTTTGCAACTCGTGCTGGGATTACAGTACCTGCCACCATGCCCAATAATTTTGTATTTTAATAGAGGCGAATACCATTTTGTGACTGGTTGGCGAACTGACCTCAAATGTTTCTCAGCCTCCCAAAGTGTGATTTACAGTAACCTGCGCCAACCATGGTTTATCTGCAAAAACGGTTCTGTAGGCCTCTAA;185;1b5fa6af-32a0-4d9f-b24a-e533e80c6e61_Basecall_1D_template;2548,24129;chr20;24135;23114076,23138665;+;13.9   chr20   23114936        chr20   23115121        DEL     185     chr20   23114936        chr20   23115121        DEL     185     ACTATATAGCAACTACTAAGCGAAACGTTACAAATATGATGAACAAGCCTAGAACACCGTGGTTCAATTGCTGCAGTGATGTAACACAACAACAAAAATGCAGCGCCTTACCAATCTCAGCACTAAAACAGCATGGGCTCTATGCCGTTCTGGTCAGTGCTGACTCCATTTGCCAGCAGGCACTGGGCTATCACGGCCTTCGCCGCAATGTTTTGCAACTCGTGCTGGGATTACAGTACCTGCCACCATGCCCAATAATTTTGTATTTTAATAGAGGCGAATACCATTTTGTGACTGGTTGGCGAACTGACCTCAAATGTTTCTCAGCCTCCCAAAGTGTGATTTACAGTAACCTGCGCCAACCATGGTTTATCTGCAAAAACGGTTCTGTAGGCCTCTAA       1b5fa6af-32a0-4d9f-b24a-e533e80c6e61_Basecall_1D_template       24135   2548,24129      chr20   23114076,23138665       +       13.9    60      2547S16M2D8M1I5M1D4M1D15M2D18M2D13M1D36M1D1M1D4M7D6M1D7M2D6M2D5M1D10M2D1M2D14M2D5M1D23M1D7M3D6M2D8M1D4M1D19M1D12M1I3M2D4M2D6M3D15M4D15M1D17M1D12M1D6M2D25M1D4M3D5M2D5M1D15M2D9M2D2M3I37M2D6M1D12M1D8M2I4M1I3M3D7M1D4M3D5M1D21M2D8M2D22M1I5M5D18M1D9M1D14M1D18M1D1M2D10M2I1M3D36M3D1M3D5M1D7M1I2M2D10M1D6M6D2M3D4M3D6M2I8M4D10M1D6M1D9M1D4M1I9M185D6M7I2M2I3M4D2M1D11M1D20M1D4M1D7M1D17M3D9M3D9M1I5M3D13M4D3M3D18M2D8M3D4M2D1M1D6M1D9M1D8M3D2M3D16M2D17M1D5M2D1M1D17M1D1M1D5M1D2M2I21M3D30M1D7M1D8M3D8M3D7M2I5M3D3M2D23M1D7M1D27M1I8M1D6M2D5M2D6M2D1M2D6M3D3M1D2M3D5M1D9M1D7M3D10M2D4M1D4M1D10M1D2M1D5M2D23M2D5M1D19M1D31M2D4M6D5M1D9M2I9M3D9M2I9M2D12M2D2M1D7M2D6M1D5M1D13M3D5M1D4M1D12M1D5M2D4M1D6M1I13M1D2M1D9M2D3M1D11M2D6M2D20M2D12M4D16M1D5M1I4M1I3M1I25M1D2M1D6M2D13M1D6M1D12M1D35M1D1M1D17M1D13M1D6M1D11M1D4M2D7M2D11M7D1M1D27M3D5M1D4M2D1M1I3M1D7M2D6M1D7M2D27M1D16M1D25M1I2M5D23M4D6M2D8M2D8M1I2M3D10M2D11M1D37M2D6M1D6M1D8M1D6M1D25M2D10M1D3M3D12M2D1M1D11M2D5M1I7M1D7M1D8M1D4M2D4M5D7M1D8M2D41M1I16M1D21M1I23M1I16M1D3M5D14M4D1M1D12M1D10M1D1M1D9M1D13M1D21M1D4M1D2
"""

print("#chr1\tBP1\tchr2\tBP2\tsupport_read\tSV_type\tread_name")

pre_BP_l = ""
BP_grp = []
SV_BP_f = open(sys.argv[1])
for line in SV_BP_f:
        line = line.replace('\n', '')
        line_l = re.split('\t', line)

        if not line_l[5] in SV_target_type_dic:
                continue

        if len(pre_BP_l) > 0 and abs(int(line_l[2]) - int(pre_BP_l[2])) > SV_range:
                if len(BP_grp) >= 1:
                        chr2_dic = {}
                        for BP_tmp in BP_grp:
                               BP_tmp_l = re.split("\t", BP_tmp)
                               chr2_dic[BP_tmp_l[3]] = 1

                        breakpoints = []
                        for chr2 in chr2_dic:
                                BP_grp_tmp = []
                                for BP_tmp in BP_grp:
                                        BP_tmp_l = re.split("\t", BP_tmp)
                                        if chr2 == BP_tmp_l[3]:
                                                BP_grp_tmp.append(BP_tmp)
                                breakpoints_tmp = get_BP(BP_grp_tmp)
                                breakpoints.extend(breakpoints_tmp)

                        for breakpoint_tmp in breakpoints:
                                breakpoint_tmp_l = re.split("\t", breakpoint_tmp)
                                breakpoint_tmp_l2 = re.split(";", breakpoint_tmp_l[3])
                                SV_type = line_l[5]
                                if line_l[5] == "FF" or line_l[5] == "RR":
                                        SV_type = "DEL"
#                                distance = int(breakpoint_tmp_l[1]) - int(breakpoint_tmp_l[0])
                                distance = get_distance(breakpoint_tmp_l[3])
                                if SV_type == "CHR":
                                        distance = "NA"

                                num_match_read = 0
                                if int(breakpoint_tmp_l[2]) >= 2:
                                        SV_name = "_".join([pre_BP_l[1], str(breakpoint_tmp_l[0]), pre_BP_l[3], str(breakpoint_tmp_l[1]), breakpoint_tmp_l[2], SV_type])
#                                        fasta_file = mafft_tmp_dir + "/" + SV_name + ".fas"
                                        num_match_read = "NA"

                                print (breakpoint_tmp_l2[0], breakpoint_tmp_l[0], breakpoint_tmp_l2[2], breakpoint_tmp_l[1], breakpoint_tmp_l[2], SV_type, str(distance), str(num_match_read), breakpoint_tmp_l[3], sep="\t")
                BP_grp = []

        BP_grp.append(line)
#        print(BP_grp)
        pre_BP_l = line_l

if len(BP_grp) >= 1:
        chr2_dic = {}
        for BP_tmp in BP_grp:
               BP_tmp_l = re.split("\t", BP_tmp)
               chr2_dic[BP_tmp_l[3]] = 1

        breakpoints = []
        for chr2 in chr2_dic:
                BP_grp_tmp = []
                for BP_tmp in BP_grp:
                        BP_tmp_l = re.split("\t", BP_tmp)
                        if chr2 == BP_tmp_l[3]:
#                               print(BP_tmp_l[3])
                               BP_grp_tmp.append(BP_tmp)
                breakpoints_tmp = get_BP(BP_grp_tmp)
                breakpoints.extend(breakpoints_tmp)

        for breakpoint_tmp in breakpoints:
                breakpoint_tmp_l = re.split("\t", breakpoint_tmp)
                breakpoint_tmp_l2 = re.split(";", breakpoint_tmp_l[3])
                BP_tmp = re.split("\t", BP_grp[0])
                SV_type = BP_tmp[5]
                if SV_type == "FF" or SV_type == "RR":
                        SV_type = "DEL"
#                distance = int(breakpoint_tmp_l[1]) - int(breakpoint_tmp_l[0])
                distance = get_distance(breakpoint_tmp_l[3])
                if SV_type == "CHR":
                        distance = "NA"

                num_match_read = 0
                if int(breakpoint_tmp_l[2]) >= 2:
                        SV_name = "_".join([pre_BP_l[1], str(breakpoint_tmp_l[0]), pre_BP_l[3], str(breakpoint_tmp_l[1]), breakpoint_tmp_l[2], SV_type])
#                        fasta_file = mafft_tmp_dir + "/" + SV_name + ".fas"
                        num_match_read = "NA"
#                        num_match_read = mafft_align(breakpoint_tmp_l[3], fasta_file, match_cutoff)

                print (breakpoint_tmp_l2[0], breakpoint_tmp_l[0], breakpoint_tmp_l2[2], breakpoint_tmp_l[1], breakpoint_tmp_l[2], SV_type, str(distance), str(num_match_read), breakpoint_tmp_l[3], sep="\t")
