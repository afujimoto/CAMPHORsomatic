import re 
import sys
import subprocess
#from Bio.Align.Applications import MafftCommandline
import numpy as np
from operator import itemgetter

#chr1   63980559        chr1    63982331        FF
#chr1   63980566        chr1    63982331        FF

def mafft_align(read_name, fasta_file):
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
	contents_l = re.split("\\\\n", contents)

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

#		if consensus_tmp == "":
#			consensus_tmp = "-"
		if consensus_tmp != "-":
			consensus += consensus_tmp

	return consensus
"""
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

		print(name, num_base, seq2)

	return num_match_read
"""
	
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
                if len(pre_BP_tmp) > 0 and abs(int(BP_tmp[2]) - int(pre_BP_tmp[2])) > SV_range:
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
                        BP_grp2 = []
#                        BP_grp2_2 = []

                BP_grp2.append("\t".join(BP_tmp))
                pre_BP_tmp = BP_tmp

        defined_BP1 = define_BP(BP_grp2, "max", 2)
        defined_BP2 = define_BP(BP_grp2, "min", 4)
        support_read = len(BP_grp2)
        support_read_name = []
        for read_tmp in BP_grp2:
                read_tmp = re.split("\t", read_tmp)
                support_read_name.append(read_tmp[0])
        breakpoint_tmp = str(defined_BP1) + "\t" + str(defined_BP2) + "\t" + str(support_read) + "\t" + ":".join(support_read_name)
        breakpoints.append(breakpoint_tmp)
        return(breakpoints)

def get_INS(BP_grp, ins_range_prop):
        length = []
#        print(BP_grp)
        for BP_tmp in BP_grp:
                BP_tmp_l = re.split("\t", BP_tmp)
#                print(BP_tmp_l[6])
                length.append(int(BP_tmp_l[6]))
        median_length = np.median(length)

#        print("===========================================")
#        print("BP_grp len = ", len(BP_grp))
        diff_BP_grp = []
        same_BP_grp = []
        for i in range(len(BP_grp)):
                BP_tmp_l = re.split("\t", BP_grp[i])
#                print(BP_tmp_l[1:7])
#                print(BP_tmp_l[6], abs(float(BP_tmp_l[6]) - float(median_length))/float(median_length), end=" ")
                if abs(float(BP_tmp_l[6]) - float(median_length))/float(median_length) > ins_range_prop:
#                        print("diff")
                        diff_BP_grp.append(BP_grp[i])
                else:
#                        print("same")
                        same_BP_grp.append(BP_grp[i])
#                        print("=========>", BP_grp[i])

        if len(same_BP_grp) == 0:
                diff_BP_grp = []
                min_distance = median_length
                min_distance_read = ""
                for i in range(len(BP_grp)):
                        BP_tmp_l = re.split("\t", BP_grp[i])
                        if abs(float(BP_tmp_l[6]) - float(median_length)) < min_distance:
                                min_distance = abs(float(BP_tmp_l[6]) - float(median_length))
                                min_distance_read = BP_grp[i]

                for i in range(len(BP_grp)):
                        if BP_grp[i] == min_distance_read:
                                same_BP_grp.append(BP_grp[i])
                        else:
                                diff_BP_grp.append(BP_grp[i])
                return min_distance_read, median_length, same_BP_grp, diff_BP_grp

        min_distance = median_length
        min_distance_read = []
        for i in range(len(BP_grp)):
               BP_tmp_l = re.split("\t", BP_grp[i])
               if abs(int(BP_tmp_l[6]) - median_length) < min_distance:
                       min_distance = abs(int(BP_tmp_l[6]) - median_length) 
                       min_distance_read = BP_grp[i]
 
        return min_distance_read, median_length, same_BP_grp, diff_BP_grp

		
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
#                print("BP_tmp========>", BP_tmp[1:5])
                if len(pre_BP_tmp) > 0 and abs(int(BP_tmp[4]) - int(pre_BP_tmp[4])) > SV_range:
                        breakpoints = get_BP2(BP_grp2)
                        breakpoints_all.extend(breakpoints)
                        BP_grp2 = []

                BP_grp2.append("\t".join(BP_tmp))
                pre_BP_tmp = BP_tmp

        breakpoints = get_BP2(BP_grp2)
        breakpoints_all.extend(breakpoints)
        return (breakpoints_all)
#BP_grp2.setdefault(read_name, []).append(["within", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, read_length, chr, str(r_start) + "," + str(r_end), str(g_start) + "," + str(g_end), bq, mq, strand, cigar])

def get_INS_info(BP_grp):
	BP_grp2 = {}
	for line in BP_grp:
		line_l = line.split("\t")
		bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_name, read_length = line_l[1], line_l[2], line_l[3], line_l[4], line_l[5], line_l[6], line_l[14], line_l[15]

		read_tmp_l = re.split(";", line_l[0])
		if "/" in read_tmp_l[8]:
			chr_l = re.split("/", read_tmp_l[8])
			bq_l = re.split("/", read_tmp_l[12])
			mq_l = re.split("/", read_tmp_l[13])
			cigar_l = re.split("/", read_tmp_l[14])
			strand_l = re.split("/", read_tmp_l[11])
			read_pos = re.split("/", read_tmp_l[7])
			genome_pos = re.split("/", read_tmp_l[10])

			seq_align_pos = []
			for i in range(len(read_pos)):
				read_pos_tmp_l = re.split(",", read_pos[i])
				genome_pos_tmp_l = re.split(",", genome_pos[i])
				seq_align_pos.append([int(read_pos_tmp_l[0]), int(read_pos_tmp_l[1]), int(genome_pos_tmp_l[0]), int(genome_pos_tmp_l[1]), chr_l[i], bq_l[i], mq_l[i], cigar_l[i], strand_l[i]])
			seq_align_pos.sort(key=itemgetter(0))

			BP_found = 0
			for i in range(len(seq_align_pos)):
				r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand = seq_align_pos[i]
				if strand == "-":
					g_start, g_end = g_end, g_start
				if (type == "DEL" or type == "INS") and chr == bp_chr1 and chr == bp_chr2 and int(g_start) <= int(bp_pos1) <= int(g_end) and int(g_start) <= int(bp_pos2) <= int(g_end):
					
					BP_grp2.setdefault(read_name, []).append(["within", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr, str(r_start) + "," + str(r_end), str(g_start) + "," + str(g_end), bq, mq, strand, cigar])
				if i < len(seq_align_pos) - 1:
					r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar2, strand2 = seq_align_pos[i + 1]
					if strand2 == "-":
						g_start2, g_end2 = g_end2, g_start2

					if (chr == bp_chr1 and chr2 == bp_chr2 and int(g_end) == int(bp_pos1) and int(g_start2) == int(bp_pos2)): 
						chr_out = chr + "/" + chr2
						r_pos_out = str(r_start) + "," +  str(r_end) + "/" + str(r_start2) + "," +  str(r_end2)
						gpos_out = str(g_start) + "," + str(g_end) + "/" + str(g_start2) + "," + str(g_end2)
						bq_out = bq + "/" + bq2
						mq_out = mq + "/" + mq2
						strand_out = strand + "/" + strand2
						cigar_out = cigar + "/" + cigar2
						BP_grp2.setdefault(read_name, []).append(["split", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr_out, r_pos_out, gpos_out, bq_out, mq_out, strand_out, cigar_out])
						BP_found = 1
					elif (chr2 == bp_chr1 and chr == bp_chr2 and int(g_start2) == int(bp_pos1) and int(g_end) == int(bp_pos2)):
						chr_out = chr2 + "/" + chr
						r_pos_out = str(r_start2) + "," +  str(r_end2) + "/" + str(r_start) + "," +  str(r_end)
						gpos_out = str(g_start2) + "," + str(g_end2) + "/" + str(g_start) + "," + str(g_end)
						bq_out = bq2 + "/" + bq
						mq_out = mq2 + "/" + mq
						strand_out = strand2 + "/" + strand
						cigar_out = cigar2 + "/" + cigar
						BP_grp2.setdefault(read_name, []).append(["split", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr_out, r_pos_out, gpos_out, bq_out, mq_out, strand_out, cigar_out])
						BP_found = 1
			if BP_found:
				BP_found = 0
				continue		
				if (type == "DEL" or type == "INS") and chr == bp_chr1 and chr == bp_chr2 and int(g_start) <= int(bp_pos1) <= int(g_end) and int(g_start) <= int(bp_pos2) <= int(g_end):
					
					BP_grp2.setdefault(read_name, []).append(["within", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr, str(r_start) + "," + str(r_end), str(g_start) + "," + str(g_end), bq, mq, strand, cigar])
			for i in range(len(seq_align_pos)):
				r_start, r_end, g_start, g_end, chr, bq, mq, cigar, strand = seq_align_pos[i]
				if strand == "-":
					g_start, g_end = g_end, g_start
				for k in range(len(seq_align_pos)):
					if i == k:
						continue
					r_start2, r_end2, g_start2, g_end2, chr2, bq2, mq2, cigar2, strand2 = seq_align_pos[k]
					if strand2 == "-":
						g_start2, g_end2 = g_end2, g_start2
					if (chr == bp_chr1 and chr2 == bp_chr2 and int(g_end) == int(bp_pos1) and int(g_start2) == int(bp_pos2)):
						chr_out = chr + "/" + chr2
						r_pos_out = str(r_start) + "," +  str(r_end) + "/" + str(r_start2) + "," +  str(r_end2)
						gpos_out = str(g_start) + "," + str(g_end) + "/" + str(g_start2) + "," + str(g_end2)
						bq_out = bq + "/" + bq2
						mq_out = mq + "/" + mq2
						strand_out = strand + "/" + strand2
						cigar_out = cigar + "/" + cigar2
						BP_grp2.setdefault(read_name, []).append(["split", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr_out, r_pos_out, gpos_out, bq_out, mq_out, strand_out, cigar_out])
					elif (chr2 == bp_chr1 and chr == bp_chr2 and int(g_start2) == int(bp_pos1) and int(g_end) == int(bp_pos2)):
						chr_out = chr + "/" + chr2
						r_pos_out = str(r_start2) + "," +  str(r_end2) + "/" + str(r_start) + "," +  str(r_end)
						gpos_out = str(g_start2) + "," + str(g_end2) + "/" + str(g_start) + "," + str(g_end)
						bq_out = bq2 + "/" + bq
						mq_out = mq2 + "/" + mq
						strand_out = strand2 + "/" + strand
						cigar_out = cigar2 + "/" + cigar
						BP_grp2.setdefault(read_name, []).append(["split", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr_out, r_pos_out, gpos_out, bq_out, mq_out, strand_out, cigar_out])
		else:
			chr = read_tmp_l[8]
			r_pos = read_tmp_l[7].split(",")
			genome_pos = read_tmp_l[10].split(",")
			bq = read_tmp_l[12]
			mq = read_tmp_l[13]
			strand = read_tmp_l[11]
			cigar = read_tmp_l[14]
			BP_grp2.setdefault(read_name, []).append(["within", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr, str(r_pos[0]) + "," + str(r_pos[1]), str(genome_pos[0]) + "," + str(genome_pos[1]), bq, mq, strand, cigar])

	return BP_grp2

def find_INS_read(BP_grp):
	INS_read_name = {}
	for BP_grp_tmp in BP_grp:
		BP_grp_tmp_l= BP_grp_tmp.split("\t")
		if BP_grp_tmp_l[5] == "INS":
			INS_read_name[BP_grp_tmp_l[14]] = 1

	return INS_read_name	

#['split', 'chr11', '88003285', 'chr5', '1845931', 'lCHR', 'NA', '11348', 'chr11/chr5', '9538,11328/8798,9522', '88003285,88005075/1845147,1845931', '17/17', '60/2', '+/+', '9537H22M2D20M2D3M1D4M3D15M2I56M2D3M1D30M1I3M1I29M1I31M2I6M2I15M1D24M2D22M5D48M4D18M2I13M1I2M5I29M2I6M1I5M1D4M1I6M1D21M3D2M1D53M1D1M1D13M1I1M2D18M3D9M3D4M1D6M2D6M2D43M2I4M1I23M4D21M1I43M1I33M1D3M3I22M1I6M1D5M1I41M1I14M1I16M1I4M1I9M2I2M2I12M2I13M2I5M1I14M1I14M1D42M1D26M1I8M1I7M1I5M2I41M5D9M2I13M3I11M1D20M1I29M1I3M1I10M1D5M1I16M1I21M1D15M2D18M1I3M1I28M2I5M2D8M2I12M1D38M1D20M1I11M1D52M2D26M1I64M1D24M1I2M2D5M2D15M1D9M1I95M1D7M2I5M2I22M2I17M1D22M1I29M20H/8797H8M1I3M1I10M1D9M3D3M1D3M3D7M1D6M1D4M1D2M1D2M1D14M4D5M1I8M1D1M1D7M1I3M1D16M2D6M1D4M1D6M2I5M1D7M1D3M1D5M1I6M1D2M2D7M1D8M1I4M1D6M1D4M2D17M3D9M1D3M2D4M1D7M1D46M1D17M1I6M1D5M1I11M1D6M1I11M3D1M1D6M1D15M1D6M1D3M1I4M1D15M1I11M1D23M2D11M1D1M1D3M3D17M2D4M4D11M1I1M1I5M1D42M1I11M1D48M1I20M1D22M1I12M1D49M1D3M1D3M1D12M1D13M1D9M1827H']

#BP_grp2.setdefault(read_name, []).append(["within", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr, str(r_start) + "," + str(r_end), str(g_start) + "," + str(g_end), bq, mq, strand, cigar])
#['within', 'chr2', '95946256', 'chr2', '95946256', 'rBR', '1154', '14029', 'chr2', '1154,14012', '95933454,95946256', '17.8', '60', '-', 

def get_INS_length(ins_info_l):
#	print("ins_info_l =>> ", ins_info_l)

	if ins_info_l[0] == "within":
		rpos_l = ins_info_l[9].split(",")
		if ins_info_l[13] == "+" and ins_info_l[5][0] == "r":
			ins_len = int(ins_info_l[7]) - int(rpos_l[1])
		elif ins_info_l[13] == "+" and ins_info_l[5][0] == "l":
			ins_len = int(rpos_l[0]) - 1
		elif ins_info_l[13] == "-" and ins_info_l[5][0] == "r":
			ins_len = int(rpos_l[0]) - 1
		elif ins_info_l[13] == "-" and ins_info_l[5][0] == "l":
			ins_len = int(ins_info_l[7]) - int(rpos_l[1])
#		print(ins_len)
		return ins_len

	elif ins_info_l[0] == "split":
		strand_l = ins_info_l[13].split("/")
#		gpos_l = ins_info_l[10].split("/")
#		gpos_l1 = gpos_l[0].split(",")
#		gpos_l2 = gpos_l[1].split(",")
#		print(ins_info_l[9])
		rpos_l = ins_info_l[9].split("/")
		rpos_l1 = rpos_l[0].split(",")
		rpos_l2 = rpos_l[1].split(",")

		max1 = int(rpos_l1[0])
		if max1 < int(rpos_l1[1]):
			max1 = int(rpos_l1[1])

		max2 = int(rpos_l2[0])
		if max2 < int(rpos_l2[1]):
			max2 = int(rpos_l2[1])

		min1 = int(rpos_l1[0])
		if min1 > int(rpos_l1[1]):
			min1 = int(rpos_l1[1])

		min2 = int(rpos_l2[0])
		if min2 > int(rpos_l2[1]):
			min2 = int(rpos_l2[1])

		if max2 > max1:
#			print("rpos", max1)
			ins_length = int(ins_info_l[7]) - max1
#			print("rpos", max1, ins_length)
			return ins_length
		elif min2 < min1:
#			print("rpos", min1)
			return min1

#		print("??????????", min1, max1, min2, max2)

def get_read_name_in_both(r_BP_grp, l_BP_grp):
	read_found_in_both = {}
	for read_name in r_BP_grp:
		if read_name in l_BP_grp:
			read_found_in_both[read_name] = 1

	return read_found_in_both

def find_close_del(cigar, mapping_start, ins_pos, ins_size, close_del_distance, close_del_prop):
#	print("mapping_start", mapping_start)
#	print("ins_pos", ins_pos)
#	print("ins_size", ins_size)

#	print(cigar)

	cigar_type1 = re.split(r'\d+', cigar)
	cigar_len1 = re.split(r'[a-zA-Z]+', cigar)
	del cigar_type1[0]
	del cigar_len1[-1]

	del_pos = {}
	target_ins_found = 0
	map_end = mapping_start - 1
	for i in range(0, len(cigar_type1)):
#		print(">>", cigar_type1[i], cigar_len1[i])
		if cigar_type1[i] == "I" or cigar_type1[i] == "S" or cigar_type1[i] == "H":
#			if cigar_type1[i] == "I":
#				print(cigar_len1[i], cigar_type1[i], map_end, ins_pos, ins_size)
			if map_end + 1 == int(ins_pos) and int(cigar_len1[i]) == ins_size:
#				print("found !!!")
				target_ins_found = 1
#			if cigar_type1[i] == "I":
#				ins_pos[map_end] = cigar_len1[i]
			continue

		if cigar_type1[i] == "D":
#			print(cigar_len1[i], cigar_type1[i], map_end)
			del_pos[map_end] = cigar_len1[i]

#		if map_end + 1 == int(ins_pos) and int(cigar_len1[i]) == ins_size:
#			target_ins_found = 1

		map_end += int(cigar_len1[i])

	close_del_l = []
	for del_pos_tmp in del_pos:
#		print("del_pos_tmp", del_pos_tmp, ins_pos, ins_pos + ins_size)
#		if 0 <= ins_pos - del_pos_tmp <= close_del_distance or 0 <= del_pos_tmp - (ins_pos + ins_size) <= close_del_distance:
		if 0 <= ins_pos - (del_pos_tmp + int(del_pos[del_pos_tmp])) <= close_del_distance or 0 <= del_pos_tmp - ins_pos <= close_del_distance:
			close_del_l.append(int(del_pos[del_pos_tmp]))

#	print("close_del_l", close_del_l)

	max_close_del = 0
	if len(close_del_l) > 0:
		max_close_del = max(close_del_l)

#	print("max_close_del", max_close_del)

	if max_close_del >= close_del_prop*float(ins_size):
#		print("!!!!!!!!!!!!!")
		return max_close_del
	elif target_ins_found == 1:
		return 0
	else:
		return -1

SV_range = int(sys.argv[2])
ins_range_prop = float(sys.argv[3])
min_mq = int(sys.argv[4])
min_supprot_length = int(sys.argv[5])
min_num_of_read_with_min_support_length = int(sys.argv[6])
close_del_distance = int(sys.argv[7])
close_del_prop = float(sys.argv[8])

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

	if len(pre_BP_l) > 0 and abs(int(line_l[2]) - int(pre_BP_l[2])) > SV_range:
		l_BP_grp = {}
		r_BP_grp = {}
		r_read_name_with_min_supprot_length = {}
		l_read_name_with_min_supprot_length = {}
		num_l_read_with_min_supprot_length = 0
		num_r_read_with_min_supprot_length = 0
		INS_support_split_read_name = {}
		INS_breakpoints = []
		if len(BP_grp) >= 1:
#			print("###############")
			chr = ""

			INS_read_name = find_INS_read(BP_grp)

			BP_grp_info = get_INS_info(BP_grp)

#0 <><><><> ['split', 'chr11', '92329738', 'chr11', '92329739', 'INS', '11633', 'chr11/chr11', '32,654/4574,11616', '92330381,92330381/92329738,92329738', '20.3/20.3', '60/60',
#["within", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, SV_length, read_length, chr, str(r_start) + "," + str(r_end), str(g_start) + "," + str(g_end), bq, mq, strand, cigar]
			INS_PASS_read_name = {}
			for BP_info_read_name in BP_grp_info:
				for i in range(len(BP_grp_info[BP_info_read_name])):
					if BP_grp_info[BP_info_read_name][i][0] == "split" and BP_grp_info[BP_info_read_name][i][5] == "INS":
						mq_l = BP_grp_info[BP_info_read_name][i][12].split("/")
						if int(mq_l[0]) >= min_mq and int(mq_l[1]) >= min_mq:
							INS_PASS_read_name[BP_info_read_name] = 1
					elif BP_grp_info[BP_info_read_name][i][0] == "within" and BP_grp_info[BP_info_read_name][i][5] == "INS" and int(BP_grp_info[BP_info_read_name][i][12]) >= min_mq:
#						print("==========================")
						g_pos_l = BP_grp_info[BP_info_read_name][i][10].split(",")
						total_ins_len = find_close_del(BP_grp_info[BP_info_read_name][i][14], int(g_pos_l[0]), int(BP_grp_info[BP_info_read_name][i][2]), int(BP_grp_info[BP_info_read_name][i][6]), close_del_distance, close_del_prop)
						if total_ins_len == 0:
							INS_PASS_read_name[BP_info_read_name] = 1

			INS_BP_grp = []
			other_BP_grp = []
			for BP_grp_tmp in BP_grp:
				BP_grp_l = re.split("\t", BP_grp_tmp)
				if BP_grp_l[5] == "INS" and BP_grp_l[14] in INS_PASS_read_name:
					INS_BP_grp.append(BP_grp_tmp)
				chr = BP_grp_l[0]

			diff_BP_grp = INS_BP_grp
			diff_BP_num = 1
			if len(diff_BP_grp):
				while diff_BP_num:
					min_distance_read, median_length, same_BP_grp, diff_BP_grp = get_INS(diff_BP_grp, ins_range_prop)
					diff_BP_num = len(diff_BP_grp)
					INS_breakpoints.append([min_distance_read, median_length, same_BP_grp])
#BP_grp2 dic
#0 ['split', 'chr11', '92329739', 'chr9', '85288869', 'lCHR', 'NA', '5105', 'chr11/chr9', '31,1624/2627,3076', '92331394,92329739/85288869,85289341', '20.6/20.6', '60/2', '-

			for BP_info_read_name in BP_grp_info:
				if BP_info_read_name in INS_read_name:
					continue

				for i in range(len(BP_grp_info[BP_info_read_name])):
					mq_l = BP_grp_info[BP_info_read_name][i][12].split("/")
					read_pos_l = BP_grp_info[BP_info_read_name][i][9].split("/")
					read_pos_l1 = read_pos_l[0].split(",")
					if int(mq_l[0]) < min_mq:
						continue
					if BP_grp_info[BP_info_read_name][i][5][0] == "l":
						l_BP_grp.setdefault(BP_info_read_name, []).append(BP_grp_info[BP_info_read_name][i])
						INS_support_split_read_name[BP_info_read_name] = 1
						if abs(int(read_pos_l1[1]) - int(read_pos_l1[0])) >= min_supprot_length:
							l_read_name_with_min_supprot_length[BP_info_read_name] = 1
							num_l_read_with_min_supprot_length += 1
					if BP_grp_info[BP_info_read_name][i][5][0] == "r":
						r_BP_grp.setdefault(BP_info_read_name, []).append(BP_grp_info[BP_info_read_name][i])
						INS_support_split_read_name[BP_info_read_name] = 1
						if abs(int(read_pos_l1[1]) - int(read_pos_l1[0])) >= min_supprot_length:
							r_read_name_with_min_supprot_length[BP_info_read_name] = 1
							num_r_read_with_min_supprot_length += 1
			
#chr11   92329738        chr11   92329739        4       INS     3946.0  NA      chr11;92329734;chr11;92329739;

		read_found_both = get_read_name_in_both(r_BP_grp, l_BP_grp)
		
		for read_name_tmp in r_read_name_with_min_supprot_length:
			if len(r_read_name_with_min_supprot_length) and read_name_tmp in read_found_both:
				num_r_read_with_min_supprot_length -= 1

		for read_name_tmp in l_read_name_with_min_supprot_length:
			if len(l_read_name_with_min_supprot_length) and read_name_tmp in read_found_both:
				num_l_read_with_min_supprot_length -= 1

		r_BP_grp2 = {}
		for read_name_tmp in r_BP_grp:
			if not read_name_tmp in read_found_both:
				r_BP_grp2[read_name_tmp] = r_BP_grp[read_name_tmp]

		l_BP_grp2 = {}
		for read_name_tmp in l_BP_grp:
			if not read_name_tmp in read_found_both:
				l_BP_grp2[read_name_tmp] = l_BP_grp[read_name_tmp]

		r_BP_grp = r_BP_grp2
		l_BP_grp = l_BP_grp2

#		print("INS_breakpoints", INS_breakpoints)

		if len(INS_breakpoints) == 0:
#			print("r", num_r_read_with_min_supprot_length, "l", num_l_read_with_min_supprot_length)
			if num_r_read_with_min_supprot_length >= min_num_of_read_with_min_support_length and num_l_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:

				read_name_l = []
				BP = []
				chr = ""
				INS_len = []
				for BP_info_read_name in l_BP_grp:
#					print(">>>>", l_BP_grp[BP_info_read_name])
					BP.append(int(l_BP_grp[BP_info_read_name][0][2]))
					for BP_grp_tmp in BP_grp:
						BP_grp_l = re.split("\t", BP_grp_tmp)
						if BP_info_read_name == BP_grp_l[14] and BP_grp_l[5][0] == "l":
							read_name_l.append(BP_grp_l[0])
							ins_len = get_INS_length(l_BP_grp[BP_info_read_name][0])
							INS_len.append(ins_len)
					chr = l_BP_grp[BP_info_read_name][0][1]

				for BP_info_read_name in r_BP_grp:
#					print(">>>>", r_BP_grp[BP_info_read_name])
					BP.append(int(r_BP_grp[BP_info_read_name][0][2]))
					for BP_grp_tmp in BP_grp:
						BP_grp_r = re.split("\t", BP_grp_tmp)
						if BP_info_read_name == BP_grp_r[14] and BP_grp_r[5][0] == "r":
							read_name_l.append(BP_grp_r[0])
							ins_len = get_INS_length(r_BP_grp[BP_info_read_name][0])
							INS_len.append(ins_len)
					chr = r_BP_grp[BP_info_read_name][0][1]

				max_INS_len = ">" + str(max(INS_len))

				median_pos = int(np.median(BP))
				
				num_support_read = len(l_BP_grp) + len(r_BP_grp)
				INS_len = 0
				BP_num = "0," + str(len(l_BP_grp)) + "," + str(len(r_BP_grp))
				
				read_name_all = ":".join(read_name_l)

				print (chr, median_pos, chr, median_pos, num_support_read, "INS", max_INS_len, BP_num, read_name_all, sep="\t") 

		for breakpoint_tmp in INS_breakpoints:
			consensus_seq = ""
			read_name_all = ""
			num_match_read = len(breakpoint_tmp[2])
			if len(breakpoint_tmp[2]) >= 2:
				read_name_l = []
				for tmp in breakpoint_tmp[2]:
					tmp_l = re.split("\t", tmp)
					read_name_l.append(tmp_l[0])

				min_distance_read_l = re.split("\t", breakpoint_tmp[0])
#chr11   92329738        chr11   92329739        4       INS     3946.0  NA      chr11;92329734;chr11;92329739;

				num_l_BP = 0
				num_r_BP = 0
				if num_l_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
					num_l_BP = len(l_BP_grp)
					for BP_grp_tmp in BP_grp:
						BP_grp_l = re.split("\t", BP_grp_tmp)
						if BP_grp_l[5][0] == "l" and BP_grp_l[14] in INS_support_split_read_name:
							read_name_l.append(BP_grp_l[0])
					
				if num_r_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
					num_r_BP = len(r_BP_grp)
					for BP_grp_tmp in BP_grp:
						BP_grp_l = re.split("\t", BP_grp_tmp)
						if BP_grp_l[5][0] == "r" and BP_grp_l[14] in INS_support_split_read_name:
							read_name_l.append(BP_grp_l[0])

				read_name_all = ":".join(read_name_l)
				BP_num = str(len(breakpoint_tmp[2])) + "," +  str(num_l_BP) + "," + str(num_r_BP)
				print (min_distance_read_l[1], min_distance_read_l[2], min_distance_read_l[3], min_distance_read_l[4], str(len(breakpoint_tmp[2]) + num_l_BP + num_r_BP), "INS", str(breakpoint_tmp[1]), BP_num, read_name_all, sep="\t")
			else:
				min_distance_read_l = re.split("\t", breakpoint_tmp[0])

				num_l_BP = 0
				num_r_BP = 0
				read_name_l = []
				for tmp in breakpoint_tmp[2]:
					tmp_l = re.split("\t", tmp)
					read_name_l.append(tmp_l[0])

				if num_l_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
					num_l_BP = len(l_BP_grp)
					for BP_grp_tmp in BP_grp:
						BP_grp_l = re.split("\t", BP_grp_tmp)
						if BP_grp_l[5][0] == "l" and BP_grp_l[14] in INS_support_split_read_name:
							read_name_l.append(BP_grp_l[0])

				if num_r_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
					num_r_BP = len(r_BP_grp)
					for BP_grp_tmp in BP_grp:
						BP_grp_l = re.split("\t", BP_grp_tmp)
						if BP_grp_l[5][0] == "r" and BP_grp_l[14] in INS_support_split_read_name:
							read_name_l.append(BP_grp_l[0])

				read_name_all = ":".join(read_name_l)
				BP_num = str(len(breakpoint_tmp[2])) + "," +  str(num_l_BP) + "," + str(num_r_BP)
				print (min_distance_read_l[1], min_distance_read_l[2], min_distance_read_l[3], min_distance_read_l[4], str(len(breakpoint_tmp[2]) + num_l_BP + num_r_BP), "INS", str(breakpoint_tmp[1]), BP_num, read_name_all, sep="\t") #min_distance_read_l[0] => read_name_all

		BP_grp = []

	BP_grp.append(line)
	pre_BP_l = line_l

INS_breakpoints = []
r_BP_grp = {}
l_BP_grp = {}
INS_support_split_read_name = {}
l_read_name_with_min_supprot_length = {}
r_read_name_with_min_supprot_length = {}
num_l_read_with_min_supprot_length = 0
num_r_read_with_min_supprot_length = 0
if len(BP_grp) >= 1:
	chr = ""

#			print(len(BP_grp))
	INS_read_name = find_INS_read(BP_grp)

	BP_grp_info = get_INS_info(BP_grp)

#0 <><><><> ['split', 'chr11', '92329738', 'chr11', '92329739', 'INS', '11633', 'chr11/chr11', '32,654/4574,11616', '92330381,92330381/92329738,92329738', '20.3/20.3', '60/60',
#["within", bp_chr1, bp_pos1, bp_chr2, bp_pos2, type, read_length, chr, str(r_start) + "," + str(r_end), str(g_start) + "," + str(g_end), bq, mq, strand, cigar]
	INS_PASS_read_name = {}
	for BP_info_read_name in BP_grp_info:
#		print("BP_info_read_name => ", BP_info_read_name)
		for i in range(len(BP_grp_info[BP_info_read_name])):
#			print(i, "<><><><>", BP_grp_info[BP_info_read_name][i])
			if BP_grp_info[BP_info_read_name][i][0] == "split" and BP_grp_info[BP_info_read_name][i][5] == "INS":
				mq_l = BP_grp_info[BP_info_read_name][i][12].split("/")
				if int(mq_l[0]) >= min_mq and int(mq_l[1]) >= min_mq:
#							print("PASS => ", BP_info_read_name)
					INS_PASS_read_name[BP_info_read_name] = 1
			elif BP_grp_info[BP_info_read_name][i][0] == "within" and BP_grp_info[BP_info_read_name][i][5] == "INS" and int(BP_grp_info[BP_info_read_name][i][12]) >= min_mq:
#				print("AAAAAAAAAAAA")
				g_pos_l = BP_grp_info[BP_info_read_name][i][10].split(",")
				total_ins_len = find_close_del(BP_grp_info[BP_info_read_name][i][14], int(g_pos_l[0]), int(BP_grp_info[BP_info_read_name][i][2]), int(BP_grp_info[BP_info_read_name][i][6]), close_del_distance, close_del_prop)
#				print("total_ins_len = ", total_ins_len)
				if total_ins_len == 0:
					INS_PASS_read_name[BP_info_read_name] = 1
						
	INS_BP_grp = []
	other_BP_grp = []
	for BP_grp_tmp in BP_grp:
		BP_grp_l = re.split("\t", BP_grp_tmp)
#				print(BP_grp_l)
		if BP_grp_l[5] == "INS" and BP_grp_l[14] in INS_PASS_read_name:
#			print("BP_grp_l =>>", BP_grp_l)
			INS_BP_grp.append(BP_grp_tmp)
		chr = BP_grp_l[0]

	diff_BP_grp = INS_BP_grp
	diff_BP_num = 1
#			INS_breakpoints = []
	if len(diff_BP_grp):
		while diff_BP_num:
#			print("1 diff_BP_num = ", len(diff_BP_grp))
			min_distance_read, median_length, same_BP_grp, diff_BP_grp = get_INS(diff_BP_grp, ins_range_prop)
			diff_BP_num = len(diff_BP_grp)
#			print("2 diff_BP_num = ", diff_BP_num)
			INS_breakpoints.append([min_distance_read, median_length, same_BP_grp])
#BP_grp2 dic
#0 ['split', 'chr11', '92329739', 'chr9', '85288869', 'lCHR', 'NA', '5105', 'chr11/chr9', '31,1624/2627,3076', '92331394,92329739/85288869,85289341', '20.6/20.6', '60/2', '-
#			print("INS_breakpoints => ", len(INS_breakpoints), INS_breakpoints)

	for BP_info_read_name in BP_grp_info:
#		print(BP_info_read_name)
		if BP_info_read_name in INS_read_name:
			continue

		for i in range(len(BP_grp_info[BP_info_read_name])):
#			print(i, BP_grp_info[BP_info_read_name][i])
			mq_l = BP_grp_info[BP_info_read_name][i][12].split("/")
			read_pos_l = BP_grp_info[BP_info_read_name][i][9].split("/")
			read_pos_l1 = read_pos_l[0].split(",")
#			print("mq_l[0]", mq_l[0], min_mq)
			if int(mq_l[0]) < min_mq:
				continue
			if BP_grp_info[BP_info_read_name][i][5][0] == "l":
				l_BP_grp.setdefault(BP_info_read_name, []).append(BP_grp_info[BP_info_read_name][i])
#				print("#", read_pos_l1)
				INS_support_split_read_name[BP_info_read_name] = 1
				if abs(int(read_pos_l1[1]) - int(read_pos_l1[0])) >= min_supprot_length:
					l_read_name_with_min_supprot_length[BP_info_read_name] = 1
					num_l_read_with_min_supprot_length += 1
			if BP_grp_info[BP_info_read_name][i][5][0] == "r":
				r_BP_grp.setdefault(BP_info_read_name, []).append(BP_grp_info[BP_info_read_name][i])
#				print("#", read_pos_l1)
				INS_support_split_read_name[BP_info_read_name] = 1
				if abs(int(read_pos_l1[1]) - int(read_pos_l1[0])) >= min_supprot_length:
					r_read_name_with_min_supprot_length[BP_info_read_name] = 1
					num_r_read_with_min_supprot_length += 1
			
#chr11   92329738        chr11   92329739        4       INS     3946.0  NA      chr11;92329734;chr11;92329739;

#		if len(INS_breakpoints) == 0:
#			if num_r_read_with_min_supprot_length >= min_num_of_read_with_min_support_length and num_l_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
#				print("Candidate !!!")

	read_found_both = get_read_name_in_both(r_BP_grp, l_BP_grp)
#	print("read_found_in_both => ", read_found_both)

#		print("r_read_name_with_min_supprot_length => ", r_read_name_with_min_supprot_length)
		
	for read_name_tmp in r_read_name_with_min_supprot_length:
		if len(r_read_name_with_min_supprot_length) and read_name_tmp in read_found_both:
			num_r_read_with_min_supprot_length -= 1

	for read_name_tmp in l_read_name_with_min_supprot_length:
		if len(l_read_name_with_min_supprot_length) and read_name_tmp in read_found_both:
			num_l_read_with_min_supprot_length -= 1

#		print("r_BP_grp => ", r_BP_grp)

	r_BP_grp2 = {}
	for read_name_tmp in r_BP_grp:
		if not read_name_tmp in read_found_both:
			r_BP_grp2[read_name_tmp] = r_BP_grp[read_name_tmp]

	l_BP_grp2 = {}
	for read_name_tmp in l_BP_grp:
		if not read_name_tmp in read_found_both:
			l_BP_grp2[read_name_tmp] = l_BP_grp[read_name_tmp]

	r_BP_grp = r_BP_grp2
	l_BP_grp = l_BP_grp2

#		print(num_r_read_with_min_supprot_length, num_l_read_with_min_supprot_length)

	if len(INS_breakpoints) == 0:
		if num_r_read_with_min_supprot_length >= min_num_of_read_with_min_support_length and num_l_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
#			print("l_BP_grp => ", len(l_BP_grp), l_BP_grp)
#			print("r_BP_grp => ", len(r_BP_grp), r_BP_grp)

			read_name_l = []
			BP = []
			chr = ""
			INS_len = []
			for BP_info_read_name in l_BP_grp:
#				print(l_BP_grp[BP_info_read_name])
				BP.append(int(l_BP_grp[BP_info_read_name][0][2]))
				for BP_grp_tmp in BP_grp:
					BP_grp_l = re.split("\t", BP_grp_tmp)
					if BP_info_read_name == BP_grp_l[14] and BP_grp_l[5][0] == "l":
						read_name_l.append(BP_grp_l[0])
						ins_len = get_INS_length(l_BP_grp[BP_info_read_name][0])
						INS_len.append(ins_len)
				chr = l_BP_grp[BP_info_read_name][0][1]

			for BP_info_read_name in r_BP_grp:
				BP.append(int(r_BP_grp[BP_info_read_name][0][2]))
				for BP_grp_tmp in BP_grp:
					BP_grp_r = re.split("\t", BP_grp_tmp)
					if BP_info_read_name == BP_grp_r[14] and BP_grp_r[5][0] == "r":
						read_name_l.append(BP_grp_r[0])
						ins_len = get_INS_length(r_BP_grp[BP_info_read_name][0])
						INS_len.append(ins_len)
#				print("INS_len => ", INS_len)
			max_INS_len = ">" + str(max(INS_len))

#				print(BP)
			median_pos = int(np.median(BP))
#				print(median_pos)
				
			num_support_read = len(l_BP_grp) + len(r_BP_grp)
			INS_len = 0
			BP_num = "0," + str(len(l_BP_grp)) + "," + str(len(r_BP_grp))
				
			read_name_all = ":".join(read_name_l)

			print (chr, median_pos, chr, median_pos, num_support_read, "INS", max_INS_len, BP_num, read_name_all, sep="\t") 

#			print("#################")
	for breakpoint_tmp in INS_breakpoints:
		consensus_seq = ""
		read_name_all = ""
		num_match_read = len(breakpoint_tmp[2])
		if len(breakpoint_tmp[2]) >= 2:
			read_name_l = []
			for tmp in breakpoint_tmp[2]:
				tmp_l = re.split("\t", tmp)
				read_name_l.append(tmp_l[0])
				min_distance_read_l = re.split("\t", breakpoint_tmp[0])
#chr11   92329738        chr11   92329739        4       INS     3946.0  NA      chr11;92329734;chr11;92329739;

			num_l_BP = 0
			num_r_BP = 0
			if num_l_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
				num_l_BP = len(l_BP_grp)
				for BP_grp_tmp in BP_grp:
					BP_grp_l = re.split("\t", BP_grp_tmp)
					if BP_grp_l[5][0] == "l" and BP_grp_l[14] in INS_support_split_read_name:
						read_name_l.append(BP_grp_l[0])
					
			if num_r_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
				num_r_BP = len(r_BP_grp)
				for BP_grp_tmp in BP_grp:
					BP_grp_l = re.split("\t", BP_grp_tmp)
					if BP_grp_l[5][0] == "r" and BP_grp_l[14] in INS_support_split_read_name:
						read_name_l.append(BP_grp_l[0])

			read_name_all = ":".join(read_name_l)
			BP_num = str(len(breakpoint_tmp[2])) + "," +  str(num_l_BP) + "," + str(num_r_BP)
			print (min_distance_read_l[1], min_distance_read_l[2], min_distance_read_l[3], min_distance_read_l[4], str(len(breakpoint_tmp[2]) + num_l_BP + num_r_BP), "INS", str(breakpoint_tmp[1]), BP_num, read_name_all, sep="\t") #min_distance_read_l[0] => read_name_all
		else:
			min_distance_read_l = re.split("\t", breakpoint_tmp[0])

			num_l_BP = 0
			num_r_BP = 0
			read_name_l = []
			for tmp in breakpoint_tmp[2]:
				tmp_l = re.split("\t", tmp)
				read_name_l.append(tmp_l[0])

			if num_l_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
				num_l_BP = len(l_BP_grp)
				for BP_grp_tmp in BP_grp:
					BP_grp_l = re.split("\t", BP_grp_tmp)

					if BP_grp_l[16] in read_found_both:
						read_found_both

					if BP_grp_l[5][0] == "l" and BP_grp_l[14] in INS_support_split_read_name:
#                                               print("l", BP_grp_l[5], BP_grp_l[0])
						read_name_l.append(BP_grp_l[0])

			if num_r_read_with_min_supprot_length >= min_num_of_read_with_min_support_length:
				num_r_BP = len(r_BP_grp)
				for BP_grp_tmp in BP_grp:
					BP_grp_l = re.split("\t", BP_grp_tmp)
					if BP_grp_l[5][0] == "r" and BP_grp_l[14] in INS_support_split_read_name:
#                                               print("r", BP_grp_l[5], BP_grp_l[0])
						read_name_l.append(BP_grp_l[0])

			read_name_all = ":".join(read_name_l)
			BP_num = str(len(breakpoint_tmp[2])) + "," +  str(num_l_BP) + "," + str(num_r_BP)
			print (min_distance_read_l[1], min_distance_read_l[2], min_distance_read_l[3], min_distance_read_l[4], str(len(breakpoint_tmp[2]) + num_l_BP + num_r_BP), "INS", str(breakpoint_tmp[1]), BP_num, read_name_all, sep="\t") #min_distance_read_l[0] => read_name_all

