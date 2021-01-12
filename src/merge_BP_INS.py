import re 
import sys
import subprocess
#from Bio.Align.Applications import MafftCommandline
from statistics import median

#chr1   63980559        chr1    63982331        FF
#chr1   63980566        chr1    63982331        FF

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
                        BP_grp2_2 = []

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

#chr2    16148441        chr2    16148441        5       INS     199.0   5,0,0   chr2;16148441;chr2;16148441;GTATTGGATGAATGTTATGTTGTATGTATTATATTATATTATATCATGTGTATACAAGAATAATGT
def define_BP_with_ovrelap(read_l0, overlap_prop):
	new_breakpoints = []

	read_l = []
	for i in range(len(read_l0)):
		read_l0_l = read_l0[i].split("\t")
		read_l.append(read_l0_l)

	read_info = {}
	for i in range(len(read_l)):
		read_info[i] = read_l[i]

	read_grp = {}
	read_grp_num = 0
	for i in range(len(read_l) - 1):
		for j in range(i + 1, len(read_l)):
#			print(read_l[i])
#			print(read_l[j])
			length1 = float(read_l[i][6].replace(">", ""))
			length2 = float(read_l[j][6].replace(">", ""))
			del_OL_length = min(length1, length2)/max(length1, length2)
#			print("del_OL_length", del_OL_length)
			if del_OL_length > 0 and del_OL_length >= overlap_prop and del_OL_length >= overlap_prop:
				found = 0
				for read_grp_num_tmp in read_grp:
					if i in read_grp[read_grp_num_tmp]:
						read_grp[read_grp_num_tmp].append(j)
						read_grp[read_grp_num_tmp] = list(set(read_grp[read_grp_num_tmp]))
						found = 1
					if j in read_grp[read_grp_num_tmp]:
						read_grp[read_grp_num_tmp].append(i)
						read_grp[read_grp_num_tmp] = list(set(read_grp[read_grp_num_tmp]))
						found = 1

				if found == 0:
					read_grp_num += 1
					read_grp[read_grp_num] = [i, j]
#	print("read_grp => ", read_grp)
#chr15;19778987;chr15;19779158;CTATAAAACGAAGCATTCTCA;171;260e3fb4-5235-44a8-99cc-e7593904140e;1456,10082/10048,13826;chr15/chr15;13844;18820120,18828991/19775274,19
#chr15;19778989;chr15;19779162;GTAGATTCTACAAAGAGTGTT;173;7031bd75-a306-432a-b3c1-36e4d8eaa9a9;4828,9876;chr15;9894;19775265,19780649;+;17;0;4827S8M2I18M1I64M1I4M1D1
#['19776992\t19777170\t1\tchr15;19776992;chr15;19777170;CGGCGTTTTCAAAGGCTAAGG;178;de24a9f0-ac13-4797-94a1-2c974fc3e3a6;5373,15967/2197,6781/18854,19	
	delete_read_grp = {}
	for read_grp_num in read_grp:
		for read_grp_num2 in read_grp:
			if read_grp_num <= read_grp_num2:
				continue
			set1 = set(read_grp[read_grp_num])
			set2 = set(read_grp[read_grp_num2])
			matched_list = list(set1 & set2)
			
			if len(matched_list):
				read_grp[read_grp_num].extend(read_grp[read_grp_num2])
				read_grp[read_grp_num] = list(set(read_grp[read_grp_num]))
				delete_read_grp[read_grp_num2] = 1

	for delete_read_grp_num in delete_read_grp:
		del read_grp[delete_read_grp_num]

#chr2    16148441        chr2    16148441        9       INS     199.0   5,0,0   chr2;16147641;chr2;16147641;TGTGTGCTGTCTCTTTTATGGAATGGGCAGAGGCTGGTTTCTCCTAAGAAGGTGACATGTTGAACA
	for read_grp_num in read_grp:
		start, end = 0, 0
		read_info_in_read_grp = []
		output_read_grp_info = []
		total_support_read = 0
		total_both_read, total_left_read, total_right_read = 0, 0, 0
		read_length_l = []
		read_length_l2 = []
		for read_num in read_grp[read_grp_num]:
			output_read_grp_info = read_info[read_num][0:8]
			total_support_read += int(read_info[read_num][4])
			
			read_num_l = read_info[read_num][7].split(",")
			total_both_read += int(read_num_l[0])
			total_left_read += int(read_num_l[1])
			total_right_read += int(read_num_l[2])
#			read_length_l.append(read_info[read_num][4])
			read_info_in_read_grp.append(read_info[read_num][8])

			if ">" in read_info[read_num][6]:
				length_tmp = read_info[read_num][6].replace(">", "")
				read_length_l2.append(float(length_tmp))
#				print(length_tmp, read_length_l2)
			else:
				read_length_l.append(float(read_info[read_num][6]))				
#				print(read_info[read_num][6], read_length_l)
				

			del read_info[read_num]

		output_read_grp_info[4] = str(total_support_read)
		output_read_grp_info[7] = str(total_both_read) + "," + str(total_left_read) + "," + str(total_right_read)

		new_del_candidate = ["\t".join(output_read_grp_info), ":".join(read_info_in_read_grp)]
		new_breakpoints.append("\t".join(new_del_candidate))

#		print("read_length_l2", read_length_l2)
#		print("read_length_l", read_length_l)		

		ins_length = 0
		if len(read_length_l2):
			ins_length = ">" + str(max(read_length_l2))
		elif len(read_length_l):
			ins_length = str(median(read_length_l))

		output_read_grp_info[6] = str(ins_length)

	for num in read_info:
		output_read_grp_info = read_info[num]
		output_read_grp_info[4] = read_info[num][4]
		new_del_candidate = ["\t".join(output_read_grp_info), read_info[num][8]]
		new_breakpoints.append("\t".join(new_del_candidate))

#	print("new_breakpoints", len(new_breakpoints), new_breakpoints)

	return new_breakpoints

def get_BP(BP_grp, overlap_prop):
        BP_grp_dic = {}
        for BP_tmp in BP_grp:
                BP_tmp_l = re.split("\t", BP_tmp)
                BP_grp_dic.setdefault(BP_tmp, {})["chr1"] = BP_tmp_l[1]
                BP_grp_dic.setdefault(BP_tmp, {})["bp1"] = int(BP_tmp_l[2])
                BP_grp_dic.setdefault(BP_tmp, {})["chr2"] = BP_tmp_l[3]
                BP_grp_dic.setdefault(BP_tmp, {})["bp2"] = int(BP_tmp_l[4])
        pre_BP_tmp = []
        BP_grp2 = []

        breakpoints_all = []
        for BP_tmp_t in sorted(BP_grp_dic.items(), key=lambda x:x[1]['bp2']):
                BP_tmp = re.split("\t", BP_tmp_t[0])
                if len(pre_BP_tmp) > 0 and abs(int(BP_tmp[4]) - int(pre_BP_tmp[4])) > SV_range:
                        breakpoints = get_BP2(BP_grp2)
                        breakpoints_all.extend(breakpoints)
                        BP_grp2 = []

                BP_grp2.append("\t".join(BP_tmp))
                pre_BP_tmp = BP_tmp

        breakpoints = get_BP2(BP_grp2)
        breakpoints_all.extend(breakpoints)

        breakpoints_all_new = define_BP_with_ovrelap(breakpoints_all, overlap_prop)

        return breakpoints_all_new

SV_range = int(sys.argv[2])
overlap_prop = float(sys.argv[3])

#chr1    143187594       chr1    143188172       2       DEL     578     NA      chr1;143187590;chr1;143188210;GAATCATCGATGAATCGAATG;620;ba610b60-9235-42dc-a61a-9c71ee035cdc;25,9
#chr1    143187239       chr1    143212767       1       DEL     25528   0       chr1;143187239;chr1;143212767;GATGGATCGTCTCATGAATGGA;25528;2139beb2-056e-4a09-b049-a3bc917a8272;2
#chr1    143187377       chr1    143212784       2       DEL     25407   NA      chr1;143187239;chr1;143212767;GATGGATCGTCTCATGAATGGA;25528;2139beb2-056e-4a09-b049-a3bc917a8272;2

pre_BP_l = ""
BP_grp = []
SV_BP_f = open(sys.argv[1])
for line in SV_BP_f:
        line = line.replace('\n', '')
        line_l = re.split('\t', line)

        if line[0] == "#":
                print(line)
                continue

        if len(pre_BP_l) > 0 and abs(int(line_l[1]) - int(pre_BP_l[1])) > SV_range:
#                print("BP_grp => ", len(BP_grp), BP_grp)
                if len(BP_grp) == 1:
                        print(BP_grp[0])
                if len(BP_grp) > 1:
                        newbreakpoints = define_BP_with_ovrelap(BP_grp, overlap_prop)
                        for BP in newbreakpoints:                
                                print(BP)

                BP_grp = []

        BP_grp.append(line)
        pre_BP_l = line_l

if len(BP_grp) >= 1:
	if len(BP_grp) == 1:
		print(BP_grp[0])
	if len(BP_grp) > 1:
		newbreakpoints = define_BP_with_ovrelap(BP_grp, overlap_prop)
		for BP in newbreakpoints:
			print(BP)
