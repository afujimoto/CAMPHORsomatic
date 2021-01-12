import sys 
import re
from operator import itemgetter

blat_f = open(sys.argv[1])
minimum_length_prop = float(sys.argv[2])
min_del_dis = int(sys.argv[3])
min_mq = int(sys.argv[4])
realign_range = int(sys.argv[5])
min_ins_ref_read_prop = float(sys.argv[6])
min_unmapped_region_length = int(sys.argv[7])
min_ins_range = int(sys.argv[8])

def get_flanking_seq(start1, end1, start2, end2, seq):
	flanking_seq = seq[start1-1:end1] + seq[start2-1:end2]
	return flanking_seq

def rev_con(seq):
	seq = seq.upper()
	seq2 = ""
	for n in seq:
		if n == "A":
			seq2 += "T"
		elif n == "T":
			seq2 += "A"
		elif n == "G":
			seq2 += "C"
		elif n == "C":
			seq2 += "C"

	return seq2[::-1]

#a37d4eb0-6e52-46a1-b1fc-5e7bb041b4f4_Basecall_1D_template/9.1   740     454,663;81,197;713,738;56,76;278,316;670,687;36,55;257,281;355,444;694,710;203,223;318,334;225,243;244,255      chr5;chr4;chr12;chr13;chr5;chr5;chr4;chr3;chr2;chr20;chr1;chr9;chr19;chr22      55702433,55702670;104928643,105050027;28108444,28108478;58749091,58749111;92781344,92781734;91028880,91028897;30471053,30471072;70056626,70056651;132342071,132342217;58743330,58743346;191223210,191223231;29994460,29994476;28850379,28850397;36832582,36832593       -;+;+;+;+;-;+;-;-;-;+;+;+;-
#2e36a1af-66ef-4fbc-9228-ceeb9750838b_Basecall_1D_template       2736    523,2724;40,540 chr1;chr1       150719288,150721717;150718459,150718973 +;+     15.35;16.26
#     60;60   522S8M1D19M1I5M1I21M1D44M1D5M3D10M1D18M3I19M1D7M1I9M2D16M1D5M1I31M5D19M3D32M3D31M1D6M1D6M1D16M2D27M1D9M1D22M1D20M1I24M1I5M8D27M2D4M1D10M1D3M1D7M2D28M1I6M1I15M1I2M3D14M1D9M1D7M1D4M1D3M1D5M2D5M2D10M1D15M1D10M1D17M1D3M2D9M1D7M2D2M1D29M1D54M4D6M1I10M1D9M1D10M2I8M2D5M2D6M2D6M1D5M1D13M1D12M1D29M4D5M1D11M4I38M1D13M1D50M1D4M4D2M5D8M1D21M1I22M1D7M1D10M1D26M2D1M2D21M4D13M4D29M4D2M1D4M1D6M2D9M2D16M1D7M1D17M1D14M1D37M1D3M1D5M1D8M1D16M2D6M3D6M2D16M5D10M1D9M4D13M4D6M2D12M1D7M2D8M1D12M1D1M1D12M2D4M1I10M1D1M1D5M1D14M1I27M1D6M1D14M1I6M1D3M2D10M1I21M1D36M4D3M5D21M2D27M1D15M1D8M17D12M1D1M1D9M3D10M2D19M2D2M1D10M1D14M3D6M4D29M3D2M1D2M3D14M1D7M1D11M2D18M1D17M1D4M3D22M1I22M2D5M1D8M1D4M1D1M1D1M1D49M2I15M1D27M2D27M1D56M1D68M3D51M12S;39H11M2D5M1D48M2D18M1D14M1D35M1D10M2D7M1I19M2D14M1D12M1D46M3D58M1I18M1D10M3I15M1I28M1D10M1D47M1I8M1D16M1D33M2I10M2196H    MD:Z:2A1A3^A10A0A0T1C4T25^A6G1A8G26^C5^CTT10^A0T15A20^A7C8^AA16^A3C6A6A18^TAAAG2T13A2^TTG0T1A29^TTT31^A6^A6^T14T1^TG27^A9^A22^T13C0T34^AATGGGTT27^TA4^G10^C3^T7^CC28G1A0C3T15^AAA0A0A0T11^C9^T7^T4^T3^T5^TT5^TT10^G15^G0C0C8^G17^A3^CA2G6^C1G5^CG2^T6G13A8^A46T0T0T0T0T3^GACG0G15^C9^A18^GG1T3^CG6^CC6^G5^A13^C4G7^A27C1^CCCA5^T15A0G2A2G26^C2G0A0C0C7^C4C0G0C0C0T1G1C16T10A0A0G0C6^C4^CCGC2^ATTTT8^T11A4G2T1G0C16G1C1^C7^C10^A1C3C0C16A2^AA1^TT21^CCGG2C0C0A1T6^TTTC0T5A1A20^ACTA2^G4^A3T0T1^AT9^TA1T8A0C1T2^A0A0A0T4^C17^C1C0C0G2C7^C37^T0A2^T5^C1G6^T14C1^TA4A1^TTA6^GG16^AAGAG10^A9^CACA13^ATAC6^AT12^A7^GG8^C12^G1^G4A0C0G5^GC12C1^A1^G0C4^G9T1A0A8C4G0G13^T6^C14A5^C3^CC0A10A1G1C2A10C1^G2A11A0G0C0G2A6A0G5G1^AAAA3^AGACA3G0G16^TG27^T15^T8^TTTTTTTTTTTTTTTTT6A1A3^A1^T9^GTC10^AA19^GC2^A10^T14^AGC6^CCTG5G23^GGC0G1^C2^CCA9G4^T4G2^T11^AC1G16^C0C16^C4^CCT20G23^GG0C4^G8^C4^C1^G1^C64^G27^CT27^G13G42^G7A41A18^ATC51;MD:Z:11^AA5^T48^AT18^A14^G35^A10^TG5C20^TT12G1^A12^C9A26G0C1A6^CCT9A37C11G6T0G0C1C1G3^C0A52^A1A8^G0C52G1^A0A15^G43

for line in blat_f:
	line = line.replace('\n', '')
	line_l = re.split('\t', line)

	bq = (re.split(";", line_l[6]))[0]

	chr_l = re.split(";", line_l[3])
	gpos_l = re.split(";", line_l[4])
	align_region = re.split(";", line_l[2])
	total_len = int(line_l[1])
	strand_l = re.split(";", line_l[5])
	cigar_l = re.split(";", line_l[8])
	mq_l = re.split(";", line_l[7])
	seq = line_l[9]
	
	seq_align_pos = []
	for i in range(0, len(align_region)):
		tmp_l = re.split(",", align_region[i])
		if float((int(tmp_l[1]) - int(tmp_l[0]))/total_len) < minimum_length_prop:
			continue
		if int(mq_l[i]) < min_mq:
			continue
		gpos_tmp_l = re.split(",", gpos_l[i])
#		seq_align_pos.append([int(tmp_l[0]), int(tmp_l[1]), int(total_len), chr_l[i], int(gpos_tmp_l[0]), int(gpos_tmp_l[1]), strand_l[i], float(bq), line_l[0], mq_l[i]])
		seq_align_pos.append([int(tmp_l[0]), int(tmp_l[1]), int(total_len), chr_l[i], int(gpos_tmp_l[0]), int(gpos_tmp_l[1]), strand_l[i], bq, line_l[0], mq_l[i]])

	if len(seq_align_pos) == 0:
		continue

	seq_align_pos.sort(key=itemgetter(0))

	mapping_start_on_read = seq_align_pos[0][0]
	mapping_end_on_read = seq_align_pos[-1][1]

#	for i in range(len(seq_align_pos)):
#		print(i, seq_align_pos[i])

#indel_region_seq = seq[realign_start - 1:realign_end]
	if mapping_start_on_read > min_unmapped_region_length:
		if seq_align_pos[0][6] == "+":
			print("NA", "NA", seq_align_pos[0][3], str(seq_align_pos[0][4]), "BR", str(mapping_start_on_read), seq[0:mapping_start_on_read + 1], line, sep="\t")
		else:
			print(seq_align_pos[0][3], str(seq_align_pos[0][5]), "NA", "NA", "BR", str(mapping_start_on_read), seq[0:mapping_start_on_read + 1], line, sep="\t")
	if total_len - mapping_end_on_read > min_unmapped_region_length:
		if seq_align_pos[-1][6] == "+":
			print( seq_align_pos[-1][3], str(seq_align_pos[-1][5]), "NA", "NA", "BR", str(total_len - mapping_end_on_read), seq[total_len - mapping_end_on_read - 1:len(seq)], line, sep="\t")
		else:
			print("NA", "NA", seq_align_pos[-1][3], str(seq_align_pos[-1][4]), "BR", str(total_len - mapping_end_on_read), seq[total_len - mapping_end_on_read - 1:len(seq)], line, sep="\t")

#>> [8, 26, 990, 'chr21', '37884344,37884362', '-', '8.2']
#>> [60, 78, 990, 'chr3', '178867520,178867538', '-', '8.2']
	for i in range(0, len(cigar_l)):
		tmp_l = re.split(",", align_region[i])
		start_l = []
		end_l = []
		gpos_tmp_l = re.split(",", gpos_l[i])

		if int(mq_l[i]) < min_mq:
			continue

		cigar_type1 = re.split(r'\d+', cigar_l[i])
		cigar_len1 = re.split(r'[a-zA-Z]+', cigar_l[i])
		del cigar_type1[0]
		del cigar_len1[-1]

		indel_start_pos = int(gpos_tmp_l[0])
		read_pos = 0
		for j in range(0, len(cigar_type1)):
			if cigar_type1[j] == "S" or cigar_type1[j] == "H":
				if cigar_type1[j] == "S":
					read_pos += int(cigar_len1[j])
				continue
			if cigar_type1[j] == "M":
				indel_start_pos += int(cigar_len1[j])
				read_pos += int(cigar_len1[j])

			if cigar_type1[j] == "D":
				if int(cigar_len1[j]) >= min_del_dis:
					realign_start = read_pos - realign_range
					realign_end = read_pos + realign_range
					indel_region_seq = seq[realign_start - 1:realign_end]
					indel_end_pos = indel_start_pos + int(cigar_len1[j])
					print(chr_l[i], str(indel_start_pos), chr_l[i], str(indel_end_pos), "DEL", cigar_len1[j], indel_region_seq, line, sep="\t")
				indel_start_pos += int(cigar_len1[j])

			if cigar_type1[j] == "I":
				if int(cigar_len1[j]) >= min_del_dis:
#					realign_start = read_pos - realign_range
#					realign_end = read_pos + realign_range
					realign_start = read_pos
					realign_end = read_pos + int(cigar_len1[j]) + 1
					indel_region_seq = seq[realign_start - 1:realign_end]
					print(chr_l[i], str(indel_start_pos), chr_l[i], str(indel_start_pos), "INS", cigar_len1[j], indel_region_seq, line, sep="\t")
				read_pos += int(cigar_len1[j])

#[2, 1580, 6692, 'chr8', 44878726, 44880519, '+', 15.1, '000a7166-808a-43d5-aad6-905b8f29ba67_Basecall_1D_template', '2']
	for i in range(0, len(seq_align_pos)):
		ins_found = 0
		for j in range(i + 2, len(seq_align_pos)):
#			print(i, seq_align_pos[i], "vs", j, seq_align_pos[j])
			if seq_align_pos[i][3] == seq_align_pos[j][3] and seq_align_pos[i][6] == seq_align_pos[j][6]:
				if seq_align_pos[i][6] == "+" and abs(seq_align_pos[i][5] - seq_align_pos[j][4]) < min_ins_range and seq_align_pos[j][0] - seq_align_pos[i][1] - 1 >= min_del_dis:
					ins_found = 1
#					print("INS found 1")
#					print("INS", i, seq_align_pos[i], "vs", j, seq_align_pos[j])
					flanking_seq = seq[seq_align_pos[i][1] + 1:seq_align_pos[j][0]]
#					print(seq_align_pos[i][1] + 1, seq_align_pos[j][0])
					if seq_align_pos[i][5] < seq_align_pos[j][4]:
						[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][5], seq_align_pos[j][3], seq_align_pos[j][4], "INS", seq_align_pos[j][0] - seq_align_pos[i][1] - 1, flanking_seq]
					else:
						[chr1, breakpoint2, chr2, breakpoint1, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][5], seq_align_pos[j][3], seq_align_pos[j][4], "INS", seq_align_pos[j][0] - seq_align_pos[i][1] - 1, flanking_seq]
				elif seq_align_pos[i][6] == "-" and abs(seq_align_pos[i][4] - seq_align_pos[j][5]) < min_ins_range and seq_align_pos[j][0] - seq_align_pos[i][1] - 1 >= min_del_dis:
					ins_found = 1
#					print("INS found 2")
#					print("INS", i, seq_align_pos[i], "vs", j, seq_align_pos[j])
#					print(seq_align_pos[i][1] + 1, seq_align_pos[j][0])
					flanking_seq = seq[seq_align_pos[i][1] + 1:seq_align_pos[j][0]]
					if seq_align_pos[i][4] < seq_align_pos[j][5]:
						[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][4], seq_align_pos[j][3], seq_align_pos[j][5], "INS", seq_align_pos[j][0] - seq_align_pos[i][1] - 1, flanking_seq]
					else:
						[chr2, breakpoint2, chr1, breakpoint1, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][4], seq_align_pos[j][3], seq_align_pos[j][5], "INS", seq_align_pos[j][0] - seq_align_pos[i][1] - 1, flanking_seq]
			if ins_found:
				print(chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq, line, sep="\t")
				ins_found = 0

	for i in range(0, len(seq_align_pos) - 1):
#		print("seq_align_pos", i, seq_align_pos[i])
#		print("seq_align_pos", i + 1, seq_align_pos[i + 1])
		[chr1, breakpoint1, chr2, breakpoint2, SV_type] = ["", "", "", "", ""]
		if seq_align_pos[i][3] != seq_align_pos[i + 1][3]: # different chr
			flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][1], seq_align_pos[i][1] + realign_range, seq)

			BP1, BP2 = 0, 0
			if seq_align_pos[i][6] == "+":
				BP1 = seq_align_pos[i][5]
			else:
				BP1 = seq_align_pos[i][4]
			if seq_align_pos[i + 1][6] == "+":
				BP2 = seq_align_pos[i + 1][4]
			else:
				BP2 = seq_align_pos[i + 1][5]

			chr1 = seq_align_pos[i][3]
			chr2 = seq_align_pos[i + 1][3]
			chr1_2, chr2_2 = sorted([chr1, chr2])
			if chr1_2 == chr1:
				[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], BP1, seq_align_pos[i + 1][3], BP2, "CHR", "NA", flanking_seq]
			else:
				[chr2, breakpoint2, chr1, breakpoint1, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], BP1, seq_align_pos[i + 1][3], BP2, "CHR", "NA", flanking_seq]
		elif seq_align_pos[i][6] == seq_align_pos[i + 1][6]: # same strand
#			print(seq_align_pos[i])
#			print(seq_align_pos[i + 1])

			if (seq_align_pos[i][6] == "+" and seq_align_pos[i + 1][0] - seq_align_pos[i][1] < min_del_dis and abs(seq_align_pos[i + 1][4] - seq_align_pos[i][5]) < min_del_dis) or (seq_align_pos[i][6] == "-" and seq_align_pos[i + 1][0] - seq_align_pos[i][1] < min_del_dis and abs(seq_align_pos[i][4] - seq_align_pos[i + 1][5]) < min_del_dis): #180709
#				print("skip1", seq_align_pos[i])
#				print("skip2", seq_align_pos[i + 1])
				continue
#				a = 1 

			if (seq_align_pos[i][6] == "+" and seq_align_pos[i + 1][0] - seq_align_pos[i][1] >= min_del_dis and float(seq_align_pos[i + 1][0] - seq_align_pos[i][1])*min_ins_ref_read_prop > abs(seq_align_pos[i + 1][4] - seq_align_pos[i][5])): # + strand ins
#				print("RRRRRRRRRRRRRRRR", float(seq_align_pos[i + 1][0] - seq_align_pos[i][1])*min_ins_ref_read_prop)
				flanking_seq = seq[seq_align_pos[i][1] + 1:seq_align_pos[i+1][0]]
#				flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][0], seq_align_pos[i+1][0] + realign_range, seq)
				if seq_align_pos[i][5] < seq_align_pos[i + 1][4]:
					[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][5], seq_align_pos[i + 1][3], seq_align_pos[i + 1][4], "INS", seq_align_pos[i + 1][0] - seq_align_pos[i][1] - 1, flanking_seq]
				else:
					[chr1, breakpoint2, chr2, breakpoint1, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][5], seq_align_pos[i + 1][3], seq_align_pos[i + 1][4], "INS", seq_align_pos[i + 1][0] - seq_align_pos[i][1] - 1, flanking_seq]
			elif (seq_align_pos[i][6] == "-" and seq_align_pos[i + 1][0] - seq_align_pos[i][1] >= min_del_dis and float(seq_align_pos[i + 1][0] - seq_align_pos[i][1])*min_ins_ref_read_prop > abs(seq_align_pos[i][4] - seq_align_pos[i+1][5])): # - strand ins
#				print("BBBBBBBBBBBBBBBBBBBBB", float(seq_align_pos[i + 1][0] - seq_align_pos[i][1])*min_ins_ref_read_prop)
				flanking_seq = seq[seq_align_pos[i][1] + 1:seq_align_pos[i+1][0]]
#				flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][0], seq_align_pos[i+1][0] + realign_range, seq)
				if seq_align_pos[i][4] < seq_align_pos[i + 1][5]:
					[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][4], seq_align_pos[i + 1][3], seq_align_pos[i + 1][5], "INS", seq_align_pos[i + 1][0] - seq_align_pos[i][1] - 1, flanking_seq]
				else:
					[chr2, breakpoint2, chr1, breakpoint1, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][4], seq_align_pos[i + 1][3], seq_align_pos[i + 1][5], "INS", seq_align_pos[i + 1][0] - seq_align_pos[i][1] - 1, flanking_seq]
			elif (seq_align_pos[i][6] == "+" and seq_align_pos[i + 1][4] - seq_align_pos[i][5] >= min_del_dis): # + strand del
				flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][0], seq_align_pos[i+1][0] + realign_range, seq)
				[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][5], seq_align_pos[i + 1][3], seq_align_pos[i + 1][4], "DEL", int(seq_align_pos[i + 1][4]) - int(seq_align_pos[i][5]), flanking_seq]
			elif (seq_align_pos[i][6] == "-" and seq_align_pos[i][4] - seq_align_pos[i+1][5] >= min_del_dis): # - strand del
				flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][0], seq_align_pos[i+1][0] + realign_range, seq)
				[chr2, breakpoint2, chr1, breakpoint1, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][4], seq_align_pos[i + 1][3], seq_align_pos[i + 1][5], "DEL", int(seq_align_pos[i][4]) - int(seq_align_pos[i+1][5]), flanking_seq]
			elif seq_align_pos[i][5] > seq_align_pos[i + 1][5] > seq_align_pos[i][4] or seq_align_pos[i + 1][5] > seq_align_pos[i][4] > seq_align_pos[i + 1][4]:
#				print("###################")
				flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][1], seq_align_pos[i+1][1] + realign_range, seq)

#				[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][5], seq_align_pos[i + 1][3], seq_align_pos[i + 1][4], "OVL", int(seq_align_pos[i + 1][4]) - int(seq_align_pos[i][5]), flanking_seq]
				if seq_align_pos[i][6] == "+":
					[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i + 1][3], seq_align_pos[i + 1][4], seq_align_pos[i][3], seq_align_pos[i][5], "OVL", int(seq_align_pos[i][5]) - int(seq_align_pos[i + 1][4]), flanking_seq]
				else:
#					print("AAAAAAAAAAAAAAAA!!!")
					[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][4], seq_align_pos[i + 1][3], seq_align_pos[i + 1][5], "OVL", int(seq_align_pos[i + 1][5]) - int(seq_align_pos[i][4]), flanking_seq]
			else:
#				if seq_align_pos[i][6] == "+" and int(seq_align_pos[i + 1][4]) < int(seq_align_pos[i][5]):
				if seq_align_pos[i][6] == "+":
					flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][1], seq_align_pos[i+1][1] + realign_range, seq)
					if int(seq_align_pos[i + 1][4]) < int(seq_align_pos[i][5]):
#						distance = int(seq_align_pos[i][5]) - int(seq_align_pos[i + 1][4])
#						print("1")
#						print(seq_align_pos[i])
#						print(seq_align_pos[i + 1])
						[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i+1][4], seq_align_pos[i + 1][3], seq_align_pos[i][5], "TRS", int(seq_align_pos[i][5]) - int(seq_align_pos[i + 1][4]), flanking_seq]
					else:
#						distance = int(seq_align_pos[i + 1][4]) - int(seq_align_pos[i][5])
						[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][5], seq_align_pos[i + 1][3], seq_align_pos[i + 1][4], "TRS", int(seq_align_pos[i + 1][4]) - int(seq_align_pos[i][5]), flanking_seq]
#				elif seq_align_pos[i][6] == "-" and int(seq_align_pos[i + 1][5]) > int(seq_align_pos[i][4]):
				elif seq_align_pos[i][6] == "-":
					flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][1], seq_align_pos[i+1][1] + realign_range, seq)
					if int(seq_align_pos[i + 1][5]) > int(seq_align_pos[i][4]):
#						distance = int(seq_align_pos[i + 1][5]) - int(seq_align_pos[i][4])
#						print("2")
#						print(seq_align_pos[i])
#						print(seq_align_pos[i + 1])
						[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i][4], seq_align_pos[i + 1][3], seq_align_pos[i + 1][5], "TRS", seq_align_pos[i + 1][5] - seq_align_pos[i][4], flanking_seq]
					else:
#						distance = int(seq_align_pos[i][4]) - int(seq_align_pos[i + 5][5])
#						print(seq_align_pos[i])
#						print(seq_align_pos[i + 1])
						[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i + 1][5], seq_align_pos[i + 1][3], seq_align_pos[i][4], "TRS", int(seq_align_pos[i][4]) - int(seq_align_pos[i + 1][5]), flanking_seq]
#				else:
#					print("STRANGE", line)

#				else:
#					flanking_seq = get_flanking_seq(seq_align_pos[i][1] - realign_range, seq_align_pos[i][1], seq_align_pos[i+1][1], seq_align_pos[i+1][1] + realign_range, seq)
#					distance = int(seq_align_pos[i][5]) - int(seq_align_pos[i + 1][4])
#					print("3")
#					print(seq_align_pos[i])
#					print(seq_align_pos[i + 1])
#					[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], seq_align_pos[i + 1][4], seq_align_pos[i + 1][3], seq_align_pos[i][5], "TRS", seq_align_pos[i + 1][4] - seq_align_pos[i][5], flanking_seq]
		elif seq_align_pos[i][6] != seq_align_pos[i + 1][6]:
#[2, 1580, 6692, 'chr8', 44878726, 44880519, '+', 15.1, '000a7166-808a-43d5-aad6-905b8f29ba67_Basecall_1D_template', '2']
			BP_gpos1 = 0
			BP_gpos2 = 0
			if seq_align_pos[i][6] == "-" and seq_align_pos[i + 1][6] == "+":
				flanking_seq1 = rev_con(seq[seq_align_pos[i][1] - realign_range: seq_align_pos[i][1]])
				flanking_seq2 = seq[seq_align_pos[i+1][1]: seq_align_pos[i+1][1] + realign_range]
				flanking_seq = flanking_seq1 + flanking_seq2

				BP_gpos1 = int(seq_align_pos[i][4]) 
				BP_gpos2 = int(seq_align_pos[i + 1][4])
			elif seq_align_pos[i][6] == "+" and seq_align_pos[i + 1][6] == "-":
				flanking_seq1 = seq[seq_align_pos[i][1] - realign_range: seq_align_pos[i][1]]
				flanking_seq2 = rev_con(seq[seq_align_pos[i+1][1]: seq_align_pos[i+1][1] + realign_range])
				flanking_seq = flanking_seq1 + flanking_seq2

				BP_gpos1 = int(seq_align_pos[i][5])
				BP_gpos2 = int(seq_align_pos[i + 1][5])
			if BP_gpos2 - BP_gpos1 > 0:
				distance = BP_gpos2 - BP_gpos1
				[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], BP_gpos1, seq_align_pos[i + 1][3], BP_gpos2, "DIF", distance, flanking_seq]
			else:
				distance = BP_gpos1 - BP_gpos2
				[chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq] = [seq_align_pos[i][3], BP_gpos2, seq_align_pos[i + 1][3], BP_gpos1, "DIF", distance, flanking_seq] #change 171128
		print(chr1, breakpoint1, chr2, breakpoint2, SV_type, distance, flanking_seq, line, sep="\t")
#618ccd46-7719-4f97-a8f0-22ba84b36202_Basecall_1D_template       chr10,chr10     63198043,63198925;63198786,63198904     840,114 935,935 -,+     16,2048 74S,21S;815H,6H 60,20   74
	
