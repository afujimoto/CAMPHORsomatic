import re
import sys
from operator import itemgetter

#a96a29f8-a632-4260-9ae6-5ad9b18b9815_Basecall_1D_template	2064	chr3	85286853	0	6918H12M2I6M1I29M1D6M1D19M1597H	*	0	0	TAATAGAAAGTATTATTTGATATTATTCTGATTTTATAATATATACTTCATTTTCATATATTTCGCTCATTTGAT	$(%"#((&*++-.)*)-,&+..30)*)()+2.2-&*,00-+&.)'%++&().20%&%%'&())'(.&''',(&*+NM:i:21	MD:Z:3A3G2C2A1A0C2G1T23T1^G2A3^T1C0A2A2A0G1A7	AS:i:31	XS:i:30	SA:Z:chrX,8237011,+,3787S9M1I26M1D12M2D2M1D13M1I15M1I9M1I9M1D7M2D8M1D2M3D12M1I4M1D4M1I23M1I12M2I16M4611S,2,71;chr4,132898266,+,8289S12M1I16M2D6M2I1M1I23M2D7M1D4M3I5M1D4M1D7M1I10M4I15M1D26M1D2M2D3M1D8M1D5M1D1M1D4M2D7M2D8M1I10M104S,0,81;	XA:Z:chr11,-127123786,6730S10M1I2M1I4M1I1M1I4M1D15M3D3M2D6M1I15M1I6M1I2M1D15M1D3M1D4M2D5M1D1M1D4M2D9M1D1M1D19M1I9M3D3M1I9M1D13M2I5M3I6M3I6M1I7M3I24M1I11M1I4M1D18M1593S,114;chr13,+35485820,1637S29M6924S,1;

def get_read_length(read_info_tmp_l):
	if read_info_tmp_l[2] == "*":
		read_length = len(read_info_tmp_l[9])
		return(read_length)

	cigar_type1 = re.split(r'\d+', read_info_tmp_l[5])
	cigar_len1 = re.split(r'[a-zA-Z]+', read_info_tmp_l[5])
	del cigar_type1[0]
	del cigar_len1[-1]

	read_length = 0
	for i in range(0, len(cigar_type1)):
		if cigar_type1[i] != "D":
			read_length += int(cigar_len1[i])

	return(read_length)

def get_start_end_on_read(read_tmp):
	read_tmp_l = re.split('\t', read_tmp)
	cigar_type1 = re.split(r'\d+', read_tmp_l[5])
	cigar_len1 = re.split(r'[a-zA-Z]+', read_tmp_l[5])
	del cigar_type1[0]
	del cigar_len1[-1]

	if int(read_tmp_l[1]) & 0x10:
		cigar_type1.reverse()
		cigar_len1.reverse()

	map_start = 0
	map_end = 0
	for i in range(0, len(cigar_type1)):
		if cigar_type1[i] == "D":
			continue
		if i == 0 and (cigar_type1[i] == "S" or cigar_type1[i] == "H"):
			map_start = int(cigar_len1[i]) + 1
			map_end += int(cigar_len1[i])
		elif cigar_type1[i] != "S" and cigar_type1[i] != "H":
			map_end += int(cigar_len1[i])
	return(map_start, map_end)

def get_start_end_on_genome(read_tmp):
	read_tmp_l = re.split('\t', read_tmp)
	cigar_type1 = re.split(r'\d+', read_tmp_l[5])
	cigar_len1 = re.split(r'[a-zA-Z]+', read_tmp_l[5])
	del cigar_type1[0]
	del cigar_len1[-1]

	map_start = int(read_tmp_l[3])
	map_end = 0
	for i in range(0, len(cigar_type1)):
		if cigar_type1[i] == "D" or cigar_type1[i] == "M":
			map_end += int(cigar_len1[i])
	map_end += map_start
	return(map_start, map_end)

def get_max_cover_alignment(start, end, read, ol_prop):
	max_cover_rate = 0
	region_length = end - start + 1

	max_aligned_read = ""
	max_match = 0
	for read_tmp in read:
		(map_start, map_end) = get_start_end_on_read(read_tmp)

def get_max_cover_alignment(start, end, read, ol_prop):
	max_cover_rate = 0
	region_length = end - start + 1

	max_aligned_read = ""
	max_match = 0
	for read_tmp in read:
#		print("read_tmp => ", read_tmp)
		(map_start, map_end) = get_start_end_on_read(read_tmp)
#		print("map_start map_end start end", map_start, map_end, start, end)
		if end < int(map_start) or int(map_end) < start:
			continue
		if start - ol_prop*region_length <= map_start and map_end <= end + ol_prop*region_length:
			if map_end - map_start + 1 > max_match:
				max_match = map_end - map_start + 1
				max_aligned_read = read_tmp
	return(max_aligned_read)

def get_uncovered_region(covered_start, covered_end, read_length):
		uncoverd_region_start = []
		uncoverd_region_end = []

		zero = 0
		for tmp in covered_start:
			if tmp == "0":
				zero += 1
		if zero == 0:
			region_boundary = [0]

		end_eq_read_length = 0
		for tmp in covered_end:
			if tmp == read_length:
				end_eq_read_length += 1
		if end_eq_read_length == 0:
			region_boundary.append(read_length)

		region_boundary.extend(covered_start)
		region_boundary.extend(covered_end)
		region_boundary.sort()
		
		for i in range(len(region_boundary) - 1):
			covered = 0
			for j in range(len(covered_start)):
				if int(covered_start[j]) <= int(region_boundary[i]) and int(region_boundary[i + 1]) <= int(covered_end[j]):
					covered = 1

			if covered == 0:
				if region_boundary[i] == 0:
					uncoverd_region_start.append(region_boundary[i])
				else:
					uncoverd_region_start.append(region_boundary[i] + 1)
				if region_boundary[i + 1] == read_length:
					uncoverd_region_end.append(region_boundary[i + 1])
				else:
					uncoverd_region_end.append(region_boundary[i + 1] - 1)
	
		return(uncoverd_region_start, uncoverd_region_end)

sam_f = open(sys.argv[1])
min_indel_size = int(sys.argv[2])
max_overlap_prop = float(sys.argv[3])
pre_name = ""
read_info = []

for line in sam_f:
	if line[0] == "@":
		continue

	line = line.replace('\n', '')
	line_l = re.split('\t', line)

	if line_l[2] == "*":
		continue

#	print(line)

#	print(read_info)

	if pre_name != line_l[0] and len(read_info) == 1:
		read_info_tmp_l = re.split("\t", read_info[0])
		cigar_type1 = re.split(r'\d+', read_info_tmp_l[5])
		cigar_len1 = re.split(r'[a-zA-Z]+', read_info_tmp_l[5])
		del cigar_type1[0]
		del cigar_len1[-1]

		for i in range(0, len(cigar_type1)):
			if (cigar_type1[i] == "D" or cigar_type1[i] == "I") and int(cigar_len1[i]) >= min_indel_size:
				read_name_l = re.split("/", read_info_tmp_l[0])
				read_name = read_name_l[0]
#				bq = read_name_l[1]
				bq = "-"
				read_length = len(read_info_tmp_l[9])
				(covered_start, covered_end) = get_start_end_on_read(read_info[0])
				(g_map_start, g_map_end) = get_start_end_on_genome(read_info[0])	
				strand = "+"
				if int(read_info_tmp_l[1]) & 0x10:
					strand = "-"
				
				aligner = ""
				for tmp in read_info_tmp_l:
					if "ALN:Z" in tmp:
						aligner = tmp
						aligner = aligner.replace("ALN:Z:", "")

				if len(aligner) == 0:
					aligner = "-"
				
				out = [read_name, str(read_length), ",".join((str(covered_start), str(covered_end))), read_info_tmp_l[2], ",".join((str(g_map_start), str(g_map_end))), strand, bq, read_info_tmp_l[4], read_info_tmp_l[5], read_info_tmp_l[9], aligner, "NA"]
				print ("\t".join(out))
				break
		read_info = []
	
	seq_align_pos = []
	cover_read_chr = []
	cover_read_chr_start = []
	cover_read_chr_end = []
	cover_read_chr_strand = []
	cover_read_start = []
	cover_read_end = []
	cover_read_alignment_cigar = []
	cover_read_alignment_seq = []
	cover_read_alignment_MDZ = []
	cover_read_mq = []
	cover_read_bq = []
	cover_read_mapper = []
	if pre_name != line_l[0] and len(read_info) > 1:
		read_length = 0
		seq_align_pos = []
		read_name = ""
		seq = ""
		for read_info_tmp in read_info:
			read_info_tmp_l = re.split("\t", read_info_tmp)
			read_name = read_info_tmp_l[0]
			cigar_type1 = re.split(r'\d+', read_info_tmp_l[5])
			del cigar_type1[0]
#			print(read_info_tmp_l[9])
			if cigar_type1[0] != "H" and cigar_type1[-1] != "H" and read_info_tmp_l[9] != "*":
				seq = read_info_tmp_l[9]
			read_length = get_read_length(read_info_tmp_l)

			(map_start, map_end) = get_start_end_on_read(read_info_tmp)
			map_length = map_end - map_start + 1
			seq_align_pos.append([map_length, read_info_tmp])

#		print(read_length)

		seq_align_pos.sort(key=itemgetter(0), reverse = True)
		read_info2 = []
		for read in seq_align_pos:
#			print("sort=>", read[1])
			read_info2.append(read[1])

		tested_uncovered_region = {}
		covered_proportion = 0
		coverd_region_start = []
		coverd_region_end = []
		uncoverd_region_start = [0]
		uncoverd_region_end = [read_length]
		pre_uncoverd_region_start = []
		pre_uncoverd_region_end = []
		while covered_proportion < 0.9:
			if pre_uncoverd_region_start == uncoverd_region_start and pre_uncoverd_region_end == uncoverd_region_end:
				break

			pre_uncoverd_region_start = uncoverd_region_start
			pre_uncoverd_region_end = uncoverd_region_end
			for i in range(len(uncoverd_region_start)):
				uncoverd_region_tmp = str(uncoverd_region_start[i]) + "\t" + str(uncoverd_region_end[i])
				if uncoverd_region_tmp in tested_uncovered_region:
					continue
	
				tested_uncovered_region[uncoverd_region_tmp] = 1
#				print(read_info2)
				max_aligned_read = get_max_cover_alignment(uncoverd_region_start[i], uncoverd_region_end[i], read_info2, max_overlap_prop)
#				print("max_aligned_read=>", max_aligned_read)
				if len(max_aligned_read) == 0:
					continue

				max_aligned_read_l = re.split('\t', max_aligned_read)
				(covered_start, covered_end) = get_start_end_on_read(max_aligned_read)
#				print("max_aligned_read_l=>", max_aligned_read_l)
#a96a29f8-a632-4260-9ae6-5ad9b18b9815_Basecall_1D_template	2064	chr3	85286853	0	6918H12M2I6M1I29M1D6M1D19M1597H	*	0	0	TAATAGAAAGTATTATTTGATATTATTCTGATTTTATAATATATACTTCATTTTCATATATTTCGCTCATTTGAT	$(%"#((&*++-.)*)-,&+..30)*)()+2.2-&*,00-+&.)'%++&().20%&%%'&())'(.&''',(&*+NM:i:21	MD:Z:3A3G2C2A1A0C2G1T23T1^G2A3^T1C0A2A2A0G1A7	AS:i:31	XS:i:30	SA:Z:chrX,8237011,+,3787S9M1I26M1D12M2D2M1D13M1I15M1I9M1I9M1D7M2D8M1D2M3D12M1I4M1D4M1I23M1I12M2I16M4611S,2,71;chr4,132898266,+,8289S12M1I16M2D6M2I1M1I23M2D7M1D4M3I5M1D4M1D7M1I10M4I15M1D26M1D2M2D3M1D8M1D5M1D1M1D4M2D7M2D8M1I10M104S,0,81;	XA:Z:chr11,-127123786,6730S10M1I2M1I4M1I1M1I4M1D15M3D3M2D6M1I15M1I6M1I2M1D15M1D3M1D4M2D5M1D1M1D4M2D9M1D1M1D19M1I9M3D3M1I9M1D13M2I5M3I6M3I6M1I7M3I24M1I11M1I4M1D18M1593S,114;chr13,+35485820,1637S29M6924S,1;
				(g_map_start, g_map_end) = get_start_end_on_genome(max_aligned_read)

				strand = "+"
				if int(max_aligned_read_l[1]) & 0x10:
					strand = "-"				
				
#				MDZ = ""
#				for tmp in max_aligned_read_l:
#					if re.search("^MD:Z:", tmp):
#						MDZ = tmp


				read_name_l = re.split("/", max_aligned_read_l[0])
				read_name = read_name_l[0]
#				ave_bq = read_name_l[1]
				ave_bq = "-"
				
#				bq = 0
#				for i in range(len(max_aligned_read_l[10])):
#					bq += ord(max_aligned_read_l[10][i]) - 33
					
#				ave_bq = str(round(bq/len(max_aligned_read_l[10]), 2))

				cover_read_chr.append(max_aligned_read_l[2])
				cover_read_chr_start.append(g_map_start)
				cover_read_chr_end.append(g_map_end)
				cover_read_chr_strand.append(strand)				
				cover_read_start.append(covered_start)
				cover_read_end.append(covered_end)
				cover_read_alignment_cigar.append(max_aligned_read_l[5])
#				cover_read_alignment_seq.append(max_aligned_read_l[9])
				cover_read_mq.append(max_aligned_read_l[4])
				cover_read_bq.append(ave_bq)

				aligner = ""
				for tmp in max_aligned_read_l:
					if "ALN:Z" in tmp:
						aligner = tmp
						aligner = aligner.replace("ALN:Z:", "")

				if len(aligner):
					cover_read_mapper.append(aligner)
				else:
					cover_read_mapper.append("-")
					
				coverd_region_start.append(covered_start)
				coverd_region_end.append(covered_end)
			(uncoverd_region_start, uncoverd_region_end) = get_uncovered_region(coverd_region_start, coverd_region_end, read_length)

		read_info = []

		cover_read_chr_pos_all = []
		for i in range(0, len(cover_read_chr_start)):
			tmp = str(cover_read_chr_start[i]) + "," + str(cover_read_chr_end[i])
			cover_read_chr_pos_all.append(tmp)

		read_pos_on_mapping = []
		for i in range(0, len(cover_read_start)):
			tmp = str(cover_read_start[i]) + "," + str(cover_read_end[i])
			read_pos_on_mapping.append(tmp)

#		out = [read_name, str(read_length), ";".join(read_pos_on_mapping), ";".join(cover_read_chr), ";".join(cover_read_chr_pos_all), ";".join(cover_read_chr_strand), ";".join(cover_read_bq), ";".join(cover_read_mq), ";".join(cover_read_alignment_cigar), ";".join(cover_read_alignment_MDZ)]
		out = [read_name, str(read_length), ";".join(read_pos_on_mapping), ";".join(cover_read_chr), ";".join(cover_read_chr_pos_all), ";".join(cover_read_chr_strand), ";".join(cover_read_bq), ";".join(cover_read_mq), ";".join(cover_read_alignment_cigar), seq, ";".join(cover_read_mapper), "NA"]
		print ("\t".join(out))

	pre_name = line_l[0]
	read_info.append(line)

if len(read_info) == 1:
	read_info_tmp_l = re.split("\t", read_info[0])
	cigar_type1 = re.split(r'\d+', read_info_tmp_l[5])
	cigar_len1 = re.split(r'[a-zA-Z]+', read_info_tmp_l[5])
	del cigar_type1[0]
	del cigar_len1[-1]

	for i in range(0, len(cigar_type1)):
		if (cigar_type1[i] == "D" or cigar_type1[i] == "I") and int(cigar_len1[i]) >= min_indel_size:
			read_name_l = re.split("/", read_info_tmp_l[0])
			read_name = read_name_l[0]
#			bq = read_name_l[1]
			bq = "-"
			read_length = len(read_info_tmp_l[9])
			(covered_start, covered_end) = get_start_end_on_read(read_info[0])
			(g_map_start, g_map_end) = get_start_end_on_genome(read_info[0])	
			strand = "+"
			if int(read_info_tmp_l[1]) & 0x10:
				strand = "-"		
			
			aligner = ""
			for tmp in read_info_tmp_l:
				if "ALN:Z" in tmp:
					aligner = tmp
					aligner = aligner.replace("ALN:Z:", "")

			if len(aligner) == 0:
				aligner = "-"

			out = [read_name, str(read_length), ",".join((str(covered_start), str(covered_end))), read_info_tmp_l[2], ",".join((str(g_map_start), str(g_map_end))), strand, bq, read_info_tmp_l[4], read_info_tmp_l[5], read_info_tmp_l[9], aligner, "NA"]
			print ("\t".join(out))
			break
	
seq_align_pos = []
cover_read_chr = []
cover_read_chr_start = []
cover_read_chr_end = []
cover_read_chr_strand = []
cover_read_start = []
cover_read_end = []
cover_read_alignment_cigar = []
cover_read_alignment_seq = []
cover_read_alignment_MDZ = []
cover_read_mq = []
cover_read_bq = []
cover_read_mapper = []
#print(read_info)
if len(read_info) > 1:
	read_length = 0
	seq_align_pos = []
	read_name = ""
	seq = ""
	for read_info_tmp in read_info:
#		print(read_info_tmp)
		read_info_tmp_l = re.split("\t", read_info_tmp)
		read_name = read_info_tmp_l[0]
		cigar_type1 = re.split(r'\d+', read_info_tmp_l[5])
		del cigar_type1[0]
#		print(read_info_tmp_l[9])
		if cigar_type1[0] != "H" and cigar_type1[-1] != "H" and read_info_tmp_l[9] != "*":
			seq = read_info_tmp_l[9]
		read_length = get_read_length(read_info_tmp_l)
#		print(read_length)

		(map_start, map_end) = get_start_end_on_read(read_info_tmp)
#		print(map_start, map_end)
		map_length = map_end - map_start + 1
		seq_align_pos.append([map_length, read_info_tmp])

#	print(read_length)

	seq_align_pos.sort(key=itemgetter(0), reverse = True)
#	print(seq_align_pos)
	read_info2 = []
	for read in seq_align_pos:
#		print("sort=>", read[0], read[1])
		read_info2.append(read[1])

	tested_uncovered_region = {}
	covered_proportion = 0
	coverd_region_start = []
	coverd_region_end = []
	uncoverd_region_start = [0]
	uncoverd_region_end = [read_length]
	pre_uncoverd_region_start = []
	pre_uncoverd_region_end = []
	while covered_proportion < 0.9:
#		print("while", pre_uncoverd_region_start, uncoverd_region_start)
		if pre_uncoverd_region_start == uncoverd_region_start and pre_uncoverd_region_end == uncoverd_region_end:
			break
		pre_uncoverd_region_start = uncoverd_region_start
		pre_uncoverd_region_end = uncoverd_region_end
		for i in range(len(uncoverd_region_start)):
#			print("=================", i)
			uncoverd_region_tmp = str(uncoverd_region_start[i]) + "\t" + str(uncoverd_region_end[i])
			if uncoverd_region_tmp in tested_uncovered_region:
				continue
	
			tested_uncovered_region[uncoverd_region_tmp] = 1
#			print("reads!!!!", read_info2)
			max_aligned_read = get_max_cover_alignment(uncoverd_region_start[i], uncoverd_region_end[i], read_info2, max_overlap_prop)
#			print("max_aligned_read=>", max_aligned_read)
			if len(max_aligned_read) == 0:
				continue

			max_aligned_read_l = re.split('\t', max_aligned_read)
			(covered_start, covered_end) = get_start_end_on_read(max_aligned_read)
#			print("max_aligned_read_l=>", max_aligned_read_l)
#a96a29f8-a632-4260-9ae6-5ad9b18b9815_Basecall_1D_template	2064	chr3	85286853	0	6918H12M2I6M1I29M1D6M1D19M1597H	*	0	0	TAATAGAAAGTATTATTTGATATTATTCTGATTTTATAATATATACTTCATTTTCATATATTTCGCTCATTTGAT	$(%"#((&*++-.)*)-,&+..30)*)()+2.2-&*,00-+&.)'%++&().20%&%%'&())'(.&''',(&*+NM:i:21	MD:Z:3A3G2C2A1A0C2G1T23T1^G2A3^T1C0A2A2A0G1A7	AS:i:31	XS:i:30	SA:Z:chrX,8237011,+,3787S9M1I26M1D12M2D2M1D13M1I15M1I9M1I9M1D7M2D8M1D2M3D12M1I4M1D4M1I23M1I12M2I16M4611S,2,71;chr4,132898266,+,8289S12M1I16M2D6M2I1M1I23M2D7M1D4M3I5M1D4M1D7M1I10M4I15M1D26M1D2M2D3M1D8M1D5M1D1M1D4M2D7M2D8M1I10M104S,0,81;	XA:Z:chr11,-127123786,6730S10M1I2M1I4M1I1M1I4M1D15M3D3M2D6M1I15M1I6M1I2M1D15M1D3M1D4M2D5M1D1M1D4M2D9M1D1M1D19M1I9M3D3M1I9M1D13M2I5M3I6M3I6M1I7M3I24M1I11M1I4M1D18M1593S,114;chr13,+35485820,1637S29M6924S,1;
			(g_map_start, g_map_end) = get_start_end_on_genome(max_aligned_read)

			strand = "+"
			if int(max_aligned_read_l[1]) & 0x10:
				strand = "-"				
				
#			MDZ = ""
#			for tmp in max_aligned_read_l:
#				if re.search("^MD:Z:", tmp):
#					MDZ = tmp


			read_name_l = re.split("/", max_aligned_read_l[0])
			read_name = read_name_l[0]
#			ave_bq = read_name_l[1]
			ave_bq = "-"
				
#			bq = 0
#			for i in range(len(max_aligned_read_l[10])):
#				bq += ord(max_aligned_read_l[10][i]) - 33
				
#			ave_bq = str(round(bq/len(max_aligned_read_l[10]), 2))

			cover_read_chr.append(max_aligned_read_l[2])
			cover_read_chr_start.append(g_map_start)
			cover_read_chr_end.append(g_map_end)
			cover_read_chr_strand.append(strand)				
			cover_read_start.append(covered_start)
			cover_read_end.append(covered_end)
			cover_read_alignment_cigar.append(max_aligned_read_l[5])
#			cover_read_alignment_seq.append(max_aligned_read_l[9])
			cover_read_mq.append(max_aligned_read_l[4])
			cover_read_bq.append(ave_bq)

			aligner = ""
			for tmp in max_aligned_read_l:
				if "ALN:Z" in tmp:
					aligner = tmp
					aligner = aligner.replace("ALN:Z:", "")

			if len(aligner):
				cover_read_mapper.append(aligner)
			else:
				cover_read_mapper.append("-")

			coverd_region_start.append(covered_start)
			coverd_region_end.append(covered_end)
		(uncoverd_region_start, uncoverd_region_end) = get_uncovered_region(coverd_region_start, coverd_region_end, read_length)
#		print(">>>>", uncoverd_region_start, uncoverd_region_end)
	read_info = []

	cover_read_chr_pos_all = []
	for i in range(0, len(cover_read_chr_start)):
		tmp = str(cover_read_chr_start[i]) + "," + str(cover_read_chr_end[i])
		cover_read_chr_pos_all.append(tmp)

	read_pos_on_mapping = []
	for i in range(0, len(cover_read_start)):
		tmp = str(cover_read_start[i]) + "," + str(cover_read_end[i])
		read_pos_on_mapping.append(tmp)

#	out = [read_name, str(read_length), ";".join(read_pos_on_mapping), ";".join(cover_read_chr), ";".join(cover_read_chr_pos_all), ";".join(cover_read_chr_strand), ";".join(cover_read_bq), ";".join(cover_read_mq), ";".join(cover_read_alignment_cigar), ";".join(cover_read_alignment_MDZ)]
	out = [read_name, str(read_length), ";".join(read_pos_on_mapping), ";".join(cover_read_chr), ";".join(cover_read_chr_pos_all), ";".join(cover_read_chr_strand), ";".join(cover_read_bq), ";".join(cover_read_mq), ";".join(cover_read_alignment_cigar), seq, ";".join(cover_read_mapper), "NA"]
	print ("\t".join(out))
