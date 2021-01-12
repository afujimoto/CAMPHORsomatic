import sys
import re

min_mq = int(sys.argv[2])
min_distance_prop = float(sys.argv[3])
min_distance = float(sys.argv[4])
r_dis_g_dis_ratio = float(sys.argv[5])
close_ins_len_prop = sys.argv[6]
min_prop_distance_from_edge = float(sys.argv[7])
min_segment_length = int(sys.argv[8])

segment_col_num = int(sys.argv[9])

ins_col_num = int(sys.argv[10])

min_read_num = int(sys.argv[11])
min_read_with_segment_num = int(sys.argv[12])
min_minimap2_read_num = int(sys.argv[13])

close_ins_len_prop_l = close_ins_len_prop.split(",")
#print(close_ins_len_prop_l, len(close_ins_len_prop_l))
if len(close_ins_len_prop_l) != 1 and len(close_ins_len_prop_l) != 3:
	print("close_ins_len_prop !!!", close_ins_len_prop, "<close_ins_len_prop> or <close_ins_len_prop>,<deletion length for close_ins_len_prop>,<close_ins_len_prop2>")
	exit()

#8797632c-7a45-4e8c-ace2-2928a60cfb55;within;chr10;74960;chr10;75215;7726;25;5975;70087;76506;12.9;60;+;-;-;178
#breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper]), r_dis_g_dis_ratio, "-"])
#breakpoint_info.append([read_tmp_l[6], "within", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], r_start, r_end, g_start, g_end, bq, mq, strand, mapper, "-", total_ins_len])

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	read_l = line_l[segment_col_num].split('|')	
	
	min_distance_cutoff = min_distance_prop*int(line_l[6])
	if min_distance > min_distance_cutoff:
		min_distance_cutoff = min_distance
	
	pass_reads_num = 0
	num_long_reads = 0
	minimap2_reads = 0
	for read_tmp in read_l:
		read_tmp_l = read_tmp.split(";")
		
		if read_tmp_l[1] == "split":
			mq_l = read_tmp_l[12].split(",")
			start_l = read_tmp_l[9].split(",")
			end_l = read_tmp_l[10].split(",")
			mapper_l = read_tmp_l[14].split(",")

			if int(mq_l[0]) < min_mq or int(mq_l[1]) < min_mq:
				continue
			if abs(int(line_l[1]) - int(read_tmp_l[3])) > min_distance_cutoff  or abs(int(line_l[3]) - int(read_tmp_l[5])) > min_distance_cutoff:
				continue
			if float(read_tmp_l[15]) > r_dis_g_dis_ratio:
				continue

			pass_reads_num += 1
		
			if abs(int(end_l[0]) - int(start_l[0])) >= min_segment_length and abs(int(end_l[1]) - int(start_l[1])) >= min_segment_length:
				num_long_reads += 1

			if mapper_l[0] == "-" and mapper_l[1] == "-":
				minimap2_reads += 1

#8797632c-7a45-4e8c-ace2-2928a60cfb55;within;chr10;74960;chr10;75215;7726;25;5975;70087;76506;12.9;60;+;-;-;178
		if read_tmp_l[1] == "within":
			if int(read_tmp_l[12]) < min_mq:
				continue
#			print(read_tmp_l)
			if abs(int(line_l[1]) - int(read_tmp_l[3])) > min_distance_cutoff  or abs(int(line_l[3]) - int(read_tmp_l[5])) > min_distance_cutoff:
				continue
	
			distance_from_edge1 = min(abs(int(read_tmp_l[9]) - int(line_l[1])), abs(int(read_tmp_l[9]) - int(line_l[3])))
			distance_from_edge2 = min(abs(int(read_tmp_l[10]) - int(line_l[1])), abs(int(read_tmp_l[10]) - int(line_l[3])))

			if distance_from_edge1 < min_prop_distance_from_edge*float(read_tmp_l[6]) or distance_from_edge2 < min_prop_distance_from_edge*float(read_tmp_l[6]):
				continue

			if len(close_ins_len_prop_l) == 1:
				close_ins_len_prop = float(close_ins_len_prop_l[0])
				if int(read_tmp_l[16]) >= close_ins_len_prop*float(int(read_tmp_l[5]) - int(read_tmp_l[3])):
#					print("close ins", close_ins_len_prop, read_tmp_l[16], int(read_tmp_l[5]) - int(read_tmp_l[3]))
					continue
			elif len(close_ins_len_prop_l) == 3:
				(close_ins_len_prop, close_ins_len_prop_len, close_ins_len_prop2) = (float(close_ins_len_prop_l[0]), int(close_ins_len_prop_l[1]), float(close_ins_len_prop_l[2]))
				if int(read_tmp_l[5]) - int(read_tmp_l[3]) >= close_ins_len_prop_len and int(read_tmp_l[16]) >= close_ins_len_prop2*float(int(read_tmp_l[5]) - int(read_tmp_l[3])):
#					print("close ins 2", close_ins_len_prop2, read_tmp_l[16], int(read_tmp_l[5]) - int(read_tmp_l[3]))
					continue
				elif int(read_tmp_l[16]) >= close_ins_len_prop*float(int(read_tmp_l[5]) - int(read_tmp_l[3])):
#					print("close ins 3", close_ins_len_prop, read_tmp_l[16], int(read_tmp_l[5]) - int(read_tmp_l[3]))
					continue


			pass_reads_num += 1

			if int(read_tmp_l[6]) >= min_segment_length:
				num_long_reads += 1

			if read_tmp_l[14] == "-":
				minimap2_reads += 1
	
	if pass_reads_num >= min_read_num and num_long_reads >= min_read_with_segment_num and minimap2_reads >= min_minimap2_read_num:
		read_filter = ",".join([str(pass_reads_num), str(num_long_reads), str(minimap2_reads)])	
		print("\t".join(line_l[0:ins_col_num]), read_filter, "\t".join(line_l[ins_col_num:len(line_l)]), sep="\t")
