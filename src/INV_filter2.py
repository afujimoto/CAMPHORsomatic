import sys
import re

min_mq = int(sys.argv[2])
#min_distance_prop = float(sys.argv[3])
min_distance_cutoff = float(sys.argv[3])
r_dis_g_dis_ratio = float(sys.argv[4])
#close_ins_len = float(sys.argv[6])
min_prop_distance_from_edge = float(sys.argv[5])
min_segment_length = int(sys.argv[6])
mapper_filter = int(sys.argv[7])

ins_col_num = int(sys.argv[8])

min_read_num = int(sys.argv[9])
min_read_with_segment_num = int(sys.argv[10])
min_minimap2_read_num = int(sys.argv[11])

#breakpoint_info.append([read_tmp_l[6], "split", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], ",".join([str(r_start), str(r_start2)]), ",".join([str(r_end), str(r_end2)]), ",".join([str(g_start), str(g_start2)]), ",".join([str(g_end), str(g_end2)]), ",".join([bq, bq2]), ",".join([mq, mq2]), ",".join([strand, strand2]), ",".join([mapper, mapper]), r_dis_g_dis_ratio, "-"])
#breakpoint_info.append([read_tmp_l[6], "within", chr, read_tmp_l[1], chr, read_tmp_l[3], read_tmp_l[9], r_start, r_end, g_start, g_end, bq, mq, strand, mapper, "-", total_ins_len])

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
#	print(">>>>", line)
	line_l = line.split("\t")

	read_l = line_l[7].split('|')	
	
#	min_distance_cutoff = min_distance_prop*int(line_l[6])
#	if min_distance > min_distance_cutoff:
#		min_distance_cutoff = min_distance
	
	pass_reads_num = 0
	num_long_reads = 0
	minimap2_reads = 0
	for read_tmp in read_l:
		read_tmp_l = read_tmp.split(";")
#		print("read_tmp_l => ", read_tmp_l)
		
		if read_tmp_l[1] == "split":
			mq_l = read_tmp_l[12].split(",")
			start_l = read_tmp_l[9].split(",")
			end_l = read_tmp_l[10].split(",")
			mapper_l = read_tmp_l[14].split(",")

			if int(mq_l[0]) < min_mq or int(mq_l[1]) < min_mq:
#				print("mq", mq_l)
				continue
			if abs(int(line_l[1]) - int(read_tmp_l[3])) > min_distance_cutoff  or abs(int(line_l[3]) - int(read_tmp_l[5])) > min_distance_cutoff:
#				print("min_distance", line_l[1], read_tmp_l[3], line_l[3], read_tmp_l[5])
				continue
			if float(read_tmp_l[15]) > r_dis_g_dis_ratio:
#				print("r_dis_g_dis_ratio", r_dis_g_dis_ratio)
				continue

			pass_reads_num += 1
		
			if abs(int(end_l[0]) - int(start_l[0])) >= min_segment_length and abs(int(end_l[1]) - int(start_l[1])) >= min_segment_length:
#				print("min_segment_length", start_l, end_l)
#				print(int(end_l[0]) - int(start_l[0]), int(end_l[1]) - int(start_l[1]))
				num_long_reads += 1

			if mapper_l[0] == "-" and mapper_l[1] == "-":
#				print("mapper_l", mapper_l)
				minimap2_reads += 1

		
		if read_tmp_l[1] == "within":
#			print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
			if int(read_tmp_l[12]) < min_mq:
#				print("mq", read_tmp_l[12])
				continue
			if abs(int(line_l[1]) - int(read_tmp_l[3])) > min_distance_cutoff  or abs(int(line_l[3]) - int(read_tmp_l[5])) > min_distance_cutoff:
#				print("min_distance", line_l[1], read_tmp_l[3], read_tmp_l[5])
				continue

#			print(read_tmp_l[9], line_l[1], line_l[3])
	
			distance_from_edge1 = min(abs(int(read_tmp_l[9]) - int(line_l[1])), abs(int(read_tmp_l[9]) - int(line_l[3])))
			distance_from_edge2 = min(abs(int(read_tmp_l[10]) - int(line_l[1])), abs(int(read_tmp_l[10]) - int(line_l[3])))

			if distance_from_edge1 < min_prop_distance_from_edge*float(read_tmp_l[6]) or distance_from_edge2 < min_prop_distance_from_edge*float(read_tmp_l[6]):
#				print("distance_from_edge1", distance_from_edge1, "distance_from_edge2", distance_from_edge2)
				continue
				
#			if int(read_tmp_l[16]) >= close_ins_len:
#				print("close_ins_len", read_tmp_l[16])
#				continue

			pass_reads_num += 1
#			print("pass_reads_num", pass_reads_num)

			if int(read_tmp_l[6]) >= min_segment_length:
#				print("min_segment_length", read_tmp_l[6])
				num_long_reads += 1

			if read_tmp_l[14] == "-":
#				print("mapper_l", read_tmp_l[14])
				minimap2_reads += 1
	
#	print("pass_reads_num", pass_reads_num)
#	print("num_long_reads", num_long_reads)
#	print("minimap2_reads", minimap2_reads)

#	if pass_reads_num >= 2 and num_long_reads >= 1 and minimap2_reads >= 1:
#		print("PASS!!")

	if pass_reads_num >= min_read_num and num_long_reads >= min_read_with_segment_num and minimap2_reads >= min_minimap2_read_num:
		read_filter = ",".join([str(pass_reads_num), str(num_long_reads), str(minimap2_reads)])	
		print("\t".join(line_l[0:ins_col_num]), read_filter, "\t".join(line_l[ins_col_num:len(line_l)]), sep="\t")
