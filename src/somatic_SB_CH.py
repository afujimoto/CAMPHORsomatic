import re
import sys
from operator import itemgetter
import subprocess
import warnings

def find_close_rep(chr, start, end, rep_dic):
	rep = []
	
	if chr in rep_dic:
		for rep_l in rep_dic[chr]:
			if int(rep_l[1]) <= start <= int(rep_l[2]) or int(rep_l[1]) <= end <= int(rep_l[2]) or start <= int(rep_l[2]) <= end:
				rep.append(";".join(rep_l))
			if end < int(rep_l[1]):
				break
	return rep

def check_run_ch(read_name1, read_name2, fastq_file_name):
	f = open(fastq_file_name)

	read_name1_all = ""
	read_name2_all = ""
	for line in f:
		if line[0] == "@" and re.search(read_name1, line):
			read_name1_all = line
		if line[0] == "@" and re.search(read_name2, line):
			read_name2_all = line
		
		if len(read_name1_all) and len(read_name2_all):
			return read_name1_all, read_name2_all

f = open(sys.argv[1])
repeat_range = int(sys.argv[2])
bam = sys.argv[3]
sr_file = sys.argv[4]
rmsk_file = sys.argv[5]
fastq_file = sys.argv[6]
segment_col_num = int(sys.argv[7])
insert_col_num = int(sys.argv[8])

#chr1    1955018 1955089 trf,18,3.9,GGAAGTGCGGCTCACACC
#chr1    1963456 1963483 trf,1,27,A
sr_dic = {}
sr_f = open(sr_file)
for line in sr_f:
	line = line.replace("\n", "")
	chain_l = re.split("\t", line)
	rep_l = re.split(",", chain_l[3])

	if len(rep_l[3]) > 1:
		continue

	sr_dic.setdefault(str(chain_l[0]), []).append(chain_l)

#chr1    99046   99116   (T)n    585,37,130,0,0,chr1,99046,99116,-248857306,+,(T)n,Simple_repeat,Simple_repeat,1,70,0,1
#chr1    104160  104196  (AC)n   585,37,29,0,0,chr1,104160,104196,-248852226,+,(AC)n,Simple_repeat,Simple_repeat,1,36,0,1
rmsk_dic = {}
rmsk_f = open(rmsk_file)
for line in rmsk_f:
	line = line.replace("\n", "")
	chain_l = re.split("\t", line)

	if len(chain_l[3]) - 3 > 1:
		continue

	rmsk_dic.setdefault(str(chain_l[0]), []).append(chain_l)

#8c345ead-d9c4-4827-b3e9-3ce3a8ed9aec;within;chr10;2291996;chr10;2292181;7589;26;7580;2288464;2296806;16.9;60;+;-;-;3|a04af95c-7bff-4924-ae7a-743a1d9c8a17;within;chr10;2292005;chr10;2292137;15522;32;15520;2283153;2298890;17.6;60;+;-;-;3
read_name = {}
for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)
	read_l = re.split('\|', line_l[segment_col_num])
	for j in range(len(read_l)):
		read_tmp_l = re.split(";", read_l[j])
		read_name[read_tmp_l[0]] = 1
f.close()

f = open(fastq_file)
total_num_read = len(read_name)
for line in f:
	if line[0] == "@":
		read_name_fastq = (line.split(" "))[0]
		read_name_fastq = read_name_fastq.replace("@", "")
		if read_name_fastq in read_name:
			line = line.replace("\n", "")
			read_name[read_name_fastq] = line
			total_num_read -= 1

	if total_num_read == 0:
		break
f.close()

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)

	read_l = re.split('\|', line_l[segment_col_num])

#	print("read_l", read_l)

	breakpoint_info = []
	for read_tmp in read_l:
		read_tmp_l = read_tmp.split(";")
		breakpoint_info.append(read_tmp_l)

#	print("breakpoint_info", breakpoint_info)

#8c345ead-d9c4-4827-b3e9-3ce3a8ed9aec;within;chr10;2291996;chr10;2292181;7589;26;7580;2288464;2296806;16.9;60;+;-;-;3|a04af95c-7bff-4924-ae7a-743a1d9c8a17;within;chr10;2292005;chr10;2292137;15522;32;15520;2283153;2298890;17.6;60;+;-;-;3
#d65e919b-c437-4b03-b0e8-d12d582d5834;split;chr11;4948517;chr11;4957169;8215;29,6290;6282,8213;4942250,4957169;4948517,4959103;14.5,14.5;60,60;+,+;-,-;0.0009246417013407304;-|e805a912-9215-4521-ad0f-a9ea0fbef3c5;split;chr11;4948526;chr11;4957159;9038;29,7119;7117,9032;4941142,4957159;4948526,4959144;18.8,18.8;60,60;+,+;-,-;0.00023166917641607784;-|bb866106-b159-4eed-b60f-47c144152793;split;chr11;4948527;chr11;4957177;7887;41,7245;7230,7882;4941327,4957177;4948527,4957827;19.1,19.1;60,60;+,+;-,-;0.0017341040462427746;-
	distance_btw_read = {}
	for i in range(len(breakpoint_info)):
		distance_btw_read[i] = 0

#	print("distance_btw_read", distance_btw_read)

	name_of_read_l = []
	num_read_with_long_segment = 0
	num_split_read = 0
	for i in range(len(breakpoint_info)):
#		print(breakpoint_info[i])
		if breakpoint_info[i][1] == "within":
			continue

		g_start = breakpoint_info[i][9].split(",")
		g_end = breakpoint_info[i][10].split(",")
		segment_len1_1 = abs(int(g_end[0]) - int(g_start[0]))
		segment_len1_2 = abs(int(g_end[1]) - int(g_start[1]))

		for j in range(len(breakpoint_info)):
			if breakpoint_info[j][1] == "within":
				continue
			if j <= i:
				continue

			g_start = breakpoint_info[j][9].split(",")
			g_end = breakpoint_info[j][10].split(",")
			segment_len2_1 = abs(int(g_end[0]) - int(g_start[0]))
			segment_len2_2 = abs(int(g_end[1]) - int(g_start[1]))

#			if abs(segment_len1_1 - segment_len2_1) < segment_length_diff_len or abs(segment_len1_2 - segment_len2_2) < segment_length_diff_len:
#@e22a5408-8e07-4b32-bf9a-e58c5c957654 runid=38469a38d9ea3789b624bd568739d7e40ba7c717 read=214 ch=383 start_time=2017-06-27T08:45:44Z
			read_name1_all = read_name[breakpoint_info[i][0]]
			read_name2_all = read_name[breakpoint_info[j][0]]
			read_name1_all_l = read_name1_all.split(" ")
			read_name2_all_l = read_name2_all.split(" ")
			read1 = int(read_name1_all_l[2].replace("read=", ""))
			read2 = int(read_name2_all_l[2].replace("read=", ""))

			if (read_name1_all_l[1] == read_name2_all_l[1] and read_name1_all_l[3] == read_name2_all_l[3] and abs(read1 - read2) <= 30) or (read_name1_all_l[1] == read_name2_all_l[1] and (abs(segment_len1_1 - segment_len2_1) < 30 or abs(segment_len1_2 - segment_len2_2) < 30)):
					name_of_read_l.append(",".join(read_name1_all_l))
					name_of_read_l.append(",".join(read_name2_all_l))
					distance_btw_read[i] += 1	
#					print(name_of_read_l)
#					print(distance_btw_read)

	right_rep = 0
	left_sr = find_close_rep(line_l[0], int(line_l[1]) - repeat_range, int(line_l[1]), sr_dic)
	left_rmsk = find_close_rep(line_l[0], int(line_l[1]) - repeat_range, int(line_l[1]), rmsk_dic)

	left_rep = 0
	right_sr = find_close_rep(line_l[0], int(line_l[3]), int(line_l[3]) + repeat_range, sr_dic)
	right_rmsk = find_close_rep(line_l[0], int(line_l[3]), int(line_l[3]) + repeat_range, rmsk_dic)

	close_rep = "|".join([";".join(left_sr), ";".join(left_rmsk), ";".join(right_sr), ";".join(right_sr)])

#8c345ead-d9c4-4827-b3e9-3ce3a8ed9aec;within;chr10;2291996;chr10;2292181;7589;26;7580;2288464;2296806;16.9;60;+;-;-;3|a04af95c-7bff-4924-ae7a-743a1d9c8a17;within;chr10;2292005;chr10;2292137;15522;32;15520;2283153;2298890;17.6;60;+;-;-;3
#d65e919b-c437-4b03-b0e8-d12d582d5834;split;chr11;4948517;chr11;4957169;8215;29,6290;6282,8213;4942250,4957169;4948517,4959103;14.5,14.5;60,60;+,+;-,-;0.0009246417013407304;-|e805a912-9215-4521-ad0f-a9ea0fbef3c5;split;chr11;4948526;chr11;4957159;9038;29,7119;7117,9032;4941142,4957159;4948526,4959144;18.8,18.8;60,60;+,+;-,-;0.00023166917641607784;-|bb866106-b159-4eed-b60f-47c144152793;split;chr11;4948527;chr11;4957177;7887;41,7245;7230,7882;4941327,4957177;4948527,4957827;19.1,19.1;60,60;+,+;-,-;0.0017341040462427746;-

	total_within_reads = 0
	strand_bias = 0
	positive_strand = 0
	negative_strand = 0
	for i in range(0, len(breakpoint_info)):
		if breakpoint_info[i][1] == "within":
			total_within_reads += 1
			if breakpoint_info[i][13] == "+":
				positive_strand += 1
			else:
				negative_strand += 1

	close_repeat_strand_bias = 0
	if (len(left_sr) or len(left_rmsk)) and negative_strand == 0 and positive_strand > 0:
		close_repeat_strand_bias = 1
	if (len(right_sr) or len(right_rmsk)) and positive_strand == 0 and negative_strand > 0:
		close_repeat_strand_bias = 1

	name_of_read = "-"
	if len(name_of_read_l) > 0:
		name_of_read = ";".join(name_of_read_l)

	number_of_diff_len_segment = 0
	for tmp in distance_btw_read:
#		print(tmp, distance_btw_read[tmp], number_of_diff_len_segment)
		if distance_btw_read[tmp] == 0:
			number_of_diff_len_segment += 1

	print("\t".join(line_l[0:insert_col_num]), number_of_diff_len_segment, close_repeat_strand_bias, str(positive_strand) + "/" + str(negative_strand), close_rep, "\t".join(line_l[insert_col_num:len(line_l)]), name_of_read, sep="\t")

