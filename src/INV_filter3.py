import sys
import re

read_col_num = int(sys.argv[2])

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	read_l = line_l[read_col_num].split('|')	

#62f92a83-b580-406f-b209-3022e16755d3;split;chr3;64549973;chr3;64550895;4321;3658,28;4314,3568;64549973,64547303;64549281,64550895;15.8,15.8;60,60;-,+;-,-;0.020828511918537376;-|f9003fcf-b129-42b5-b0ee-786609c9f916;split;chr3;64550105;chr3;64551014;6170;86,1773;1733,5831;64551823,64551014;64550105,64555422;14.9,14.9;60,60;-,+;-,-;0.006482982171799027;-
	
	pass_reads_num = 0
	num_long_reads = 0
	minimap2_reads = 0
	segment = []
	for read_tmp in read_l:
		read_tmp_l = read_tmp.split(";")
		
		if read_tmp_l[1] == "split":
			mq_l = read_tmp_l[12].split(",")
			start_l = read_tmp_l[9].split(",")
			end_l = read_tmp_l[10].split(",")
			mapper_l = read_tmp_l[14].split(",")
			strand_l = read_tmp_l[13].split(",")

			start1, start2, end1, end2, strand1, strand2 = 0, 0, 0, 0, "", ""
			if start_l[0] == read_tmp_l[3]:
				(start1, end1, strand1) = (start_l[0], end_l[0], strand_l[0])
			elif end_l[0] == read_tmp_l[3]:
				(start1, end1, strand1) = (end_l[0], start_l[0], strand_l[0])
			if start_l[1] == read_tmp_l[5]:
				(start2, end2, strand2) = (start_l[1], end_l[1], strand_l[1])
			elif end_l[1] == read_tmp_l[5]:
				(start2, end2, strand2) = (end_l[1], start_l[1], strand_l[1])
		
			if int(start1) > int(end1):
				start1, end1 = end1, start1

			if int(start2) > int(end2):
				start2, end2 = end2, start2

			segment.append([int(start1), int(end1), int(start2), int(end2), strand1, strand2])

	pair = []
	for i in range(len(segment)):
		for j in range(i + 1, len(segment)):
			if segment[i][0] <= segment[j][0] <= segment[i][1] or segment[i][0] <= segment[j][1] <= segment[i][1] or segment[j][0] <= segment[i][0] <= segment[j][1]:
				pair.append(i)
				pair.append(j)
			elif segment[i][2] <= segment[j][2] <= segment[i][3] or segment[i][2] <= segment[j][3] <= segment[i][3] or segment[j][2] <= segment[i][2] <= segment[j][3]:
				pair.append(i)
				pair.append(j)
 
#	print("pair", pair)
	if len(pair) > 0:
#		print(">>>>>>>>>>>>>>>>>>")
#		print("samtools view /share/amed_snt/WORK/fujimoto/nanopore/180615/RK014_C_WGS.All/RK014_C_WGS.All.fastq.minimap2.merge.sort.bam ", line_l[0] + ":" + str(int(line_l[1]) - 1000) + "-" + str(int(line_l[1]) + 1000), " > /share/amed_snt/WORK/fujimoto/nanopore/180908/RK014_C_" + line_l[0] + "_" + line_l[1] + ".sam")	
#		print("samtools view /share/amed_snt/WORK/fujimoto/nanopore/180615/RK014_C_WGS.All/RK014_C_WGS.All.fastq.minimap2.merge.sort.bam ", line_l[2] + ":" + str(int(line_l[3]) - 1000) + "-" + str(int(line_l[3]) + 1000), " > /share/amed_snt/WORK/fujimoto/nanopore/180908/RK014_C_" + line_l[2] + "_" + line_l[3] + ".sam")	
		print(line)
