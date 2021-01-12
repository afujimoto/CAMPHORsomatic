import re
import sys
from operator import itemgetter
import subprocess

f = open(sys.argv[1])
cutoff_OL_prop = float(sys.argv[2])
segment_col = int(sys.argv[3])
minimum_length_prop = 0
for line in f:
	line = line.replace("\n", "")
	line_l = re.split("\t", line)

	read_l = re.split('\|', line_l[segment_col])

	breakpoint_info = []
	for tmp in read_l:
		tmp_l = tmp.split(";")
		breakpoint_info.append(tmp_l)

#	print(breakpoint_info)
	
#d907a4ec-64df-4bbe-bf15-09557e7333a5;split;chr10;10458836;chr10;10460901;4355;31,3176;3168,4344;10462157,10460901;10458836,10462160;17.8,17.8;60,60;-,+;-,-;0.0018369690011481056;-|b2a65c76-6cd9-4e7f-a8b2-78278789a014;split;chr10;10458903;chr10;10460894;5794;32,3966;3956,5733;10463039,10460894;10458903,10462979;17.1,17.1;60,60;-,+;-,-;0.0017259233690024164;-

	num_OL_read = 0
	for i in range(len(breakpoint_info)):
		g_start_l = breakpoint_info[i][9].split(",")
		g_end_l = breakpoint_info[i][10].split(",")
		g_start1, g_start2, g_end1, g_end2 = int(g_start_l[0]), int(g_start_l[1]), int(g_end_l[0]), int(g_end_l[1])
		if g_start1 > g_end1:
			g_start1, g_end1 = g_end1, g_start1
		if g_start2 > g_end2:
			g_start2, g_end2 = g_end2, g_start2

#		print(g_start1, g_end1, g_start2, g_end2)

		OL = 0
		if g_end1 - g_start1 < g_end2 - g_start2:
#			print(g_end1 - g_start1, g_end2 - g_start2)
			if g_start2 <= g_start1 <= g_end2 and g_start2 <= g_end1 <= g_end2:
#				print("AAAAAAAAAAAAAAAAA")
				OL = float(g_end1 - g_start1)/float(g_end1 - g_start1)
			elif g_start2 <= g_start1 <= g_end2 and g_end2 < g_end1:
#				print("BBBBBBBBBBBBBBBBBBB")
				OL = float(g_end2 - g_start1)/float(g_end1 - g_start1)
			elif g_start2 <= g_end1 <= g_end2 and g_start1 < g_start2:
#				print("CCCCCCCCCCCCCC")
				OL = float(g_end1 - g_start2)/float(g_end1 - g_start1)

		if g_end2 - g_start2 < g_end1 - g_start1:
#			print(g_end2 - g_start2, g_end1 - g_start1)
			if g_start1 <= g_start2 <= g_end1 and g_start1 <= g_end2 <= g_end1:
#				print("DDDDDDDDDDDDDDDDD")
				OL = float(g_end2 - g_start2)/float(g_end2 - g_start2)
			elif g_start1 <= g_start2 <= g_end1 and g_end1 < g_end2:
#				print("EEEEEEEEEEEEEEEEEEEEE")
				OL = float(g_end1 - g_start2)/float(g_end2 - g_start2)
			elif g_start1 <= g_end2 <= g_end1 and g_start2 < g_start1:
#				print("FFFFFFFFFFFFFFFFFF")
				OL = float(g_end2 - g_start1)/float(g_end2 - g_start2)

#		print("OL = ", OL)

		if OL > cutoff_OL_prop:
#			print(OL)
			num_OL_read += 1

	total = len(breakpoint_info)
#	print("OL=", num_OL_read, total)

	if (num_OL_read and total >= 3) or (num_OL_read == 0 and total >= 2): #180915 >0 => >= 2 #180918 max(num_OL_read, total - num_OL_read) >= 2# 180918_2 (num_OL_read and total >= 3) or (num_OL_read == 0 and total >= 2)
		print(line)
