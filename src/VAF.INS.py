import sys
import subprocess
import re
import pysam

#chr10   1293159 chr10   1293297 2       DEL     139     NA      2       NS;2/DIS;2/MQ;2/CN;2/CI;2       0.8/12/176/0.69 0.0/NA/0/NA     1,1     NA     #chr10   1301507 chr10   1301684 2       DEL     178     NA      2       NS;2/DIS;2/MQ;2/CN;2/CI;2       0.0/NA/0/NA     0.27/NA/0/NA    1,1     NA      
#def get_depth(chr, pos, bam):
#	cmd = "samtools view " + bam + " " + chr + ":" + pos + "-" + pos
#	contents = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
#	contents = str(contents)
#	contents = contents.replace("b\'", "")
#	contents_l = re.split("\\\\n", contents)
#	num = len(contents_l)
#	return num

bamfile = pysam.AlignmentFile(sys.argv[2], "rb")
min_VAF = float(sys.argv[3])
min_mq = int(sys.argv[4])
prop_low_q_reads = float(sys.argv[5])
read_num_col = int(sys.argv[6])
slice_col = int(sys.argv[7])

f = open(sys.argv[1])
for line in f:
	line = line.replace("\n", "")
	line_l = line.split("\t")
	
#	depth1 = get_depth(line_l[0], str(int(line_l[1]) - 100), sys.argv[2])
#	depth2 = get_depth(line_l[2], str(int(line_l[3]) + 100), sys.argv[2])

#	print(line_l)

	depth1 = 0
	low_mq_read1 = 0
#	print(line_l[0], int(line_l[1]) - 100, int(line_l[1]) - 99)
	for read in bamfile.fetch(line_l[0], int(line_l[1]) - 100, int(line_l[1]) - 99):
#		print(read)
		if read.mapping_quality < min_mq:
#			print(read)
			low_mq_read1 += 1
		depth1 += 1

	depth2 = 0
	low_mq_read2 = 0
	for read in bamfile.fetch(line_l[2], int(line_l[3]) + 99, int(line_l[3]) + 100):
		if read.mapping_quality < min_mq:
#			print(read)
			low_mq_read2 += 1
		depth2 += 1

#	print(low_mq_read1, depth1, low_mq_read2, depth2)

#	print("depth1 = ", depth1, depth1_2)
#	print("depth2 = ", depth2, depth2_2)

	VAF = round(float(line_l[read_num_col])/(float(depth1 + depth2)/2), 3)

	prop_low_q_reads_1, prop_low_q_reads_2 = 0, 0
	if depth1 > 0:
		prop_low_q_reads_1 = float(low_mq_read1)/float(depth1)
	if depth2 > 0:
		prop_low_q_reads_2 = float(low_mq_read2)/float(depth2)

	if VAF >= min_VAF and prop_low_q_reads_1 <= prop_low_q_reads and prop_low_q_reads_2 <= prop_low_q_reads:
		low_q_raed_prop = str(low_mq_read1) + "," + str(depth1) + "/" + str(low_mq_read2) + "," + str(depth2)
		print("\t".join(line_l[0:slice_col]), str(depth1) + "," + str(depth2), str(VAF), low_q_raed_prop, "\t".join(line_l[slice_col:len(line_l)]), sep="\t")
#	else:
#		low_q_raed_prop = str(low_mq_read1) + "," + str(depth1) + "/" + str(low_mq_read2) + "," + str(depth2)
#		print("#############", low_q_raed_prop, VAF, line)
