import sys

##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">
##FILTER=<ID=STRANDBIAS,Description="Strand is biased if Strandbias_pval< 0.01.">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /archive/data/amed_snt/WORK/james/NA19240_SV_calling_results/20201203_minimap2_v2.17_alignments/NA1
#9240_all.fastq.0.minimap2.sort.bam
#chr1    66239   0       ATTATATTATATAATATATAATATAAATATAATATAA   N       .       PASS    PRECISE;SVMETHOD=Snifflesv1.0.12;CHR2=chr1;END=66276;STD_quant_start=6.164414;STD_quant_stop=5.412947;Kurtosis_quant_start=0.221510;Kurtosis_quant_stop=1.767848;SVTYPE=DEL;SUPTYPE=AL;SVLEN=-37;STRANDS=+-;STRANDS2=9,7,9,7;RE=16;REF_strand=11,13;Strandbias_pval=0.54037;AF=0.666667 GT:DR:DV        0/1:8:16
#chr1    90388   1       N       CAGTTAACCTGCTGCTTCCTGGAGAGGGACAGTCCCTCACAACCCTCTGTCTCTTGCCAAG   .       PASS    IMPRECISE;SVMETHOD=Snifflesv1.0.12;CHR2=chr1;END=90396;STD_quant_start=23.171103;STD_quant_stop=22.329353;Kurtosis_quant_start=-0.096051;Kurtosis_quant_stop=-0.945865;SVTYPE=INS;SUPTYPE=AL;SVLEN=56;STRANDS=+-;STRANDS2=7,7,7,7;RE=14;REF_strand=9,14;Strandbias_pval=0.733171;AF=0.608696    GT:DR:DV        0/1:9:14

#chr1    904500  chr1    904606  3       DEL     103     2,2,2   6,6     0.5     0,6/0,6 3       0       1/2     |||     ERR3219856.808362;within;chr1;904481;chr1;904581;23
#chr1    909110  chr1    909514  2       DEL     417     2,2,2   10,10   0.2     0,10/0,10       2       0       2/0     |||     ERR3219854.4807567;within;chr1;909101;chr1;
#chr1    934105  chr1    934942  6       DEL     1000    4,4,4   20,23   0.279   0,20/0,23       6       0       1/3     |||     ERR3219853.2509625;within;chr1;934067;chr1;

##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency.">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INVDUP,Description="InvertedDUP with unknown boundaries">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">


print('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency.">')
print('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">')
print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">')
print('##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">')
print('##INFO=<ID=ND,Number=1,Type=Integer,Description="number of reads in matched-normal sample (average number of reads at 100bp upstream and downstream positions from breakpoints)">')
print('##ALT=<ID=DEL,Description="Deletion">')
print('##ALT=<ID=INS,Description="Insertion">')
print('##ALT=<ID=INV,Description="Inversion">')
print('##ALT=<ID=TRA,Description="Translocation">')
print('##ALT=<ID=CHR,Description="Chromosomal Translocation">')
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

ID=0

DEL_f = open(sys.argv[1])
for line in DEL_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if float(line_l[13]) > 1:
		line_l[13] = "1"

	depth_l = line_l[8].split(",")
	depth = (int(depth_l[0]) + int(depth_l[1]))/2

	INFO=";".join(["CHR2="+line_l[2], "END="+line_l[3], "AF="+line_l[13], "RE="+line_l[4], "LEN="+str(int(line_l[3]) - int(line_l[1])), "ND="+str(depth)])
	print(line_l[0], line_l[1], ID, "N", "<DEL>", ".", "PASS", INFO, sep="\t")

	ID += 1

#chr1    820880  chr1    820880  13      INS     241.0   13,0,0  11,11,11        15,14   0.897   1,15/0,14       ERR3219857.2119650;within
#chr1    876001  chr1    876001  5       INS     1921.0  3,2,0   3,3,3   23,24   0.17    1,23/1,24       ERR3219854.5592732;within;chr1;87

INS_f = open(sys.argv[2])
for line in INS_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if float(line_l[10]) > 1:
		line_l[10] = "1"

	depth_l = line_l[9].split(",")
	depth = (int(depth_l[0]) + int(depth_l[1]))/2

	if ">" in line_l[6]:
		INFO=";".join(["CHR2="+line_l[2], "END="+line_l[3], "AF="+line_l[10], "RE="+line_l[4], "LEN="+line_l[6], "ND="+str(depth)])
	else:
		INFO=";".join(["CHR2="+line_l[2], "END="+line_l[3], "AF="+line_l[10], "RE="+line_l[4], "LEN="+str(int(float(line_l[6]))), "ND="+str(depth)])
	print(line_l[0], line_l[1], ID, "N", "<INS>", ".", "PASS", INFO, sep="\t")

	ID += 1

#chr20   36999458        chr20   47554552        3       DIF     10555082        2,2,2   17,21   0.158   0,17/1,21       3       ERR3219857.2782840;split;chr20;36999453;ch
#chr20   37001800        chr20   47554552        4       DIF     10552771        3,3,3   24,21   0.178   4,24/1,21       4       ERR3219853.1365735;split;chr20;37001781;ch


INV_f = open(sys.argv[3])
for line in INV_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if float(line_l[12]) > 1:
		line_l[12] = "1"

	depth_l = line_l[9].split(",")
	depth = (int(depth_l[0]) + int(depth_l[1]))/2

	INFO=";".join(["CHR2="+line_l[2], "END="+line_l[3], "AF="+line_l[2], "RE="+line_l[4], "LEN="+str(int(float(line_l[6]))), "ND="+str(depth)])
	print(line_l[0], line_l[1], ID, "N", "<INV>", ".", "PASS", INFO, sep="\t")

	ID += 1

TRA_f = open(sys.argv[4])
for line in TRA_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if float(line_l[12]) > 1:
		line_l[12] = "1"

	depth_l = line_l[9].split(",")
	depth = (int(depth_l[0]) + int(depth_l[1]))/2

	INFO=";".join(["CHR2="+line_l[2], "END="+line_l[3], "AF="+line_l[12], "RE="+line_l[4], "LEN="+str(int(float(line_l[6]))), "ND="+str(depth)])
	print(line_l[0], line_l[1], ID, "N", "<TRA>", ".", "PASS", INFO, sep="\t")

	ID += 1

#chr21   10414914        chr7    152401119       26      CHR     NA      25,25,25        23,65   0.591   3,23/2,65       26      ERR3219853.698600;split;chr21;10414914;chr
#chr22   11068963        chr7    152414188       8       CHR     NA      8,8,8   73,32   0.152   0,73/3,32       8       ERR3219855.734408;split;chr22;11068963;chr7;152414

CHR_f = open(sys.argv[5])
for line in CHR_f:
	line = line.replace("\n", "")
	line_l = line.split("\t")

	if float(line_l[11]) > 1:
		line_l[11] = "1"

	depth_l = line_l[8].split(",")
	depth = (int(depth_l[0]) + int(depth_l[1]))/2

	INFO=";".join(["CHR2="+line_l[2], "END="+line_l[3], "AF="+line_l[11], "RE="+line_l[4], "LEN=NA", "ND="+str(depth)])
	print(line_l[0], line_l[1], ID, "N", "<CHR>", ".", "PASS", INFO, sep="\t")

	ID += 1
