# CAMPHORsomatic
somatic SV caller for long-reads

Overview
1. Idetify SV candidates from bam files for cancer and matched normal samples　　
2. Compare SV candidates in cancer and matched normal samples and indetify somatic SV candidates     　　

## Requirement
python3       
pysam module of python

perl

samtools (0.1.18 or higher)

## Input file
** Two bam files (bam sorted by read name and bam sorted by genome coordinate)        
** Index file (.bai) for bam file sorted by genome coordinate         
** Fastq file of the sequence data

If samtools is not installed in the environment, the path to execution file of samtools can be specified in the config file (pram.config).  

## Output file format
vcf file of SVs (SV.vcf)

## Usage
```
cd <path to CAMPHOR>　　
sh CAMPHOR.sh <bam of normal sample (sorted by read name)> <bam of normal sample (sorted by genome coordinate)> <fastq> <output>　　
```

## Example
```
git clone https://github.com/afujimoto/CAMPHOR.git　
cd CAMPHOR/CAMPHOR　　
sh CAMPHOR.sh ./example/NA18943.chr22.sort_by_name.test.bam ./example/NA18943.chr22.sort.test.bam ./example/NA18943.chr22.sort.test.fastq test
```

## Parameter setting in configuration file
We consider the parameter set in the provided configuration apprppreate for 20x coverage WGS data.  
If you would like to use different parameters, please make changes in the parm.config file.　　   
In the current file, minimum variant alele frequency (MIN_VAF), minimu number of reads (MIN_READ_NUMBER) and minimum indel length (MIN_INDEL_LENGTH) (bp) are set ot 0.15, 2 and 100.　　  
For greater depth of coverage, larger minimu number of reads sould be apprppreate. 　　  

We developed this method with nanopore sequence data basecalled by albacore (total error rate =~ 15%), and set minimum indel length to 100bp to remove false potitives. But newer basecallers increase the accuracy and, smaller minimum indel length (50bp or smaller) and lower minimum variant alele frequency (0.1 or lower) can be used.

## Repeat filtering
If repeat infmarmaiton files are privided, our method filter SV candisates with repeat infromation (Repeat masker, Tandm repeat finder, Segmental duplication, Self-chain). This filter can increase specitficity.
Please prepare anntaiton files with the fllowing procedures.

Repeat masker　　     
Download rmsk.txt from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
```
grep Simple_repeat <path to rmsk.txt>|python .src/repeat/rmsk.py /dev/stdin > ./data/rmsk.txt
```

Tandem repeat     
Download simpleRepeat.txt from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/　　
```
python ./src/repeat/simpleRepeat.py <path to simpleRepeat.txt>|sort -k1,1 -k2,2g > ./data/simplerepeat.txt　　
```
  
Segmental duplication　　     
Download genomicSuperDups.txt from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/　　
```
python ./src/repeat/seg_dup.py <path to genomicSuperDups.txt>|sort -k1,1 -k2,2g > ./data/seg_dup.txt
```

Self-chain　　     
Download chainSelf.txt file from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/　　
```
python .src/repeat/ucsc_selfchain.py <path to chainSelf.txt> | sort -k1,1 -k2,2g > ./data/chainSelf.txt
```

## Preformance
Performance of this tool is provided in Fujimoto et al. (in rivision).

## Licence
GPL

## Contact

Akihiro Fujimoto - afujimoto@m.u-tokyo.ac.jp

http://www.humgenet.m.u-tokyo.ac.jp/index.en.html
