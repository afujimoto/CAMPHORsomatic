# CAMPHORsomatic
somatic SV caller for long-reads

Overview
1. Idetify SV candidates from bam files for cancer and matched normal samples　　
2. Compare SV candidates between cancer and matched normal samples and indetify somatic SV candidates     　　

## Requirement
python3       
pysam module of python

perl

samtools (0.1.18 or higher)

## Input file
** Two bam files (one bam sorted by read name and another sorted by genome coordinates)        
** Index file (.bai) for bam file sorted by genome coordinate         
** Fastq file of the sequence data

If samtools is not installed in the environment, the path to the execution file of samtools can be specified within the config file (pram.config).  

## Output file format
vcf file of SVs (somatic_SV.vcf)

## Usage
```
cd <path to CAMPHOR>　　
sh CAMPHOR_SVcall.sh <bam of cancer sample (sorted by read name)> <bam of cancer sample (sorted by genome coordinate)> <output directory of cancer>   
sh CAMPHOR_SVcall.sh <bam of normal sample (sorted by read name)> <bam of normal sample (sorted by genome coordinate)> <output directory of normal>   
sh CAMPHOR_comparison.sh <output directory of cancer> <output directory of normal> <bam of cancer sample (sorted by genome coordinate)> <bam of normal sample (sorted by genome coordinate)> <fastq file of cancer> <output directory>   
```

## Example
```
git clone https://github.com/afujimoto/CAMPHOR.git　
cd CAMPHOR/CAMPHOR　　
sh CAMPHOR.sh ./example/NA18943.chr22.sort_by_name.test.bam ./example/NA18943.chr22.sort.test.bam ./example/NA18943.chr22.sort.test.fastq test
```

## Parameter setting in configuration file
We consider the parameter set in the provided configuration appropriate for 20x coverage WGS data.  

We developed this method with nanopore sequence data basecalled by albacore (total error rate =~ 15%), and set minimum indel length to 100bp to remove false positives. But newer basecallers have increased accuracy and, smaller minimum indel length (50bp or smaller) can be used. For this, users can change the "MIN_INDEL_LENGTH" within pram.config file.

## Repeat filtering
Our method filters SV candisates with the provided repeat infromation (Repeat masker, Tandm repeat finder, Segmental duplication, Self-chain).     
Please prepare anntaiton files with the following procedures.       


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
Performance of this tool is provided in Fujimoto et al. (in revision).

## Licence
GPL

## Contact

Akihiro Fujimoto - afujimoto@m.u-tokyo.ac.jp

http://www.humgenet.m.u-tokyo.ac.jp/index.en.html
