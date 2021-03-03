# CAMPHORsomatic
somatic SV caller for long-reads

Overview
1.	Identify SV candidates from bam files for cancer and matched normal samples　　　　
2.	Compare SV candidates between cancer and matched normal samples and identify somatic SV candidates　　　　　　
  　　
## Requirement
python3       
pysam module of python

perl

samtools (0.1.18 or higher)

If samtools is not installed in the environment, the path to the execution file of samtools can be specified within the config file (pram.config).     

## Input file
** Two bam files (one bam sorted by read name and another sorted by genome coordinates) for cancer and matched-normal samples        
** Index file (.bai) for bam files sorted by genome coordinates         
** Fastq file of cancer


## Output file format
vcf file of SVs (somatic_SV.vcf)

## Usage
```
cd <path to CAMPHOR>　　
sh CAMPHOR_SVcall.sh <bam of cancer sample (sorted by read name)> <bam of cancer sample (sorted by genome coordinate)> <output directory of cancer>   
sh CAMPHOR_SVcall.sh <bam of normal sample (sorted by read name)> <bam of normal sample (sorted by genome coordinate)> <output directory of normal>   
sh CAMPHOR_comparison.sh <output directory of cancer> <output directory of normal> <bam of cancer sample (sorted by genome coordinate)> <bam of normal sample (sorted by genome coordinate)> <fastq file of cancer> <output directory of somatic SV>   
```

## Example
```
git clone https://github.com/afujimoto/CAMPHORsomatic　       
cd CAMPHORsomatic      　　
sh CAMPHOR_SVcall.sh ./example/sample1.sort_by_name.test.bam ./example/sample1.sort.test.bam sample1   
sh CAMPHOR_SVcall.sh ./example/sample2.sort_by_name.test.bam ./example/sample2.sort.test.bam sample2
sh CAMPHOR_comparison.sh ./sample1 ./sample2 ./example/sample1.sort.test.bam ./example/sample2.sort.test.bam ./example/sample1.sort.test.fastq SV       

```

## Parameter setting in configuration file
We consider the parameter set in the provided configuration appropriate for 20x coverage WGS data.  

We developed this method with nanopore sequence data basecalled by albacore (total error rate =~ 15%), and set minimum indel length to 100bp to remove false positives. But newer basecallers have increased accuracy and, smaller minimum indel length (50bp or smaller) can be used. For this, users can change the "MIN_INDEL_LENGTH" within the pram.config file.

## Repeat filtering
Our method filters SV candidates with the provided repeat information (Repeat masker, Tandem repeat finder, Segmental duplication, Self-chain).        
Please prepare annotation files with the following procedures.     


Make a directory for repeat files in CAMPHORsomatic directory　　
```
mkdir data
```

Repeat masker　　     
Download rmsk.txt from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
```
grep Simple_repeat <path to rmsk.txt>|python ./src/repeat/rmsk.py /dev/stdin > ./data/rmsk.txt
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
python ./src/repeat/ucsc_selfchain.py <path to chainSelf.txt> | sort -k1,1 -k2,2g > ./data/chainSelf.txt
```
## Data of noramal samples
CAMPHOR_comparison.sh compares cancer SVs and normal SV candidates, and removes germline SVs. For this comparison, <SV type>_candidate.txt0 files in <output directory of normal> are used. Users can merge these SV files of other normal samples, and save the same name in a new directory. The new directory can be used as <output directory of normal> in analysis with CAMPHOR_comparison.sh. This analysis increases power to remove germline SVs.

## Performance
False positive rate was estimated to be ~7% with PCR in a liver cancer sample set (Fujimoto et al. in revision).


## Licence
GPL

## Contact

Akihiro Fujimoto - afujimoto@m.u-tokyo.ac.jp

http://www.humgenet.m.u-tokyo.ac.jp/index.en.html
