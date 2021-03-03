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
     
## Installation and usage via Docker
Install Docker in yor computer, and run the below comamnds to install and run test.
```
git clone https://github.com/afujimoto/CAMPHORsomatic
cd <path to CAMPHORsomatic>
docker build -t camphorsomatic .
docker run --rm -it -v $PWD/sv:/CAMPHOR/SV -v $PWD/sample1:/CAMPHOR/sample1 -v $PWD/sample2:/CAMPHOR/sample2 camphorsomatic
```

If you want to run for your own data, please run the below comamnds.
```
git clone https://github.com/afujimoto/CAMPHORsomatic
cd <path to CAMPHORsomatic>
docker build -t camphorsomatic .
docker run --rm -it -v <path to directory of cancer bam>:/cancer_input -v <path to directory of normal bam>:/normal_input -v <path to output directory>:/out -v <path to output directory of cancer>:/cancer_output -v <path to output directory of normal>:/normal_output -v <path to output directory>:/output camphorsomatic sh CAMPHOR_SVcall.sh cancer_input/<name of bam file of cancer sample (sorted by read name)> cancer_input/<name of bam of cancer sample (sorted by genome coordinate)> cancer_output && sh CAMPHOR_SVcall.sh normal_input/<name bam file of normal sample (sorted by read name)> normal_input/<name bam file of normal sample (sorted by genome coordinate)> normal_output && sh CAMPHOR_comparison.sh cancer_output normal_output cancer_input/<name of bam of cancer sample (sorted by genome coordinate)> normal_input/<name bam file of normal sample (sorted by genome coordinate)> /cancer_input/<name of fastq file of cancer> out
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
    
   
Automatatic obtaining the annotation data used from UCSC  
Download and format change for repeat information can be perfoemd with the commnds below.
```
cd CAMPHORsomatic
mkdir ./data/
curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz | zcat | grep Simple_repeat | python3 ./src/repeat/rmsk.py /dev/stdin > ./data/rmsk.txt
curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz | zcat | python3 ./src/repeat/simpleRepeat.py /dev/stdin | sort -k1,1 -k2,2g > ./data/simplerepeat.txt
curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz | zcat | python3 ./src/repeat/seg_dup.py /dev/stdin | sort -k1,1 -k2,2g > ./data/seg_dup.txt
curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chainSelf.txt.gz | zcat | python3 ./src/repeat/ucsc_selfchain.py /dev/stdin | sort -k1,1 -k2,2g > ./data/chainSelf.txt
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
