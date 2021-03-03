FROM ubuntu:20.04

RUN apt-get update \
  && apt-get -y --no-install-recommends install \
    curl \
    python3 \
    python3-pip \
    samtools \
  && apt-get -y clean \
  && rm -rf /var/lib/apt/lists/*

RUN pip3 install numpy pysam

RUN curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz | zcat | grep Simple_repeat | gzip -c > /rmsk.txt.gz
RUN curl -LO http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
RUN curl -LO http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
RUN curl -LO http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chainSelf.txt.gz

RUN mkdir -p /CAMPHOR/data
COPY ./ /CAMPHOR

WORKDIR /CAMPHOR

RUN zcat /rmsk.txt.gz | python3 ./src/repeat/rmsk.py /dev/stdin > data/rmsk.txt
RUN zcat /simpleRepeat.txt.gz | python3 ./src/repeat/simpleRepeat.py /dev/stdin | sort -k1,1 -k2,2g > data/simplerepeat.txt
RUN zcat /genomicSuperDups.txt.gz | python3 ./src/repeat/seg_dup.py /dev/stdin | sort -k1,1 -k2,2g > data/seg_dup.txt
RUN zcat /chainSelf.txt.gz | zcat | python3 ./src/ucsc_selfchain.py /dev/stdin | sort -k1,1 -k2,2g > data/chainSelf.txt

CMD sh CAMPHOR_SVcall.sh ./example/sample1.sort_by_name.test.bam ./example/sample1.sort.test.bam sample1 \
  && sh CAMPHOR_SVcall.sh ./example/sample2.sort_by_name.test.bam ./example/sample2.sort.test.bam sample2 \
  && sh CAMPHOR_comparison.sh ./sample1 ./sample2 ./example/sample1.sort.test.bam ./example/sample2.sort.test.bam ./example/sample1.sort.test.fastq SV
