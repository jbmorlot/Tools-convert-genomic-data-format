#!/bin/bash

DIR=$(dirname $0)
bedtools genomecov -i $1.tagAlign -bg -g $DIR/hg19.chrom.sizes.txt > $1.bedGraph
bedSort $1.bedGraph $1.2.bedGraph
./wigToBigWig $1.2.bedGraph $DIR/hg19.chrom.sizes.txt $1.bigWig
rm -f $1.bedGraph $1.2.bedGraph
