#!/bin/bash
#
# Usage:
#   ./annotate-score.sh annotations.txt.gz target-file.tsv
#
# Example:
#   ./annotate-score.sh mpc-2.0.b38.txt.gz MIS-chr22.tsv | gzip -c > output.txt.gz
#
# Input format:
#   zcat annotations.txt.gz | head -1
#       1       69094   G       A       2.7340307082    1:69094:G:A
#
#   # only the first three columns are mandatory
#   zless target-file.tsv | head -2
#       varID	chrom	pos	ENST	column	AlphaMissense
#       22:15690708:A:G	22	15690708	ENST00000343518	1	0.7361
#

if [ $# != 2 ]; then
    echo "Usage:"
    echo "  ./annotate-score.sh annotations.txt.gz target-file.tsv"
    echo
    echo "Example:"
    echo "  ./annotate-score.sh mpc-2.0.b38.txt.gz MIS-chr22.tsv | gzip -c > output.txt.gz"
    echo
    echo "Input format:"
    echo "  $ zcat annotations.txt.gz | head -1"
    echo "  1       69094   G       A       2.7340307082    1:69094:G:A"
    echo
    echo "  # only the first three columns are mandatory"
    echo "  $ zless target-file.tsv | head -2"
    echo "  varID	chrom	pos	ENST	column	AlphaMissense"
    echo "  22:15690708:A:G	22	15690708	ENST00000343518	1	0.7361"
    echo
    exit;
fi

DIR=`mktemp -d rmme.XXXXXXXXXX`

annot-tsv --version >/dev/null 2>/dev/null

# Install annot-tsv if not available in PATH
if [ $? != 0 ]; then
    cd $DIR
    wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
    tar -jxf htslib-1.19.1.tar.bz2
    cd htslib-1.19.1
    ./configure
    make
    export PATH=$PATH:`realpath .`
    cd ../../
fi

ANNOTS=$1
INPUT=$2

#set -x

zless $INPUT | head -1 | sed 's,\r,,g' | awk '{OFS="\t"}{print $line,"score" }' > $DIR/dat.hdr.txt
zless $INPUT | tail -n +2 | gzip -c > $DIR/dat.txt.gz
CHR=(`zcat $DIR/dat.txt.gz | cut -f2 | sort | uniq`)
if [ ${#CHR[@]} != 1 ]; then
    echo "Expected single chromosome in $INPUT\n"
    exit;
fi

# number of columns in the output file
NCOL=`cat $DIR/dat.hdr.txt | sed 's,\t,\n,g' | wc -l`


(cat $DIR/dat.hdr.txt;
 zcat $ANNOTS | awk '$1=="'$CHR'"' | annot-tsv -t $DIR/dat.txt.gz -c 1,2,2:2,3,3 -m 6:1 -f 5:$NCOL)

rm -rf $DIR