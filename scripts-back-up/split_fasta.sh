#!/bin/bash
FILES=$1
DIR=$2
mkdir -p $DIR
for file in $FILES
do
echo $FILES
while read line
do
echo Line: $line
    if [[ $line == '' ]]; then
        echo "space"
    elif [[ ${line:0:1} == '>' ]]; then
        echo "header"
        outfile=${DIR}/${file}.${line#>}.fasta
        echo $line > $outfile
    else
        echo "contig"
        echo $line >> $outfile
    fi
done
done
