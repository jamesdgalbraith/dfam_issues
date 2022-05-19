#!/bin/bash

export QUERY

CONTIGS=`grep -c ">" seq/${QUERY}`

# Split fasta
echo Splitting
Rscript scripts/splitter.R -f seq/${QUERY} -p ${CONTIGS} -t DNA -o data/split

# List split files
ls data/split/${QUERY}_seq_* | sed 's|data/split/||' > ${QUERY}_split_seq_list.txt

# rpstblastn for known domains
echo Blasting for domains
mkdir -p data/self_blast
parallel --bar --jobs 128 -a ${QUERY}_split_seq_list.txt 'blastn -query data/split/{} -subject data/split/{} -out data/self_blast/{}.self_blast.out -outfmt "6 std slen" -evalue 1e-50 -task dc-megablast'

find . -type f -name "${QUERY}_seq_*.self_blast.out" -exec cat {} + > data/${QUERY}.self_blast.out

# delete temp files
find . -mindepth 2 -name "${QUERY}_seq_*" -delete
rm ${QUERY}_split_seq_list.txt
