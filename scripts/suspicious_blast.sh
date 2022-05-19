#!/bin/bash

export QUERY

CONTIGS=`grep -c ">" seq/${QUERY}`

# Split fasta
echo Splitting
Rscript scripts/splitter.R -f seq/${QUERY} -p ${CONTIGS} -t DNA -o data/split

# List split files
ls data/split/${QUERY}_seq_* | sed 's|data/split/||' > ${QUERY}_split_seq_list.txt

# rpstblastn for known domains
echo blastxing
mkdir -p data/suspicious_blast
parallel --bar --jobs 128 -a ${QUERY}_split_seq_list.txt 'blastx -query data/split/{} -db ~/databases/new_NR/nr -out data/suspicious_blast/{}.suspicious_blast.out -outfmt "6 std slen"'
> data/${QUERY}.self_blast.out

find . -type f -name "${QUERY}_seq_*.suspicious_blast.out" -exec cat {} + > data/${QUERY}.suspicious_blast.out

# delete temp files
find . -mindepth 2 -name "${QUERY}_seq_*" -delete
rm ${QUERY}_split_seq_list.txt
