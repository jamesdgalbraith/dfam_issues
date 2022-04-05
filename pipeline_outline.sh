#!/bin/bash

export QUERY

# Split fasta
echo Splitting
Rscript scripts/splitter.R -f seq/${QUERY} -p 2560 -t DNA -o data/split

# List split files
ls data/split/${QUERY}_seq_* | sed 's|data/split/||' > ${QUERY}_split_seq_list.txt

# rpstblastn for known domains
echo Blasting for domains
mkdir -p data/rps_out
parallel --bar --jobs 128 -a ${QUERY}_split_seq_list.txt 'rpstblastn -query data/split/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
cat data/rps_out/${QUERY}_seq_* > data/${QUERY}.rps.out

# dc-megablast against RepBase
echo Blasting against RepBase
mkdir -p data/repbase_out
parallel --bar --jobs 128 -a ${QUERY}_split_seq_list.txt 'blastn -query data/split/{} -db /media/projectDrive_1/databases/repbase/RepBase5May2021.fasta -out data/repbase_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs" -task dc-megablast -num_threads 1'
cat data/repbase_out/${QUERY}_seq_* > data/${QUERY}.repbase.out

# Run TRF
mkdir -p data/trf_out
parallel --bar --jobs 128 -a ${QUERY}_split_seq_list.txt 'trf data/split/{} 2 7 7 80 10 50 500 -d -h'

# convert trf to gff
parallel --bar --jobs 128 -a ${QUERY}_split_seq_list.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/trf_out/{}.trf.gff'

# compile trf
cat data/trf_out/${QUERY}_seq_*.trf.gff > data/${QUERY}.trf.gff

# remove temp files
rm  data/split/${QUERY}_seq_* ${QUERY}*.2.7.7.80.10.50.500.dat data/trf_out/${QUERY}_seq_* ${QUERY}_split_seq_list.txt data/rps_out/${QUERY}_seq_*txt

echo "all done"