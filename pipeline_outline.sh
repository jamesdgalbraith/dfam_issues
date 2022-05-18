#!/bin/bash

usage() { echo "Usage: $0 [-q query FASTA]" 1>&2; exit 1; }

while getopts 'q:T' flag; do
  case "${flag}" in
    q) QUERY=${OPTARG} ;;
    # t) THREADS=${OPTARG} ;;
    *) usage
       exit 1 ;;
  esac
done

if [ -z "${QUERY}" ]; then
    usage
    exit 1
fi

export QUERY

# Split fasta
echo Splitting
Rscript scripts/splitter.R -f seq/${QUERY} -p 1280 -t DNA -o data/split

# List split files
ls data/split/${QUERY}_seq_* | sed 's|data/split/||' > ${QUERY}_split_seq_list.txt

# rpstblastn for known domains
echo Blasting for domains
mkdir -p data/rps_out
parallel --bar --jobs 124 -a ${QUERY}_split_seq_list.txt 'rpstblastn -query data/split/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
cat data/rps_out/${QUERY}_seq_* > data/${QUERY}.rps.out

# dc-megablast against RepBase
echo Blasting against RepBase
mkdir -p data/repbase_out
parallel --bar --jobs 124 -a ${QUERY}_split_seq_list.txt 'blastn -query data/split/{} -db /media/projectDrive_1/databases/repbase/RepBase5May2021.fasta -out data/repbase_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs" -task dc-megablast -num_threads 1'
cat data/repbase_out/${QUERY}_seq_* > data/${QUERY}.repbase.out

# Run TRF
mkdir -p data/trf_out
parallel --bar --jobs 124 -a ${QUERY}_split_seq_list.txt 'trf data/split/{} 2 7 7 80 10 50 500 -d -h'

# convert trf to gff
parallel --bar --jobs 124 -a ${QUERY}_split_seq_list.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/trf_out/{}.trf.gff'

# compile trf
cat data/trf_out/${QUERY}_seq_*.trf.gff > data/${QUERY}.trf.gff

# delete temp files
find . -mindepth 2 -name "${QUERY}_seq_*" -delete
find . -name "${QUERY}*.2.7.7.80.10.50.500.dat" -delete
rm ${QUERY}_split_seq_list.txt

# run clean up Rscript
Rscript cleanup.R -file ${QUERY}

echo "all done"