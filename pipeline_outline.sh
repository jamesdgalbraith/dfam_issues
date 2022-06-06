#!/bin/bash

usage() { echo "Usage: $0 [-q query FASTA] [-t threads to use]" 1>&2; exit 1; }

threads=1
while getopts 'q:t:' flag; do
  case "${flag}" in
    q) query=${OPTARG} ;;
    t) threads=${OPTARG} ;;
    *) usage
       exit 1 ;;
  esac
done

if [ -z "${query}" ]; then
    usage
    exit 1
fi

if ! [[ "$threads" =~ ^[0-9]+$ ]] ; then
  echo "Error: threads must be a positive, whole number."
  exit 1
elif [ "$threads" -eq 0 ]; then
  echo "Error: threads must be a positive, whole number."
  exit 1
fi

query_name=$( echo $query | sed 's/.*\///' )

export query
export threads

# Split fasta
echo Splitting
Rscript scripts/splitter.R -f ${query} -p 1280 -t DNA -o data/split

# List split files
ls data/split/${query_name}_seq_* | sed 's|data/split/||' > ${query_name}_split_seq_list.txt

# rpstblastn for known domains
echo Blasting for domains
mkdir -p data/rps_out
parallel --bar --jobs $threads -a ${query_name}_split_seq_list.txt 'rpstblastn -query data/split/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
cat data/rps_out/${query_name}_seq_* > data/${query_name}.rps.out

# dc-megablast against RepBase
echo Blasting against RepBase
mkdir -p data/repbase_out
parallel --bar --jobs $threads -a ${query_name}_split_seq_list.txt 'blastn -query data/split/{} -db /media/projectDrive_1/databases/repbase/RepBase5May2021.fasta -out data/repbase_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs" -task dc-megablast -num_threads 1'
cat data/repbase_out/${query_name}_seq_* > data/${query_name}.repbase.out

# Run TRF
mkdir -p data/trf_out
parallel --bar --jobs $threads -a ${query_name}_split_seq_list.txt 'trf data/split/{} 2 7 7 80 10 50 500 -d -h'

# convert trf to gff
parallel --bar --jobs $threads -a ${query_name}_split_seq_list.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/trf_out/{}.trf.gff'

# compile trf
cat data/trf_out/${query_name}_seq_*.trf.gff > data/${query_name}.trf.gff

# delete temp files
find . -mindepth 2 -name "${query_name}_seq_*" -delete
find . -name "${query_name}*.2.7.7.80.10.50.500.dat" -delete
rm ${query_name}_split_seq_list.txt

# run clean up Rscript
Rscript cleanup.R -file ${query}

echo "all done"