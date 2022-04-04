#!/bin/bash

# split repeatmodeler
Rscript ~/projectDrive_2/DDE_Pipeline_2/splitter.R -f echisCarinatus-families.fa -p 256 -t DNA -o data/split

# list split files
ls data/split/echisCarinatus-families.fa_seq_* | sed 's|data/split/||' > echisCarinatus-families.fa_split_seq_list.txt

# make rps_out dir
mkdir -p data/rps_out

# rpstblastn for known domains
parallel --bar --jobs 128 -a split_seq_list.txt 'rpstblastn -query data/split/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
cat data/rps_out/echisCarinatus-families.fa_seq_* > data/echisCarinatus-families.fa.rps.out

# remove temp files
rm data/rps_out/echisCarinatus-families.fa_seq_* data/split/echisCarinatus-families.fa_seq_*

# run tandem repeat finder
trf echisCarinatus-families.fa 2 7 7 80 10 50 500 -d -h

# convert trf to gff
python trf2gff.py -d echisCarinatus-families.fa.2.7.7.80.10.50.500.dat -o data/echisCarinatus-families.fa.trf.gff
rm echisCarinatus-families.fa.2.7.7.80.10.50.500.dat