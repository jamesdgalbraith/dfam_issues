setwd("~/projectDrive_1/earlgrey_vipers/cleanup/")

largest_component$no <- 1:nrow(largest_component)
largest_component[largest_component$te_width>=13000,]

library(tidyverse)
library(plyranges)
library(BSgenome)
library(svglite)

# most common repeats
eg_gff <- read_gff("new_cleanup/echisCarinatus.filteredRepeats.sorted.gff")
eg_gff$ID <- tolower(gsub("\"", "", eg_gff$ID))
eg_gff$Tstart <- as.integer(gsub("\"", "", eg_gff$Tstart))
eg_gff$Tend <- as.integer(gsub("\"", "", eg_gff$Tend))

largest_component <- tibble(indiv_width = width(eg_gff), ID = eg_gff$ID) %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(total_width = sum(indiv_width)) %>%
  dplyr::ungroup() %>%
  dplyr::select(ID, total_width) %>%
  base::unique() %>%
  dplyr::arrange(-total_width)

largest_component

# read in library
combined_library <- readDNAStringSet("new_cleanup/echisCarinatus_combined_library.fasta")

# extract curated repeats
combined_library_tbl <- tibble(seqnames = names(combined_library), te_width = width(combined_library)) %>%
  mutate(ID  = tolower(sub("#.*", "", seqnames)), class = sub(" .*", "", sub(".*#", "", seqnames))) %>%
  tidyr::separate(class, into = c("class", "subclass"), sep = "/")

# read in genome
genome_seq <- readDNAStringSet(filepath = "../cleaned_up_rm_new/Echis_ycKpl.FINAL.fasta")
names(genome_seq)

largest_component <- largest_component %>%
  inner_join(combined_library_tbl)

for(i in 1:50){
  
  to_plot <- as_tibble(eg_gff[eg_gff$ID == largest_component$ID[i]]) %>%
    arrange(-score)
  
  if(nrow(to_plot) > 1000){to_plot <- to_plot[1:1000,]}
  
  ggplot() + geom_segment(data = to_plot, aes(x = Tstart, xend = Tend, y = 1:nrow(to_plot), yend = 1:nrow(to_plot))) +
    geom_segment(aes(x = c(0,combined_library_tbl[combined_library_tbl$ID == largest_component$ID[i],]$te_width),
                     xend = c(0, combined_library_tbl[combined_library_tbl$ID == largest_component$ID[i],]$te_width),
                     y = c(0, 0), yend = c(nrow(to_plot), nrow(to_plot))), colour = "red") +
    ggtitle(combined_library_tbl[combined_library_tbl$ID == largest_component$ID[i],]$seqnames) +
    scale_x_continuous(limits = c(0, combined_library_tbl[combined_library_tbl$ID == largest_component$ID[i],]$te_width))
  
  
  ggsave(paste0("new_cleanup/cleanup_plots/", largest_component$ID[i], ".svg"))
  
  
  
  to_plot_ranges <- as_granges((to_plot %>% mutate(start = start - 20, end = end + 20) %>% dplyr::select(-width)))
  to_plot_ranges <- to_plot_ranges[1:50]  
  to_plot_seq <- getSeq(genome_seq, to_plot_ranges)
  
  names(to_plot_seq) <- paste0(seqnames(to_plot_ranges), ":", ranges(to_plot_ranges), "(", strand(to_plot_ranges), ")")
  
  to_plot_seq <- c(combined_library[names(combined_library) == combined_library_tbl[combined_library_tbl$ID == largest_component$ID[i],]$seqnames],
                   to_plot_seq)
  
  writeXStringSet(to_plot_seq, paste0("new_cleanup/to_align/", largest_component$ID[i], ".fasta"))
  
  system(paste0("mafft --localpair --thread 32 new_cleanup/to_align/", largest_component$ID[i], ".fasta > new_cleanup/aligned/", largest_component$ID[i], ".fasta"))
}

# remove repeats which have been curated
# cleaned_up <- readDNAStringSet("cleaned_up_50.fasta")
# cleaned_up <- NULL
unclean_combined_library <- combined_library[!sub("#.*", "", names(combined_library)) %in% sub("#.*", "", names(cleaned_up))]

# read in trf
trf_out_ranges <- read_gff3("new_cleanup/crotalusTigris_de_novo_repeat_library_iter4.trf.gff")

# convert data types, remove manually curated, calculate length of satellites, ensure over 2.5 copies
large_period_trf_out <- as_tibble(trf_out_ranges) %>%
  mutate(seqnames = as.character(seqnames)) %>%
  # filter((grepl("rnd", seqnames) | startsWith(seqnames, "DR")),
  #        !sub("#.*", "", seqnames) %in% sub("#.*", "", names(cleaned_up))) %>%
  mutate(copies = as.double(copies),
         period = as.double(period),
         consensus_size = as.double(consensus_size),
         length = ceiling(period*copies)) %>%
  dplyr::filter(period > 20, copies >=3) %>%
  dplyr::select(seqnames, start, end, period, copies, length)

inner_join(large_period_trf_out,
           (combined_library_tbl %>% mutate(seqnames = sub(" .*", "", seqnames)))) %>%
  filter(length >= 0.25*te_width)

trf_out <- as_tibble(trf_out_ranges) %>%
  mutate(seqnames = as.character(seqnames)) %>%
  # filter((grepl("rnd", seqnames) | startsWith(seqnames, "DR")),
  #        !sub("#.*", "", seqnames) %in% sub("#.*", "", names(cleaned_up))) %>%
  mutate(copies = as.double(copies),
         period = as.double(period),
         consensus_size = as.double(consensus_size),
         length = ceiling(period*copies)) %>%
  dplyr::filter(copies > 3) %>%
  dplyr::select(seqnames, start, end, period, copies, length)

combined_library_tbl$names <- sub(" .*", "", combined_library_tbl$seqnames)

combined_library_tbl %>% mutate

likely_satellite_tbl <- trf_out %>%
  dplyr::rename(names = seqnames) %>%
  inner_join(combined_library_tbl) %>%
  as_granges() %>%
  reduce() %>%
  as_tibble() %>%
  mutate(seqnames = as.character(seqnames)) %>%
  inner_join(combined_library_tbl) %>%
  group_by(seqnames) %>%
  mutate(sum_width = sum(width), start = 1) %>%
  ungroup() %>%
  filter(sum_width >= 0.2 * te_width) %>%
  dplyr::select(names, start, te_width) %>%
  dplyr::rename(seqnames = names, end = te_width) %>%
  base::unique()

# separate out likely satellites
likely_satellite <- unclean_combined_library[sub(" .*", "", names(unclean_combined_library)) %in% likely_satellite_tbl$seqnames,]
likely_not_satellite <- unclean_combined_library[!names(unclean_combined_library) %in% names(likely_satellite),]

# read in rps_data
# rps_blast_out <- read_tsv(file = "echisCarinatus/echisCarinatus_combined_library.rps.out",
                          col_names = c("qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                        "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle")) %>%
  separate(stitle, into = c("ref", "abbrev", "full"), sep = ", ") %>%
  mutate(ID = sub("#.*", "", tolower(qseqid))) %>%
  filter(startsWith(tolower(qseqid), "rnd") | startsWith(tolower(qseqid), "dr")) %>%
  filter(ID %in% tolower(sub("#.*", "", names(likely_not_satellite))))

# read in "suitable" domains
suitable_domains <- read_tsv("~/projectDrive_1/repbase_check/suitable_domains.tsv", col_names = c("sseqid", "ref", "full"))

# identify repeats containing bad and suitable domains
definately_bad_domains <- rps_blast_out %>%
  filter(startsWith(abbrev, "PBP") | # taste receptors
           startsWith(abbrev, "7tm") | abbrev %in% c("ANF_receptor", "NCD3G") | # vomeronasal genes
           startsWith(abbrev, "ZnMc") | # venom genes
           abbrev == "KRAB" | # KRAB zinc fingers
           startsWith(tolower(full), "cytochrome") |
           (abbrev == "SCAN" & !endsWith(qseqid, "#LTR/Gypsy")) |
           (startsWith(abbrev, "Ig") & !grepl("#LTR/", qseqid))
  )

definately_bad_domains <- rps_blast_out %>%
  filter(qseqid %in% definately_bad_domains$qseqid)

good_domains <- rps_blast_out %>%
  filter(ref %in% suitable_domains$ref) %>%
  filter(!qseqid %in% definately_bad_domains$qseqid)

good_domains <- rps_blast_out %>%
  filter(qseqid %in% good_domains$qseqid)

potentially_bad_domains <- rps_blast_out %>%
  filter(!ref %in% suitable_domains$ref) %>%
  filter(!qseqid %in% good_domains$qseqid,
         !qseqid %in% definately_bad_domains$qseqid)

potentially_bad_domains <- rps_blast_out %>%
  filter(qseqid %in% potentially_bad_domains$qseqid)

View(as_tibble(as.data.frame(table(potentially_bad_domains$full))) %>% arrange(-Freq))

likely_classified_correctly <- likely_not_satellite[!sub("#.*", "", names(likely_not_satellite)) %in% c(potentially_bad_domains$ID, definately_bad_domains$ID),]
likely_classified_incorrectly <- likely_not_satellite[sub("#.*", "", names(likely_not_satellite)) %in% c(potentially_bad_domains$ID, definately_bad_domains$ID),]

suspicious_protein_absent_LINEs <-
  likely_classified_correctly[grepl("#LINE", names(likely_classified_correctly)) &
                                !sub("#.*", "", names(likely_classified_correctly)) %in% sub("#.*", "", rps_blast_out$qseqid) &
                                (startsWith(names(likely_classified_correctly), "DR") | startsWith(names(likely_classified_correctly), "rnd"))]

suspicious_protein_absent_LTRs <-
  likely_classified_correctly[grepl("#LTR", names(likely_classified_correctly)) &
                                !sub("#.*", "", names(likely_classified_correctly)) %in% sub("#.*", "", rps_blast_out$qseqid) &
                                width(likely_classified_correctly) > 500 &
                                (startsWith(names(likely_classified_correctly), "DR") | startsWith(names(likely_classified_correctly), "rnd"))]

suspicious_protein_absent <- c(suspicious_protein_absent_LINEs, suspicious_protein_absent_LTRs)

# identify and rename suspicious LTRs and LINEs
writeXStringSet(x = suspicious_protein_absent, filepath = "suspicious_protein_absent_seq.fasta")

system("/home/james/bin/ncbi-blast/bin/blastn -task dc-megablast -query suspicious_protein_absent_seq.fasta -db ../../repbase_check/RepBase5May2021.fasta -outfmt \"6 qseqid qstart qend qlen pident length sseqid\" -num_threads 64 -out suspicious_protein_absent_check.out")

suspicious_blast <- read_tsv("suspicious_protein_absent_check.out", col_names = c("seqnames", "qstart", "qend", "qlen", "pident", "length", "sseqid")) %>%
  filter(length >= 0.5* qlen) %>%
  group_by(seqnames) %>%
  dplyr::slice(1) %>%
  ungroup()

not_suspicious_protein_absent <- suspicious_protein_absent[sub("#.*", "", names(suspicious_protein_absent)) %in% sub("#.*", "",suspicious_blast$seqnames)]

rather_suspicious_protein_absent <- suspicious_protein_absent[!sub("#.*", "", names(suspicious_protein_absent)) %in% sub("#.*", "",suspicious_blast$seqnames)]

likely_classified_correctly <- likely_classified_correctly[!names(likely_classified_correctly) %in% names(rather_suspicious_protein_absent)]

# compile manually curates and those cleaned for satellite and domains
ready_for_export <- c(cleaned_up, likely_classified_correctly)

# remove sequences which are over 50% N
too_much_N <- as_tibble(alphabetFrequency(ready_for_export)) %>%
  mutate(seqnames = names(ready_for_export), repeat_width = width(ready_for_export)) %>%
  filter(N >= 0.50*repeat_width)
ready_for_export <- ready_for_export[!names(ready_for_export) %in% too_much_N$seqnames]

# write new library to file
writeXStringSet(ready_for_export, "crotalusTigris_iter4_library_Mar9.fasta")
