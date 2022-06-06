# Import librarires
library(optparse, quietly = T)

# Get optparse
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-d", "--dfam"), type="logical", default=FALSE,
              help="Set to true if data analysed contains vertebrate sequences from DFAM", metavar="character")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# # for testing
# opt <- list(file = "seq/cleaned_Echis_RepeatModeler.fasta")

if(is.null(opt$file)){
  stop("Input file must be provided.")
}

library(tidyverse, quietly = T, warn.conflicts = F, verbose = F)
library(plyranges, quietly = T, warn.conflicts = F, verbose = F)
library(BSgenome, quietly = T, warn.conflicts = F, verbose = F)


if (grepl("seq/", opt$file)){
  opt$file <- sub("seq/", "", opt$file)
}

# Read in and preprocess all data
combined_library <- suppressMessages(readDNAStringSet(paste0("seq/", opt$file)))
repbase_library <- suppressMessages(readDNAStringSet("/media/projectDrive_1/databases/repbase/RepBase5May2021.fasta"))
repbase_classes <- read_tsv("data/repbase_classes.tsv", show_col_types = F)
trf_out_ranges <- read_gff3(paste0("data/", opt$file, ".trf.gff"))
rps_blast_out <- read_tsv(file = paste0("data/", opt$file, ".rps.out"),
                          col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                        "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle"),
                          show_col_types = F) %>%
  separate(stitle, into = c("ref", "abbrev", "full"), sep = ", ")
repbase_blast_out <- read_tsv(file = paste0("data/", opt$file, ".repbase.out"),
                              col_names = c("seqnames", "qstart", "qend", "qlen", "rb_ref", "sstart", "send", "slen",
                                            "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs"),
                              show_col_types = F) %>%
  dplyr::filter(evalue <= 1e-50)
suitable_domains <- read_tsv("data/suitable_domains.tsv", col_names = c("sseqid", "ref", "abbrev", "full"),
                             show_col_types = F)

# Tabulate repbase
repbase_library_tbl <- as_tibble(as.data.frame(names(repbase_library))) %>%
  dplyr::rename(seqnames = `names(repbase_library)`) %>%
  tidyr::separate(seqnames, into = c("rb_ref", "rb_subclass", "rb_species"), sep = "\t") %>%
  mutate(te_width = width(repbase_library)) %>%
  inner_join(repbase_classes) %>%
  mutate(rb_subclass = gsub(" ", "_", rb_subclass))


# Tabulate input data
combined_library_tbl <- tibble(seqnames = sub(" .*", "", names(combined_library)),
                               te_width = width(combined_library)) %>%
  mutate(class = sub(" .*", "", sub(".*#", "", seqnames))) %>%
  tidyr::separate(class, into = c("class", "subclass"), sep = "/")

# DFAM specific filter
if(opt$dfam == T){
  
  dfam_list <- read_tsv("data/Dfam.3.2.seq_list.txt", col_names = "seqnames", show_col_types = F)
  combined_library <- combined_library[sub("#.*", "", names(combined_library)) %in% dfam_list$seqnames]
  
  # determine species
  combined_library_tbl$species = sub(".*@", "", names(combined_library))
  
  # remove curated
  combined_library_tbl <- combined_library_tbl[startsWith(combined_library_tbl$seqnames, "DR"),]
  trf_out_ranges <- trf_out_ranges[startsWith(x = as.character(seqnames(trf_out_ranges)), "DR")]
  rps_blast_out <- rps_blast_out[startsWith(rps_blast_out$seqnames, "DR"),]
  repbase_blast_out <- repbase_blast_out[startsWith(repbase_blast_out$seqnames, "DR"),]
  
}

# Tabulate TRF data, convert data types, calculate length of satellites, ensure over 2.1 copies
trf_out_tbl <- as_tibble(trf_out_ranges) %>%
  mutate(seqnames = as.character(seqnames),
         copies = as.double(copies),
         period = as.double(period),
         consensus_size = as.double(consensus_size),
         length = ceiling(period*copies)) %>%
  dplyr::filter(copies >=2.1) %>%
  dplyr::select(seqnames, start, end, period, copies, length)

# Calculate proportion of TE which is tandem repeat
combined_library_tbl_tr <- plyranges::as_granges(trf_out_tbl) %>%
  IRanges::reduce() %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(seqnames = as.character(seqnames)) %>%
  dplyr::select(seqnames, width) %>%
  dplyr::group_by(seqnames) %>%
  dplyr::mutate(tr_width = base::sum(width)) %>%
  dplyr::ungroup() %>%
  dplyr::select(seqnames, tr_width) %>%
  base::unique() %>%
  full_join(combined_library_tbl) %>%
  dplyr::mutate(tr_width = ifelse(is.na(tr_width), 0, tr_width),
                tr_prop = tr_width/te_width)

# Determine best repbase hit
best_repbase_blast_out <- repbase_blast_out %>%
  dplyr::group_by(seqnames) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  inner_join(repbase_library_tbl) %>%
  dplyr::select(-te_width, -evalue) %>%
  inner_join(combined_library_tbl_tr) %>%
  dplyr::filter(qcovs > 50, tr_prop < 0.1) %>%
  dplyr::filter(!grepl("sat", class, ignore.case = T)) # removes satellites

# Identify repeats with no repbase or CDD hits
oddly_unknown <- combined_library_tbl_tr %>% 
  filter(!seqnames %in% repbase_blast_out$seqnames, !seqnames %in% rps_blast_out$seqnames)

# Identify repeats with no CDD hits
no_domains <- combined_library_tbl_tr[!combined_library_tbl_tr$seqnames %in% rps_blast_out$seqnames,]

# Filter repeats containing suitable domains
rps_blast_out_suitable <- rps_blast_out %>%
  filter(sseqid %in% suitable_domains$sseqid)

# Filter repeats containing suspicious domains
rps_blast_out_suspicious <- rps_blast_out %>%
  filter(!sseqid %in% suitable_domains$sseqid)

# Filter repeats containing only suitable domains
rps_blast_out_only_suitable <- rps_blast_out %>%
  filter(seqnames %in% rps_blast_out_suitable$seqnames, !seqnames %in%  rps_blast_out_suspicious$seqnames)

# Filter repeats containing only suspicious domains
rps_blast_out_only_suspicious <- rps_blast_out %>%
  filter(!seqnames %in% rps_blast_out_suitable$seqnames)

# Filter repeats with both suitable and suspicious domains
rps_blast_out_chimeric <- rps_blast_out %>%
  filter(seqnames %in% rps_blast_out_suitable$seqnames, seqnames %in% rps_blast_out_suspicious$seqnames)

# Filter repeats with both suitable and suspicious domains
definately_bad_domains <- rps_blast_out %>%
  filter(startsWith(abbrev, "PBP") | # taste receptors
           startsWith(abbrev, "7tm") | abbrev %in% c("ANF_receptor", "NCD3G") | # vomeronasal genes
           startsWith(abbrev, "ZnMc") | # venom genes
           abbrev == "KRAB" | # KRAB zinc fingers
           startsWith(tolower(full), "cytochrome") | # cytochrome p450
           (startsWith(abbrev, "Ig") & !grepl("#LTR/", seqnames)) # immunoglobin
  )

# Determine families only contain un suitable domains
only_bad_domains <- rps_blast_out %>%
  filter(seqnames %in% definately_bad_domains$seqnames, !seqnames %in% rps_blast_out_suitable$seqnames)

# identify families over 50% "N"
too_much_nothingness <- as_tibble(alphabetFrequency(combined_library)) %>%
  mutate(total = rowSums(alphabetFrequency(combined_library)),
         seqnames = names(combined_library)) %>%
  filter(N >= 0.5*total)

likely_trash <- base::unique(c(definately_bad_domains$seqnames,
                               rps_blast_out_only_suspicious$seqnames, 
                               sub(" .*", "", too_much_nothingness$seqnames),
                               combined_library_tbl_tr[combined_library_tbl_tr$tr_prop > 0.5 &
                                                         !grepl(pattern = "Satellite", x = combined_library_tbl_tr$seqnames),]$seqnames))

# good seqs a) do not contain only suspicious domains b) do not contain definately bad domains c) are less than 20% "N" and are less than 50% tandem repeat
good_seqs <- combined_library[!sub(" .*", "", names(combined_library)) %in% likely_trash]
trash_seq <- combined_library[sub(" .*", "", names(combined_library)) %in% c(definately_bad_domains$seqnames,
                                                                             rps_blast_out_only_suspicious$seqnames)]

# write to file
writeXStringSet(good_seqs, paste0("out/cleaned_", opt$file))
writeXStringSet(trash_seq, paste0("out/trash_", opt$file))
