library(tidyverse)
library(plyranges)
library(BSgenome)
library(svglite)
library(wordcloud)

# set analysis set
to_analyse <- "seq/dfam_Hemiptera.fasta"

if (grepl("seq/", to_analyse)){
  to_analyse <- sub("seq/", "", to_analyse)
}

source_dfam <- F

# most common repeats
# read in libraries
combined_library <- readDNAStringSet(paste0("seq/", to_analyse))
repbase_library <- readDNAStringSet("/media/projectDrive_1/databases/repbase/RepBase5May2021.fasta")
repbase_classes <- read_tsv("data/repbase_classes.tsv")
dfam_list <- read_tsv("data/Dfam.3.2.seq_list.txt", col_names = "seqnames")

combined_library[!startsWith(names(combined_library), "DR")]
combined_library[startsWith(names(combined_library), "DR")]

# tablise repbase
repbase_library_tbl <- as_tibble(as.data.frame(names(repbase_library))) %>%
  dplyr::rename(seqnames = `names(repbase_library)`) %>%
  tidyr::separate(seqnames, into = c("rb_ref", "rb_subclass", "rb_species"), sep = "\t") %>%
  mutate(te_width = width(repbase_library)) %>%
  inner_join(repbase_classes) %>%
  mutate(rb_subclass = gsub(" ", "_", rb_subclass))

if(source_dfam == T){

  combined_library <- combined_library[sub("#.*", "", names(combined_library)) %in% dfam_list$seqnames]

}

# extract curated repeats
combined_library_tbl <- tibble(seqnames = sub(" .*", "", names(combined_library)),
                               te_width = width(combined_library)) %>%
  mutate(class = sub(" .*", "", sub(".*#", "", seqnames))) %>%
  tidyr::separate(class, into = c("class", "subclass"), sep = "/")

# 
if(source_dfam == T){
  
  # add species
  combined_library_tbl$species = sub(".*@", "", names(combined_library))
  
  # remove curated
  combined_library_tbl <- combined_library_tbl[startsWith(combined_library_tbl$seqnames, "DR"),]
  
}

# read in trf
trf_out_ranges <- read_gff3(paste0("data/", to_analyse, ".trf.gff"))

if(source_dfam == T){
  
  # remove curated
  trf_out_ranges <- trf_out_ranges[startsWith(x = as.character(seqnames(trf_out_ranges)), "DR")]
  
}

# convert data types, calculate length of satellites, ensure over 2.5 copies
trf_out_tbl <- as_tibble(trf_out_ranges) %>%
  mutate(seqnames = as.character(seqnames),
         copies = as.double(copies),
         period = as.double(period),
         consensus_size = as.double(consensus_size),
         length = ceiling(period*copies)) %>%
  dplyr::filter(copies >=2.1) %>%
  dplyr::select(seqnames, start, end, period, copies, length)

# calculate proportion which is tandem repeat
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

# class_to_plot <- "LINE"
# 
# ggplot() +
#   geom_histogram(mapping = aes(combined_library_tbl_tr[combined_library_tbl_tr$class==class_to_plot,]$tr_prop), binwidth = 0.025) + theme_bw() +
#   scale_x_continuous(name = paste0("Proportion of TE width identified as simple repeat (n = ",
#                                    nrow(combined_library_tbl_tr[combined_library_tbl_tr$class==class_to_plot,]), ")")) +
#   ggtitle(label = class_to_plot)
# 
# ggplot() +
#   geom_histogram(mapping = aes(combined_library_tbl_tr[combined_library_tbl_tr$class==class_to_plot & combined_library_tbl_tr$subclass != "Penelope",]$te_width), binwidth = 500) + theme_bw() +
#   scale_x_continuous(name = "Length (bp)") + 
#   ggtitle(label = paste0("Dfam ", class_to_plot, "s"))

# read in rps_data
rps_blast_out <- read_tsv(file = paste0("data/", to_analyse, ".rps.out"),
                          col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                        "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle")) %>%
  separate(stitle, into = c("ref", "abbrev", "full"), sep = ", ")

if (source_dfam == T) {
  rps_blast_out <- rps_blast_out[startsWith(rps_blast_out$seqnames, "DR"),]
}

# read in repbase blast data
repbase_blast_out <- read_tsv(file = paste0("data/", to_analyse, ".repbase.out"),
                          col_names = c("seqnames", "qstart", "qend", "qlen", "rb_ref", "sstart", "send", "slen",
                                        "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs")) %>%
  dplyr::filter(evalue <= 1e-50)

if (source_dfam == T) {
  repbase_blast_out <- repbase_blast_out[startsWith(repbase_blast_out$seqnames, "DR"),]
  }

# joining repbase data with other data
best_repbase_blast_out <- repbase_blast_out %>%
  dplyr::group_by(seqnames) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  inner_join(repbase_library_tbl) %>%
  dplyr::select(-te_width, -evalue) %>%
  inner_join(combined_library_tbl_tr) %>%
  dplyr::filter(qcovs > 50, tr_prop < 0.1) %>%
  dplyr::filter(!grepl("sat", class, ignore.case = T))

combined_library_tbl_tr
correctos <- best_repbase_blast_out[best_repbase_blast_out$class == best_repbase_blast_out$rb_class,]
best_repbase_blast_out[best_repbase_blast_out$class != best_repbase_blast_out$rb_class,] %>%
  dplyr::select(seqnames, class, subclass, qlen, rb_ref, rb_class, rb_subclass, slen, qcovs, tr_prop)
# 
# as_tibble(as.data.frame(table(rps_blast_out_only_suitable[grepl("#LINE", rps_blast_out_only_suitable$seqnames)
#                             & !rps_blast_out_only_suitable$seqnames %in% correctos$seqnames,]$abbrev)))
# 
# repbase_blast_out %>%
#   dplyr::group_by(seqnames) %>%
#   dplyr::slice(1) %>%
#   dplyr::ungroup() %>%
#   inner_join(repbase_library_tbl) %>%
#   dplyr::select(-te_width, -evalue) %>%
#   inner_join(combined_library_tbl_tr) %>%
#   dplyr::filter(qcovs <= 50) %>%
#   dplyr::filter(!grepl("sat", class, ignore.case = T)) %>%
#   dplyr::select(seqnames, class, subclass, qlen, rb_ref, rb_class, rb_subclass, slen, qcovs, tr_prop) %>%
#   dplyr::arrange(-tr_prop)
# 
# best_repbase_blast_out %>%
#   filter(rb_class == class)
# 
# 
# best_repbase_blast_out %>%
#   dplyr::select(seqnames, class, subclass, qlen, rb_ref, rb_class, rb_subclass, slen, qcovs, tr_prop) %>%
#   dplyr::filter(class != rb_class, class != "Unknown")

# identify sequences with no homology to RepBase or CDD
oddly_unknown <- combined_library_tbl_tr %>% 
  filter(!seqnames %in% repbase_blast_out$seqnames, !seqnames %in% rps_blast_out$seqnames)

no_domains <- combined_library_tbl_tr[!combined_library_tbl_tr$seqnames %in% rps_blast_out$seqnames,]

# read in "suitable" domains
suitable_domains <- read_tsv("data/suitable_domains.tsv", col_names = c("sseqid", "ref", "abbrev", "full"))

# Filter repeats containing suitable domains
rps_blast_out_suitable <- rps_blast_out %>%
  filter(sseqid %in% suitable_domains$sseqid)

rps_blast_out_suspicious <- rps_blast_out %>%
  filter(!sseqid %in% suitable_domains$sseqid)

rps_blast_out_only_suitable <- rps_blast_out %>%
  filter(seqnames %in% rps_blast_out_suitable$seqnames, !seqnames %in%  rps_blast_out_suspicious$seqnames)

rps_blast_out_only_suspicious <- rps_blast_out %>%
  filter(!seqnames %in% rps_blast_out_suitable$seqnames)

writeXStringSet(filepath = paste0("out/only_suspicious_", to_analyse),
                x = combined_library[sub(" .*", "", names(combined_library)) %in% rps_blast_out_only_suspicious$seqnames])

rps_blast_out_chimeric <- rps_blast_out %>%
  filter(seqnames %in% rps_blast_out_suitable$seqnames, seqnames %in% rps_blast_out_suspicious$seqnames)

definately_bad_domains <- rps_blast_out %>%
  filter(startsWith(abbrev, "PBP") | # taste receptors
           startsWith(abbrev, "7tm") | abbrev %in% c("ANF_receptor", "NCD3G") | # vomeronasal genes
           startsWith(abbrev, "ZnMc") | # venom genes
           abbrev == "KRAB" | # KRAB zinc fingers
           startsWith(tolower(full), "cytochrome") | # cytochrome p450
           (startsWith(abbrev, "Ig") & !grepl("#LTR/", seqnames)) # immunoglobin
  )

only_bad_domains <- rps_blast_out %>%
  filter(seqnames %in% definately_bad_domains$seqnames, !seqnames %in% rps_blast_out_suitable$seqnames)

# View(rps_blast_out_only_suspicious)
# length(combined_library_tbl_tr$seqnames)
# length(no_domains$seqnames)
# length(base::unique(rps_blast_out_only_suitable$seqnames))
# length(base::unique(rps_blast_out_only_suspicious$seqnames))
# length(base::unique(rps_blast_out_only_suspicious[rps_blast_out_only_suspicious$seqnames %in% definately_bad_domains$seqnames,]$seqnames))
# length(base::unique(rps_blast_out_chimeric$seqnames))
# length(base::unique(rps_blast_out_chimeric[rps_blast_out_chimeric$seqnames %in% definately_bad_domains$seqnames,]$seqnames))
# as_tibble(as.data.frame(table(rps_blast_out_only_suspicious$abbrev))) %>% arrange(-Freq)
# 
# length(oddly_unknown$seqnames)
# as.data.frame(table(oddly_unknown$class))
# 
# trf_out_ranges[seqnames(trf_out_ranges) %in% combined_library_tbl_tr[combined_library_tbl_tr$tr_prop > 0.5,]$seqnames]

# identify families over 20% "N"
too_much_nothingness <- as_tibble(alphabetFrequency(combined_library)) %>%
  mutate(total = rowSums(alphabetFrequency(combined_library)),
         seqnames = names(combined_library)) %>%
  filter(N >= 0.5*total)

rps_blast_out_chimeric %>% filter(seqnames %in% too_much_nothingness$seqnames)

# good seqs a) do not contain only suspicious domains b) do not contain definately bad domains c) are less than 20% "N" and are less than 50% tandem repeat
good_seqs <- combined_library[!names(combined_library) %in% c(definately_bad_domains$seqnames,
                                                 too_much_nothingness$seqnames,
                                                 combined_library_tbl_tr[combined_library_tbl_tr$tr_prop > 0.5 &&
                                                                           !grepl("Satellite", names(combined_library_tbl_tr)),]$seqnames)]

writeXStringSet(good_seqs, paste0("out/cleaned_", to_analyse))


# self_blast_out <- read_tsv(paste0("data/", to_analyse, ".self_blast.out"),
#          col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen",
#                        "start", "end", "sstart", "send", "evalue", "bitscore", "slen")) %>%
#   dplyr::select(-sseqid)
# 
# for_seg_plotting <- as_tibble(reduce_ranges(as_granges(self_blast_out %>% dplyr::filter(start != sstart, end != send, length < 0.9*slen)), slen = mean(slen))) %>%
#   group_by(seqnames) %>%
#   mutate(sum_width = sum(width), seqnames = as.character(seqnames)) %>%
#   dplyr::slice(1) %>%
#   ungroup() %>%
#   dplyr::filter(sum_width >= 0.5*slen) %>%
#   arrange(-slen) %>%
#   inner_join(oddly_unknown)
# 
# table(self_blast_out[self_blast_out$start != self_blast_out$sstart &
#                self_blast_out$end != self_blast_out$send &
#                self_blast_out$length < 0.9*self_blast_out$slen,] $seqnames) %>%
#   base::as.data.frame() %>% 
#   dplyr::as_tibble() %>% 
#   dplyr::mutate(seqnames = as.character(Var1)) %>%
#   dplyr::filter(!grepl("Satellite", seqnames))
# 
# ggplot(self_blast_out[self_blast_out$seqnames == for_seg_plotting$seqnames[i],]) +
#   geom_segment(aes(x = start, xend = end, y = sstart, yend = send, colour = send-sstart > 0), show.legend = F) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   ggtitle(for_seg_plotting$seqnames[i])
# 
# 
# ggplot(oddly_unknown[!oddly_unknown$class %in% c("Unknown", "Simple_repeat", "Satellite") & 
#                        !grepl("RNA", oddly_unknown$class),]) + geom_point(aes(x = tr_width, y = tr_prop, colour = class))
# 
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$seqnames %in% rps_blast_out_only_suitable$seqnames,]) + geom_point(aes(x = tr_width, y = tr_prop, colour = class))
# 
# very_odd <- combined_library_tbl_tr[combined_library_tbl_tr$seqnames %in% rps_blast_out_only_suitable$seqnames,] %>%
#   filter(tr_prop > 0.2) %>%
#   arrange(-tr_prop)
# 
# best_repbase_blast_out %>%
#   dplyr::select(seqnames, length, qlen, pident, rb_ref, rb_class, class, tr_prop, qcovs, species) %>%
#   filter(class != rb_class, class != "Unknown")
#   
# best_repbase_blast_out %>%
#   dplyr::filter(rb_class == "Satellite", class != "Satellite") %>%
#   dplyr::select(seqnames, length, qlen, pident, rb_ref, rb_class, class, tr_prop, qcovs, species)
# 
# 
# 
# plottable <- very_odd$seqnames
# 
# i=2
# 
# trf_out_tbl[trf_out_tbl$seqnames == plottable[i],] %>%
#   as_granges()
# 
# for_plotting <- rps_blast_out %>%
#   dplyr::filter(seqnames == plottable[i],
#                 send - sstart + 1 >= 0.2*slen)
# 
# for_plotting <- for_plotting %>%
#   dplyr::mutate(n = (nrow(for_plotting):1)/nrow(for_plotting),
#                 start = ifelse(qstart < qend, qstart, qend),
#                 end = ifelse(qstart > qend, qstart, qend)) %>%
#   dplyr::select(seqnames, start, end, qlen, abbrev, n, bitscore) %>%
#   arrange(start, -bitscore) %>%
#   mutate(overlap = lag(end) > start,
#          overlap = ifelse(is.na(overlap), FALSE, overlap)) %>%
#   filter(overlap == FALSE)
# 
# for_plotting <- for_plotting %>%
#   mutate(n = (nrow(for_plotting):1)/nrow(for_plotting),
#          y_coord = (2 * n - (1/nrow(for_plotting)))/2)
# 
# ggplot(data = for_plotting) +
#   geom_rect(aes(xmin = start, xmax = end, ymin = n-(1/nrow(for_plotting)), ymax = n, fill = abbrev), show.legend = F, colour = "white", size = 0.5) +
#   geom_label(aes(label=str_wrap(abbrev,12), x=(start + end)/2, y = y_coord), size=3.5) +
#   theme_bw() + scale_y_continuous(labels = NULL, breaks = NULL, name = NULL) + ggtitle(for_plotting$seqnames[1]) +
#   scale_x_continuous(limits=c(0,max(c(for_plotting$qlen))), name = "Position (bp)")
# 
# 
# for_cloud <- as_tibble(as.data.frame(table(sub("_.*", "", rps_blast_out_suspicious$abbrev)))) %>%
#   mutate(Var1 = as.character(Var1)) %>%
#   arrange(-Freq)
# 
# wordcloud(words = substring(for_cloud$Var1, 1, 25), freq = for_cloud$Freq, random.order = F, colors = brewer.pal(8, "Dark2"),
#           rot.per = 0.35)
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$class == "Unknown",]) + geom_histogram(mapping = aes(tr_prop), binwidth = 0.025) +
#   scale_x_continuous(limits = c(-0.025, 1.025))
# 
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$class == "LINE" & combined_library_tbl_tr$subclass != "Penelope",]) +
#   geom_histogram(aes(te_width), binwidth = 250) +
#   scale_y_continuous(trans = "log10")
# 
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$class == "LINE" & combined_library_tbl_tr$subclass == "Penelope",]) +
#   geom_histogram(aes(te_width), binwidth = 250)+
#   scale_y_continuous(trans = "log10")
# 
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$class == "SINE",]) +
#   geom_histogram(aes(te_width), binwidth = 250)+
#   scale_y_continuous(trans = "log10")
# 
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$class == "LTR",]) +
#   geom_histogram(aes(te_width), binwidth = 250)
# 
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$class == "DNA",]) +
#   geom_histogram(aes(te_width, fill = sub("-.*", "" ,subclass)), binwidth = 250)
# 
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$class == "RC",]) +
#   geom_histogram(aes(te_width), binwidth = 250)
# 
# ggplot(combined_library_tbl_tr[combined_library_tbl_tr$class == "LINE" & combined_library_tbl_tr$subclass != "Penelope",]) +
#   geom_point(aes(x = te_width, y = tr_prop))
# 
# combined_library_tbl_tr[combined_library_tbl_tr$tr_prop > 0.1 &
#                           combined_library_tbl_tr$class == "LINE" &
#                           combined_library_tbl_tr$subclass != "Penelope",]
# 
# rps_blast_out[startsWith(rps_blast_out$abbrev, "PBP1_CaSR"),]$full
#   }
# 
# combined_library_tbl_tr[combined_library_tbl_tr$te_width > 6000 &
#                           combined_library_tbl_tr$class == "LINE" &
#                           combined_library_tbl_tr$subclass != "Penelope",] %>%
#   arrange(-tr_prop)
