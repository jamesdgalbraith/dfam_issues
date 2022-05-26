library(tidyverse)

# script fishes out eukaryote proteins and removes hypothetical proteins

odb10v1_eukaryotes <- read_tsv("data/odb10v1_level2species.tab", col_names = c("top", "species", "X3", "X4")) %>%
  dplyr::select(top, species) %>%
  dplyr::filter(top == 2759)

odb10v1_genes_eukaryotes <- read_tsv("data/odb10v1_genes.tab", col_names = F) %>%
  dplyr::filter(X2 %in% odb10v1_eukaryotes$species)

odb10v1_genes_eukaryotes_output <- odb10v1_genes_eukaryotes %>%
  dplyr::select(X1, X2, X4, X8)

odb10v1_genes_eukaryotes_output <- odb10v1_genes_eukaryotes_output %>%
  filter(!startsWith(x = X8, prefix = "hypothetical protein")) %>%
  filter(!startsWith(x = X8, prefix = "Hypothetical protein"))

write_tsv(odb10v1_genes_eukaryotes_output, "data/eukaryotes_all_data.tsv", col_names = F)