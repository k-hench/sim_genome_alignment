#!/usr/bin/env Rscript
library(valr)
library(tidyverse)
library(prismatic)
library(Biostrings)
library(glue)
set.seed(42)

dir.create("data/genomes/", showWarnings = FALSE)
dir.create("img/", showWarnings = FALSE)
dir.create("results/bed", showWarnings = FALSE, recursive = TRUE)

n_inv <- 2
l_inv <- 3e3

base_dna <- readDNAStringSet("data/pre/final_hap_1_f.fa")
base_size <- length(base_dna$hap_1)
base_break <- rbinom(1, size = base_size, prob = .66)

cut_coordinates <- \(seed = 42, buff_size = 150){
  set.seed(seed)
  tibble(chrom = c("chr1", "chr2"),
         start = c(0 +  rbinom(1, size = buff_size, prob = .66),
                   base_break +  rbinom(1, size = buff_size, prob = .66)),
         end = c(base_break - rbinom(1, size = buff_size, prob = .66),
                 base_size - rbinom(1, size = buff_size, prob = .66))) |>
    mutate(size = end - start)
}

export_bed <- \(seed = 42, ...){
  cut_coordinates(seed = seed, ...) |>
    create_invs(n = n_inv, seed = seed, length = l_inv)|>
    select(chrom, start, end) |>
    mutate(chrom = glue("g{seed}_{chrom}")) |>
    arrange(chrom, start) |>
    write_tsv(glue("results/bed/inversions_{seed}.bed"))
}

create_invs <- \(genome_in, 
                 seed = 42,
                 n = n_inv, 
                 length = l_inv){
  set.seed(seed)
  bed_random(genome_in,
             n = n,
             length = length,
             seed = seed) |> 
    left_join(genome_in |>
                dplyr::select(chrom, ostart = start)) |> 
    mutate(gstart = ostart + start,
           gend = ostart + end)
}

# cut_coordinates(2, 250) |>
#   create_invs(n = 5)

p0 <- ggplot() +
  geom_vline(xintercept = c(0,base_break, base_size),
             lty = 3, linewidth =.4) +
  geom_linerange(data =  tibble(chrom = c("chr1", "chr2"),
                                start = c(0, base_break),
                                end = c(base_break, base_size)),
                 aes(y = 0, xmin = start, xmax = end),
                 color = "gray90",
                 linewidth = 3) +
  map(1:3, \(idx){geom_linerange(data = cut_coordinates(idx, 250),
                                 aes(y = idx, xmin = start, xmax = end,
                                     color = glue("g{idx}")),
                                 linewidth = 3)}) +
  map(1:3, \(idx){geom_linerange(data = cut_coordinates(idx, 250) |>
                                   create_invs(n = n_inv, seed = idx, length = l_inv),
                                 aes(y = idx, xmin = gstart, xmax = gend),
                                 color = rgb(0,0,0,.3),
                                 linewidth = 4.5)}) +
  scale_color_manual("genome",
                     values = c("gray30",
                                "#f46d43",
                                "#66c2a5") |>
                       clr_desaturate(.1) |>
                       set_names(nm = c("g1", "g2", "g3"))) +
  theme_minimal()

ggsave(plot = p0, 
       "img/base_genomes.svg",
       width = 7, height = 2)

modify_genome <- \(idx, opt_args = list()){
  set.seed(idx)
  template_dna <- readDNAStringSet(glue("data/pre/final_hap_{idx}_f.fa"))
  proto_genome <- cut_coordinates(idx, buff_size = 250)
  inversions <- proto_genome |>
    create_invs(n = n_inv, seed = idx)
  
  genome_out <- proto_genome |> 
    mutate(seq = map2(start, end, \(s,e){template_dna[[1]][s:e]}))
  
  for(k in nrow(inversions)){
    idx_tmp <- which(genome_out$chrom == inversions$chrom[[k]])
    seq_tmp <- (genome_out[idx_tmp,]$seq)[[1]]
    genome_out[idx_tmp,]$seq[[1]] <- c(seq_tmp[1:(inversions$start[[k]]- 1)],
                                       reverseComplement(seq_tmp[inversions$start[[k]]:(inversions$end[[k]] - 1)]),
                                  seq_tmp[inversions$end[[k]]:genome_out$size[[idx_tmp]]])
    rm(idx_tmp, seq_tmp)
  }
  genome_out
}

export_genome <- \(idx){
  writeXStringSet(DNAStringSet(modify_genome(idx = idx)$seq |> 
                                 set_names(nm = c("chr1", "chr2"))),
                  filepath = glue("data/genomes/genome_{idx}.fa"))
}

1:3 |> walk(export_genome)
1:3 |> walk(export_bed)
