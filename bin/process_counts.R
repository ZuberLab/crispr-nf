#!/usr/bin/env Rscript

################################################################################
# process featureCounts output to be compatible with MAGeCK
# part of CRISPR / shRNA screen pre-processing pipeline
# 
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2017/08/31
################################################################################

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------
# library
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# command arguments
# 1) first argument should be flag whether to combine counts by name
# 2) second argument should be path to original library text file

args         <- commandArgs(trailingOnly = TRUE)
count_file   <- args[1]
library_file <- args[2]

# ------------------------------------------------------------------------------
# import
# ------------------------------------------------------------------------------
feature_counts <- read_tsv(count_file, comment = "#")
library        <- read_tsv(library_file)

# ------------------------------------------------------------------------------
# process
# ------------------------------------------------------------------------------
lane <- str_replace(basename(count_file), "_fc.txt", "")

id2gene <- library %>%
  select(id, gene)

counts <- feature_counts %>%
  select(-Chr, -Start, -End, -Strand, -Length) %>%
  rename(id = Geneid) %>%
  gather(sample_name, count, -id) %>%
  mutate(sample_name = str_replace_all(sample_name, paste0(lane, "_|.bam"), "")) %>%
  spread(sample_name, count) %>%
  inner_join(id2gene, by = "id") %>%
  select(id, gene, everything())

# ------------------------------------------------------------------------------
# export
# ------------------------------------------------------------------------------
counts %>%
  write_tsv(paste0(lane, ".txt"))