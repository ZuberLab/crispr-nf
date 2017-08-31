#!/usr/bin/env Rscript

################################################################################
# combine processed counts by summing counts by sample name across lanes
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
args <- commandArgs(trailingOnly = TRUE)

# ------------------------------------------------------------------------------
# import & process
# ------------------------------------------------------------------------------
files <- args
names(files) <- str_replace(basename(files), ".txt", "")
counts <- lapply(files, read_tsv) %>%
  lapply(gather, sample_name, count, -id, -gene) %>%
  bind_rows(.id = "lane") %>%
  group_by(id, gene, sample_name) %>%
  summarize(count = sum(count)) %>%
  ungroup %>%
  spread(sample_name, count)

# ------------------------------------------------------------------------------
# export
# ------------------------------------------------------------------------------
counts %>%
  write_tsv("counts.txt")