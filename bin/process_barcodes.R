#!/usr/bin/env Rscript

################################################################################
# process sample annotation to obtain input demultiplexing
# part of CRISPR / shRNA screen pre-processing pipeline
# 
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2017/09/20
################################################################################

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------
# library
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

args       <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

# input_file <- "barcodes.txt"

# ------------------------------------------------------------------------------
# import
# ------------------------------------------------------------------------------
# barcode file should be tab-separated text file with three columns:
# 1) sample_name
# 2) lane
# 3) barcode
# 4) ... any description of samples (e.g. genotypic information)

raw <- read_tsv(input_file)

# ------------------------------------------------------------------------------
# process
# ------------------------------------------------------------------------------
raw %>%
  select(lane, sample_name, barcode) %>%
  arrange(sample_name) %>%
  nest(-lane) %>%
  walk2(.x = .$data, .y = .$lane, .f = ~ write_tsv(.x, paste0(.y, ".txt"), col_names = FALSE))
