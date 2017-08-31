#!/usr/bin/env Rscript

################################################################################
# process library text file
# part of CRISPR / shRNA screen pre-processing pipeline
# 
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2017/08/30
################################################################################

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------
# library
library(Biostrings)
library(readr)
library(stringr)
library(purrr)
library(magrittr)
library(tibble)

# command arguments
# 1) first argument should be path to library text file
# 2) second argument should be 
args       <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
padding    <- args[2]

# ------------------------------------------------------------------------------
# import
# ------------------------------------------------------------------------------
# library file should be tab-separated text file with three columns:
# 1) id      : unique id of sgRNA / shRNA
# 2) sequence: sequence of sgRNA / shRNA as it appears in the sequencing reads
# 3) gene    : gene targeted by sgRNA / shRNA (grouping variable)

raw <- read_tsv(input_file)

# check for id and sequence duplication
stopifnot(!any(duplicated(raw$id)))
stopifnot(!any(duplicated(raw$sequence)))

# ------------------------------------------------------------------------------
# process
# ------------------------------------------------------------------------------
# generate FASTA file for bowtie2 index
# TODO: this has to be corrected for strandedness
fasta <- raw$sequence %>%
  toupper %>%
  str_pad(pad = toupper(padding), width = 20, side = "right") %>%
  set_names(raw$id) %>%
  DNAStringSet 

# generate SAF annotation file for featureCount (subread package)
saf <- tibble(GeneID = raw$id, 
              Chr    = raw$id, 
              Start  = 1, 
              End    = 20, 
              Strand = "*")

# ------------------------------------------------------------------------------
# export
# ------------------------------------------------------------------------------
fasta %>%
  writeXStringSet("library.fasta", format = "fasta")

saf %>%
  write_tsv("library.saf")