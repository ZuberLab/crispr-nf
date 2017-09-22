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
# 1) path to library text file
# 2) strandedness of library (forward, reverse)
# 3) base (in forward orientation) to make guides length 20
args         <- commandArgs(trailingOnly = TRUE)
input_file   <- args[1]
strandedness <- args[2]
padding_base <- toupper(args[3])

stopifnot(strandedness %in% c("forward", "reverse"))
stopifnot(padding_base %in% c("A", "T", "C", "G"))

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
strand_side <- ifelse(strandedness == "forward", "right", "left")
padding     <- ifelse(strandedness == "forward", padding_base, chartr("ATGC", "TACG", padding_base))

fasta <- raw$sequence %>%
  toupper %>%
  str_pad(pad = padding, width = 20, side = strand_side) %>%
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
