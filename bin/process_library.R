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
# 3) base (in forward orientation) to fill guide/shrna to maximal length (typically 20 or 22)
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
# 2) group   : group targeted by sgRNA / shRNA (e.g gene id or domain)
# 3) sequence: sequence of sgRNA / shRNA as it appears in the sequencing reads


raw <- read_tsv(input_file)

# check for id and sequence duplication
stopifnot(!any(duplicated(raw$id)))
stopifnot(!any(duplicated(raw$sequence)))

# ------------------------------------------------------------------------------
# process
# ------------------------------------------------------------------------------
# generate FASTA file for bowtie2 index
strand_side <- ifelse(strandedness == "forward", "left", "right")
padding     <- ifelse(strandedness == "forward", chartr("ATGC", "TACG", padding_base), padding_base)

seq_length <- max(nchar(raw$sequence))

fasta <- raw$sequence %>%
  toupper %>%
  str_pad(pad = padding, width = seq_length, side = strand_side) %>%
  set_names(raw$id) %>%
  DNAStringSet 

# generate SAF annotation file for featureCount (subread package)
saf <- tibble(GeneID = raw$id, 
              Chr    = raw$id, 
              Start  = 1, 
              End    = seq_length, 
              Strand = "*")

# ------------------------------------------------------------------------------
# export
# ------------------------------------------------------------------------------
fasta %>%
  writeXStringSet("library.fasta", format = "fasta")

saf %>%
  write_tsv("library.saf")
