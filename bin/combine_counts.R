#!/usr/bin/env Rscript

################################################################################
# combine processed counts by summing counts by sample name across lanes
# counts are stord as annotated expression set and MAGeCK input file
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
library(Biobase)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)

# command arguments
args         <- commandArgs(trailingOnly = TRUE)
library_file <- args[1]
barcode_file <- args[2]
count_files  <- args[3:length(args)]

# library_file <- "library.txt"
# barcode_file <- "barcodes.txt"
# count_files <- c("results/CAF21ANXX_6_20170127B_20170128/counts/CAF21ANXX_6_20170127B_20170128.txt",
#                  "results/CAF21ANXX_6_20170127B_20170130/counts/CAF21ANXX_6_20170127B_20170130.txt")

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------
read_featurecounts <- function(path) {
  read_tsv(path, comment = "#") %>%
    select(-Chr, -Start, -End, -Strand, -Length) %>%
    rename(id = Geneid)
}

# ------------------------------------------------------------------------------
# combine counts
# ------------------------------------------------------------------------------
names(count_files) <- str_replace(basename(count_files), ".txt", "")
pattern <- paste(c(paste0(names(count_files), "_"), ".bam"), collapse = "|")

counts <- lapply(count_files, read_featurecounts) %>%
  lapply(gather, sample_name, count, -id) %>%
  bind_rows %>%
  mutate(sample_name = str_replace_all(sample_name, pattern, "")) %>%
  group_by(id, sample_name) %>%
  summarize(count = sum(count)) %>%
  ungroup %>%
  spread(sample_name, count) %>%
  data.frame(row.names = 1) %>%
  as.matrix

# ------------------------------------------------------------------------------
# build expression set
# ------------------------------------------------------------------------------
features <- read_tsv(library_file) %>%
  arrange(id) %>%
  data.frame %>%
  AnnotatedDataFrame
rownames(features) <- features$id

pheno <- read_tsv(barcode_file) %>%
  select(-lane, -barcode) %>%
  arrange(sample_name) %>%
  distinct %>%
  data.frame %>%
  AnnotatedDataFrame
rownames(pheno) <- pheno$sample_name

eset <- ExpressionSet(assayData = counts,
                      featureData = features,
                      phenoData = pheno)

eset %>%
  write_rds("counts.rds")

# ------------------------------------------------------------------------------
# MAGeCK output
# ------------------------------------------------------------------------------
eset_to_mageck <- function(eset) {
  id_name <- colnames(fData(eset))[1]
  id2gene <- fData(eset)[, 1:2]
  counts <- exprs(eset) %>%
    data.frame %>%
    rownames_to_column(var = id_name)
  inner_join(id2gene, counts, by = id_name)
}

eset_to_mageck(eset) %>%
  write_tsv("counts_mageck.txt")