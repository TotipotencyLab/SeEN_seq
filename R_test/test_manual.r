library(yaml)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(rtracklayer)
library(ShortRead)
library(Rsamtools)
library(QuasR)

# Set up variables for testing
wd <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/dev/lab_pipelines/SeEN_seq"
setwd(wd)
config <- yaml::read_yaml("config.yaml")
str(config)

# /////////////////////////////////////////////////////////////////////////////////////////////////
# region Test fastq search
# /////////////////////////////////////////////////////////////////////////////////////////////////
source(paste0(wd, "/scripts/find_fastq.r"))
sample_df_path <- config$sample_table_path
sample_df <- read.table(sample_df_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Locating fastq with prefix specified in sample_df
sample_df <- find_fastq(sample_df, config$fastq_dir)
# Annotating the library layout (i.e., single-end or paired-end) of each sample
sample_df <- sample_df_add_lib_layout(sample_df)

#region make QuasR input
source(paste0(wd, "/scripts/make_QuasR_input.r"))
# Transform the sample_df to QuasR input table
QuasR_df_list <- sample_df_to_QuasR_table(sample_df, output_prefix = "./results/QuasR/QuasR_sample", Return = TRUE, verbose = TRUE)


# region Run QuasR

source(paste0(wd, "/scripts/QuasR_wrapper.r"))

# Testing the base function
# Align and count reads
SE_count_mat <- QuasR_wrapper_base(
  sample_path = QuasR_df_list$input_SE,
  genome = config$ref_seq_path,
  config = config,
  verbose = TRUE
)

# Run QuasR for paired-end samples
PE_count_mat <- QuasR_wrapper_base(
  sample_path = QuasR_df_list$input_PE,
  genome = config$ref_seq_path,
  config = config,
  verbose = TRUE
)

# Testing the main wrapper function
count_mat <- QuasR_wrapper(config)

