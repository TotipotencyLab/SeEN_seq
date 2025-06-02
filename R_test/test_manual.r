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

# Plot libraries
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)

# Set up variables for testing
wd <- "/fs/pool/pool-toti-bioinfo/bioinfo/siwat_chad/dev/lab_pipelines/SeEN_seq"
wd <- "/Users/chad/Lab/dev/lab_pipelines/SeEN_seq"
wd <- getwd()
setwd(wd)
config <- yaml::read_yaml("config.yaml")
str(config)

sample_df <- read.table(config$sample_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ref_df <- read.table(config$ref_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

fastq_dir <- config$fastq_dir
fastq_ext_regex <- "((\\.fq)|(\\.fastq))(\\.gz)?$"
verbose=TRUE

ref_fa <- Biostrings::readDNAStringSet(config$ref_seq_path)
ref_gr <- GRanges(
  seqnames = names(ref_fa),
  ranges = IRanges(start = 1, end = width(ref_fa)),
)
# Because we count per each entry of fasta, there shouldn't be any duplicated names
names(ref_gr) <- as.character(seqnames(ref_gr))

# /////////////////////////////////////////////////////////////////////////////////////////////////
# region Test fastq search
# /////////////////////////////////////////////////////////////////////////////////////////////////
source(paste0(wd, "/scripts/find_fastq.r"))
# Locating fastq with prefix specified in sample_df
sample_df <- find_fastq(sample_df, config$fastq_dir)
# Annotating the library layout (i.e., single-end or paired-end) of each sample
sample_df <- sample_df_add_lib_layout(sample_df)

#region make QuasR input
source(paste0(wd, "/scripts/make_QuasR_input.r"))
# Transform the sample_df to QuasR input table

QuasR_dir <- paste0(config$result_dir, "/", config$project_name, "/QuasR")
QuasR_output_prefix <- paste0(QuasR_dir, "/", config$project_name)
QuasR_count_mat_path <- paste0(QuasR_output_prefix, "_count_matrix.txt")

QuasR_df_list <- sample_df_to_QuasR_table(sample_df, output_prefix = QuasR_output_prefix, Return = TRUE, verbose = TRUE)


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
# write.table(count_mat, file = paste0(wd, "/results/count_matrix.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
# Heatmap(count_mat, cluster_rows = FALSE, cluster_columns = FALSE)



# Test constructing SE object
source(paste0(wd, "/scripts/Summarized_SeEN_Experiment.r"))
valid_SumExp_inputs(count_mat = count_mat, sample_table=sample_df, ref_table=ref_df, behavior = "message")
SeEN_se <- Summarized_SeEN_Experiment(count_matrix=count_mat, sample_table=sample_df, ref_table=ref_df, config=config)

SeEN_se <- SeEN_se.library_size_norm(SeEN_se=SeEN_se, pseudo_count=1)
SeEN_se <- SeEN_se.enrichment(SeEN_se=SeEN_se, comparisons=config$compare_fraction)


# SeEN_se to Xlsx file ----
source(paste0(wd, "/scripts/write_SeEN_xlsx.r"))
SeEN_se.write_xlsx(SeEN_se, file=paste0(wd, "/results/", config$project_name, "/", config$project_name, "_SeEN_data.xlsx"))

source(paste0(wd, "/scripts/plot_SeEN_enrichment.r"))
colData(SeEN_se)
SeEN_se.plot_enrichment(SeEN_se)
SeEN_se.plot_enrichment(SeEN_se, color_by=NULL, group_by=NULL, facet_formula=NULL)
SeEN_se.plot_enrichment(SeEN_se, group_by="sample", color_by="condition")
SeEN_se.plot_enrichment(SeEN_se, group_by="sample", color_by="condition", facet_formula="condition~.")
SeEN_se.plot_enrichment(SeEN_se, group_by="sample", color_by="condition", average_group=TRUE)

SeEN_se.plot_enrichment(SeEN_se, group_by="condition", color_by="condition", average_group=TRUE)
SeEN_se.plot_enrichment(SeEN_se, group_by="condition", color_by="condition", average_group=TRUE, 
                        facet_formula="condition~.")
SeEN_se.plot_enrichment(SeEN_se, group_by="condition", color_by="condition", average_group=TRUE, 
                        NA_replace.x=-1, facet_formula="condition~.")
SeEN_se.plot_enrichment(SeEN_se, group_by="condition", color_by="condition", average_group=TRUE, 
                        NA_replace.x=-1, facet_formula="condition~.", 
                        show_plot_as="line and point", show_group_avg_as="errorbar and ribbon")
SeEN_se.plot_enrichment(SeEN_se, group_by="condition", color_by="condition", average_group=TRUE, 
                        NA_replace.x=-1, facet_formula="condition~.", 
                        show_plot_as="line", show_group_avg_as="errorbar")
SeEN_se.plot_enrichment(SeEN_se, plot_data=TRUE)
SeEN_se.plot_enrichment(SeEN_se, average_group=TRUE, plot_data=TRUE)


