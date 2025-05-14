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
sample_df <- find_fastq(sample_df, config$fastq_dir)
sample_df <- sample_df_add_lib_layout(sample_df)