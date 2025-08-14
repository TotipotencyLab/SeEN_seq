# Making input table for QuasR
sample_df_to_QuasR_table <- function(sample_df, output_prefix="./QuasR_sample", Return=FALSE, verbose=FALSE) {
  # Dependencies:
  #  - find_fastq.r: sample_df_add_lib_layout function
  
  # Corner cases:
  # 1) sample_df has mixed of paired and single-end library layout
  # 2) sample_df has samples without fastq files
  
  # Check library layout
  # This should add new column: `lib_layout` annotating whether each sample is single-end ("single") or paired-end ("paired")
  # Sample without fastq files will be NA
  sample_df <- sample_df_add_lib_layout(sample_df)

  # lib_layout column is in a way
  if (any(is.na(sample_df$lib_layout))) {
    tmp_df <- sample_df[is.na(sample_df$lib_layout), , drop = FALSE]
    sample_df <- sample_df[!is.na(sample_df$lib_layout), , drop = FALSE]
    warning(
      "The input contains samples without fastq files. These will be ignored",
      paste0("Samples without fastq files:\n", paste0(tmp_df$sample[is.na(tmp_df$lib_layout)], collapse = "\n"))
    )
  }
  
  # split single-end and paired-end samples
  sample_SE_df <- sample_df[sample_df$lib_layout == "single", ]
  sample_PE_df <- sample_df[sample_df$lib_layout == "paired", ]
  if (nrow(sample_SE_df) > 0 && nrow(sample_PE_df) > 0) {
    message("This project contains both samples with single-end and paired-end fastq files. Making separated QuasR input tables for each library layout.")
  }
  
  if (nrow(sample_SE_df) > 0 || nrow(sample_PE_df) > 0) {
    out_dir <- dirname(output_prefix)
    if(!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
  }

  # Placeholders for outputs
  QuasR_input_SE_file <- character(0)
  QuasR_input_PE_file <- character(0)
  QuasR_SE_df <- data.frame()
  QuasR_PE_df <- data.frame()


  # NOTE: QuasR doesn't seems to understand the symbolic link, make sre we feed in the absolute path here.
  if(nrow(sample_SE_df) > 0){
    # Input table for SE required two column: FileName and SampleName
    QuasR_SE_df <- sample_SE_df %>%
      select(FileName = fastq_SE, SampleName = unique_name) %>%
      rowwise() %>%
      mutate(FileName = tools::file_path_as_absolute(FileName)) %>%
      ungroup()
    QuasR_input_SE_file <- paste0(output_prefix, "_sample_SE.txt")
    if(verbose) message("Writing QuasR input table for single-end samples to ", QuasR_input_SE_file)
    write.table(QuasR_SE_df, file = QuasR_input_SE_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  if(nrow(sample_PE_df) > 0){
    # Input table for PE required three column: FileName1, FileName2 and SampleName
    QuasR_PE_df <- sample_PE_df %>%
      select(FileName1 = fastq_R1, FileName2 = fastq_R2, SampleName = unique_name) %>%
      rowwise() %>%
      mutate(
        FileName1 = tools::file_path_as_absolute(FileName1),
        FileName2 = tools::file_path_as_absolute(FileName2)
      ) %>%
      ungroup()
    QuasR_input_PE_file <- paste0(output_prefix, "_sample_PE.txt")
    if(verbose) message("Writing QuasR input table for paired-end samples to ", QuasR_input_PE_file)
    write.table(QuasR_PE_df, file = QuasR_input_PE_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  if (Return) {
    return(list(
      input_SE = QuasR_input_SE_file,
      input_PE = QuasR_input_PE_file,
      SE = QuasR_SE_df,
      PE = QuasR_PE_df
    ))
  }
}
