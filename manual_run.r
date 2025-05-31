# Trying to run this manually without sourcing other script

# //////////////////////////////////////////////////////////////////////////////
# Set up ====
# //////////////////////////////////////////////////////////////////////////////

## Load Library ---------------------------------------------------------------
# General
library(yaml)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(parallel)

# Genomics
require(GenomicRanges)
require(Biostrings)
library(ShortRead)
library(Rsamtools)
library(QuasR)
# library(GenomicFeatures)
# library(rtracklayer)

# Plot libraries
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)


## key variable ---------------------------------------------------------------
wd <- getwd()
config <- yaml::read_yaml("config.yaml")
sample_df_path <- config$sample_table_path
ref_table_path <- config$ref_table_path

sample_df <- read.table(sample_df_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ref_df <- read.table(ref_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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


# ..........................................................................----
# //////////////////////////////////////////////////////////////////////////////
#region QuasR ====
# //////////////////////////////////////////////////////////////////////////////
# GOAL: getting the `count_mat` variable 
#       -- count matrix of reads from each sample mapped to reference construct.

# Set up expected QuasR output file paths
QuasR_dir <- paste0(config$result_dir, "/", config$project_name, "/QuasR")
QuasR_output_prefix <- paste0(QuasR_dir, "/", config$project_name)
QuasR_count_mat_path <- paste0(QuasR_output_prefix, "_count_matrix.txt")

if(file.exists(QuasR_count_mat_path)){
  message("QuasR count matrix already exists at: ", QuasR_count_mat_path, "\nRead count matrix from this file instad of re-running QuasR.")
  count_mat <- read.table(QuasR_count_mat_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
}else{
  # Running QuasR
  
  ## Find FastQ ===============================================================
  try({sample_df <- as_tibble(sample_df)})
  sample_df <- mutate(sample_df, fastq_R1 = NA, fastq_R2 = NA, fastq_SE = NA, fastq_unknown = NA)
  
  f <- list.files(fastq_dir, recursive = TRUE, full.name = TRUE, pattern = fastq_ext_regex)
  # placeholder: keep row index of strange sample-fastq mappping
  unfound_sample <- integer()
  unknown_sample <- integer()
  for (i in seq_along(sample_df$fastq_prefix)) {
    # TODO: implement both startsWith and grepl (i.e., fastq_prefix is regex)
    f_hit <- base::startsWith(basename(f), sample_df$fastq_prefix[i])
    if (any(f_hit)) {
      # assign found fastq files to the sample_df
      fq_hit <- f[f_hit]
      if (length(f[f_hit]) == 1) {
        sample_df$fastq_SE[i] <- fq_hit
      } else {
        # Check whether it has R1 and R2 pattern
        R1_hit_idx <- grep(basename(fq_hit), pattern = paste0("[_.]R?1", fastq_ext_regex))
        R2_hit_idx <- grep(basename(fq_hit), pattern = paste0("[_.]R?2", fastq_ext_regex))
        # if they are not overlap, and each have only one hit
        if (length(intersect(R1_hit_idx, R2_hit_idx)) == 0 && length(R1_hit_idx) == 1 && length(R2_hit_idx) == 1) {
          sample_df$fastq_R1[i] <- fq_hit[R1_hit_idx]
          sample_df$fastq_R2[i] <- fq_hit[R2_hit_idx]
          fq_hit <- fq_hit[-unique(c(R1_hit_idx, R2_hit_idx))] # Removed from the current hit
        }
        
        # If there are still some left, assign them to fastq_unknown
        if (length(fq_hit) > 0) {
          sample_df$fastq_unknown[i] <- paste0(fq_hit, collapse = "; ")
          unknown_sample <- c(unknown_sample, i)
        }
      }
    } else {
      unfound_sample <- c(unfound_sample, i)
    }
  }
  
  # Warning of strange sample-fastq mapping
  if(length(unfound_sample) == nrow(sample_df)) {
    warning(
      "Couldn't find any fastq files for all samples. Please check if the directory for the fastq files is correct.\n",
      "Received fastq directory: ", fastq_dir, "\n"
    )
  } else if (length(unfound_sample) > 0) {
    tmp_df <- sample_df[unfound_sample, , drop = FALSE]
    warning("Some samples are not found in fastq directory:\n", 
            paste0(tmp_df$sample, "  -  ", tmp_df$fastq_prefix, collapse = "\n"))
  }
  if (length(unknown_sample) > 0) {
    tmp_df <- sample_df[unknown_sample, , drop = FALSE]
    warning("Some samples have unknown fastq files associated to them:\n", 
            paste0(tmp_df$sample, "  -  ", tmp_df$fastq_prefix, collapse = "\n"))
  }
  
  
  # Remove newly created columns (i.e., by this function) with all NA values
  all_na_col_idx <- sapply(sample_df, function(x) all(is.na(x)))
  target_col_idx <- colnames(sample_df) %in% c("fastq_R1", "fastq_R2", "fastq_SE", "fastq_unknown")
  sample_df <- sample_df[, !(all_na_col_idx & target_col_idx), drop = FALSE]
  # sample_df <- sample_df %>% select(where(~ any(!is.na(.))))
  
  
  ## Annotate library layout ==================================================
  # Assuming that sample_df has mixed of paired and single-end fastq files
  # Annotating library layout (single- or paired-end (SE, or PE)) for each sample
  proj_lib_layout <- c(single = FALSE, paired = FALSE)
  SE_cols <- "fastq_SE"
  PE_cols <- c("fastq_R1", "fastq_R2")
  
  PE_row_idx <- integer()
  if (all(PE_cols %in% colnames(sample_df))) {
    PE_row_idx <- which(rowSums(is.na(sample_df[, PE_cols])) == 0)
    if (length(PE_row_idx) > 0) proj_lib_layout["paired"] <- TRUE
  }
  SE_row_idx <- integer()
  if (all(SE_cols %in% colnames(sample_df))) {
    SE_row_idx <- which(rowSums(is.na(sample_df[, SE_cols])) == 0)
    if (length(SE_row_idx) > 0) proj_lib_layout["single"] <- TRUE
  }
  
  if(length(intersect(SE_row_idx, PE_row_idx)) > 0) {
    mixed_layout_idx <- intersect(SE_row_idx, PE_row_idx)
    tmp_df <- sample_df[mixed_layout_idx, , drop=FALSE]
    warning(
      "Some samples have both paired (PE) and single-end (SE) fastq files associated to them. The PE will take priority here.\nAffected samples:\n",
      paste0(tmp_df$sample, "  -  ", tmp_df$fastq_prefix, collapse = "\n")
    )
  }
  
  sample_df$lib_layout <- NA
  sample_df$lib_layout[SE_row_idx] <- "single"
  sample_df$lib_layout[PE_row_idx] <- "paired"
  
  
  if (any(is.na(sample_df$lib_layout))) {
    tmp_df <- sample_df[is.na(sample_df$lib_layout), , drop = FALSE]
    sample_df <- sample_df[!is.na(sample_df$lib_layout), , drop = FALSE]
    warning(
      "The input contains samples without fastq files. These will be ignored",
      paste0("Samples without fastq files:\n", paste0(tmp_df$sample[is.na(tmp_df$lib_layout)], collapse = "\n"))
    )
  }
  
  
  ## Run QuasR ================================================================
  ### Prepare inputs ----------------------------------------------------------
  # split single-end and paired-end samples
  sample_SE_df <- sample_df[sample_df$lib_layout == "single", ]
  sample_PE_df <- sample_df[sample_df$lib_layout == "paired", ]
  if (nrow(sample_SE_df) > 0 && nrow(sample_PE_df) > 0) {
    message("This project contains both samples with single-end and paired-end fastq files. Making separated QuasR input tables for each library layout.")
  }
  
  if(!dir.exists(QuasR_dir)) {
    message("Creating output directory for QuasR results: ", out_dir)
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Make the input table for QuasR
  QuasR_input_SE_file <- character(0)
  QuasR_input_PE_file <- character(0)
  
  # NOTE: QuasR doesn't seems to understand the symbolic link, make sre we feed in the absolute path here.
  if(nrow(sample_SE_df) > 0){
    # Input table for SE required two column: FileName and SampleName
    QuasR_SE_df <- sample_SE_df %>%
      select(FileName = fastq_SE, SampleName = unique_name) %>%
      rowwise() %>%
      mutate(FileName = tools::file_path_as_absolute(FileName)) %>%
      ungroup()
    QuasR_input_SE_file <- paste0(QuasR_output_prefix, "_sample_SE.txt")
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
    QuasR_input_PE_file <- paste0(QuasR_output_prefix, "_sample_PE.txt")
    if(verbose) message("Writing QuasR input table for paired-end samples to ", QuasR_input_PE_file)
    write.table(QuasR_PE_df, file = QuasR_input_PE_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  
  # config$QuasR_params
  QuasR_params = config$QuasR_params
  # Placeholder 
  QuasR_aln_dir <- QuasR_dir
  QuasR_cache_dir <- QuasR_dir
  
  # make temp dir
  if (!is.null(QuasR_params$aln_dir)) {
    QuasR_aln_dir <- paste0(QuasR_dir, "/", QuasR_params$aln_dir)
    if (!dir.exists(QuasR_aln_dir)) {
      dir.create(QuasR_aln_dir, recursive = TRUE)
    }
  }
  
  if (!is.null(QuasR_params$cache_dir)) {
    QuasR_cache_dir <- paste0(QuasR_dir, "/", QuasR_params$cache_dir)
    if (!dir.exists(QuasR_cache_dir)) {
      dir.create(QuasR_cache_dir, recursive = TRUE)
    }
  }
  
  ### Run QuasR ---------------------------------------------------------------
  # cluster parameters
  cl <- NULL # QuasR default
  # request_core <- QuasR_params$n_core
  if (is.null(QuasR_params$n_core)) QuasR_params$n_core <- 1
  avail_core <- parallel::detectCores()
  use_core <- min(QuasR_params$n_core, avail_core)
  if (use_core > 1) {
    # NOTE: "FORK" is only available on Unix-like systems
    # This may cause bug for Windows users
    if (verbose)  message("Using ", use_core, " cores for QuasR mapping.")
    cl <- makeCluster(use_core, type = "FORK") 
  }
  
  # Running this as a loop to keep the code structure for SE and PE to be the same
  QuasR_inputs <- c(SE = QuasR_input_SE_file, PE = QuasR_input_PE_file)
  # placeholders for results
  count_mat_list <- list()
  for(i in seq_along(QuasR_inputs)){
    if(!file.exists(QuasR_inputs[i])) next
    
    if(verbose) message("Ruuning QuasR alignment for ", names(QuasR_inputs)[i], " samples using input file: ", QuasR_inputs[i])
    
    # Alignment
    aln_obj <- qAlign(
      # sampleFile = QuasR_df_list$input_PE,
      sampleFile = QuasR_inputs[i],
      genome = config$ref_seq_path,
      aligner = QuasR_params$aligner,
      projectName = config$project_name,
      alignmentsDir = QuasR_aln_dir,
      alignmentParameter = QuasR_params$aln_params,
      splicedAlignment = QuasR_params$spliced_aln,
      cacheDir = QuasR_cache_dir,
      clObj = cl
    )
    m <- qCount(proj = aln_obj, query = ref_gr, clObj = cl)
    m <- m[ , !colnames(m) %in% c("width"), drop=FALSE] # Remove the default "width" column
    count_mat_list[[ names(QuasR_inputs)[i] ]] <- m
  }
  
  # Getting count matrix
  if(length(count_mat_list) > 1){
    # Combine matrix
    # Asumption:
    # 1) no duplicated sample names
    # 2) All rownames are matches between all results (should be the case because they used the same reference)
    count_mat <- NULL
    for(i in seq_along(count_mat_list)){
      if(i == 1){
        count_mat <- count_mat_list[[i]]
      }else{
        count_mat <- cbind(count_mat, count_mat_list[[i]])
      }
    }
    
  }else if(length(count_mat_list) == 1){
    count_mat <- count_mat_list[[1]]
    
  }else{
    stop("No count matrix produced from the QuasR alignment. Please check the input files and parameters.")
  }
  
  write.table(count_mat, file = QuasR_count_mat_path, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  # Stop cluster
  if (inherits(cl, "cluster")) {
    stopCluster(cl)
  }
} # End if(file.exists(QuasR_count_mat_path)) -- check if QuasR count matrix already exists



# .........................................................................----
# /////////////////////////////////////////////////////////////////////////////
#region SummarizedExperiment ====
# /////////////////////////////////////////////////////////////////////////////
# GOAL: Gathering count, sample_table, ref_table to SummarizedExperiment (SE) object

## Check inputs ----
# Check whether row/column names of count matrix matches names of reference and sample
unfound_ref <- rownames(count_mat)[!rownames(count_mat) %in% ref_df$ref_name]
if(length(unfound_ref) > 0){
  if(length(unfound_ref) == nrow(count_mat)){
    stop("None of the row names in the count matrix matches with the reference table. ")
  }else{
    unfound_txt <- paste0(unfound_ref, collapse = ", ")
    if(length(unfound_ref) > 5){
      unfound_txt <- paste0(paste0(head(unfound_ref, 5), collapse=", "), " ... ")
    }
    stop("Some of the row names in the count matrix do not match with the reference table:\n",
         "Missing ref (", length(unfound_ref), "/", nrow(count_mat), "):\n\t", unfound_txt)
  }
}

unfound_sample <- colnames(count_mat)[!colnames(count_mat) %in% sample_df$unique_name]
if(length(unfound_sample) > 0){
  if(length(unfound_sample) == ncol(count_mat)){
    stop("None of the column names in the count matrix matches with the sample table. ")
  }else{
    unfound_txt <- paste0(unfound_sample, collapse = ", ")
    if(length(unfound_sample) > 5){
      unfound_txt <- paste0(paste0(head(unfound_sample, 5), collapse=", "), " ... ")
    }
    stop("Some of the column names in the count matrix do not match with the sample table:\n",
         "Missing sample (", length(unfound_sample), "/", ncol(count_mat), "):\n\t", unfound_txt)
  }
}

# Making SE object ------------------------------------------------------------
## Row data = reference ----
# Make sure row/column names are matches with those of ref/sample tables
row_df <- ref_df %>% 
  tibble::column_to_rownames(var="ref_name")
row_df <- row_df[rownames(count_mat), , drop = FALSE]

## Col data = samples ----
tmp_col_df <- sample_df %>% 
  tibble::column_to_rownames(var="unique_name")
tmp_col_df <- tmp_col_df[colnames(count_mat), , drop = FALSE]

col_df <- tmp_col_df %>% 
  dplyr::select(-fraction, -starts_with("fastq_")) %>% 
  distinct() %>% 
  `rownames<-`(NULL) %>% 
  tibble::column_to_rownames(var="sample")

## Split counts matrix ----
fraction_type <- unique(tmp_col_df$fraction)
unq_sample <- unique(tmp_col_df$sample)

# Set up a placeholder for any count matrix extracted from each fraction type
zeros_mat <- matrix(0, nrow=nrow(count_mat), ncol=length(unq_sample), 
                  dimnames=list(rownames(count_mat), unq_sample))

count_mat_fraction <- list()
for(i in seq_along(fraction_type)){
  cur_name <- paste0("counts_", fraction_type[i])
  name_map <- tmp_col_df %>% 
    dplyr::filter(fraction == fraction_type[i]) %>% 
    tibble::rownames_to_column(var="unique_name") %>% 
    dplyr::pull(sample, name=unique_name)
  
  if(any(duplicated(name_map))){
    warning("Duplicated sample name found for the current fraction: ", fraction_type[i], 
            "\nThis may cause issues when extracting count matrix for this fraction.")
  }
  
  m <- count_mat[ , names(name_map), drop=FALSE]
  colnames(m) <- name_map[colnames(m)]
  missing_samples <- setdiff(unq_sample, colnames(m))
  if(length(missing_samples) > 0){
    # Add missing samples with zero counts
    m <- cbind(m, zeros_mat[ , missing_samples, drop=FALSE])
  }
  count_mat_fraction[[cur_name]] <- m[ , unq_sample, drop=FALSE] # reorder the columns to match with unq_sample
}
# str(count_mat_fraction, max.level=1)

SeEN_se <- SummarizedExperiment(
  assays = count_mat_fraction, 
  colData = col_df, 
  rowData = row_df
)


# .........................................................................----
# /////////////////////////////////////////////////////////////////////////////
#region Analysis ====
# /////////////////////////////////////////////////////////////////////////////

