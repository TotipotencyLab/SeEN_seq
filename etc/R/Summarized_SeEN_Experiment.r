# Gathering count, sample_table, ref_table to SummarizedExperiment (SE) object

valid_SumExp_inputs <- function(count_matrix, ref_table, sample_table, behavior="stop") {
  # Define the behavior when inputs are not matching
  if (!behavior %in% c("stop", "warn", "message", "ignore")) {
    warning("Invalid behavior argument. It should be one of 'stop', 'warn', 'message', or 'ignore'. Switch to 'stop'.")
    behavior <- "stop"
  }
  response_fn = switch(behavior,
    stop = stop,
    warn = warning,
    message = message,
    ignore = function(...) NULL
  )

  # count_matrix vs ref_table ----
  # Check whether row/column names of count matrix matches names of reference and sample
  unfound_ref <- rownames(count_matrix)[!rownames(count_matrix) %in% ref_table$ref_name]
  if (length(unfound_ref) > 0) {
    if (length(unfound_ref) == nrow(count_matrix)) {
      response_fn("None of the row names in the count matrix matches with the reference table. ")
      return(FALSE)
    } else {
      unfound_txt <- paste0(unfound_ref, collapse = ", ")
      if (length(unfound_ref) > 5) {
        unfound_txt <- paste0(paste0(head(unfound_ref, 5), collapse = ", "), " ... ")
      }
      response_fn(
        "Some of the row names in the count matrix do not match with the reference table:\n",
        "Missing ref (", length(unfound_ref), "/", nrow(count_matrix), "):\n\t", unfound_txt
      )
      return(FALSE)
    }
  }

  # count_matrix vs sample_table ----
  unfound_sample <- colnames(count_matrix)[!colnames(count_matrix) %in% sample_table$unique_name]
  if (length(unfound_sample) > 0) {
    if (length(unfound_sample) == ncol(count_matrix)) {
      response_fn("None of the column names in the count matrix matches with the sample table. ")
      return(FALSE)
    } else {
      unfound_txt <- paste0(unfound_sample, collapse = ", ")
      if (length(unfound_sample) > 5) {
        unfound_txt <- paste0(paste0(head(unfound_sample, 5), collapse = ", "), " ... ")
      }
      response_fn(
        "Some of the column names in the count matrix do not match with the sample table:\n",
        "Missing sample (", length(unfound_sample), "/", ncol(count_matrix), "):\n\t", unfound_txt
      )
      return(FALSE)
    }
  }

  return(TRUE)
}


# .........................................................................----
# /////////////////////////////////////////////////////////////////////////////
#region Main function ----
# /////////////////////////////////////////////////////////////////////////////

Summarized_SeEN_Experiment <- function(count_matrix, sample_table, ref_table, config=NULL){
  require(SummarizedExperiment)
  require(dplyr)
  require(tibble)

  # Check input validity
  input_pass <- valid_SumExp_inputs(count_matrix=count_matrix, sample_table=sample_table, ref_table=ref_table, behavior = "message")
  if (!input_pass) {
    stop("Invalid inputs for SummarizedExperiment. Please check the count matrix, sample table, and reference table.")
  }

  # Row data ------------------------------------------------------------------
  # Row data = reference sequences
  row_df <- ref_table %>%
    tibble::column_to_rownames(var = "ref_name")
  row_df <- row_df[rownames(count_matrix), , drop = FALSE]

  # Col data -----------------------------------------------------------------
  # NOTE: 
  # Sample table is expected to have unique_name and sample columns
  # The `unique_name` is corresponding to the individual fraction of the SeEN-seq experiment, hence individual fastq files.
  #     This column should match the column names of the count matrix.
  # The `sample` column is the actual sample names (i.e., each lane of the SeEN-seq experiment).
  #     Thus it can be corresponded to multiple fractions.
  #     The sample name will be used for the column names of the SummarizedExperiment object.
  
  tmp_col_df <- sample_table %>%
    tibble::column_to_rownames(var = "unique_name")
  tmp_col_df <- tmp_col_df[colnames(count_matrix), , drop = FALSE]
  unq_sample <- unique(tmp_col_df$sample)

  # TODO: Improve this chunk to dynamically detect the columns that 
  #       didn't cause duplicated `sample` columns after `dplyr::distinct()` function
  
  # Get number of unique value in each other column
  col_unq_val_count <- data.frame()
  for(i in seq_along(unq_sample)){
    df <- tmp_col_df[tmp_col_df$sample == unq_sample[i], , drop=FALSE] %>% 
      lapply(FUN=function(x){length(unique(x))}) %>% 
      as.data.frame()
    rownames(df) <- unq_sample[i]
    col_unq_val_count <- rbind(col_unq_val_count, df)
  }
  valid_col_flag <- sapply(col_unq_val_count, FUN=function(x){all(x == 1)})
  valid_col <- names(valid_col_flag)[valid_col_flag]
  
  if(any(!valid_col_flag)){
    message(
      "The following columns will not be included in the sample table of SummarizedExperiment object\n\t", 
      paste0(names(valid_col_flag)[!valid_col_flag], collapse=", ")
    )
  }
  
  col_df <- tmp_col_df %>%
    dplyr::select(all_of(valid_col)) %>% # 'sample' should be one of them
    distinct() %>%
    `rownames<-`(NULL) %>%
    tibble::column_to_rownames(var = "sample")

  # counts matrix -------------------------------------------------------------
  # NOTE:
  # Count matrix is now a mix of different fractions of the SeEN-seq experiment (individual fastq file).
  # For each fraction name, extract the count matrix and store it as a list (assays for SummarizedExperiment).
  # In case certain fraction doesn't have certain samples, the columns of those samples will be filled with zeros.

  fraction_names <- unique(tmp_col_df$fraction)
  # for filling empty samples in the future
  zeros_mat <- matrix(0,
    nrow = nrow(count_matrix), ncol = length(unq_sample),
    dimnames = list(rownames(count_matrix), unq_sample)
  )

  count_matrix_fraction <- list() # placeholder for assays of SummarizedExperiment
  for (i in seq_along(fraction_names)) {
    cur_name <- paste0("counts_", fraction_names[i])
    
    # Dictionary of mapping unique name to sample (new column names of count matrix)
    name_map <- tmp_col_df %>%
      dplyr::filter(fraction == fraction_names[i]) %>%
      tibble::rownames_to_column(var = "unique_name") %>%
      dplyr::pull(sample, name = unique_name)

    if (any(duplicated(name_map))) {
      warning(
        "Duplicated sample name found for the current fraction: ", fraction_names[i],
        "\nThis may cause issues when extracting count matrix for this fraction."
      )
    }

    # Extract the count matrix for the current fraction
    m <- count_matrix[, names(name_map), drop = FALSE]
    colnames(m) <- name_map[colnames(m)]
    
    # Finds missing samples for this fraction and fills with zero counts
    missing_samples <- setdiff(unq_sample, colnames(m))
    if (length(missing_samples) > 0) {
      m <- cbind(m, zeros_mat[, missing_samples, drop = FALSE])
    }
    # if (!is.null(config$pseudo_count)) {
    #   m <- m + config$pseudo_count
    # }

    count_matrix_fraction[[cur_name]] <- m[, unq_sample, drop = FALSE] # reorder the columns to match with unq_sample
  }

  # Metadata ------------------------------------------------------------------
  meta_list <- list(fraction_names = fraction_names)
  if (!is.null(config)) {
    meta_list$project_name <- config$project_name
    meta_list$ref_seq_path <- config$ref_seq_path
    meta_list$ref_table_path <- config$ref_table_path
    meta_list$sample_table_path <- config$sample_table_path
  }
  
  # SumExp --------------------------------------------------------------------
  SeEN_se <- SummarizedExperiment(
    assays = count_matrix_fraction,
    colData = col_df,
    rowData = row_df,
    metadata = meta_list
  )

  return(SeEN_se)
}

# .........................................................................----
# /////////////////////////////////////////////////////////////////////////////
#region Process Counts ----
# /////////////////////////////////////////////////////////////////////////////

SeEN_se.library_size_norm <- function(SeEN_se, pseudo_count=1, output_assay_prefix="libnorm_counts", quiet=FALSE){
  if(output_assay_prefix == ""){
    warning("The output assay prefix is empty, revert to default 'libnorm_counts'.")
    output_assay_prefix <- "libnorm_counts"
  }
  fraction_names <- metadata(SeEN_se)$fraction_names
  
  new_assay_names <- c()
  for (i in seq_along(fraction_names)) {
    use_assay_name <- paste0("counts_", fraction_names[i])
    if(!use_assay_name %in% assayNames(SeEN_se)){
      warning("Couldn't find the assay ", use_assay_name, ". Skip normalizing this fraction.")
      next
    }
    m <- assay(SeEN_se, paste0("counts_", fraction_names[i]))
    m <- m + pseudo_count # Adding pseudo count
    lib_size_mat <- matrix(colSums(m), nrow = nrow(m), ncol = ncol(m), byrow = TRUE)
    libnorm_m <- m / lib_size_mat
    cur_assay_name <- paste0(output_assay_prefix, "_", fraction_names[i])
    new_assay_names <- c(new_assay_names, cur_assay_name)
    assay(SeEN_se, cur_assay_name) <- libnorm_m
  }
  
  if(length(new_assay_names) > 0){
    if(!quiet){
      cat("Adding these new assays to the input: ", paste0(new_assay_names, collapse=", "))
    }
  }else{
    warning("No new assays added to the input")
  }

  return(SeEN_se)
}

SeEN_se.enrichment <- function(SeEN_se, comparisons, pseudo_count = 0,
                               from_assay_perfix = "libnorm_counts", output_assay_prefix = "enrichment") {
  # pseudo_count is expected to be added already from the library size normalization step.
  #     Adding this argument just in case when we want to change when pseudo count is added in the pipeline.
  
  # Validate inputs -----------------------------------------------------------
  # NOTE: on comparisons:
  # `comparisons` is a list that should be corresponded to `config$compare_fraction`
  # Expect comparisons to be something like this:
  # comparisons <- list(
  #   "test_fraction_1" = c("control_fraction_1", "control_fraction_2"),
  #   "test_fraction_2" = c("control_fraction_1")
  # )
  # test_fraction_1, test_fraction_2, control_fraction_1, and control_fraction_2 are part of the fraction_names

  if (!is.list(comparisons)) {
    stop("The comparisons should be a list of named vectors, where each name is the test fraction and the vector contains control fractions.")
  }

  requested_fractions <- unique(c(names(comparisons), unlist(comparisons)))
  fraction_names <- metadata(SeEN_se)$fraction_names
  if (!all(requested_fractions %in% fraction_names)) {
    unfound_fraction <- requested_fractions[!requested_fractions %in% fraction_names]
    warning(
      "Some of the requested fractions are not found in the available fractions of the SeEN object. Please check the fraction names.",
      "\nAvaliable fractions: ", paste0(fraction_names, collapse = ", "),
      "\nRequested fractions: ", paste0(requested_fractions, collapse = ", "),
      "\nUnfound fractions:   ", paste0(unfound_fraction, collapse = ", ")
    )
  }

  if (output_assay_prefix == "") {
    warning("The output assay prefix is empty, revert to default 'enrichment'.")
    output_assay_prefix <- "enrichment"
  }

  # Enrichment ----------------------------------------------------------------
  inf_warning <- FALSE
  for (i in seq_along(comparisons)) {
    test_frac <- names(config$compare_fraction)[i]
    for (j in seq_along(config$compare_fraction[[i]])) {
      ctrl_frac <- config$compare_fraction[[i]][j]
      cur_assay_name <- paste0(output_assay_prefix, "_", test_frac, "_vs_", ctrl_frac)
      if (!test_frac %in% fraction_names || !ctrl_frac %in% fraction_names) {
        warning("Fraction ", test_frac, " or ", ctrl_frac, " not found in the fractions of the input SeEN_se object. Skipping enrichment calculation for this pair.")
        next
      }
      test_m <- assay(SeEN_se, paste0(from_assay_perfix, "_", test_frac))
      ctrl_m <- assay(SeEN_se, paste0(from_assay_perfix, "_", ctrl_frac))
      # Calculating Enrichment
      enrichment_m <- log2((test_m + pseudo_count) / (ctrl_m + pseudo_count))
      if (any(is.infinite(enrichment_m))) inf_warning = TRUE
      assay(SeEN_se, cur_assay_name) <- enrichment_m
    }
  }

  if (inf_warning) {
    msg <- "Infinite values found in the enrichment matrix."
    if (pseudo_count == 0) {
      msg <- paste0(msg, " Please consider using `pseudo_count` argument to avoid this issue.")
    }
    warning(msg)
  }

  return(SeEN_se)
}
