# /////////////////////////////////////////////////////////////////////////////////////////////////
#region QuasR_wrapper_base
# /////////////////////////////////////////////////////////////////////////////////////////////////

QuasR_wrapper_base <- function(sample_path, genome, config, verbose=TRUE){
  if(length(sample_path) == 0){
    # do this because we want the wrapper to be able to handle both single-end and paired-end samples,
    #   which can be empty if no fastq files are found for that library layout type
    return(NULL)
  }
  # Base version of QuasR wrapper
  require(QuasR)
  # require(ShortRead)
  require(GenomicRanges)
  require(Biostrings)
  
  # Set up ----------------------------------------------------------------------------------------
  params <- config$QuasR_params

  # make temp dir
  if (!is.null(params$aln_dir)) {
    if (!dir.exists(params$aln_dir)) {
      dir.create(params$aln_dir, recursive = TRUE)
    }
  }

  if (!is.null(params$cache_dir)) {
    if (!dir.exists(params$cache_dir)) {
      dir.create(params$cache_dir, recursive = TRUE)
    }
  }

  # cluster parameters
  cl <- NULL # QuasR default
  request_core <- params$n_core
  if (is.null(request_core)) request_core <- 1
  avail_core <- parallel::detectCores()
  use_core <- min(request_core, avail_core)
  if (use_core > 1) {
    # NOTE: "FORK" is only available on Unix-like systems
    # This may cause bug for Windows users
    if (verbose)  message("Using ", use_core, " cores for QuasR mapping.")
    cl <- makeCluster(use_core, type = "FORK") 
  }

  # Run QuasR Alignment ---------------------------------------------------------------------------
  aln_obj <- qAlign(
    # sampleFile = QuasR_df_list$input_PE,
    sampleFile = sample_path,
    genome = config$ref_seq_path,
    aligner = params$aligner,
    projectName = params$project_name,
    alignmentsDir = params$aln_dir,
    alignmentParameter = params$aln_params,
    splicedAlignment = params$spliced_aln,
    cacheDir = params$cache_dir,
    clObj = cl
  )

  # Count aligned reads per construct (entry of FastA) --------------------------------------------
  # ref_fa <- readFasta(config$ref_seq_path)
  ref_fa <- Biostrings::readDNAStringSet(config$ref_seq_path)
  ref_gr <- GRanges(
    seqnames = names(ref_fa),
    ranges = IRanges(start = 1, end = width(ref_fa)),
  )
  # Because we count per each entry of fasta, there shouldn't be any duplicated names
  names(ref_gr) <- as.character(seqnames(ref_gr))

  count_mat <- qCount(proj = aln_obj, query = ref_gr, clObj = cl)

  if (inherits(cl, "cluster")) {
    stopCluster(cl)
  }

  return(count_mat)
}



# /////////////////////////////////////////////////////////////////////////////////////////////////
#region QuasR_wrapper
# /////////////////////////////////////////////////////////////////////////////////////////////////

QuasR_wrapper <- function(config){
  # Top level wrapper for QuasR, only need config list as input

  if(!is.list(config)){
    # Attempt to read yaml file
    require(yaml)
    config <- yaml::read_yaml(config)
  }

  # Get sample input
  sample_df_path <- config$sample_table_path
  sample_df <- read.table(sample_df_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Locating fastq with prefix specified in sample_df
  sample_df <- find_fastq(sample_df, config$fastq_dir)
  # Annotating the library layout (i.e., single-end or paired-end) of each sample
  sample_df <- sample_df_add_lib_layout(sample_df)
  
  # Extract input for QuasR for each library layout
  QuasR_df_list <- sample_df_to_QuasR_table(sample_df, output_prefix = "./results/QuasR/QuasR_sample", Return = TRUE, verbose = TRUE)

  # Alignment and count reads ---------------------------------------------------------------------
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

  # Summarize count matrix ------------------------------------------------------------------------
  count_mat <- NULL
  if (!is.null(SE_count_mat)) {
    count_mat <- SE_count_mat
  }

  if (!is.null(PE_count_mat)) {
    if (is.null(count_mat)) {
      count_mat <- PE_count_mat
    } else {
      count_mat <- cbind(count_mat, PE_count_mat)
    }
  }

  return(count_mat)
}