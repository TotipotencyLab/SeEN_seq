
## Chunk 0 ----
t0 <- Sys.time()
required_pkg <- c("dplyr", "tibble", "tidyr", "parallel", "QuasR", "SummarizedExperiment", "ggplot2")
sapply(required_pkg, library, character.only=TRUE)

## Chunk 1 ----
# General analysis project setting
config <- list(
  project_name = "Test_SeEN_seq",
  output_dir = "./results/Test_SeEN_seq",
  sample_table_path = "example/sample_metadata_PE.tsv",
  ref_table_path = "./example/ref_metadata.tsv",
  ref_fasta_path = "./example/ref_sequence.fasta",
  n_core=8,
  enrichment=list(c("Bound", "Unbound")), # test vs background
  log2_pseudo_count=1
)


# Make sure output folder exist
QuasR_dir <- paste0(config$output_dir, "/QuasR")
plot_dir <- paste0(config$output_dir, "/plot")
dir.create(config$output_dir, showWarnings=FALSE)
dir.create(QuasR_dir, showWarnings=FALSE)
dir.create(paste0(QuasR_dir, "/cache"), showWarnings=FALSE)
dir.create(plot_dir, showWarnings=FALSE)


## Chunk 2 ----
sample_df <- read.table(config$sample_table_path, sep="\t", header=TRUE,stringsAsFactors=FALSE)
# Extract relavant columns for QuasR
QuasR_fq_col <- list(SE="FileName", PE=c("FileName1", "FileName2"))
QuasR_fq_col_flag <- sapply(QuasR_fq_col, FUN=function(x){all(x %in% colnames(sample_df))})
QuasR_fq_col_found <- QuasR_fq_col[QuasR_fq_col_flag]
# Report wrong input column
if(length(QuasR_fq_col_found) != 1){stop("Incorrect sample table column names")}
# Extract input columns for QuasR and save
QuasR_input_path <- paste0(config$output_dir, "/", config$project_name, "_QuasR_input.tsv")
QuasR_df <- dplyr::select(sample_df, all_of(QuasR_fq_col_found[[1]]), SampleName=sample)
# Get absolute path of FASTQ file
for(fq_col in unlist(QuasR_fq_col)){
  if(fq_col %in% colnames(QuasR_df)){
    QuasR_df[[fq_col]] <- sapply(QuasR_df[[fq_col]],
                                 tools::file_path_as_absolute)
  }
}
write.table(QuasR_df, file=QuasR_input_path, sep="\t",
            row.names=FALSE, col.names=TRUE, quote=FALSE)

## Chunk 3 ----
# Setup CPU cluster for parallel run
cl <- NULL # placeholder for non-parallel run
avail_cores <- parallel::detectCores()
if((avail_cores > 1) && (config$n_core > 1)){
  use_core <- min(config$n_core, avail_cores, na.rm=TRUE)
  cl <- parallel::makeCluster(use_core, type = "FORK")
}
# Alignment with QuasR
QuasR_proj <- QuasR::qAlign(
  sampleFile = QuasR_input_path,
  genome = config$ref_fasta_path,
  aligner = "Rbowtie",
  projectName = config$project_name,
  alignmentsDir = QuasR_dir,
  alignmentParameter = NULL,
  splicedAlignment = FALSE,
  cacheDir = paste0(QuasR_dir, "/cache"),
  clObj = cl)
# Making reference for read count
ref_fa <- Biostrings::readDNAStringSet(config$ref_fasta_path)
ref_gr <- GRanges(seqnames = names(ref_fa), ranges = IRanges(start =
                                                               1, end = width(ref_fa)))
names(ref_gr) <- as.character(seqnames(ref_gr))
# Count read mapped to construct
count_matrix <- QuasR::qCount(proj = QuasR_proj, query = ref_gr,
                              clObj = cl)
count_matrix <- count_matrix[,which(colnames(count_matrix)!="width")]
if(!is.null(cl) && inherits(cl, "cluster")) parallel::stopCluster(cl)


## Chunk 4 ----
count_matrix_path <- paste0(config$output_dir, "/",
                            config$project_name, "_count_matrix.txt")
write.table(count_matrix, file=count_matrix_path, sep="\t",
            row.names=TRUE, col.names=TRUE, quote=FALSE)


## Tmp chunk ----
t1 <- Sys.time()
# count_matrix <- count_matrix[rownames(count_matrix)!="WT601_bp000" , ]


## Chunk 5 ----
# Count unique value of each column per 'group'
group_col_unq_vals <- sapply(
  split(sample_df, f=sample_df$group),
  FUN=function(df){
    sapply(df, function(x){length(unique(x))})
  })
invalid_group_df_colnames <-
  rownames(group_col_unq_vals)[apply(group_col_unq_vals > 1, MARGIN=1,
                                     FUN=any)]
if(length(invalid_group_df_colnames) > 0){message("These columns will
be removed from the group metadata table:\n ", paste0(invalid_group_df_colnames, collapse=", "))}
# Construct group metadata table
group_df <- sample_df %>%
  dplyr::select(-all_of(invalid_group_df_colnames)) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames(var="group")


## Chunk 6 ----
fraction_names <- unique(sample_df$fraction)
group_name <- unique(sample_df$group)
# placeholder matrix for empty count
zeros_mat <- matrix(0, nrow=nrow(count_matrix),
                    ncol=length(group_name), dimnames=list(rownames(count_matrix),
                                                           group_name))
# Collecting reads
count_mat_fraction <- list()
for(i in seq_along(fraction_names)){
  cur_name <- paste0("counts_", fraction_names[i])
  sample_group_map <- sample_df %>%
    dplyr::filter(fraction == fraction_names[i]) %>%
    dplyr::pull(group, name=sample)
  m <- count_matrix[ , names(sample_group_map), drop=FALSE]
  colnames(m) <- sample_group_map[colnames(m)]
  missing_samples <- setdiff(group_name, colnames(m))
  if(length(missing_samples) > 0){
    # Add missing samples with zero counts
    m <- cbind(m, zeros_mat[ , missing_samples, drop=FALSE])
  }
  m <- m + config$log2_pseudo_count
  # reorder the columns to match with group_name
  count_mat_fraction[[cur_name]] <- m[ , group_name, drop=FALSE]
}


## Chunk 7 ----
ref_df <- read.table(config$ref_table_path, sep="\t", header=TRUE,
                     stringsAsFactors=FALSE, row.names="ref_name")
ref_df <- ref_df[rownames(count_matrix), ]
SeEN_se <- SummarizedExperiment(assays=count_mat_fraction, colData=group_df[group_name , ], rowData=ref_df[rownames(count_matrix), ])


# Chunk 8 ----
for(i in seq_along(fraction_names)){
  m <- assay(SeEN_se, paste0("counts_", fraction_names[i]))
  lib_size_mat <- matrix(colSums(m), nrow=nrow(m), ncol=ncol(m),
                         byrow=TRUE)
  cur_assay <- paste0("libnorm_counts_", fraction_names[i])
  assay(SeEN_se, cur_assay) <- m/lib_size_mat
}


## Chunk 9 ----
for(i in seq_along(seq_along(config$enrichment))){
  frac_test <- config$enrichment[[i]][1]
  frac_ctrl <- config$enrichment[[i]][2]
  cur_assay<- paste0("enrichment_", frac_test, "_vs_", frac_ctrl)
  test_m <- assay(SeEN_se, paste0("libnorm_counts_", frac_test))
  ctrl_m <- assay(SeEN_se, paste0("libnorm_counts_", frac_ctrl))
  assay(SeEN_se, cur_assay) <- log2(test_m/ctrl_m)
}

## Chunk 10 ----
# Extract enrichment matrix and transform to tidy format
tidy_enrich_df <- assay(SeEN_se, "enrichment_Bound_vs_Unbound") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="ref_name") %>%
  tidyr::gather(key="group", value="enrichment", -ref_name)
# Adding metadata columns from group and reference table
tidy_enrich_df2 <- tidy_enrich_df %>%
  dplyr::left_join(tibble::rownames_to_column(group_df, var="group"),
                   by="group") %>%
  dplyr::left_join(tibble::rownames_to_column(ref_df,var="ref_name"),
                   by="ref_name")

# Calculate average enrichment score across replicates
enrich_df <- tidy_enrich_df2 %>%
  group_by(position, condition) %>%
  dplyr::reframe(
    avg_enrichment = mean(enrichment),
    sd_enrichment = sd(enrichment),
    sd_enrichment = replace(sd_enrichment, is.na(sd_enrichment), 0)
  )


## Chunk 11 ----
x_breaks <- c(1, seq(0, max(enrich_df$position, na.rm=T)+10, by=10)[-1])
p <- enrich_df %>%
  dplyr::filter(position > 0) %>% # ignore the construct with no motif
  ggplot(aes(x=position, y=avg_enrichment, group=condition)) +
  geom_hline(yintercept=0, color="black") +
  geom_ribbon(aes(ymin=avg_enrichment-sd_enrichment,
                  ymax=avg_enrichment+sd_enrichment,
                  fill=condition), alpha=0.3) +
  geom_line(aes(color=condition)) +
  geom_point(aes(color=condition), size=1) +
  scale_x_continuous(breaks=x_breaks) +
  scale_y_continuous(limits=c(NA, 7.5), breaks=seq(-0, 7.5, by=2.5)) +
  facet_grid(condition~.) +
  labs(x="Position", y="Enrichment Score") +
  theme_classic() +
  theme(legend.position="bottom", 
        panel.grid.major.x=element_line(color="gray", linewidth=0.5),
        strip.text.y=element_text(size=12, angle=0),
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background.y=element_blank())
ggsave(plot=p, filename=paste0(plot_dir, "/", config$project_name, "_enrichment_plot.pdf"), width=8, height=5, dpi=300, units="in")


## Chunk 12 ----
# RDS object for future reanalysis
saveRDS(SeEN_se, file=paste0(config$output_dir, "/", config$project_name, "_SE.rds"))
# Save to Excel file
result_list <- c(
  list(raw_count =
         tibble::rownames_to_column(as.data.frame(count_matrix), var="sample")),
  lapply(assays(SeEN_se),
         FUN=function(m){tibble::rownames_to_column(as.data.frame(m), var="sample")}),
  list(
    ref_metadata = rownames_to_column(as.data.frame(rowData(SeEN_se)), var="ref_name"),
    sampl_group_metadata = rownames_to_column(as.data.frame(colData(SeEN_se)), var="ref_name"), 
    sample_metadata = sample_df
  ))
writexl::write_xlsx(result_list, path=paste0(config$output_dir, "/",config$project_name, "_results.xlsx"))

t2 <- Sys.time()

