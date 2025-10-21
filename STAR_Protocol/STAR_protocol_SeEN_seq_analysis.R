# From text
## MARK: c1
# In R console
library("tools")
library("tidyverse")
library("parallel")
library("QuasR")
library("Biostrings")
library("GenomicRanges")
library("SummarizedExperiment")
library("writexl")

# required_pkg <- c("tidyverse", "parallel", "QuasR", "SummarizedExperiment", "ggplot2")
# sapply(required_pkg, library, character.only=TRUE)

## MARK:  C2: setup
# General analysis project setting
wd <- getwd() # assume we are in the GitHub directory
config <- list(
  project_name = "Test_SeEN_seq",
  output_dir = "./results/Test_SeEN_seq",
  sample_table_path = "./example/sample_metadata_PE.tsv",
  ref_table_path = "./example/ref_metadata.tsv",
  ref_fasta_path = "./example/ref_sequence.fasta",
  n_core=8,
  enrichment=list(c("Bound", "Unbound")), # test vs background
  log2_pseudo_count=1
)
# Make sure output folder exist
dir.create(config$output_dir, showWarnings=FALSE, recursive=TRUE)
# Make sure the paths in config are absolute
for(i in seq_along(config)){
  if(grepl("(_path$)|(_dir$)", names(config)[i])){
    config[[i]] <- file_path_as_absolute(config[[i]])
  }
}
QuasR_dir <- paste0(config$output_dir, "/QuasR")
plot_dir <- paste0(config$output_dir, "/plot")
dir.create(QuasR_dir, showWarnings=FALSE)
dir.create(paste0(QuasR_dir, "/cache"), showWarnings=FALSE)
dir.create(plot_dir, showWarnings=FALSE)




## MARK: C3
setwd(wd)
sample_df <- read.table(config$sample_table_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
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
write.table(QuasR_df, file=QuasR_input_path, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



## MARK: C3 -- SE
# sample_df <- read.table(str_replace(config$sample_table_path, "_PE", "_SE"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# # Extract relavant columns for QuasR
# QuasR_fq_col <- list(SE = "FileName", PE = c("FileName1", "FileName2"))
# QuasR_fq_col_flag <- sapply(QuasR_fq_col, FUN = function(x) {
#   all(x %in% colnames(sample_df))
# })
# QuasR_fq_col_found <- QuasR_fq_col[QuasR_fq_col_flag]
# # Report wrong input column
# if (length(QuasR_fq_col_found) != 1) {
#   stop("Incorrect sample table column names")
# }
# # Extract input columns for QuasR and save
# QuasR_input_path <- paste0(config$output_dir, "/", config$project_name, "_QuasR_input.tsv")
# QuasR_df <- dplyr::select(sample_df, all_of(QuasR_fq_col_found[[1]]), SampleName = sample)
# # Get absolute path of FASTQ file
# for (fq_col in unlist(QuasR_fq_col)) {
#   if (fq_col %in% colnames(QuasR_df)) {
#     QuasR_df[[fq_col]] <- sapply(
#       QuasR_df[[fq_col]],
#       tools::file_path_as_absolute
#     )
#   }
# }
# write.table(QuasR_df,
#   file = QuasR_input_path, sep = "\t",
#   row.names = FALSE, col.names = TRUE, quote = FALSE
# )


## MARK: C4
# Setup CPU cluster for parallel run
cl <- NULL # placeholder for non-parallel run
avail_cores <- parallel::detectCores()
if((avail_cores > 1) && (config$n_core > 1)){
  use_core <- min(config$n_core, avail_cores, na.rm=TRUE)
  cl <- parallel::makeCluster(use_core, type = "FORK") 
}
# Alignment with QuasR
setwd(QuasR_dir)
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
ref_gr <- GRanges(seqnames = names(ref_fa), ranges = IRanges(start = 1, end = width(ref_fa)))
names(ref_gr) <- as.character(seqnames(ref_gr))
# Count read mapped to construct
count_matrix <- QuasR::qCount(proj = QuasR_proj, query = ref_gr, clObj = cl)
count_matrix <- count_matrix[,which(colnames(count_matrix)!="width")] 
if(!is.null(cl) && inherits(cl, "cluster")) parallel::stopCluster(cl)


## MARK: C5
setwd(wd)
count_matrix_path <- paste0(config$output_dir, "/", config$project_name, "_count_matrix.txt")
write.table(count_matrix, file=count_matrix_path, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)


## MARK: C6
# Count unique value of each column per 'lane'
lane_col_unq_vals <- sapply(
  split(sample_df, f=sample_df$lane), 
  FUN=function(df){
    sapply(df, function(x){length(unique(x))})
  })
invalid_lane_df_colnames <- rownames(lane_col_unq_vals)[apply(lane_col_unq_vals > 1, MARGIN=1, FUN=any)]
if(length(invalid_lane_df_colnames) > 0){message("These columns will be removed from the lane metadata table:\n  ", paste0(invalid_lane_df_colnames, collapse=", "))}
# Construct lane metadata table
lane_df <- sample_df %>% 
  dplyr::select(-all_of(invalid_lane_df_colnames)) %>% 
  dplyr::distinct() %>% 
  tibble::column_to_rownames(var="lane")


## MARK: C7
fraction_names <- unique(sample_df$fraction)
lane_name <- unique(sample_df$lane)
# placeholder matrix for empty count
zeros_mat <- matrix(0, nrow=nrow(count_matrix), ncol=length(lane_name), dimnames=list(rownames(count_matrix), lane_name))
# Collecting reads
count_mat_fraction <- list()
for(i in seq_along(fraction_names)){
  cur_name <- paste0("counts_", fraction_names[i])
  sample_lane_map <- sample_df %>% 
    dplyr::filter(fraction == fraction_names[i]) %>% 
    dplyr::pull(lane, name=sample)
  m <- count_matrix[ , names(sample_lane_map), drop=FALSE]
  colnames(m) <- sample_lane_map[colnames(m)]
  missing_samples <- setdiff(lane_name, colnames(m))
  if(length(missing_samples) > 0){
    # Add missing samples with zero counts
    m <- cbind(m, zeros_mat[ , missing_samples, drop=FALSE])
  }
  m <- m + config$log2_pseudo_count
  # reorder the columns to match with lane_name
  count_mat_fraction[[cur_name]] <- m[ , lane_name, drop=FALSE] 
}


## MARK: C8
ref_df <- read.table(config$ref_table_path, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names="ref_name")
ref_df <- ref_df[rownames(count_matrix), ]
SeEN_se <- SummarizedExperiment(assays=count_mat_fraction, colData=lane_df[lane_name , ], rowData=ref_df[rownames(count_matrix), ])



## MARK: C9
for(i in seq_along(fraction_names)){
  m <- assay(SeEN_se, paste0("counts_", fraction_names[i]))
  lib_size_mat <- matrix(colSums(m), nrow=nrow(m), ncol=ncol(m), byrow=TRUE)
  cur_assay <- paste0("libnorm_counts_", fraction_names[i])
  assay(SeEN_se, cur_assay) <- m/lib_size_mat
}


## MARK: C10
for(i in seq_along(seq_along(config$enrichment))){
  frac_test <- config$enrichment[[i]][1]
  frac_ctrl <- config$enrichment[[i]][2]
  cur_assay<- paste0("enrichment_", frac_test, "_vs_", frac_ctrl)
  test_m <- assay(SeEN_se, paste0("libnorm_counts_", frac_test))
  ctrl_m <- assay(SeEN_se, paste0("libnorm_counts_", frac_ctrl))
  assay(SeEN_se, cur_assay) <- log2(test_m/ctrl_m)
}

## MARK: C11
SeEN_se_2_tidy <- function(se, assay_name, value_colname=NULL){
  if (is.null(value_colname)) value_colname=assay_name
  tidy_df <- as.data.frame(assay(se, assay_name)) %>% 
    rownames_to_column(var="ref_name") %>% 
    gather(key="lane", value=!!sym(value_colname), -ref_name) %>%
    # Adding metadata columns from lane and reference table
    left_join(rownames_to_column(as.data.frame(colData(se)), var="lane"), by="lane") %>%
    left_join(rownames_to_column(as.data.frame(rowData(se)), var="ref_name"), by="ref_name")
  return(tidy_df)
}

## MARK: C12: 
# Extract raw counts from unbound fraction of lane with no TF added
unbound_df <- SeEN_se_2_tidy(SeEN_se[, colData(SeEN_se)$Protein_conc_uM == 0], assay_name="counts_Unbound", value_colname="count")
# Detecting outlier
unbound_df <- unbound_df %>%
  dplyr::group_by(lane) %>%
  dplyr::mutate(
    outlier = count %in% boxplot.stats(count)$out,
    label = dplyr::if_else(outlier, true="Outlier", false="Normal"),
    outlier_label = dplyr::if_else(outlier, true=position, false=NA)
  ) %>% dplyr::ungroup()
# Show under/over-represented construct for each lane
print(dplyr::filter(unbound_df, outlier)[, c("lane", "position", "count")])

## MARK: C13
x_breaks <- c(1, seq(0, max(unbound_df$position, na.rm = T) + 10, by = 10)[-1])
p_unbound <- unbound_df %>%
  ggplot(aes(x = position, y = count, group = position)) +
  geom_hline(yintercept = 0 + config$log2_pseudo_count) +
  geom_col(aes(fill = label), alpha = 0.5) +
  geom_point(aes(color = label)) +
  geom_text(aes(label = outlier_label), vjust = -1) +
  scale_y_log10() + scale_x_continuous(breaks=x_breaks) +
  facet_grid(lane ~ ., scales = "free_y") +
  labs(title = "Raw reads from unbound fraction without TF") +
  theme_bw()
ggsave(plot = p_unbound, filename = paste0(plot_dir, "/", config$project_name, "_unbound_count_plot.pdf"), width = 8, height = 5, dpi = 300, units = "in")


## MARK: C14
# Extract enrichment matrix and transform to tidy format
tidy_enrich_df <- SeEN_se_2_tidy(SeEN_se, assay_name = "enrichment_Bound_vs_Unbound", value_colname = "enrichment")
# Calculate average enrichment score across replicates
enrich_df <- tidy_enrich_df %>%
  group_by(position, condition) %>% 
  dplyr::reframe(
    avg_enrichment = mean(enrichment),
    sd_enrichment = sd(enrichment),
    sd_enrichment = replace(sd_enrichment, is.na(sd_enrichment), 0)
  )


## MARK: C15
p <- enrich_df %>%
  dplyr::filter(position > 0) %>% # ignore the construct with no motif
  ggplot(aes(x=position, y=avg_enrichment, group=condition)) +
  geom_hline(yintercept=0, color="black") +
  geom_ribbon(aes(ymin=avg_enrichment-sd_enrichment,
                  ymax=avg_enrichment+sd_enrichment,
                  fill=condition), alpha=0.3) +
  geom_line(aes(color=condition)) +
  geom_point(aes(color=condition)) +
  scale_x_continuous(breaks=x_breaks) +
  scale_y_continuous(limits=c(NA, 7.5), breaks=seq(-0, 7.5, by=2.5)) +
  facet_grid(condition~.) +
  labs(x="Position", y="Enrichment Score") +
  theme_classic() +
  theme(legend.position="bottom", 
        strip.text.y=element_text(size=12, angle=0),
        strip.background.y=element_blank())
ggsave(plot=p, filename=paste0(plot_dir, "/", config$project_name, "_enrichment_plot.pdf"), width=8, height=5, dpi=300, units="in")


## MARK: C16
# RDS object for future reanalysis
saveRDS(SeEN_se, file=paste0(config$output_dir, "/", config$project_name, "_SE.rds"))
# Save to Excel file
result_list <- c(
  list(raw_count = tibble::rownames_to_column(as.data.frame(count_matrix), var="sample")),
  lapply(assays(SeEN_se), FUN=function(m){tibble::rownames_to_column(as.data.frame(m), var="sample")}),
  list(
    ref_metadata = rownames_to_column(as.data.frame(rowData(SeEN_se)), var="ref_name"),
    sampl_lane_metadata = rownames_to_column(as.data.frame(colData(SeEN_se)), var="ref_name"),
    sample_metadata = sample_df
))
writexl::write_xlsx(result_list, path=paste0(config$output_dir, "/", config$project_name, "_results.xlsx"))
