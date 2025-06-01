# Write SeEN-seq associated data to Excel files

SeEN_se.write_xlsx <- function(SeEN_se, file){
  require(dplyr)
  require(tibble)
  require(writexl)
  require(SummarizedExperiment)
  
  out_list <- list()
  r_df <- as.data.frame(rowData(SeEN_se)) %>% 
    rownames_to_column(var="ref_name")
  c_df <- as.data.frame(colData(SeEN_se)) %>% 
    rownames_to_column(var="sample")
  out_list[["reference"]] <- r_df
  out_list[["sample"]] <- c_df
  
  assay_list <- as.list(assays(SeEN_se))
  assay_list <- lapply(assay_list, function(m){
    as.data.frame(m) %>% 
      rownames_to_column(var="ref_name") # Make this consistent with row data
  })
  
  out_list <- c(out_list, assay_list)
  # str(out_list, max.level=1)
  
  write_xlsx(out_list, path=file)
}
