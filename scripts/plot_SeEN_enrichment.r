SeEN_se.plot_enrichment <- function(SeEN_se, x="motif_position", y=NULL, y_regex="^enrichment_", group_by="sample", color_by=NULL, 
                                    NA_replace.x = -1, NA_replace.y = NULL,
                                    show_plot_as = c("line"), show_group_avg_as = c("ribbon"),
                                    average_group=FALSE, facet_formula=NULL, 
                                    plot_data=FALSE){
  # Inputs:
  # SeEN_se <SummarizedExperiment>: SummarizedExperiment object with enrichment data
  # x <character>: one of the column names of the row data (reference info) of SeEN_se
  # y <character>: An explicit way to specify which assay to use for plotting (see `assayNames(SeEN_se)` for available assays)
  # y_regex <character-regular expression>: will be used if `y` is NULL; guessing which assay to be used
  # group_by <character>: one of the column names of the col data (sample info) of SeEN_se
  # color_by <character>: one of the column names of the col data (sample info) of SeEN_se
  # average_group <logical>: whether the value should be average within the group
  
  # Extract row and column data from SE object
  r_df <- as.data.frame(rowData(SeEN_se)) %>% 
    rownames_to_column(var="ref_name")
  c_df <- as.data.frame(colData(SeEN_se)) %>% 
    rownames_to_column(var="sample")
  
  # Process aesthetic variables ===============================================
  if(is.null(group_by)) group_by = character(0)
  if(is.null(color_by)) color_by = character(0)
  
  # Guess what assay to use if user didn't provide us this info
  if(is.null(y)){
    message("The input for `y` is not given, guessing what assay to use ... ")
    y_candidates <- assayNames(SeEN_se)[grepl(assayNames(SeEN_se), pattern=y_regex)]
    if(length(y_candidates) == 0){
      stop("No suitable assays for plotting. Please specify `y` as of of the `assayNames(SeEN_se)`")
    }else if(length(y_candidates) > 1){
      warning("Multiple candidate assays found! Only using the first one.", 
              "\n\tUnused assays: ", paste0(y_candidates[-1], collapse=", "))
    }
    y = y_candidates[1]
    message("Using assay '", y, "' as y-axis")
  }else{
    if(!y %in% assayNames(SeEN_se)){
      stop("The given y (", y, ") is not part of the available assays in `SeEN_se` input. Please make sure it's part of assayNames(SeEN_se)")
    }
  }
  
  # Extract data ==============================================================
  expr_df <- assay(SeEN_se, y) %>% 
    as.data.frame() %>% 
    rownames_to_column(var="ref_name")
  plot_df <- expr_df %>% 
    gather(key="sample", value="value", -ref_name) %>% 
    left_join(c_df, by="sample") %>%
    left_join(r_df, by="ref_name") %>% 
    # replace_na(list(position=-1)) %>%  # For row data
    as_tibble()
  
  # Replace NA value of x and y axes 
  # NOTE: This was added because the WT construct that has no motif position is marked as NA in motif_position (x-axis) of the row data
  na_replace_list <- list()
  if(!is.null(NA_replace.x)){
    na_replace_list[[x]] <- NA_replace.x
  }
  
  if(!is.null(NA_replace.y)){
    na_replace_list$value <- NA_replace.y
  }
  
  if(length(na_replace_list) > 0){
    require(tidyr)
    plot_df <- replace_na(plot_df, replace=na_replace_list)
  }
  
  # Average over group
  if(average_group){
    # require(glue)
    group_colnames <- unique(c(group_by, x, color_by))
    cat("Average data over these columns:", paste0(group_colnames, collapse=", "), "\n")
    plot_df <- plot_df %>% 
      # replace_na(setNames(list(-1), nm=x)) %>% 
      group_by(across(all_of(group_colnames))) %>% 
      reframe(
        sd = ifelse(length(value) > 1, sd(value, na.rm=TRUE), 0),
        value = mean(value, na.rm=TRUE)
      )
  }
  
  # Check duplicated X values
  dup_x <- plot_df %>% 
    group_by(across(all_of(unique(c(x, group_by, color_by))))) %>% 
    reframe(dup_x=n() > 1) %>% 
    dplyr::pull(dup_x) %>% 
    any()
  
  if(dup_x){
    warning("The plot may contain duplicated data point on X-axis within the same group of data. Please consider checking the group_by and color_by arguments.")
  }
  
  if(plot_data){
    return(plot_df)
  }
  
  # Making a plot =============================================================
  ## Plot type inputs ---------------------------------------------------------
  # NOTE: This looks complicated, but hopefully it would allow us to add new feature 
  #       more easily by just editing in default_plot_type and avail_plot_type variables.
  
  check_plot_option <- function(input, options){
    sapply(options, function(x) any(grepl(input, pattern=x, ignore.case=TRUE)))
  }
  
  # Main data plot type
  show_plot_as <- tolower(show_plot_as)
  default_plot_type <- c("line")
  avail_plot_type <- c(point="point", line="line") # The value will be used as regex
  show_plot_setting <- check_plot_option(input=show_plot_as, options=avail_plot_type)
  
  if(sum(show_plot_setting) == 0){
    warning("The input 'show_plot_as' doesn't matched any of the accepted options: ", paste0(avail_plot_type, collapse=", "), "\nSwtiching to default setting")
    show_plot_as <- default_plot_type
    show_plot_setting <- check_plot_option(input=show_plot_as, options=avail_plot_type)
  }
  
  # Ribbon (i.e., group average)
  show_group_avg_as <- tolower(show_group_avg_as)
  default_group_avg_type <- c("ribbon")
  avail_group_avg_type <- c(ribbon="ribbon", errorbar="error_?bar") # The value will be used as regex
  show_group_avg_setting <- check_plot_option(input=show_group_avg_as, options=avail_group_avg_type)
  
  if(sum(show_group_avg_setting) == 0){
    warning("The input 'show_plot_as' doesn't matched any of the accepted options: ", paste0(avail_plot_type, collapse=", "), "\nSwtiching to default setting")
    show_group_avg_as <- default_plot_type
    show_group_avg_setting <- check_plot_option(input=show_group_avg_as, options=avail_group_avg_type)
  }
  
  
  ## Plot Asthetic ------------------------------------------------------------
  ## Dynamically build plot aesthetic mapping
  # global aes
  p_aes <- aes(x=.data[[x]], y=value)
  # Element aes
  # line_aes <- aes()
  sd_aes <- aes(ymin=value-sd, ymax=value+sd)
  if(length(group_by) > 0) p_aes$group <- substitute(.data[[group_by]])
  if(length(color_by) > 0){
    p_aes$color <- substitute(.data[[color_by]])
    # p_aes$fill <- substitute(.data[[color_by]])
    # line_aes$color <- substitute(.data[[color_by]])
    sd_aes$fill <- substitute(.data[[color_by]]) # For ribbon fill
  }
  
  # Plotting ------------------------------------------------------------------
  p <- plot_df %>% 
    ggplot(mapping=p_aes) + 
    geom_hline(yintercept = 0, linetype="dashed", color="grey25") +
    labs(y=y) +
    theme_bw()
  
  if(show_plot_setting["line"]){
    p <- p + geom_line(alpha=0.7)
  }
  if(show_plot_setting["point"]){
    p <- p + geom_point(size=0.75, alpha=0.7)
  }
  
  if(average_group){
    if(show_group_avg_setting["ribbon"]){
      p <- p + geom_ribbon(mapping=sd_aes, color=NA, alpha=0.5)
    }
    if(show_group_avg_setting["errorbar"]){
      p <- p + geom_errorbar(mapping=sd_aes, alpha=0.5)
    }
  }
  
  if(is.character(facet_formula)){
    try({facet_formula <- as.formula(facet_formula)})
  }
  if(is(facet_formula, "formula")){
    p <- p +
      facet_grid(facet_formula) +
      theme(strip.text.y=element_text(angle=0), 
            strip.background.y=element_rect(fill=NA, color=NA))
  }
  
  return(p)
}