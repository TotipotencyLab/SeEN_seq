

make_SeENseq_templates <- function(template, insert, step_size, insert_start_nt=1, insert_end_nt=Inf, 
                                   adaptor_5="", adaptor_3="", 
                                   warn_if=character(0), 
                                   show_output=FALSE, outfile=NULL, quiet=FALSE){
  # Input
  #  template <character, length=1>
  #    DNA template for the sequence to be inserted
  #  insert <character, length=1>
  #    Inserting sequence (e.g., transcription factor binding motif sequence)
  #    NB: Consider making this accept the PWM motif, or multiple sequences in the future?
  #  step_size <integer/numeric, length=1>
  #    Step size (in base pair) of motif insertion sliding
  #  adaptor_5, adaptor_3 <character, length=1> [Optional]
  #    Adaptor sequence to be added to the 5' and 3' end of each designed oligo, respectively.
  #  watch_list <character vector> [Optional]
  #    List of sequences to watch for. For example, restriction sites of restriction enzymes, motifs of other transcription factors, etc.
  #  quiet <logical, length=1>
  
  # ASCII Color Code for console printing of each sequence types
  # ref: https://gist.github.com/JBlond/2fea43a3049b38287e5e9cefc87b2124
  fmt=list(
    reset="\033[0m",
    adpt = "\033[45m\033[1;30m", # purple background, black foreground
    tmpl = "\033[46m\033[1;30m", # cyan background, black foreground
    inst = "\033[43m\033[1;30m",  # yellow background, black foreground (bold)
    warn = "\033[41m\033[1;30m" # Red background, Black foreground
  )
  
  format_print_txt <- function(txt, ascii_format=NULL){
    txt <- paste0(txt[!is.na(txt)], collapse="")
    if(is.null(ascii_format)) return(txt)
    if(length(txt) > 0 && txt != ""){
      txt <- paste0(ascii_format, txt, "\033[0m") # Reset format at the end
    }
    return(txt)
  }
  
  require(stringr)
  # require(dplyr)
  
  template_nt <- base::strsplit(template, split="")[[1]]
  template_length <- length(template_nt)
  insert_nt <- base::strsplit(insert, split="")[[1]]
  insert_length <- length(insert_nt)
  insert_start <- seq(insert_start_nt, min(insert_end_nt, (length(template_nt)-insert_length+1)), by=step_size)
  
  RC_map <- c(A="T", C="G", G="C", `T`="A") # Reverse complement dictionary
  
  out_seq <- c()
  out_print <- c()
  
  seq_df <- data.frame()
  for(i in seq_along(insert_start)){
    # Get the sequence before and after the insert sites.
    tmpl_1 <- template_nt[0:(insert_start[i]-1)]
    insert_end <- insert_start[i] + insert_length - 1
    if(insert_end < template_length){
      tmpl_2 <- template_nt[(insert_end+1):template_length]
    }else if(insert_end == template_length){
      tmpl_2 <- character(0)
    }else{
      tmpl_2 <- character(0)
      warning("Insertion generate sequence longer than the template.")
    }
    # tmpl_2 <- template_nt[insert_end:min(insert_end, template_length)]
    
    seq_body <- toupper(paste0(c(tmpl_1, insert, tmpl_2), collapse=""))
    warn=FALSE # placeholder
    for(w in seq_along(warn_if)){
      cur_warn_name <- names(warn_if)[w]
      cur_warn_seq <- toupper(warn_if[w])
      # Consider using BioString for this?
      cur_warn_seq_rc <- paste0(RC_map[strsplit(cur_warn_seq, split="")[[1]]], collapse="") # Only works with ATCG nucleotide
      if(grepl(seq_body, pattern=cur_warn_seq)){
        warn_pos <- stringr::str_locate_all(toupper(seq_body), pattern=cur_warn_seq)[[1]][ , "start"]
        message("Sequence from watch list (", w, " - ", cur_warn_seq, " - ", cur_warn_name, ") is detected in insert position ", i, "\n  ",
                "-- found position: ", paste0(warn_pos, collapse=", "))
        warn=TRUE
        if(show_output) Sys.sleep(1)
      }
      # Reverse complement
      if(grepl(seq_body, pattern=cur_warn_seq_rc)){
        warn_pos <- stringr::str_locate_all(seq_body, pattern=cur_warn_seq_rc)[[1]][ , "start"]
        message("Sequence from watch list (", w, " - ", cur_warn_seq_rc, " - ", cur_warn_name, " [RevComp]) is detected in insert position ", i, "\n  ",
                "-- found position: ", paste0(warn_pos, collapse=", "))
        warn=TRUE
        if(show_output) Sys.sleep(1)
      }
    }
    
    out_seq[i] <- paste0(c(adaptor_5, seq_body, adaptor_3), collapse="")
    print_nt <- c(insert_start[i], " - ",
      format_print_txt(adaptor_5, ascii_format=fmt$adpt),
      format_print_txt(tmpl_1, ascii_format=fmt$tmpl),
      format_print_txt(insert, ascii_format=fmt$inst),
      format_print_txt(tmpl_2, ascii_format=fmt$tmpl),
      format_print_txt(adaptor_3, ascii_format=fmt$adpt),
      character(0)
    )
    out_print[i] <- paste0(print_nt, collapse="")
    if(show_output) cat(out_print[i], "\n", sep="")
    
    df <- data.frame(insert_start = insert_start[i], sequence=out_seq[i], warning=warn)
    seq_df <- rbind(seq_df, df)
  }
  
  # Formatting output
  out_ext <- ifelse(!is.null(outfile), yes=tolower(tools::file_ext(outfile)), no="")
  if(out_ext %in% c("txt", "tsv")){
    write.table(seq_df, file=outfile, sep="\t", row.names=FALSE, col.names=TRUE,quote=FALSE)
  }else if(out_ext %in% c("csv")){
    write.csv(seq_df, file=outfile, row.names=FALSE, col.names=TRUE,quote=FALSE)
  }else{
    return(seq_df)
  }
}



# Example ============================================================================================

if(F){
  W601_seq = "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT"
  ESRRB_motif = "TCAAGGCCA"
  step_size = 5
  adaptor_5 = "tatctcttgtagccgatgtgctaatggccatgaacgagcgtgtggAAGCATGGACTTACAGCTATTGCGTGATATC"
  adaptor_3 = "GATATCCATGGTTACCAACCCTATGGTCACTtattcgtcacgacaacaggacatccaatgttctgcacaaggctttc"
  watch_list = c(EcoRI="GAATTC")
  show_output=FALSE
  quiet=FALSE
  
  seq_df <- make_SeENseq_templates(template=W601_seq, insert=ESRRB_motif, step_size=1, adaptor_5=adaptor_5, adaptor_3=adaptor_3, warn_if=watch_list)
  seq_df <- make_SeENseq_templates(template=W601_seq, insert=ESRRB_motif, step_size=1, adaptor_5=adaptor_5, adaptor_3=adaptor_3, warn_if=watch_list, show_output=TRUE)
  make_SeENseq_templates(template=W601_seq, insert=ESRRB_motif, step_size=10, adaptor_5=adaptor_5, adaptor_3=adaptor_3, warn_if=watch_list, show_output=TRUE, 
                         outfile=paste0("./tests/SeEN_seq_template_seq.txt"))
  
  # Manual inspect
  str_view_all(seq_df$sequence[10], pattern=watch_list)
  stringr::str_locate_all(toupper(seq_df$sequence[10]), pattern=toupper(watch_list))[[1]][ , "start"]
}

