
#' Quick function to clean the right hand side variables dataset
#' 
#' @param rhs_dt `data.frame` object for the right hand side variables
#' 
#' @export
#' 
clean_rhs <- function(rhs_dt){
  
  ### let us reconstruct the column names 
  # first select the relevant columns
  relcols <- colnames(rhs_dt)[grepl(pattern = "^X\\d+|^\\.\\.\\.\\d+", 
                                    x = colnames(rhs_dt))]
  
  current_label <- NA
  
  new_labels <- character(length(relcols))
  
  for (i in seq_along(relcols)){
    
    if (grepl("^X\\d+$", relcols[i])) {
      current_label <- relcols[i]      # update current X-label
      new_labels[i] <- current_label
    } else if (grepl("^\\.\\.\\.\\d+$", relcols[i])) {
      new_labels[i] <- current_label  # assign to last seen X-label
    } else {
      new_labels[i] <- cols[i]  # for ID columns etc.
    }
    
  }
  
  relcols <- new_labels
  
  rowtwo_list <- 
    rhs_dt[1,] |> 
    as.list() |> 
    unlist() |> 
    unname() |>
    str_replace_na(replacement = "")
  
  new_labels <- c(colnames(rhs_dt)[!grepl(pattern = "^X\\d+|^\\.\\.\\.\\d+", 
                                          x = colnames(rhs_dt))],
                  new_labels)
  
  colnames(rhs_dt) <- 
    tibble(rowtwo_list,
           new_labels) |>
    mutate(y = ifelse(grepl(pattern = "^X\\d+", x = new_labels), "y", ""),
           varnames = paste0(new_labels, y, rowtwo_list)) |>
    dplyr::select(varnames) |>
    unlist() |>
    unname()
  
  rhs_dt <- rhs_dt[-1,]
  
  return(rhs_dt)
  
} 
