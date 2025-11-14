extract_param_cols <- function(colnames_vec, param_names) {
  out <- list()
  
  for (nm in param_names) {
    # pattern: start with 'nm' and then optional "[" or end of string
    pattern <- paste0("^", nm, "(\\[|$)")
    out[[nm]] <- grep(pattern, colnames_vec)
  }
  
  out
}