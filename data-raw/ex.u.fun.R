ex.u.fun <- function(O, DATA){
  RESULT <- NULL
  RESULT$Q <- -O
  RESULT$E <- O
  return(RESULT)
}
usethis::use_data(ex.u.fun, overwrite = TRUE)
