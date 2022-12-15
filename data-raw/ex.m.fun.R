ex.m.fun <- function(K, DATA){
  # Total number of erroneous contacts in the new set
  errCont <- sum(
    # Contacts not satisfying the threshold value of the gap between their two regions
    ( (K[,"j"]-K[,"i"]) <= DATA$gaplimit ) |
      # Contacts also present in the original set
      ( K[,"j"]==DATA$IJ_ORIG[,"j"]        ) |
      # Duplicated contacts in the new set
      ( duplicated(K)                      )
  )
  # Percentage of erroneous contacts
  O <- (errCont/DATA$numContacts)*100
  return(O)
}
usethis::use_data(ex.m.fun, overwrite = TRUE)
